#!/usr/bin/env python
"""
Repair pdb2gmx -missing side-chain H atoms left at (0, 0, 0).

Those atoms make minimization diverge (Fmax=inf) and grompp abort at NVT.

NOTE: this module is imported by the topology generator that runs inside the
GROMACS singularity image, whose interpreter is Python 2.7. Keep it strictly
Python 2/3 compatible (no f-strings, no type annotations, no PEP 585 generics).
"""

from collections import defaultdict

_ANCHOR_NAMES = ("CB", "CG", "CA", "N", "C", "O")
_OFFSETS = (
    (0.10, 0.00, 0.00),
    (-0.10, 0.00, 0.00),
    (0.00, 0.10, 0.00),
    (0.00, -0.10, 0.00),
    (0.00, 0.00, 0.10),
    (0.00, 0.00, -0.10),
)


def _parse_atom_line(line):
    return {
        "resnr": int(line[0:5]),
        "resname": line[5:10].strip(),
        "atomname": line[10:15].strip(),
        "atomnr": int(line[15:20]),
        "x": float(line[20:28]),
        "y": float(line[28:36]),
        "z": float(line[36:44]),
    }


def _format_atom_line(atom):
    return (
        "{0:5d}{1:>5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}".format(
            atom["resnr"],
            atom["resname"],
            atom["atomname"],
            atom["atomnr"],
            atom["x"],
            atom["y"],
            atom["z"],
        )
    )


def _is_zero(x, y, z, tol):
    return abs(x) + abs(y) + abs(z) < tol


def fix_gro_zero_coords(gro_path, tol=1e-3):
    """
    Move atoms at the origin onto their residue heavy-atom anchor.

    Returns:
        Number of atoms repaired.
    """
    with open(gro_path) as handle:
        lines = handle.read().splitlines()

    if len(lines) < 3:
        return 0

    natoms = int(lines[1].strip())
    header = lines[0]
    atom_lines = lines[2 : 2 + natoms]
    box_lines = lines[2 + natoms :]

    atoms = [_parse_atom_line(line) for line in atom_lines]
    by_res = defaultdict(list)
    for idx, atom in enumerate(atoms):
        by_res[atom["resnr"]].append((idx, atom))

    fixed = 0
    offset_idx = defaultdict(int)

    for idx, atom in enumerate(atoms):
        if not _is_zero(atom["x"], atom["y"], atom["z"], tol):
            continue

        anchor = None
        for name in _ANCHOR_NAMES:
            for _, candidate in by_res[atom["resnr"]]:
                if candidate["atomname"] == name and not _is_zero(
                    candidate["x"], candidate["y"], candidate["z"], tol
                ):
                    anchor = candidate
                    break
            if anchor is not None:
                break

        if anchor is None:
            for _, candidate in by_res[atom["resnr"]]:
                if candidate is atom:
                    continue
                if not _is_zero(candidate["x"], candidate["y"], candidate["z"], tol):
                    anchor = candidate
                    break

        if anchor is None:
            raise ValueError(
                "No anchor atom for zero-coordinate entry "
                "{0} {1} (atom {2}) in {3}".format(
                    atom["resname"], atom["atomname"], atom["atomnr"], gro_path
                )
            )

        k = offset_idx[atom["resnr"]]
        dx, dy, dz = _OFFSETS[k % len(_OFFSETS)]
        offset_idx[atom["resnr"]] += 1

        atom["x"] = anchor["x"] + dx
        atom["y"] = anchor["y"] + dy
        atom["z"] = anchor["z"] + dz
        fixed += 1

    if fixed == 0:
        return 0

    out_lines = [header, str(natoms)]
    out_lines.extend(_format_atom_line(atom) for atom in atoms)
    out_lines.extend(box_lines)
    with open(gro_path, "w") as handle:
        handle.write("\n".join(out_lines) + "\n")

    return fixed


def find_zero_coord_atoms(gro_path, tol=1e-3):
    with open(gro_path) as handle:
        lines = handle.read().splitlines()
    natoms = int(lines[1].strip())
    zeros = []
    for i, line in enumerate(lines[2 : 2 + natoms], start=1):
        atom = _parse_atom_line(line)
        if _is_zero(atom["x"], atom["y"], atom["z"], tol):
            zeros.append((i, "{0} {1}".format(atom["resname"], atom["atomname"])))
    return zeros
