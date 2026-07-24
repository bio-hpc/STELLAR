#!/usr/bin/env bash
# Debug: confirmar que este script se está ejecutando
echo "DEBUG create_pdb.sh: Script iniciado con $# argumentos"
echo "DEBUG create_pdb.sh: Argumentos recibidos: $@"

if [ "$#" -ne 5 ]; then
    echo "ERROR: The following parameters are needed:"
    echo "1º xtc"
    echo "2º gro"
    echo "3º out_pdb"
    echo "4 num step"
    echo "5º prefix_gromacs"
    exit
fi

xtc=$1
gro=$2
out_pdb=$3
step=$4
prefix_gromacs=$5

echo "DEBUG create_pdb.sh: step='$step', xtc='$xtc'"

#variable steps no se usa
if [ "$step" == "-1" ];then 
    # Valor por defecto
    DEFAULT_DUMP=99999999
    
    # Intentar extraer el número de frames del XTC usando múltiples métodos
    # Método 1: Buscar "Last frame" o "Last step"
    detected_dump=$(${prefix_gromacs} check -f ${xtc} 2>&1 | grep -iE "last.*frame|last.*step" | grep -oE "[0-9]+" | tail -1)
    
    # Método 2: Si no funciona, buscar "Step" y extraer el último número
    if [ -z "$detected_dump" ] || [ "$detected_dump" == "" ]; then
        detected_dump=$(${prefix_gromacs} check -f ${xtc} 2>&1 | grep -i "Step" | grep -oE "[0-9]+\.[0-9]+|[0-9]+" | tail -1 | cut -d. -f1)
    fi
    
    # Método 3: Buscar "frames" en la salida
    if [ -z "$detected_dump" ] || [ "$detected_dump" == "" ]; then
        detected_dump=$(${prefix_gromacs} check -f ${xtc} 2>&1 | grep -i "frames" | grep -oE "[0-9]+" | head -1)
    fi
    
    # Método 4: Buscar cualquier número grande después de "Step" o "frame"
    if [ -z "$detected_dump" ] || [ "$detected_dump" == "" ]; then
        check_output=$(${prefix_gromacs} check -f ${xtc} 2>&1)
        detected_dump=$(echo "$check_output" | grep -iE "step|frame" | grep -oE "[0-9]{4,}" | tail -1)
    fi
    
    # Verificar si mpi está definida, si no, usar cadena vacía
    if [ -z "${mpi}" ]; then
        mpi_cmd=""
    else
        mpi_cmd="${mpi}"
    fi
    
    # Determinar si usar -dump o no
    if [ ! -z "$detected_dump" ] && [ "$detected_dump" != "" ] && [[ "$detected_dump" =~ ^[0-9]+$ ]]; then
        # Se detectó un valor válido, usar -dump
        echo "INFO: Número de frames detectado: $detected_dump"
        echo 0 | ${prefix_gromacs} trjconv${mpi_cmd} -f ${xtc} -s ${gro} -o ${out_pdb} -dump ${detected_dump}
    else
        # No se detectó, ejecutar sin -dump (trjconv extraerá el último frame por defecto)
        echo "WARNING: No se pudo extraer el número de frames del XTC. Ejecutando sin -dump (extraerá último frame)"
        echo 0 | ${prefix_gromacs} trjconv${mpi_cmd} -f ${xtc} -s ${gro} -o ${out_pdb}
    fi
else
    echo "echo 0 |${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -dump $step"
    echo 0 |${prefix_gromacs} trjconv${mpi} -f ${xtc} -s ${gro} -o ${out_pdb} -dump $step
fi





#dump=$(echo $step_md*0.002 | bc) #0.002 es el paso
#execute "echo 0 |${prefix_gromacs} trjconv${mpi} -f ${out_molec}_complex_md.xtc -s ${out_molec}_complex_md.gro -o ${out_molec}_complex_md_no_center.pdb -dump $dump"
#execute "echo 0 |${prefix_gromacs} trjconv${mpi} -f ${out_molec}_complex_md_center.xtc -s ${out_molec}_complex_md.gro -o ${out_molec}_complex_md.pdb -dump $dump"