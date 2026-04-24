#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOWNLOAD_DIR="${ROOT_DIR}/download"
TMP_DIR="${DOWNLOAD_DIR}/tmp"
ASSET_ID="1yjLWm5uBZqOT_6pRgnd_wXpSFjZicmvt"
ASSET_NAME="STELLAR.tar"
ASSET_PATH="${DOWNLOAD_DIR}/${ASSET_NAME}"

mkdir -p "${DOWNLOAD_DIR}" "${TMP_DIR}"

download_with_gdown() {
  if command -v gdown >/dev/null 2>&1; then
    echo "Downloading ${ASSET_NAME} from Google Drive with gdown..."
    gdown --id "${ASSET_ID}" --output "${ASSET_PATH}"
    return 0
  fi
  return 1
}

download_with_python_fallback() {
  echo "gdown not found. Falling back to Python downloader..."
  python3 - <<'PY'
import os
import re
import sys
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import HTTPCookieProcessor, Request, build_opener

asset_id = "1yjLWm5uBZqOT_6pRgnd_wXpSFjZicmvt"
out_path = Path(os.environ["ASSET_PATH"])

base = "https://drive.google.com/uc"
params = urlencode({"export": "download", "id": asset_id})
url = f"{base}?{params}"

opener = build_opener(HTTPCookieProcessor())

def fetch(u: str):
    req = Request(u, headers={"User-Agent": "Mozilla/5.0"})
    return opener.open(req)

resp = fetch(url)
data = resp.read()

text_probe = data[:200000].decode("utf-8", errors="ignore")
match = re.search(r'confirm=([0-9A-Za-z_]+)', text_probe)

if match:
    confirm = match.group(1)
    url = f"{base}?{urlencode({'export': 'download', 'confirm': confirm, 'id': asset_id})}"
    resp = fetch(url)
    data = resp.read()

if not data:
    print("Download failed: empty response from Google Drive.", file=sys.stderr)
    sys.exit(1)

out_path.write_bytes(data)
print(f"Saved {out_path}")
PY
}

safe_extract() {
  local tar_file="$1"
  local destination="$2"

  while IFS= read -r entry; do
    if [[ "${entry}" == /* ]] || [[ "${entry}" == *".."* ]]; then
      echo "Unsafe path detected inside ${tar_file}: ${entry}" >&2
      exit 1
    fi
  done < <(tar -tf "${tar_file}")

  tar -xf "${tar_file}" -C "${destination}"
}

if [[ ! -f "${ASSET_PATH}" ]]; then
  if ! download_with_gdown; then
    export ASSET_PATH
    download_with_python_fallback
  fi
else
  echo "${ASSET_NAME} already exists at ${ASSET_PATH}; skipping download."
fi

rm -rf "${TMP_DIR}"
mkdir -p "${TMP_DIR}"

echo "Extracting ${ASSET_NAME}..."
safe_extract "${ASSET_PATH}" "${TMP_DIR}"

echo "Expanding nested .tar files found in ${ASSET_NAME}..."
while IFS= read -r nested_tar; do
  safe_extract "${nested_tar}" "${ROOT_DIR}"
done < <(find "${TMP_DIR}" -type f -name "*.tar" | sort)

echo "Copying non-tar asset content to repository root..."
shopt -s dotglob nullglob
for item in "${TMP_DIR}"/*; do
  base_name="$(basename "${item}")"
  if [[ "${base_name}" == *.tar ]]; then
    continue
  fi
  cp -a "${item}" "${ROOT_DIR}/"
done
shopt -u dotglob nullglob

echo "Assets deployed successfully."
