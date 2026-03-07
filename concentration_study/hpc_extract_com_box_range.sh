#!/usr/bin/env bash
set -u

# Non-interactive batch extraction of per-molecule COM + box dimensions.
# Runs across folders like 8_mmc ... 150_mmc on HPC.
#
# Usage:
#   bash hpc_extract_com_box_range.sh
#   bash hpc_extract_com_box_range.sh 8 150
#   bash hpc_extract_com_box_range.sh 8 150 /path/to/concentration_study

MIN_MMC="${1:-8}"
MAX_MMC="${2:-150}"
ROOT_DIR="${3:-$(pwd)}"
STEP_MMC="${4:-4}"

LOG_FILE="${ROOT_DIR}/batch_extract_$(date +%Y%m%d_%H%M%S).log"

echo "Root directory: ${ROOT_DIR}" | tee -a "${LOG_FILE}"
echo "Range: ${MIN_MMC}_mmc to ${MAX_MMC}_mmc (step ${STEP_MMC})" | tee -a "${LOG_FILE}"

if ! command -v gmx >/dev/null 2>&1; then
  echo "ERROR: gmx not found in PATH" | tee -a "${LOG_FILE}"
  exit 1
fi

process_dir() {
  local dir="$1"
  local n="$2"

  echo "" | tee -a "${LOG_FILE}"
  echo "===== Processing ${dir} (${n} molecules) =====" | tee -a "${LOG_FILE}"

  if [ ! -d "${dir}" ]; then
    echo "SKIP: directory not found" | tee -a "${LOG_FILE}"
    return 0
  fi

  cd "${dir}" || { echo "ERROR: cannot cd ${dir}" | tee -a "${LOG_FILE}"; return 1; }

  local tpr="mmc_${n}_325K.tpr"
  local xtc="mmc_${n}_325K.xtc"

  if [ ! -f "${tpr}" ]; then
    tpr=$(ls -1 *.tpr 2>/dev/null | head -n 1 || true)
  fi
  if [ ! -f "${xtc}" ]; then
    xtc=$(ls -1 *.xtc 2>/dev/null | head -n 1 || true)
  fi

  if [ -z "${tpr:-}" ] || [ ! -f "${tpr}" ]; then
    echo "FAIL: no .tpr found" | tee -a "${LOG_FILE}"
    cd "${ROOT_DIR}" || true
    return 1
  fi
  if [ -z "${xtc:-}" ] || [ ! -f "${xtc}" ]; then
    echo "FAIL: no .xtc found" | tee -a "${LOG_FILE}"
    cd "${ROOT_DIR}" || true
    return 1
  fi

  echo "Using TPR: ${tpr}" | tee -a "${LOG_FILE}"
  echo "Using XTC: ${xtc}" | tee -a "${LOG_FILE}"

  # Build a default index first, then auto-detect the X2IT group id.
  if ! printf "q\n" | gmx make_ndx -f "${tpr}" -o index_default.ndx >> "${LOG_FILE}" 2>&1; then
    echo "FAIL: gmx make_ndx (default index)" | tee -a "${LOG_FILE}"
    cd "${ROOT_DIR}" || true
    return 1
  fi

  local x2it_group
  x2it_group=$(awk '
    /^\[/ {
      name=$0
      gsub(/\[/, "", name)
      gsub(/\]/, "", name)
      gsub(/^ +| +$/, "", name)
      groups[++count]=name
    }
    END {
      for (i=1; i<=count; i++) {
        if (groups[i] == "X2IT") {
          print i-1
          exit
        }
      }
    }
  ' index_default.ndx)

  if [ -z "${x2it_group}" ]; then
    echo "FAIL: could not find X2IT group in default index" | tee -a "${LOG_FILE}"
    cd "${ROOT_DIR}" || true
    return 1
  fi

  echo "Detected X2IT group id: ${x2it_group}" | tee -a "${LOG_FILE}"

  cat > make_ndx_auto.in <<EOF
keep ${x2it_group}
splitres 0
q
EOF

  if ! gmx make_ndx -f "${tpr}" -o mmc_molecules.ndx < make_ndx_auto.in >> "${LOG_FILE}" 2>&1; then
    echo "FAIL: gmx make_ndx" | tee -a "${LOG_FILE}"
    cd "${ROOT_DIR}" || true
    return 1
  fi

  # Extract COM for groups 1-N (group 0 is combined, groups 1-N are individual split molecules).
  if ! seq 1 "${n}" | gmx traj -f "${xtc}" -s "${tpr}" -n mmc_molecules.ndx -ox com_all.xvg -com -ng "${n}" >> "${LOG_FILE}" 2>&1; then
    echo "FAIL: gmx traj -com" | tee -a "${LOG_FILE}"
    cd "${ROOT_DIR}" || true
    return 1
  fi

  # Extract time-varying box dimensions (select System automatically).
  if ! printf "0\n" | gmx traj -f "${xtc}" -s "${tpr}" -ob box.xvg >> "${LOG_FILE}" 2>&1; then
    echo "FAIL: gmx traj -ob" | tee -a "${LOG_FILE}"
    cd "${ROOT_DIR}" || true
    return 1
  fi

  # Quick sanity check for expected COM columns: 1 + 3*N
  local first_data
  first_data=$(awk '!/^[@#]/{print; exit}' com_all.xvg)
  local cols
  cols=$(echo "${first_data}" | awk '{print NF}')
  local expected=$((1 + 3 * n))
  if [ "${cols}" -lt "${expected}" ]; then
    echo "WARN: columns=${cols}, expected>=${expected}" | tee -a "${LOG_FILE}"
  else
    echo "OK: COM columns=${cols}, expected=${expected}" | tee -a "${LOG_FILE}"
  fi

  echo "DONE: ${dir}" | tee -a "${LOG_FILE}"
  cd "${ROOT_DIR}" || true
  return 0
}

cd "${ROOT_DIR}" || exit 1

FAIL_COUNT=0
SUCCESS_COUNT=0
DISCOVERED_COUNT=0
SELECTED_COUNT=0

mapfile -t CANDIDATE_DIRS < <(find "${ROOT_DIR}" -maxdepth 2 -type d -name '*_mmc' | sort -V)
DISCOVERED_COUNT=${#CANDIDATE_DIRS[@]}
echo "Discovered directories: ${DISCOVERED_COUNT}" | tee -a "${LOG_FILE}"
if [ "${DISCOVERED_COUNT}" -eq 0 ]; then
  echo "ERROR: no *_mmc directories found under ${ROOT_DIR}" | tee -a "${LOG_FILE}"
  exit 1
fi

for dir in "${CANDIDATE_DIRS[@]}"; do
  dir_base=$(basename "${dir}")
  if [[ ! "${dir_base}" =~ ^([0-9]+)_mmc$ ]]; then
    continue
  fi
  n="${BASH_REMATCH[1]}"

  if [ "${n}" -lt "${MIN_MMC}" ] || [ "${n}" -gt "${MAX_MMC}" ]; then
    continue
  fi

  if [ $(( (n - MIN_MMC) % STEP_MMC )) -ne 0 ]; then
    continue
  fi

  SELECTED_COUNT=$((SELECTED_COUNT + 1))
  if process_dir "${dir}" "${n}"; then
    SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
  else
    FAIL_COUNT=$((FAIL_COUNT + 1))
  fi
done

if [ "${SELECTED_COUNT}" -eq 0 ]; then
  echo "WARNING: found *_mmc directories, but none matched range/step (${MIN_MMC}-${MAX_MMC}, step ${STEP_MMC})" | tee -a "${LOG_FILE}"
fi

echo "" | tee -a "${LOG_FILE}"
echo "Discovered=${DISCOVERED_COUNT}, selected=${SELECTED_COUNT}" | tee -a "${LOG_FILE}"
echo "Finished. success=${SUCCESS_COUNT}, fail=${FAIL_COUNT}" | tee -a "${LOG_FILE}"
echo "Log: ${LOG_FILE}" | tee -a "${LOG_FILE}"

exit 0
