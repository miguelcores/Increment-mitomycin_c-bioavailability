#!/usr/bin/env bash
#SBATCH --job-name=mmc_com_extract
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.log
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --ntasks=1

set -x  # Print all commands executed

MIN_MMC="${1:-8}"
MAX_MMC="${2:-150}"
ROOT_DIR="${3:-.}"
STEP_MMC="${4:-4}"

LOG_FILE="${ROOT_DIR}/batch_extract_${SLURM_JOB_ID}.log"
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

{
  echo "========================================"
  echo "SLURM COM+Box Extraction Started"
  echo "========================================"
  echo "Timestamp: ${TIMESTAMP}"
  echo "Job ID: ${SLURM_JOB_ID}"
  echo "Host: $(hostname)"
  echo "PWD: $(pwd)"
  echo "Root directory: ${ROOT_DIR}"
  echo "Range: ${MIN_MMC}_mmc to ${MAX_MMC}_mmc (step ${STEP_MMC})"
  echo ""
  echo "Environment:"
  echo "  PATH=${PATH}"
  echo "  User: $(whoami)"
  echo ""

  # Load GROMACS module if needed (adjust for your HPC)
  if command -v module >/dev/null 2>&1; then
    echo "Loading GROMACS module..."
    module load GROMACS/2021.5-foss-2021b 2>&1 || {
      echo "WARNING: Could not load module, trying gmx from PATH"
    }
  fi

  echo "GROMACS version:"
  gmx -version 2>&1 | head -3

  if ! command -v gmx >/dev/null 2>&1; then
    echo "ERROR: gmx not found in PATH" >&2
    exit 1
  fi

  echo ""
  echo "========================================"
  echo "Starting extraction loop"
  echo "========================================"
  echo ""

  cd "${ROOT_DIR}" || {
    echo "ERROR: cannot cd to ${ROOT_DIR}" >&2
    exit 1
  }

  SUCCESS_COUNT=0
  FAIL_COUNT=0
  DISCOVERED_COUNT=0
  SELECTED_COUNT=0

  echo "Discovering candidate directories (*_mmc) under ${ROOT_DIR} (maxdepth=2)..."
  mapfile -t CANDIDATE_DIRS < <(find "${ROOT_DIR}" -maxdepth 2 -type d -name '*_mmc' | sort -V)
  DISCOVERED_COUNT=${#CANDIDATE_DIRS[@]}
  echo "Discovered directories: ${DISCOVERED_COUNT}"
  if [ "${DISCOVERED_COUNT}" -eq 0 ]; then
    echo "ERROR: no *_mmc directories found under ${ROOT_DIR}" >&2
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

    echo ""
    echo "===== Processing ${dir_base} (${n} molecules) ====="
    echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"

    cd "${dir}" || {
      echo "ERROR: cannot cd ${dir}" >&2
      FAIL_COUNT=$((FAIL_COUNT + 1))
      continue
    }

    pwd
    ls -lh *.tpr *.xtc 2>&1 | head -5

    # Find TPR/XTC files (prefer production files)
    tpr="mmc_${n}_325K.tpr"
    xtc="mmc_${n}_325K.xtc"
    if [ ! -f "${tpr}" ]; then
      tpr=$(ls -1 *.tpr 2>/dev/null | head -n 1 || true)
    fi
    if [ ! -f "${xtc}" ]; then
      xtc=$(ls -1 *.xtc 2>/dev/null | head -n 1 || true)
    fi

    if [ -z "${tpr}" ] || [ ! -f "${tpr}" ]; then
      echo "ERROR: no .tpr file" >&2
      FAIL_COUNT=$((FAIL_COUNT + 1))
      cd "${ROOT_DIR}" || true
      continue
    fi

    if [ -z "${xtc}" ] || [ ! -f "${xtc}" ]; then
      echo "ERROR: no .xtc file" >&2
      FAIL_COUNT=$((FAIL_COUNT + 1))
      cd "${ROOT_DIR}" || true
      continue
    fi

    echo "TPR: ${tpr}"
    echo "XTC: ${xtc}"

    # Step 1: Create default index to detect X2IT group
    echo ""
    echo "Step 1: Creating default index..."
    if ! printf "q\n" | gmx make_ndx -f "${tpr}" -o index_default.ndx 2>&1 | tail -20; then
      echo "ERROR: gmx make_ndx failed" >&2
      FAIL_COUNT=$((FAIL_COUNT + 1))
      cd "${ROOT_DIR}" || true
      continue
    fi

    # Step 2: Find X2IT group id
    echo ""
    echo "Step 2: Detecting X2IT group..."
    x2it_group=""
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
      echo "ERROR: could not find X2IT group" >&2
      echo "Available groups in index_default.ndx:"
      grep "^\[" index_default.ndx 2>&1
      FAIL_COUNT=$((FAIL_COUNT + 1))
      cd "${ROOT_DIR}" || true
      continue
    fi

    echo "Found X2IT at group id: ${x2it_group}"

    # Step 3: Create selection script for make_ndx
    echo ""
    echo "Step 3: Creating split index..."
    cat > make_ndx_auto.in <<EOF
keep ${x2it_group}
splitres 0
q
EOF
    cat make_ndx_auto.in

    if ! gmx make_ndx -f "${tpr}" -o mmc_molecules.ndx < make_ndx_auto.in 2>&1 | tail -20; then
      echo "ERROR: gmx make_ndx (split) failed" >&2
      FAIL_COUNT=$((FAIL_COUNT + 1))
      cd "${ROOT_DIR}" || true
      continue
    fi

    echo "Index file created. Groups:"
    grep "^\[" mmc_molecules.ndx 2>&1 | head -10

    # Step 4: Extract COM for all molecules (groups 1-N are split molecules)
    echo ""
    echo "Step 4: Extracting COM coordinates..."
    echo "Selecting groups: $(seq 1 ${n} | tr '\n' ' ')"
    if ! seq 1 "${n}" | gmx traj -f "${xtc}" -s "${tpr}" -n mmc_molecules.ndx \
         -ox com_all.xvg -com -ng "${n}" 2>&1 | tail -30; then
      echo "ERROR: gmx traj -com failed" >&2
      FAIL_COUNT=$((FAIL_COUNT + 1))
      cd "${ROOT_DIR}" || true
      continue
    fi

    # Verify output
    if [ ! -f com_all.xvg ]; then
      echo "ERROR: com_all.xvg not created" >&2
      FAIL_COUNT=$((FAIL_COUNT + 1))
      cd "${ROOT_DIR}" || true
      continue
    fi

    first_data=""
    cols=""
    expected=""
    first_data=$(awk '!/^[@#]/{print; exit}' com_all.xvg)
    cols=$(echo "${first_data}" | awk '{print NF}')
    expected=$((1 + 3 * n))

    echo "COM file sanity check:"
    echo "  Columns: ${cols} (expected: ${expected})"
    echo "  Frames: $(grep -c -v '^[@#]' com_all.xvg || echo 0)"
    echo "  Sample data: ${first_data:0:100}..."

    # Step 5: Extract box dimensions
    echo ""
    echo "Step 5: Extracting box dimensions..."
    if ! printf "0\n" | gmx traj -f "${xtc}" -s "${tpr}" -ob box.xvg 2>&1 | tail -20; then
      echo "ERROR: gmx traj -ob failed" >&2
      FAIL_COUNT=$((FAIL_COUNT + 1))
      cd "${ROOT_DIR}" || true
      continue
    fi

    # Verify output
    if [ ! -f box.xvg ]; then
      echo "ERROR: box.xvg not created" >&2
      FAIL_COUNT=$((FAIL_COUNT + 1))
      cd "${ROOT_DIR}" || true
      continue
    fi

    echo "Box file created:"
    echo "  Frames: $(grep -c -v '^[@#]' box.xvg || echo 0)"
    echo "  Sample: $(awk '!/^[@#]/{print; exit}' box.xvg)"

    echo ""
    echo "SUCCESS: ${dir_base}"
    SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    cd "${ROOT_DIR}" || true
  done

  if [ "${SELECTED_COUNT}" -eq 0 ]; then
    echo "WARNING: found *_mmc directories, but none matched range/step (${MIN_MMC}-${MAX_MMC}, step ${STEP_MMC})"
  fi

  echo ""
  echo "========================================"
  echo "Extraction complete"
  echo "========================================"
  echo "Discovered: ${DISCOVERED_COUNT}"
  echo "Selected: ${SELECTED_COUNT}"
  echo "Successful: ${SUCCESS_COUNT}"
  echo "Failed: ${FAIL_COUNT}"
  echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
  echo ""

} 2>&1 | tee "${LOG_FILE}"

echo "Full log written to: ${LOG_FILE}"
exit 0
