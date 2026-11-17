#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# Wait for the user-created READY sentinel, then run the real pipeline once.
while true; do
  if [ -f 03_Raw_Data/READY ]; then
    bash run_all.sh config/pipeline_config.yaml >> Logs/watch_and_run.log 2>&1
    exit $?
  fi
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] waiting for 03_Raw_Data/READY" >> Logs/watch_and_run.log
  sleep 600
done
