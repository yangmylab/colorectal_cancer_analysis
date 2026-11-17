#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

CONFIG="${1:-config/pipeline_config.yaml}"
shift || true

mkdir -p Logs
LOG="Logs/run_all_$(date +%Y%m%d_%H%M%S).log"
Rscript 06_Analysis/1_R/1_Rscripts/run_pipeline.R --config "$CONFIG" "$@" 2>&1 | tee "$LOG"
