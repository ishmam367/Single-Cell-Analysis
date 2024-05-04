#!/bin/bash

# Find the project directory
PROJECT_DIR=$(dirname "$(realpath $0)")
datasets=$(find "${PROJECT_DIR}/dataset" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;)

for dataset_id in $datasets; do
    echo "---------------------------- Pipelines started for $dataset_id ----------------------------"

    python3 "${PROJECT_DIR}/pipeline.py" \
        --dataset_id "${dataset_id}"

    echo "---------------------------- Pipelines ended for $dataset_id ----------------------------"

    echo ""
done
