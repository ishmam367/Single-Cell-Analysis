#!/bin/bash

# Find the project directory
PROJECT_DIR=$(dirname "$(realpath $0)")
dataset_id="diabetesII_h1" #possible Id's = diabetesII_h1, pbmc3k, mpn_HI108


echo "---------------------------- Pipelines started for $dataset_id ----------------------------"

python3 "${PROJECT_DIR}/pipeline.py" \
    --dataset_id "${dataset_id}"

echo "---------------------------- Pipelines ended for $dataset_id ----------------------------"