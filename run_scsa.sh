#!/bin/bash

# Parse input arguments
while getopts ":i:o:" opt; do
  case ${opt} in
    i ) # Input file argument
      input=$OPTARG
      ;;
    o ) # Output file argument
      output=$OPTARG
      ;;
    \? ) # Invalid option
      echo "Usage: $0 -i <input_file> -o <output_file>"
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

# Check if the input and output files are provided
if [ -z "$input" ] || [ -z "$output" ]; then
  echo "Input or output file not provided. Exiting."
  exit 1
fi

PROJECT_DIR=$(dirname "$(realpath $0)")
SCSA_DIR="$PROJECT_DIR/SCSA"

echo "Project Directory: $PROJECT_DIR"
echo "SCSA Directory: $SCSA_DIR"

# Run SCSA script with the provided input and output file names
python3 "$SCSA_DIR/SCSA.py" \
    --input "$input" \
    --output "$output" \
    --db "$SCSA_DIR/whole_v2.db" \
    --source scanpy \
    --foldchange 1.5 \
    --pvalue 0.01 \
    --species Human \
    --outfmt txt \
    --Gensymbol \
    --noprint





















# #!/bin/bash

# # Parse command line arguments
# while [[ "$#" -gt 0 ]]; do
#     case $1 in
#         -o|--output) output="$2"; shift ;;
#         *) echo "Unknown parameter passed: $1"; exit 1 ;;
#     esac
#     shift
# done

# # Check if the output file is provided
# if [ -z "$output" ]; then
#     echo "Output file not provided. Exiting."
#     exit 1
# fi

# PROJECT_DIR=$(dirname "$(realpath $0)")
# SCSA_DIR="$PROJECT_DIR/SCSA"

# echo "Project Directory: $PROJECT_DIR"
# echo "SCSA Directory: $SCSA_DIR"

# # Run SCSA script with the provided output file name
# python3 "$SCSA_DIR/SCSA.py" \
#     --input "$PROJECT_DIR/diff_exp_result.csv" \
#     --output "$output" \
#     --db "$SCSA_DIR/whole_v2.db" \
#     --source scanpy \
#     --foldchange 1.5 \
#     --pvalue 0.01 \
#     --species Human \
#     --outfmt txt \
#     --Gensymbol \
#     --noprint




# 'PROJECT_DIR=$(dirname "$(realpath $0)")
# SCSA_DIR="$PROJECT_DIR/SCSA"

# echo "Project Directory: $PROJECT_DIR"
# echo "SCSA Directory: $SCSA_DIR"

# python3 "$SCSA_DIR/SCSA.py" \
#     --input "$PROJECT_DIR/diff_exp_result.csv" \
#     --output "$PROJECT_DIR/scsa_result.txt" \
#     --db "$SCSA_DIR/whole_v2.db" \
#     --source scanpy \
#     --foldchange 1.5 \
#     --pvalue 0.01 \
#     --species Human \
#     --outfmt txt \
#     --Gensymbol \
#     --noprint'