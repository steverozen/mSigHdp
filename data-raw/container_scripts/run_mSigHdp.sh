#!/bin/bash

PROJECT_DIR="."
SEED="519"
K_GUESS="26"
OUTPUT_DIR="."
singularity exec -e $PROJECT_DIR/msighdp_2.1.0.3.sif nohup Rscript --vanilla $PROJECT_DIR/test_mSigHdp.R $PROJECT_DIR $SEED $K_GUESS $OUTPUT_DIR > $PROJECT_DIR/msighdp.log 2>&1 &
