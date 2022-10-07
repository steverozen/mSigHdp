#!/bin/bash

PROJECT_DIR="/home/e0012078/container/msighdp"
SEED="519"
K_GUESS="26"
OUTPUT_DIR="/home/e0012078/output/msighdp"

singularity exec -e $PROJECT_DIR/msighdp.sif nohup Rscript --vanilla $PROJECT_DIR/test_mSigHdp.R $PROJECT_DIR $SEED $K_GUESS $OUTPUT_DIR > $PROJECT_DIR/msighdp.log 2>&1 &