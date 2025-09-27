#!/bin/bash


OUTS_DIR="cellranger_test2/outs"


SUMMARY_DIR="cellranger_summary"


mkdir -p "$SUMMARY_DIR"


cp "$OUTS_DIR/web_summary.html" "$SUMMARY_DIR/"
cp "$OUTS_DIR/metrics_summary.csv" "$SUMMARY_DIR/"


# cp "$OUTS_DIR/cloupe.cloupe" "$SUMMARY_DIR/"
# cp -r "$OUTS_DIR/analysis" "$SUMMARY_DIR/"

echo "Summary files are collected in: $SUMMARY_DIR"
