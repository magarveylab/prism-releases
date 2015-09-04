#!/bin/bash

MODEL=$1
SEQUENCES=$2

echo $MODEL
echo $SEQUENCES

# run hmmsearch
hmmsearch $MODEL $SEQUENCES
