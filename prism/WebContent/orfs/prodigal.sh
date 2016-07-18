#!/bin/bash

MODE=$1
OUTPUT=$2
INPUT=$3
PROTEIN=$4
NUCLEOTIDE=$5

# run prodigal
prodigal -p $MODE -o $OUTPUT -i $INPUT -a $PROTEIN -d $NUCLEOTIDE
