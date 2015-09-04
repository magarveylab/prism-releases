#!/bin/bash

INPUT=$1
PROTEIN=$2
NUCLEOTIDE=$3

# run prodigal
prodigal -i $INPUT -a $PROTEIN -d $NUCLEOTIDE
