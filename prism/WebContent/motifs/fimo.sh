#!/bin/bash

MOTIF=$1
QUERY=$2

echo $MODEL
echo $SEQUENCES

# run fimo
fimo --text -verbosity 1 $MOTIF $QUERY