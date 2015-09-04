#!/bin/bash

for file in *.fasta; do
	basename=$(basename $file)
	filename=${basename%.*}	
	makeblastdb -in "$file" -dbtype prot -out "$filename"
done