#!/bin/bash

DATABASE=$1
QUERY=$2

if [ "$#" -ne 3 ]; then
	# run blastp
	blastp -db $DATABASE -query $QUERY
else 
	EVALUE=$3
	# run blastp
	blastp -db $DATABASE -query $QUERY -evalue $EVALUE
fi


