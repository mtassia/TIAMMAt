#!/bin/bash

for FASTA in *.renamed_headers.fasta
do
	NEW_NAME=`echo "${FASTA}" | sed 's/.renamed_headers//'`
	printf "Old name was ${FASTA}, new name is ${NEW_NAME}\n"
	mv ${FASTA} ${NEW_NAME}
done
