#!/bin/bash
# $1 input is list of pfam accessions (e.g., PF00065)
# Stockholm2Fasta.py has to be in the same directory as this script

S2F=`echo "$0" | sed 's/Grab_models.sh$/Stockholm2Fasta.py/'`

while read PFAM_ACCESSION
do
	echo "...Downloading ${PFAM_ACCESSION} profile HMM"
	wget "https://pfam.xfam.org/family/${PFAM_ACCESSION}/hmm"		
	NAME=`head hmm | grep 'NAME' | awk '{print $2}'`
	mv hmm ${PFAM_ACCESSION}_${NAME}.hmm

	echo "...Downloading ${PFAM_ACCESSION} seed alignment"
	wget "https://pfam.xfam.org/family/${PFAM_ACCESSION}/alignment/seed"
	${S2F} seed ${PFAM_ACCESSION}_${NAME}.fasta
	rm seed
done < $1
