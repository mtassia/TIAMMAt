#!/bin/bash
# $1 input is list of pfam accessions (e.g., PF00065)
# $2 input is Pfam-A.hmm
# $3 input is Pfam-A.seed
# Stockholm2Fasta.py has to be in the same directory as this script

S2F=`echo "$0" | sed 's/Grab_models.sh$/Stockholm2Fasta.py/'`

if [[ ! -s ${1} ]];then
	echo "Input does cannot be found. Exiting..."
	exit 0
fi

while read PFAM_ACCESSION
do
	echo "...Obtaining ${PFAM_ACCESSION} profile HMM"
	#wget "https://pfam.xfam.org/family/${PFAM_ACCESSION}/hmm"		
	hmmfetch $2 ${PFAM_ACCESSION} > hmm
	NAME=`head hmm | grep 'NAME' | awk '{print $2}'`
	SHORT_ACCESSION=`head hmm | grep 'ACC' | awk '{print $2}' | awk -F"." '{print $1}'`
	mv hmm ${SHORT_ACCESSION}_${NAME}.hmm

	echo "...Obtaining ${PFAM_ACCESSION} seed alignment"
	#wget "https://pfam.xfam.org/family/${PFAM_ACCESSION}/alignment/seed"
	esl-afetch $3 ${PFAM_ACCESSION} > seed
	${S2F} seed ${SHORT_ACCESSION}_${NAME}.fasta
	rm seed
done < $1
