#!/usr/bin/env bash
# Stockholm2Fasta.py has to be in the same directory as this script

#If no arguments are passed to Grab_models.sh, print help output
if [[ "$#" == 0 ]]; then
	printf "\t-i [pfam_accession_list.txt]\tSpecify path to file containing a new-line delimited list of pfam accessions (e.g., PF01582.23)"
	printf "\n\t-p [Pfam-A.hmm]\t\t\tSpecify path to uncompressed Pfam database"
	printf "\n\t-s [Pfam-A.seed]\t\tSpecify path to uncompressed and indexed Pfam seed dataset\n"
	exit 1
fi

#Program command variables & execution check
S2F=`echo "$0" | sed 's/Grab_models.sh$/Stockholm2Fasta.py/'`
HMMFETCH=`which hmmfetch 2>/dev/null`
ESLAFETCH=`which esl-afetch 2>/dev/null`

if ! [[ -x $HMMFETCH ]];then
	EXECUTABLE_ERROR=`echo "TRUE"`
	printf "ERROR: Cannot find hmmfetch executable in \$PATH\n"
fi

if ! [[ -x $ESLAFETCH ]];then
	EXECUTABLE_ERROR=`echo "TRUE"`
	printf "ERROR: Cannot find esl-afetch executable in \$PATH\n"
fi

if [[ ${EXECUTABLE_ERROR} == "TRUE" ]]; then #If any program cannot be executed, exit
        printf "Exiting...\n"
        exit 1
fi

#User-input variable definitions
while getopts ":i:p:s:" opt; do
	case $opt in
		i) # Set path to user-specified Pfam accession list
			ACCESSION_LIST=`echo "${OPTARG}"`
			;;
		p) # Set path to Pfam database
			PFAM_DB_PATH=`echo "${OPTARG}"`
			;;
		s) # Set path to Pfam seed database
			PFAM_SEEDS=`echo "${OPTARG}"`
			;;
		\?) # Error call if invalid flag provided
			echo "Invalid option -${OPTARG}"
			exit 1
			;;
		:) # Error call if flag reguires an argument and none provided
			echo "Flag -$OPTARG requires an argument"
			exit 1
			;;
	esac
done

#Input check
ERROR=`echo "FALSE"`

if [[ ! -s ${ACCESSION_LIST} ]];then
	echo "List of Pfam accessions cannot be found."
	ERROR=`echo "TRUE"`
fi
if [[ ! -s ${PFAM_DB_PATH} ]];then
	echo "Pfam database cannot be found."
	ERROR=`echo "TRUE"`
fi
if [[ ! -s ${PFAM_SEEDS} ]];then
	echo "Pfam seed database cannot be found."
	ERROR=`echo "TRUE"`
fi
if [[ ! -s ${PFAM_SEEDS}.ssi ]];then
	echo "Pfam seed database is not indexed. Index using esl-afetch --index [Pfam-A.seeds]."
	ERROR=`echo "TRUE"`
fi
if [[ ${ERROR} == "TRUE" ]]; then #If any errors above, exit.
	printf "\nPlease check the ERROR reports above and try again.\nExiting...\n"
	exit 1
fi

#Grab models and seeds
while read PFAM_ACCESSION
do
	if [[ `grep -c ${PFAM_ACCESSION} $PFAM_SEEDS | xargs` == 0 ]];then
		echo "ERROR: ${PFAM_ACCESSION} not found, skipping..."
		continue
	fi

	if [[ ` echo "${PFAM_ACCESSION}" | grep -c "\." | xargs` != 1 ]];then
		echo "ERROR: ${PFAM_ACCESSION} does not possess version number pull from seed databse; skipping..."
		continue
	fi

	echo "...Obtaining ${PFAM_ACCESSION} profile HMM"
	$HMMFETCH ${PFAM_DB_PATH} ${PFAM_ACCESSION} > hmm
	NAME=`head hmm | grep 'NAME' | awk '{print $2}'`
	SHORT_ACCESSION=`head hmm | grep 'ACC' | awk '{print $2}' | awk -F"." '{print $1}'`
	mv hmm ${SHORT_ACCESSION}_${NAME}.hmm

	echo "...Obtaining ${PFAM_ACCESSION} seed alignment"
	$ESLAFETCH $PFAM_SEEDS ${PFAM_ACCESSION} > seed
	$S2F seed ${SHORT_ACCESSION}_${NAME}.fasta
	rm seed
done < ${ACCESSION_LIST}
