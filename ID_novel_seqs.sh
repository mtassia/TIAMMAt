#!/bin/bash

#Execute in SECOND_SEARCH directory from IRDS

for NOVEL_LIST in *.sequences_only_identified_after_revision.txt
do
	TAXON_ID=`echo "${NOVEL_LIST}" | awk -F"." '{print $1}'`
	DOMAIN_ID=`echo "${NOVEL_LIST}" | awk -F"." '{print $2}'`
	FASTAFILE_ID=`echo "${TAXON_ID}.${DOMAIN_ID}_base_and_revised_present.fasta"`
	
	if [[ -s ${NOVEL_LIST} ]]; then #Run two different loops depending on whether novel seqs were identfied
		sed "/>/s/>.*/&_${TAXON_ID}_${DOMAIN_ID}/" ${FASTAFILE_ID} > ${FASTAFILE_ID}.renaming_tmp #Sed to add domain IDs to the end of each fasta header and make a new file for renamed seqs

		while read NOVEL_ACCESSION #Read the newly identified seq list
		do
			sed -i "/>/s/${NOVEL_ACCESSION}.*/&_NOV/" ${FASTAFILE_ID}.renaming_tmp #Append "NOV" to headers of sequences which weren't identfied before domain revision
		done < ${NOVEL_LIST}
		mv ${FASTAFILE_ID}.renaming_tmp ${TAXON_ID}.${DOMAIN_ID}_base_and_revised_present.IDs_renamed.fasta #Rename the temporary renaming file
	else
		sed "/>/s/>.*/&_${TAXON_ID}_${DOMAIN_ID}/" ${FASTAFILE_ID} > ${TAXON_ID}.${DOMAIN_ID}_base_and_revised_present.IDs_renamed.fasta #For fasta files without newly identfied sequences, append the domain ID to header.
	fi
done
