# RelaxedDomainSearch
---

Code associated with domain HMM model revision for use in non-model species protein evolution studies.

Description: [More to come]

---
### PIPELINE:
![RelaxedDomainSearch Pipeline](https://github.com/mtassia/RelaxedDomainSearch/blob/master/Program_diagram.png)

---
### DEPENDENCIES:
- `HMMer` (available at http://hmmer.org/)
- `Pull_coordinates.py` (Included in this repo)
- `Select_contigs.pl` (Included in this repo)
- `Best_fit_domains.py` (Included in this repo)
Each of these programs must be installed into `$PATH`

---
### USAGE:
It's recommended that for all path arguments, absolute paths be used. 

```
-d [directory]  Specify directory containing datasets [DEFAULT: pwd]
-h              Print help
-m [directory]  Specify directory containing pfam models and seeds [MANDATORY; cannot be the same directory as datasets]
-o [string]     Set name for output directory [default: RDS_output_{DATE}]
-p [path]       Specify path to Pfam database [MANDATORY]
-t [integer]    Specify number of threads for hmmscan and hmmsearch
```

**Example Command:** 
```
ImprovedRelaxedDomainSearch -d /path/to/Proteomes/ -m /path/to/Target_Pfams/ -p /path/to/Pfam-A.hmm -t 4
```

---
### OUTPUTS:
**For the following inputs:**
- `-d /path/to/Proteomes/` contains:
  >- ProteomeA.fasta
  >- ProteomeB.fasta
- `-m /path/to/Target_Pfams/` contains:
  >- PF00069_Pkinase.fasta
  >- PF00069_Pkinase.hmm
  >- PF00240_Ubiquitin.fasta
  >- PF00240_Ubiquitin.hmm
- `-o [default]`

**Outputs the following:**
```bash
RDS_output_[WkDay]_[Month]_[Year]/ #Default output directory created by program
  1. IDENTIFICATION_STATISTICS.txt #TSV containing domain counts before & after revision
  2. SCRIPT_LOG.txt #Data log reporting operations and immediate results during program execution
  
RDS_output_[WkDay]_[Month]_[Year]/MODELS/ #Directory containing compressed models for initial HMMsearch step
  3. PF00069_Pkinase.hmm.h3* #Four files associated with initial Pkinase domain compression from HMMpress
  4. [Repeat (3) for Ubiquitin domain]
  
RDS_output_[WkDay]_[Month]_[Year]/Pkinase_PF00069/ #Output directory for Pkinase domain revision
  5. ProteomeA.hmmsearch_for_Pkinase.domtblout #HMMsearch domain table output for potential Pkinase domains in ProteomeA.fasta
  6. ProteomeA.Pkinase_present.fasta #Fasta file containing sequences with potential Pkinase domains found in ProteomeA.fasta
  7. ProteomeA.Pkinase_present.hmmscan_against_pfam.domtblout #HMMscan output (domain table) of above sequences annotated with all Pfam domains
  8. ProteomeA.Pkinase_present.hmmscan_against_pfam.domtblout.besthits.tsv #Above table filtered for domains meeting HMMer inclusion threshold and overlap
  9. ProteomeA.Pkinase_present.hmmscan_against_pfam.Pkinase_coordinates.tsv #Coordinate TSV for all best-fit Pkinase domains in ProteomeA.fasta
  10. ProteomeA.Pkinase_present.hmmscan_against_pfam.target_list.txt #List of sequences with best-fit Pkinase domains in ProteomA.fasta
  11. [Repeat (5-10) for ProteomeB]
  
RDS_output_[WkDay]_[Month]_[Year]/Pkinase_PF00069/MODEL_REVISION/ #Directory used for Pkinase HMM profile revision; final revised domain model in RDS_output_[WkDay]_[Month]_[Year]/REVISED_MODELS/
  12. ProteomeA.Pkinase_regions_extracted.fasta #Pkinase domain sequences from ProteomeA.fasta
  13. [Repeat (12) for ProteomeB]
  14. PF00069_Pkinase.fasta #Pkinase base domain seed
  15. Pkinase_from_all_datasets.fasta #Pkinase domain seed + domains extracted from each dataset
  16. Pkinase_from_all_datasets.msa.trimmed.stockholm #Stockholm alignment of Pkinase domain with non-homologous ends trimmed

RDS_output_[WkDay]_[Month]_[Year]/Pkinase_PF00069/SECOND_SEARCH/ #Directory contains search for Pkinase domains following revision
  17. Pfam_with_Pkinase_revisions.hmm #Pfam database appended with revised domain profile
  18. Pfam_with_Pkinase_revisions.hmm.h3* #Four files associated with Pfam_with_Pkinase_revisions.hmm compression
  19. ProteomeA.hmmsearch_for_Pkinase_base_and_revised.domtblout #HMMsearch domain table output for potential Pkinase domains (either base or revised) in ProteomeA.fasta
  20. ProteomeA.Pkinase_base_and_revised_present.fasta #Fasta file containing sequences with potential Pkinase domains (base or revised) found in ProteomeA.fasta
  21. ProteomeA.Pkinase_base_and_revised_present.hmmscan_against_Pkinase_REV_appended_Pfam.domtblout #HMMscan output (domain table) of above sequences annotated with all Pfam domains (including revised Pkinase domain)
  22. ProteomeA.Pkinase_base_and_revised_present.hmmscan_against_Pkinase_REV_appended_Pfam.domtblout.besthits.tsv #Best-hits file associated with above file
  23. ProteomeA.Pkinase_base_and_revised_present.hmmscan_against_Pkinase_REV_appended_Pfam.target_list.txt #List of sequences with best-fit Pkinase domains (base or revised) in ProteomA.fasta
  24. ProteomeA.Pkinase.sequences_only_identified_after_revision.txt #Sequences in ProteomeA.fasta with Pkinase only found after revision
  25. [Repeat (19-24) for ProteomeB]
  
RDS_output_[WkDay]_[Month]_[Year]/Ubiquitin_PF00240/ #Output directory for Ubiquitin domain revision
  26. [Repeat (5-25) for Ubiquitin domain]
  
RDS_output_[WkDay]_[Month]_[Year]/REVISED_MODELS/ #Directory containing all revised models
  27. Pkinase_REVISION.hmm #Revised Pkinase domain HMM profile
  28. Pkinase_base_and_revision.hmm #Concatonation of base and revised Pkinase domain HMM profiles
  29. Pkinase_base_and_revision.hmm.h3* #Four files associated with Pkinase_base_and_revision.hmm compression
  30. [Repeat (27-29) for Ubiquitin domain]
  
RDS_output_[WkDay]_[Month]_[Year]/FINAL_HMMSCAN/ #Directory containing domain annotations following all domain revisions
  31. Pfam_with_all_revisions.hmm #Pfam database appanded with all revised domain profiles
  32. Pfam_with_all_revisions.hmm.h3* #Four files associated with Pfam_with_all_revisions.hmm compression
  33. ProteomeA.Pkinase_present_after_revision.fasta #Sequences associated with ProteomeA.Pkinase_base_and_revised_present.hmmscan_against_Pkinase_REV_appended_Pfam.target_list.txt above
  34. ProteomeA.Pkinase_present_after_revision.hmmscan_vs_revised_Pfam.domtblout #HMMscan output (domain table) of above sequences annotated with all Pfam domains (including all revised domains)
  35. ProteomeA.Pkinase_present_after_revision.hmmscan_vs_revised_Pfam.domtblout.besthits.tsv #Best-hits file associated with the above file
  36. [Repeat (33-35) for each ProteomeA+Ubiquitin, ProteomeB+Pkinase, and ProteomeB+Ubiquitin]
```
