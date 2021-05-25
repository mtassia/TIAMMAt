# TIAMMAt: Taxon-Informed Adjustment of Markov-Model Attributes

#### Revising domain profile HMMs using empirical datasets to capture sequence variation in non-model species and improve homologous sequence identification.
#### TIAMMAt is currently under manuscript revision and is undergoing updates; please contact me for further detail if you would like to use the tool before publication.
---
### DESCRIPTION:
_Preamble: This pipeline and it's application to immune-gene evolution is in prep for publication_

In its current state, Pfam-A profile HMMs are derived of representative seed alignments encompassing curated sequences from select taxa. Due to taxonomic bias, domain seeds appear to reflect a heavy biomedical species bias; and, as such, standard Pfam-A domain models appear to underestimate the number of homologuos domains within non-model species transcriptome/genome datasets. TIAMMAt aims to improve species/sequence diversity intrinsically captured within individual Pfam domain profile HMMs.

##### The program is organized into three main blocks (detailed diagram can be seen below):
1) Loop among amino acid datasets for best-hit motifs to the target domain(s)
2) Loop through each domain profile HMM for revision to capture sequence variation of homologous domains within target datasets
3) Final scan among amino acid datasets for all Pfam-A entries included all revised domains

---
### PIPELINE:
![RelaxedDomainSearch Pipeline](https://github.com/mtassia/RelaxedDomainSearch/blob/master/Program_diagram.png)

---
### DEPENDENCIES:
- `Pull_coordinates.py` (Included in this repo; requires `Python3` + `BioPython`; `BioPython` can be installed with `pip install biopython`).
- `Select_contigs.pl` (Written by J.D. White; Included in this repo; requires `Perl`)
- `Best_fit_domains.py` (Included in this repo; requires `Python3` + `re` and `sys` modules)

- `HMMER` (Available at http://hmmer.org/; users should follow the installation instructions provided by HMMER prior to running TIAMMAt)

TIAMMAt has been tested and is compatible with `HMMER` versions 3.1 and 3.3.2

---
### USAGE:
It's recommended that for all path arguments, absolute paths be used.

```
-d [directory]  Specify directory containing protein datasets [DEFAULT: pwd]
-h              Print help
-m [directory]  Specify directory containing pfam models and seeds [MANDATORY; cannot be the same directory as datasets]
-o [string]     Set name for output directory [default: TIAMMAt_output_{DATE}]
-p [path]       Specify path to Pfam database [MANDATORY]
-t [integer]    Specify number of threads for hmmscan and hmmsearch
```

**Example Command:**
```
tiammat -d /path/to/Proteomes/ -m /path/to/Target_Pfams/ -p /path/to/Pfam-A.hmm -t 4
```

---
### INPUTS:
For each domain of interest, the seed (an unaligned fasta file obtainable from the "**Alignments**" section for any domain in Pfam) and the model (raw HMM obtained from the "**Curation & model**" section for any given domain in Pfam) must be downloaded from the Pfam server (http://pfam.xfam.org/). Within the model directory (`-m`), the prefix naming convention must be identical for each domain (e.g., *PF00069_Pkinase*.fasta & *PF00069_Pkinase*.hmm). `Grab_models.sh` can be used to generate the files required in the model directory (see **Support Scripts** section below). Additionally, ensure the domain models and the Pfam database input (`-p`) are the same version before running TIAMMAt.

**Example input structure:**
- `-d /path/to/Proteomes/` contains:
  >- ProteomeA.fasta
  >- ProteomeB.fasta
- `-m /path/to/Target_Pfams/` contains:
  >- PF00069_Pkinase.fasta
  >- PF00069_Pkinase.hmm
  >- PF00240_Ubiquitin.fasta
  >- PF00240_Ubiquitin.hmm
- `[-o default]`
- `-p /path/to/Pfam.hmm`

Absolute paths for `-d`, `-m`, and `-p` are recommended.

`tiammat` *should not* be run with nucleotide input - cannot currently detect the sequence alphabet of the input.

See *Known Issues* section below.

### OUTPUTS:
**Output structure:**
```bash
TIAMMAt_output_[WkDay]_[Month]_[Year]/ #Default output directory created by program
  1. IDENTIFICATION_STATISTICS.txt #TSV containing domain counts before & after revision
  2. SCRIPT_LOG.txt #Data log reporting operations and immediate results during program execution

TIAMMAt_output_[WkDay]_[Month]_[Year]/MODELS/ #Directory containing compressed models for initial HMMsearch step
  3. PF00069_Pkinase.hmm.h3* #Four files associated with initial Pkinase domain compression from HMMpress
  4. [Repeat (3) for Ubiquitin domain]

TIAMMAt_output_[WkDay]_[Month]_[Year]/Pkinase_PF00069/ #Output directory for Pkinase domain revision
  5. ProteomeA.hmmsearch_for_Pkinase.domtblout #HMMsearch domain table output for potential Pkinase domains in ProteomeA.fasta
  6. ProteomeA.hmmsearch.Pkinase_present.fasta #Fasta file containing sequences with potential Pkinase domains found in ProteomeA.fasta from HMMsearch
  7. ProteomeA.Pkinase_present.hmmscan_against_pfam.domtblout #HMMscan output (domain table) of above sequences annotated with all Pfam domains
  8. ProteomeA.Pkinase_present.hmmscan_against_pfam.domtblout.besthits.tsv #Above table filtered for domains meeting HMMer inclusion threshold and overlap
  9. ProteomeA.Pkinase_present.hmmscan_against_pfam.Pkinase_coordinates.tsv #Coordinate TSV for all best-fit Pkinase domains in ProteomeA.fasta
  10. ProteomeA.Pkinase_present.hmmscan_against_pfam.target_list.txt #List of sequences with best-fit Pkinase domains in ProteomA.fasta
  11. [Repeat (5-10) for ProteomeB]

TIAMMAt_output_[WkDay]_[Month]_[Year]/Pkinase_PF00069/MODEL_REVISION/ #Directory used for Pkinase HMM profile revision; final revised domain model in TIAMMAt_output_[WkDay]_[Month]_[Year]/REVISED_MODELS/
  12. ProteomeA.Pkinase_regions_extracted.fasta #Pkinase domain sequences from ProteomeA.fasta
  13. [Repeat (12) for ProteomeB]
  14. PF00069_Pkinase.fasta #Pkinase base domain seed
  15. Pkinase_from_all_datasets.fasta #Pkinase domain seed + domains extracted from each dataset
  16. Pkinase_from_all_datasets.msa.trimmed.stockholm #Stockholm alignment of Pkinase domain with non-homologous ends trimmed

TIAMMAt_output_[WkDay]_[Month]_[Year]/Pkinase_PF00069/SECOND_SEARCH/ #Directory contains search for Pkinase domains following revision
  17. Pfam_with_Pkinase_revisions.hmm #Pfam database appended with revised domain profile
  18. Pfam_with_Pkinase_revisions.hmm.h3* #Four files associated with Pfam_with_Pkinase_revisions.hmm compression
  19. ProteomeA.hmmsearch_for_Pkinase_base_and_revised.domtblout #HMMsearch domain table output for potential Pkinase domains (either base or revised) in ProteomeA.fasta
  20. ProteomeA.hmmsearch.Pkinase_base_and_revised_present.fasta #Fasta file containing sequences with potential Pkinase domains (base or revised) found in ProteomeA.fasta from HMMsearch
  21. ProteomeA.Pkinase_base_and_revised_present.hmmscan_against_Pkinase_REV_appended_Pfam.domtblout #HMMscan output (domain table) of above sequences annotated with all Pfam domains (including revised Pkinase domain). IMPORTANT: Sequences annotated in this file do not necessarily contain a best-hit target domain.
  22. ProteomeA.Pkinase_base_and_revised_present.hmmscan_against_Pkinase_REV_appended_Pfam.domtblout.besthits.tsv #Best-hits file associated with above file
  23. ProteomeA.Pkinase_base_and_revised_present.hmmscan_against_Pkinase_REV_appended_Pfam.target_list.txt #List of sequences with best-fit Pkinase domains (base or revised) in ProteomA.fasta
  24. ProteomeA.Pkinase.sequences_only_identified_after_revision.txt #Sequences in ProteomeA.fasta with Pkinase only found after revision
  25. [Repeat (19-24) for ProteomeB]

TIAMMAt_output_[WkDay]_[Month]_[Year]/Ubiquitin_PF00240/ #Output directory for Ubiquitin domain revision
  26. [Repeat (5-25) for Ubiquitin domain]

TIAMMAt_output_[WkDay]_[Month]_[Year]/REVISED_MODELS/ #Directory containing all revised models
  27. Pkinase_REVISION.hmm #Revised Pkinase domain HMM profile
  28. Pkinase_base_and_revision.hmm #Concatonation of base and revised Pkinase domain HMM profiles
  29. Pkinase_base_and_revision.hmm.h3* #Four files associated with Pkinase_base_and_revision.hmm compression
  30. [Repeat (27-29) for Ubiquitin domain]

TIAMMAt_output_[WkDay]_[Month]_[Year]/FINAL_HMMSCAN/ #Directory containing domain annotations following all domain revisions
  31. Pfam_with_all_revisions.hmm #Pfam database appanded with all revised domain profiles
  32. Pfam_with_all_revisions.hmm.h3* #Four files associated with Pfam_with_all_revisions.hmm compression
  33. ProteomeA.Pkinase_present_after_revision.fasta #Sequences associated with ProteomeA.Pkinase_base_and_revised_present.hmmscan_against_Pkinase_REV_appended_Pfam.target_list.txt above -- i.e, sequences with best-hits to either base or revised domain
  34. ProteomeA.Pkinase_present_after_revision.hmmscan_vs_revised_Pfam.domtblout #HMMscan output (domain table) of above sequences annotated with all Pfam domains (including all revised domains)
  35. ProteomeA.Pkinase_present_after_revision.hmmscan_vs_revised_Pfam.domtblout.besthits.tsv #Best-hits file associated with the above file
  36. [Repeat (33-35) for each ProteomeA+Ubiquitin, ProteomeB+Pkinase, and ProteomeB+Ubiquitin]
```
**Major Output Explanation:**
Summary statistics in `IDENTIFICATION_STATISTICS.txt` are a useful metric for checking the results of TIAMMAt. For easy viewing, try `column -t IDENTIFICATION_STATISTICS.txt`. The column headers of `IDENTIFICATION_STATISTICS.txt` are as follows:
* 1 = Dataset ID
* 2 = Domain ID
* 3 = Number of sequences from pre-revision dataset with best-fit target domain
* 4 = Number of domain occurances within column (3)
* 5 = Number of sequences in post-revision dataset with best-fit target domain (base or revised)
* 6 = Number of base domain occurances within column (5)
* 7 = Number of revised domain occurances within column (6)

Users primarily interested in the revised domain profile-HMMs can find their revised models in the `REVISED_MODELS` directory, and those interested in the number of domain-containing proteins per dataset can find those data in the `FINAL_HMMSCAN` directory.

To check input, the top lines of `SCRIPT_LOG.txt` report the input files.

---
### KNOWN ISSUES:
* `-d` argument cannot take current directory shortcut, `./` (i.e., `./[directory_name]` does not work); `-d` works with absolute path and assumes `pwd` as relative path (i.e., `[directory_name]` is `./[directory_name]`)
* If base domain accession used for input does not match the accession present in Pfam (e.g., if TIR domain accession PF01582.22 is a target domain (`-m`), but not present in the local Pfam database (`-p`)), TIAMMAt will complete its run with no findings (even if homologous domains are present). See **INPUTS** section above.

---
### UPCOMING CHANGES:
* Incorporating a compression loop into `tiammat` to tarball intermediate model revision directories to reduce final output size.
* Input QC testing for input amino acid file, not nucleotide.
* Improvements to error reporting.

---
### SUPPORT SCRIPTS:
These programs must be manually executed - they are not run by TIAMMAt. Organized alphabetically below.

**`Domain_svgwrite.py`:** Uses the `[name].besthits.tsv` input from TIAMMAt and generates a domain diagram object per annotated sequence into a single editable svg canvas.
* *DEPENDENCIES:* `Python3`, `re`, `sys`, `svgwrite`
* *USAGE:* `python3 Domain_svgwrite.py [name].besthit.tsv`
* *OUTPUT:* `[name].besthits.tsv.domain_diagram.svg`

![Domain_diagram](https://github.com/mtassia/RelaxedDomainSearch/blob/master/Domain_SVGwrite_example.PNG)

**`ID_novel_seqs.sh`:** Can be executed in the `SECOND_SEARCH` directory to rename headers to incorperate the revised domain and flag newly identified sequences with a 'NOV' tag *[!DEPRECATED!]*.
* *USAGE:* Execute command in `SECOND_SEARCH/` directory
* *OUTPUT:* `[Taxon_ID].[Domain_ID]_base_and_revised_present.IDs_renamed.fasta`

**`Grab_models.sh`:** Reads a list of pfam accessions and automatically downloads both the input `*.hmm` and seed fasta files into `pwd`. Automatically runs `Stockholm2fasta.py` below. This script is particularly useful when revising multiple Pfam models.
* *DEPENDENCIES:* `Python3`, `BioPython` (both required for `Stockholm2fasta.py` which is run by `Grab_models.sh`)
* *USAGE:* `Grab_models.sh [pfam_accession_list.txt]`
* *INPUT FILE EXAMPLE:*
```
PF00554
PF00605
PF00619
PF01582
PF05729
PF10401
PF11648
PF12721
PF13676
```
* *OUTPUT:* `[Pfam_accession]_[Pfam_name].hmm` & `[Pfam_accession]_[Pfam_name].fasta` for every accession in `[pfam_accession_list.txt]`; example:
```
PF00554_RHD_DNA_bind.fasta
PF00554_RHD_DNA_bind.hmm
PF00605_IRF.fasta
PF00605_IRF.hmm
PF00619_CARD.fasta
PF00619_CARD.hmm
PF01582_TIR.fasta
PF01582_TIR.hmm
PF05729_NACHT.fasta
PF05729_NACHT.hmm
PF10401_IRF-3.fasta
PF10401_IRF-3.hmm
PF11648_RIG-I_C-RD.fasta
PF11648_RIG-I_C-RD.hmm
PF12721_RHIM.fasta
PF12721_RHIM.hmm
PF13676_TIR_2.fasta
PF13676_TIR_2.hmm
```

**`Stockholm2fasta.py`:** Reads input stockholm alignment (in this case, the seed alignment downloaded from Pfam) and converts it to an unaligned fasta file. This can be used independently of `Grab_models.sh`
* *DEPENDENCIES:* `Python3`, `BioPython`
* *USAGE:* `Stockholm2fasta.py [input_stockholm] [output_fasta]`
