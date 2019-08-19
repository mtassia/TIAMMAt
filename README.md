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
  IDENTIFICATION_STATISTICS.txt #TSV containing domain counts before & after revision
  SCRIPT_LOG.txt #Data log reporting operations and immediate results during program execution
  
RDS_output_[WkDay]_[Month]_[Year]/MODELS/ #Directory containing compressed models for initial HMMsearch step
  PF00069_Pkinase.hmm.* #Four files associated with initial domain compression from HMMpress
  PF00420_Ubiquitin.hmm.* #Four files associated with domain initial compression from HMMpress
  
RDS_output_[WkDay]_[Month]_[Year]/Pkinase_PF00069/ #Output directory for Pkinase domain revision
  ProteomeA.hmmsearch_for_Pkinase.domtblout #HMMsearch domain table output for potential Pkinase domains in ProteomeA.fasta
  ProteomeA.Pkinase_present.fasta #Fasta file containing sequences with potential Pkinase domains found in ProteomeA.fasta
  ProteomeA.Pkinase_present.hmmscan_against_pfam.domtblout #HMMscan output (domain table) of above sequences annotated with all Pfam domains
  ProteomeA.Pkinase_present.hmmscan_against_pfam.domtblout.besthits.tsv #Above table filtered for domains meeting HMMer inclusion threshold and overlap
  ProteomeA.Pkinase_present.hmmscan_against_pfam.Pkinase_coordinates.tsv #Coordinate TSV for all best-fit Pkinase domains in ProteomeA.fasta
  ProteomeA.Pkinase_present.hmmscan_against_pfam.target_list.txt #List of sequences with best-fit Pkinase domains in ProteomA.fasta
  [Repeat above files for ProteomeB]
```
