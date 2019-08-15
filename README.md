# RelaxedDomainSearch
---

Code associated with domain HMM model revision for use in non-model species protein evolution studies.

Description: [More to come]

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

### OUTPUTS:

