# GSE158634 — audit of the GEO-to-publication sample crosswalk

**Status: candidate mapping only; not author-confirmed.**  
**Dataset:** [GSE158634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158634)  
**Publication:** Zhang et al., *Cell Death & Disease* (2021), PMID [33664254](https://pubmed.ncbi.nlm.nih.gov/33664254/), DOI [10.1038/s41419-021-03525-x](https://doi.org/10.1038/s41419-021-03525-x).

## What the sources establish

The GEO series lists 20 GSM samples in four title blocks:

- 8 normal-skin samples: `NS_001`, `NS_002`, `NS_004`, `NS_005`, `NS_006`–`NS_009`;
- 8 actinic-keratosis samples: `AK_010`–`AK_017`;
- 4 cutaneous SCC samples: `cSCC_018`–`cSCC_021`.

The paper's Table 2 describes eight SKH-1 mice: four UVF controls (one normal-skin specimen per animal) and four UVR animals (each contributing normal skin, AK-1, AK-2 and cSCC). This yields exactly the same tissue counts as the GEO series.

**Important limitation:** Table 2 does not print GSM accessions, and the GEO SOFT/BioSample records do not contain an explicit animal identifier for each sample. The correspondence below is therefore inferred from the order and numbering of GEO titles relative to the order of specimens in Table 2. It is useful for review and sensitivity planning, but it is not an independently verified accession-to-animal key.

## Candidate correspondence

| GEO title | GSM | Candidate paper subject | Candidate paper sample | Candidate role |
|---|---|---:|---|---|
| NS_001 | GSM4804794 | 1 | NC-1 | UVF normal |
| NS_002 | GSM4804795 | 2 | NC-2 | UVF normal |
| NS_004 | GSM4804796 | 3 | NC-3 | UVF normal |
| NS_005 | GSM4804797 | 4 | NC-4 | UVF normal |
| NS_006 | GSM4804798 | 5 | 1-NC | UVR normal |
| NS_007 | GSM4804799 | 6 | 2-NC | UVR normal |
| NS_008 | GSM4804800 | 7 | 3-NC | UVR normal |
| NS_009 | GSM4804801 | 8 | 4-NC | UVR normal |
| AK_010 | GSM4804802 | 5 | 1-AK-1 | UVR AK-1 |
| AK_011 | GSM4804803 | 5 | 1-AK-2 | UVR AK-2 |
| AK_012 | GSM4804804 | 6 | 2-AK-1 | UVR AK-1 |
| AK_013 | GSM4804805 | 6 | 2-AK-2 | UVR AK-2 |
| AK_014 | GSM4804806 | 7 | 3-AK-1 | UVR AK-1 |
| AK_015 | GSM4804807 | 7 | 3-AK-2 | UVR AK-2 |
| AK_016 | GSM4804808 | 8 | 4-AK-1 | UVR AK-1 |
| AK_017 | GSM4804809 | 8 | 4-AK-2 | UVR AK-2 |
| cSCC_018 | GSM4804810 | 5 | 1-cSCC | UVR cSCC |
| cSCC_019 | GSM4804811 | 6 | 2-cSCC | UVR cSCC |
| cSCC_020 | GSM4804812 | 7 | 3-cSCC | UVR cSCC |
| cSCC_021 | GSM4804813 | 8 | 4-cSCC | UVR cSCC |

The candidate animal labels `UVR_1`–`UVR_4` in the generated crosswalk correspond to paper subjects 5–8, respectively. They are local analysis labels, not identifiers reported by the authors.

## Consequences for statistical design

If the candidate correspondence is later confirmed:

1. The UVF group has four independent animal-level normal-skin observations.
2. The UVR group has four animals, each with four tissue specimens (normal, AK-1, AK-2 and cSCC).
3. AK-1 and AK-2 from the same animal are repeated lesions, not independent biological replicates.
4. A model must account for within-animal dependence. Treating all 20 samples as independent would be pseudoreplication.
5. The current pipeline must continue to keep candidate animal IDs in separate `publication_candidate_*` columns. They must not populate verified `animal_id` or `exposure_group`, and must not unlock inferential differential-expression analysis.

## Remaining validation requests

A definitive mapping requires at least one of:

- a GSM-to-subject/sample key supplied by the corresponding author;
- a GEO/SRA/BioSample update containing explicit animal/sample identifiers;
- original sample-sheet or sequencing-core metadata that links library IDs to the specimens in Table 2.

The article and GEO also disagree on sequencing instrument (article methods: NextSeq 500; GEO platform: HiSeq 1500; SOFT processing text: HiSeq 2500). The processing text also mentions both `hg19` and `mm9`. These are separate provenance issues and are not resolved by the candidate sample crosswalk.

## Pipeline implementation

The current implementation in `src/akstats/transcriptomics.py` writes the candidate crosswalk to `geo_publication_candidate_crosswalk.csv` and the candidate repeated-measures audit to `geo_publication_candidate_design.csv`. The status remains explicitly provisional, and the verified experimental-design audit remains blocked when explicit animal IDs are absent.
