# NCBI integration

AKStats retrieves GSE158634 from the official NCBI GEO distribution.

## Retrieval architecture

1. `NCBIGEOClient.search_geo_series()` queries the Entrez GEO DataSets database (`db=gds`) with ESearch.
2. `summarize_geo_series()` retrieves Entrez summaries with ESummary when a GEO DataSet UID is available.
3. `fetch_series_soft()` downloads `GSE158634_family.soft.gz` from the GEO FTP directory. This file is the source of sample titles, characteristics, exposure labels, animal identifiers (when deposited), platform and processing notes.
4. `fetch_count_matrix()` downloads `GSE158634_3.count.csv.gz`, the supplementary count matrix.
5. Linked SRA experiments (when present in GEO sample relations) are retrieved through EFetch (`db=sra`) to inspect instrument models, SRA sample accessions (`SRS`) and external BioSample accessions (`SAMN`). These identifier types are stored separately.
6. Linked BioSample accessions are enriched through EFetch (`db=biosample`) to recover explicit animal/treatment attributes.
7. `transcriptomics.run_transcriptomics()` reconciles the sample accessions, checks non-negative integer counts, exports sample-level QC and builds an animal-level design audit.

E-utilities provide Entrez metadata; the complete SOFT record and large count matrix are downloaded from GEO's official FTP endpoints, following NCBI's programmatic-access guidance.

## Run

The regular `akstats` command calls this integration as part of the reproducible pipeline:

```bash
pip install -e ".[dev]"
akstats
```

Or call the client directly:

```python
from pathlib import Path
from akstats.ncbi import NCBIGEOClient
from akstats.transcriptomics import run_transcriptomics

client = NCBIGEOClient()
resources = client.fetch_study(Path("data/raw"))
print(resources["entrez_uids"])

status = run_transcriptomics(Path("outputs"), client=client)
print(status["status"], status["inference_status"])
```

## Optional configuration

- `NCBI_EMAIL`: contact email included in Entrez requests.
- `NCBI_API_KEY`: optional NCBI API key. It is not required for normal low-volume use.

Do not commit API keys. Set them as environment variables or repository secrets. NCBI documents a default limit of 3 requests/second without a key and up to 10 requests/second with a key. The client enforces the corresponding minimum interval between sequential E-utilities requests, caches bulk GEO files, and retries transient HTTP errors.

## Scientific safeguards

- Sample groups are parsed from the family SOFT metadata; the expected GSM list is used only as a completeness check.
- Normal UVF and UVR samples are kept distinct when exposure is explicitly recorded.
- Verified animal IDs are extracted only from explicit animal/mouse/subject characteristics or linked BioSample attributes. The numeric suffix in a sample title is not treated as a verified animal identifier.
- A provisional crosswalk aligns GEO title numbering with Table 2 of PMID 33664254 (four UVF controls; four UVR animals with matched NS, AK-1, AK-2 and cSCC samples). It is exported separately for manual review and never used for inference. See the [detailed crosswalk evidence audit](GSE158634_SAMPLE_CROSSWALK_AUDIT.md) for the candidate GSM-to-paper-sample table and its assumptions.
- Multiple AK lesions from one animal are repeated observations, not independent biological replicates.
- No sample is removed solely because its library-size QC flag is raised.
- Gene-symbol mapping is exact and case-insensitive. Ensembl identifiers are not guessed.
- Genome assembly and annotation release are audited across source statements. The current SOFT processing text mentions both hg19 and mm9; this conflict is surfaced and must not be silently resolved.
- GEO/SRA report HiSeq 1500, the SOFT processing note mentions HiSeq 2500, and the article methods report NextSeq 500 (PMID 33664254). Verify this three-way instrument discrepancy before reuse.
- The current implementation reports descriptive log2-CPM summaries only. Inferential differential-expression modeling remains blocked until the animal-level design is reconstructed and reviewed.

## NCBI references

- [GEO programmatic access](https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html)
- [E-utilities documentation](https://www.ncbi.nlm.nih.gov/books/NBK25499/)
- [GSE158634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158634)
