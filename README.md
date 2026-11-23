# AKStats

**Reproducible exploratory bioinformatics analysis of the ARID2–MAPK axis across actinic keratosis (AK) and cutaneous squamous cell carcinoma (cSCC).**

AKStats converts the original Colab analysis into a modular Python package and produces a reproducible bundle of research outputs: CSV tables, high-resolution figures, an Excel workbook, and Markdown/HTML reports.

> **Research-use notice:** this repository is not a diagnostic tool, clinical decision-support system, or validated prognostic model. Results are exploratory and require interpretation alongside the source publication, study design, assay limitations, and independent biological evidence.

## Research portal

The AKStats portal is organized into two connected pages, published through GitHub Pages:

| Page | What you'll find | Link |
|---|---|---|
| **Project presentation** | Scientific question, research scope, analytical strategy, team, and project context. | [Open the landing page](https://iamdatagabs.github.io/akstats/) |
| **Interactive dashboard** | Aggregate indicators, filters by tissue stage and pathway module, interactive charts, underlying records, and CSV export. | [Open the dashboard](https://iamdatagabs.github.io/akstats/dashboard.html) |

The pages are connected through their navigation and calls to action: use **Explore the dashboard** to move from the presentation to the data explorer, or **Presentation / Know the project** to return to the research overview.

- **Repository:** https://github.com/iamdatagabs/akstats
- **GitHub Pages deployment:** https://github.com/iamdatagabs/akstats/actions

The dashboard is generated from the research pipeline and keeps analytical levels separate: group-level MAF summaries (ARID2 vs MAPK), gene-level MAF summaries (gene × stage, with module annotation), and module-level ΔMAF. The module selector filters gene-level records and the module chart; it does not relabel or filter group-level MAF summaries. The portal also presents panel size, source and filtered variant counts, patient counts, patient-level MAF summaries, prevalence above the operational 2% threshold with Wilson 95% confidence intervals, paired stage tests, and patient-level co-occurrence. The tissue stages (normal skin, AK, cSCC) are clinicopathological descriptors of the source cohort; they are not patient outcomes.

The current human mutation dataset does not provide individual-level age, sex, comorbidities, treatment response, recurrence, survival, or other longitudinal clinical outcomes. The independent GSE158634 mouse RNA-seq layer is retrieved directly from NCBI: Entrez E-utilities are used for GEO accession discovery/summary and linked BioSample enrichment, while the official GEO FTP distribution provides the family SOFT sample characteristics and supplementary integer count matrix. The live retrieval reconciled all 20 samples and 22,090 count features. Current SOFT metadata does not explicitly expose animal IDs or UVF/UVR assignments for each sample, so inferential differential-expression analysis remains blocked; sample-title suffixes are never used to invent animal links. The source processing metadata also contains a genome-build conflict (hg19 and mm9) and a sequencing-platform/instrument discrepancy (GPL18480 / HiSeq 1500 versus a HiSeq 2500 processing note), which must be resolved before inferential analysis. The dashboard therefore does not label MAF or ΔMAF as clinical endpoints and must not be used for diagnosis, prognosis, or treatment decisions.

## Research team

- **Gabriel Passos de Jesus, MSc.** — orientador; Departamento de Medicina, Universidade Estadual de Maringá (UEM), Brasil. Research coordination, computational engineering, bioinformatics workflow, statistical analysis, and reproducible software. [Currículo Lattes](http://lattes.cnpq.br/9928924222649468).
- **Giovanna Lyra** — pesquisadora; curso de graduação em Medicina, Universidade Santo Amaro (UNISA), Brasil. [Currículo Lattes](http://lattes.cnpq.br/8496468420504604).

## Research scope

- Functional variants in ARID2 and a predefined set of MAPK pathway genes.
- MAF across normal skin, actinic keratosis, and invasive cSCC.
- Variant-level change from AK to cSCC (ΔMAF).
- Patient-level co-occurrence at a configurable 2% MAF threshold.
- Patient-level MAF summaries, Wilson confidence intervals for threshold prevalence, and paired stage tests.
- Exploratory non-parametric statistics and pathway-module summaries.
- NCBI GEO integration for GSE158634 retrieves the family SOFT metadata and integer count matrix. The workflow reconciles all 20 GEO sample accessions, preserves normal UVF and UVR groups separately when exposure is available, audits repeated AK lesions by explicit animal ID, and exports sample-level QC. It does not infer animal IDs from title suffixes. Gene-symbol mapping is reported only when identifiers match exactly; genome assembly/annotation remain explicitly unresolved unless present in the source metadata. No differential-expression model is run unless the animal-level design can be reconstructed and reviewed.

## Data provenance

The mutation workflow downloads supplementary Excel files associated with DOI 10.1038/s41467-025-65687-y:

- [Supplementary panel (MOESM8)](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-65687-y/MediaObjects/41467_2025_65687_MOESM8_ESM.xlsx)
- [Supplementary mutation table (MOESM9)](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-65687-y/MediaObjects/41467_2025_65687_MOESM9_ESM.xlsx)
- [NCBI GEO GSE158634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158634)
- [GEO family SOFT metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158634/soft/GSE158634_family.soft.gz)
- [GEO supplementary gene-count matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE158nnn/GSE158634/suppl/GSE158634_3.count.csv.gz)
- [NCBI GEO programmatic access](https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html)

Original MOESM8/MOESM9 workbooks are versioned under `data/raw/`; processed outputs are generated separately. GSE158634 is a distinct mouse RNA-seq experiment (20 samples; GPL18480, Illumina HiSeq 1500), so it can provide independent biological context but cannot be joined to the human patient records. The NCBI client uses E-utilities for Entrez metadata and the official GEO FTP endpoints for the full SOFT record and count matrix. The FTP files are cached under `data/raw/`. Consult the source article for cohort definitions, experimental methods, and data-use terms.

## Quick start

Python 3.10 or newer is required.

```bash
git clone https://github.com/iamdatagabs/akstats.git
cd akstats
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\\Scripts\\Activate.ps1
pip install -e ".[dev]"
akstats
```

The first run requires internet access to retrieve the source files. Later runs reuse the local cache. Optional NCBI settings are `NCBI_EMAIL` and `NCBI_API_KEY`; the API key is not required for low-volume use. E-utilities requests identify the application as `akstats` and respect NCBI rate limits. Tests use synthetic fixtures and do not download research data:

```bash
pytest
ruff check src tests
```

## Final products

Running `akstats` creates:

| Artifact | Purpose |
|---|---|
| `outputs/AKStats_research_results.xlsx` | Consolidated workbook with filtered variants, descriptive MAF, modules, co-occurrence, Fisher tests, and study overview |
| `outputs/research_report.md` | Research report with dataset counts, tables, figures, provenance, and limitations |
| `outputs/research_report.html` | Portable HTML version of the report |
| `outputs/run_manifest.json` | Machine-readable record of source, run date, cohort counts, and outputs |
| `outputs/maf_descriptive.csv` | Descriptive MAF statistics by group and stage |
| `outputs/module_summary.csv` | MAF and ΔMAF summaries by pathway module |
| `outputs/patient_cooccurrence.csv` | Patient-level ARID2/MAPK presence matrix |
| `outputs/patient_level_summary.csv` | Patient-level MAF summaries and prevalence with Wilson 95% confidence intervals |
| `outputs/patient_level_tests.csv` | Paired Friedman/Wilcoxon tests and exact McNemar tests |
| `outputs/geo_transcriptomics_status.json` | NCBI retrieval, sample reconciliation, design and validation status for GSE158634 RNA-seq |
| `outputs/geo_sample_metadata.csv` | Sample-level metadata parsed from the official GEO family SOFT file |
| `outputs/geo_experimental_design.csv` | Animal-level repeated-measures audit; empty when animal IDs are not explicitly available |
| `outputs/geo_sample_qc.csv` | Library size and detected-feature QC for the 20 GEO samples |
| `outputs/geo_target_expression.csv` | Descriptive log2-CPM summaries for mapped ARID2–MAPK targets and Rora |
| `outputs/geo_target_contrasts.csv` | Differences between stage medians; not differential-expression tests |
| `outputs/fisher_tests.csv` | Fisher exact tests by stage |
| `outputs/maf_by_stage.png` | MAF distribution figure (300 dpi) |
| `outputs/delta_maf_by_module.png` | ΔMAF distribution by pathway module (300 dpi) |
| `data/processed/arid2_mapk_variants.csv` | Filtered and annotated variant records |

### Download from GitHub Actions

The **Build research products** workflow can be started manually from the Actions tab. It executes the pipeline and publishes generated files as a downloadable artifact named `AKStats-research-products` (30-day retention). The **Publish research portal** workflow regenerates the aggregate dashboard dataset and deploys the static site to GitHub Pages whenever relevant code, data, or site files change.

## Architecture

    akstats/
    ├── .github/workflows/
    │   ├── ci.yml
    │   └── research-products.yml
    ├── data/raw/                 # Download cache; not versioned
    ├── data/processed/           # Derived datasets; not versioned
    ├── outputs/                  # Tables, figures and reports; not versioned
    ├── src/akstats/
    │   ├── analysis.py            # Transformations and statistical summaries
    │   ├── cli.py                 # Command-line interface
    │   ├── config.py              # Sources, thresholds, pathway definitions
    │   ├── data.py                # Human mutation-source download and validation
    │   ├── ncbi.py                # NCBI E-utilities and GEO FTP client
    │   ├── transcriptomics.py     # GEO SOFT parsing, design audit and RNA-seq QC
    │   ├── pipeline.py            # End-to-end orchestration
    │   ├── plots.py               # Figure generation
    │   └── report.py              # Excel, Markdown, HTML and manifest exports
    ├── tests/
    ├── Fonte.ipynb                # Original exploratory notebook
    └── pyproject.toml

## Statistical interpretation and limitations

1. **Experimental unit:** multiple variants may belong to the same patient. Variant rows are not independent biological replicates. The legacy row-wise Friedman and Wilcoxon calculations are exploratory only. New patient-level summaries aggregate the maximum observed variant MAF per patient, group and stage; this is a descriptive operational measure, not total tumor burden. Paired patient-level tests use patients as the unit of analysis.
2. **MAF is not clone size:** ΔMAF is a descriptive change in observed variant allele fraction. Copy-number changes, purity, sampling, and assay characteristics can affect it.
3. **Threshold choice:** the 2% MAF threshold is an operational parameter, not a universal biological cutoff.
4. **Co-occurrence is not causation:** Fisher's exact test uses a small patient cohort. Report cell counts and uncertainty alongside p-values.
5. **Functional annotation:** consequence classes were inherited from the exploratory notebook and require review against the study's annotation protocol.
6. **GEO sample design:** normal skin from UVF and UVR conditions must remain distinct. Multiple AK lesions from the same animal are repeated observations, not independent biological replicates. Animal identifiers must come from explicit sample metadata; title suffixes are not animal IDs. No inferential RNA-seq model is run until the experimental unit and matched design are verified.
7. **RNA-seq identifiers and annotation:** the supplementary file is treated as integer counts only after numeric/integer validation. Gene-symbol matching is exact and case-insensitive; Ensembl or other identifiers are not guessed. Genome build and annotation release must be documented before gene-level biological interpretation.

## Original notebook

The initial exploratory workflow remains preserved in `Fonte.ipynb`. The Python package is the maintainable implementation; notebook outputs should not be treated as the authoritative reproducibility record.

## Citation

Cite this repository and the underlying source article when using the software. Add the final article citation and a software release DOI when a release is archived.

## License

MIT. This software license does not replace the terms governing third-party datasets.
