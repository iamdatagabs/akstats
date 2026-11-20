# AKStats

**Reproducible exploratory bioinformatics analysis of the ARID2–MAPK axis across actinic keratosis (AK) and cutaneous squamous cell carcinoma (cSCC).**

AKStats refactors the original Google Colab workflow into a modular Python package. It retrieves public supplementary mutation tables, validates source schemas, annotates functional variants, summarizes mutant allele fractions (MAF), evaluates patient-level co-occurrence, and exports tables and figures.

> **Research-use notice:** this repository is not a diagnostic tool, clinical decision-support system, or validated prognostic model. Results are exploratory and require interpretation alongside the source publication, study design, assay limitations, and independent biological evidence.

## Research scope

- Functional variants in ARID2 and a predefined set of MAPK pathway genes.
- MAF across normal skin, actinic keratosis, and invasive cSCC.
- Variant-level change from AK to cSCC (ΔMAF).
- Patient-level co-occurrence at a configurable 2% MAF threshold.
- Exploratory non-parametric statistics and pathway-module summaries.
- A GEO series-matrix parser for GSE158634. Probe IDs are retained as probes; gene-level interpretation requires documented platform annotation.

## Data provenance

The mutation workflow downloads supplementary Excel files associated with DOI 10.1038/s41467-025-65687-y:

- [Supplementary panel (MOESM8)](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-65687-y/MediaObjects/41467_2025_65687_MOESM8_ESM.xlsx)
- [Supplementary mutation table (MOESM9)](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-65687-y/MediaObjects/41467_2025_65687_MOESM9_ESM.xlsx)
- [NCBI GEO GSE158634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158634)

Raw files are cached locally and are not committed. Consult the source article for cohort definitions, experimental methods, and data-use terms.

## Architecture

    akstats/
    ├── .github/workflows/ci.yml
    ├── data/raw/                 # Download cache; not versioned
    ├── data/processed/            # Derived datasets; not versioned
    ├── outputs/                   # Tables and figures; not versioned
    ├── src/akstats/
    │   ├── analysis.py            # Transformations and statistical summaries
    │   ├── cli.py                 # Command-line interface
    │   ├── config.py              # Sources, thresholds, pathway definitions
    │   ├── data.py                # Download, parsing, schema validation
    │   ├── pipeline.py            # End-to-end orchestration
    │   └── plots.py               # Figure generation
    ├── tests/
    ├── Fonte.ipynb                # Original exploratory notebook
    └── pyproject.toml

## Quick start

Python 3.10 or newer is required.

    git clone https://github.com/iamdatagabs/akstats.git
    cd akstats
    python -m venv .venv

Activate the environment, then install and run:

    # Linux/macOS
    source .venv/bin/activate
    # Windows PowerShell
    .venv\Scripts\Activate.ps1

    pip install -e ".[dev]"
    akstats

The first run requires internet access to retrieve the source files. Later runs reuse the local cache. Tests use synthetic fixtures and do not download research data:

    pytest
    ruff check src tests

## Outputs

| Artifact | Description |
|---|---|
| data/processed/arid2_mapk_variants.csv | Filtered and annotated variant records |
| outputs/maf_descriptive.csv | Descriptive MAF statistics by group and stage |
| outputs/module_summary.csv | MAF and ΔMAF summaries by pathway module |
| outputs/patient_cooccurrence.csv | Patient-level ARID2/MAPK presence matrix |
| outputs/fisher_tests.csv | Fisher exact tests by stage |
| outputs/exploratory_tests.json | Exploratory Friedman and Wilcoxon results |
| outputs/maf_by_stage.png | MAF distribution figure |
| outputs/delta_maf_by_module.png | ΔMAF distribution by pathway module |

## Statistical interpretation and limitations

1. **Experimental unit:** multiple variants may belong to the same patient. Variant rows are not independent biological replicates. The row-wise Friedman and Wilcoxon calculations retained from the notebook are labelled exploratory and must not be interpreted as patient-level inference.
2. **MAF is not clone size:** ΔMAF is a descriptive change in observed variant allele fraction. Copy-number changes, purity, sampling, and assay characteristics can affect it.
3. **Threshold choice:** the 2% MAF threshold is an operational parameter, not a universal biological cutoff.
4. **Co-occurrence is not causation:** Fisher's exact test uses a small patient cohort. Report cell counts and uncertainty alongside p-values.
5. **Functional annotation:** consequence classes were inherited from the exploratory notebook and require review against the study's annotation protocol.
6. **GEO probes are not genes:** expression data must be mapped using the correct GPL platform annotation before making claims about ARID2 or MAPK gene expression. Selecting the most variable probes does not establish gene-specific validation.

## Original notebook

The initial exploratory workflow remains preserved in Fonte.ipynb. The Python package is the maintainable implementation; notebook outputs should not be treated as the authoritative reproducibility record.

## Citation

Cite this repository and the underlying source article when using the software. Add the final article citation and a software release DOI when a release is archived.

## License

MIT. This software license does not replace the terms governing third-party datasets.
