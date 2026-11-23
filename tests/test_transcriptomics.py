
import pandas as pd

from akstats.transcriptomics import (
    _metadata_consistency,
    EXPECTED_SAMPLES,
    experimental_design,
    load_geo_counts,
    merge_biosample_metadata,
    merge_sra_metadata,
    parse_geo_family_soft,
    summarize_geo_counts,
)


def _soft_fixture(include_animals=True):
    lines = [
        "^SERIES = GSE158634",
        "!Series_title = UV-induced cSCC progression",
        "!Series_overall_design = Normal skin, AK and cSCC samples",
    ]
    samples = []
    accessions = sorted(EXPECTED_SAMPLES)
    for index, accession in enumerate(accessions):
        if index < 4:
            title, stage, exposure, animal = f"NS_UVF_{index + 1:02}", "normal skin", "UVF", None
        elif index < 8:
            animal_num = index - 3
            title, stage, exposure, animal = f"NS_UVR_{animal_num:02}", "normal skin", "UVR", animal_num
        elif index < 16:
            animal_num = (index - 8) // 2 + 1
            lesion_num = (index - 8) % 2 + 1
            title, stage, exposure, animal = f"AK_{animal_num:02}_{lesion_num}", "actinic keratosis", "UVR", animal_num
        else:
            animal_num = index - 15
            title, stage, exposure, animal = f"cSCC_{animal_num:02}", "cutaneous squamous cell carcinoma", "UVR", animal_num
        samples.append((accession, title, stage, exposure, animal))
    for accession, title, stage, exposure, animal in samples:
        lines.extend([
            f"^SAMPLE = {accession}",
            f"!Sample_title = {title}",
            f"!Sample_source_name_ch1 = {stage}",
            f"!Sample_characteristics_ch1 = tissue: {stage}",
            f"!Sample_characteristics_ch1 = exposure: {exposure}",
        ])
        if include_animals and animal is not None:
            lines.append(f"!Sample_characteristics_ch1 = animal_id: Mouse {animal:02}")
        lines.append(f"!Sample_relation = BioSample: SAMN000000{len(lines):02}")
        lines.append(f"!Sample_relation = SRA: SRX{len(lines):07}")
        lines.append("!Sample_platform_id = GPL18480")
    return "\n".join(lines) + "\n"


def test_soft_parser_extracts_sample_metadata_without_guessing_animal_ids():
    metadata, series = parse_geo_family_soft(_soft_fixture(include_animals=True))
    assert series["accession"] == "GSE158634"
    assert len(metadata) == 20
    assert metadata["sample_accession"].nunique() == 20
    assert metadata["stage"].value_counts().to_dict() == {
        "Normal skin": 8, "Actinic keratosis": 8, "Cutaneous SCC": 4
    }
    assert metadata["analysis_group"].value_counts()["Normal skin (UVF)"] == 4
    assert metadata["analysis_group"].value_counts()["Normal skin (UVR)"] == 4
    assert metadata.loc[metadata["stage"].eq("Actinic keratosis"), "animal_id"].notna().all()
    assert metadata.loc[metadata["stage"].eq("Normal skin"), "animal_id"].isna().sum() == 4


def test_design_audit_recognizes_repeated_ak_lesions_within_four_animals():
    metadata, _ = parse_geo_family_soft(_soft_fixture(include_animals=True))
    design, status = experimental_design(metadata)
    assert len(design) == 4
    assert design["n_ak_lesions"].eq(2).all()
    assert design["complete_matched_set"].all()
    assert status["design_ready_for_paired_model"] is True


def test_missing_animal_ids_block_inference_instead_of_inventing_pairs():
    metadata, _ = parse_geo_family_soft(_soft_fixture(include_animals=False))
    design, status = experimental_design(metadata)
    assert design.empty
    assert status["animal_id_available"] is False
    assert status["inference_status"] == "blocked_missing_animal_id"


def test_counts_join_to_geo_accessions_and_keep_uvf_uvr_normals_separate(tmp_path):
    metadata, _ = parse_geo_family_soft(_soft_fixture(include_animals=True))
    data = {"feature_id": ["Arid2", "Rora", "Mapk1"]}
    for accession in metadata["sample_accession"]:
        data[accession] = [100, 50, 25]
    data[metadata["sample_accession"].iloc[0]] = [1, 1, 1]
    path = tmp_path / "counts.csv.gz"
    pd.DataFrame(data).to_csv(path, index=False, compression="gzip")

    counts, joined = load_geo_counts(path, metadata)
    assert counts.shape == (3, 20)
    assert set(counts.columns) == EXPECTED_SAMPLES
    qc, summary, contrasts = summarize_geo_counts(counts, joined)
    assert len(qc) == 20
    flagged = qc.loc[qc["sample_accession"].eq(metadata["sample_accession"].iloc[0]), "qc_flag"].iloc[0]
    assert flagged == "review_low_library_size"
    assert {"Arid2", "Rora", "Mapk1"} == set(summary["Gene"])
    assert {"Normal skin (UVF)", "Normal skin (UVR)"} <= set(summary["Analysis_group"])
    assert contrasts["interpretation"].str.contains("not differential-expression").all()


def test_identifier_mismatch_does_not_invent_gene_mapping(tmp_path):
    metadata, _ = parse_geo_family_soft(_soft_fixture(include_animals=True))
    data = {"feature_id": ["ENSMUSG00000000001", "ENSMUSG00000000002"]}
    for accession in metadata["sample_accession"]:
        data[accession] = [100, 50]
    path = tmp_path / "counts.csv.gz"
    pd.DataFrame(data).to_csv(path, index=False, compression="gzip")
    counts, joined = load_geo_counts(path, metadata)
    qc, summary, contrasts = summarize_geo_counts(counts, joined)
    assert len(qc) == 20
    assert summary.empty
    assert contrasts.empty



def test_biosample_metadata_enriches_only_explicit_exposure_and_animal_fields():
    metadata, _ = parse_geo_family_soft(_soft_fixture(include_animals=False))
    accession = metadata.loc[0, "sample_accession"]
    biosample_accession = metadata.loc[0, "biosample_accession"]
    enriched = merge_biosample_metadata(metadata, {
        biosample_accession: {
            "accession": biosample_accession,
            "attributes": {"animal_id": ["Mouse 07"], "treatment": ["UVR"]},
        }
    })
    row = enriched.loc[enriched["sample_accession"].eq(accession)].iloc[0]
    assert row["animal_id"] == "07"
    assert row["exposure_group"] == "UVR"
    assert row["analysis_group"] == "Normal skin (UVR)"


def test_metadata_consistency_flags_conflicting_genome_and_instrument_records():
    metadata = pd.DataFrame({
        "data_processing": [
            "mapped to human genome (hg19); Genome_build: mm9; HiSeq 2500 instrument"
        ],
        "platform_id": ["GPL18480"],
    })
    status = _metadata_consistency(metadata)
    assert status["genome_builds_detected"] == ["hg19", "mm9"]
    assert status["genome_build_consistency"] == "conflicting_source_metadata"
    assert status["platform_instrument_warning"] is not None



def test_sra_metadata_enriches_only_explicit_experiment_relations():
    metadata, _ = parse_geo_family_soft(_soft_fixture(include_animals=False))
    first = metadata.iloc[0]
    enriched = merge_sra_metadata(metadata, {
        first["sra_experiment_accession"]: {
            "accession": first["sra_experiment_accession"],
            "instrument_model": "Illumina HiSeq 2500",
            "biosample_accession": "SAMN12345678",
        }
    })
    row = enriched.loc[enriched["sample_accession"].eq(first["sample_accession"])].iloc[0]
    assert row["sra_instrument_model"] == "Illumina HiSeq 2500"
    assert row["biosample_accession"] == "SAMN12345678"
