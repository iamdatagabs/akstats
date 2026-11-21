import pandas as pd

from akstats.transcriptomics import SAMPLES, load_geo_counts, summarize_geo_counts


def test_geo_count_matrix_maps_all_samples_and_keeps_cohorts_separate(tmp_path):
    data = {"feature_id": ["Arid2", "Rora", "Mapk1"]}
    for accession in SAMPLES:
        data[accession] = [100, 50, 25]
    path = tmp_path / "counts.csv.gz"
    pd.DataFrame(data).to_csv(path, index=False, compression="gzip")
    counts, metadata = load_geo_counts(path)
    assert counts.shape == (3, 20)
    assert metadata["animal_id"].isna().all()
    assert metadata["stage"].value_counts().to_dict() == {
        "Normal skin": 8, "Actinic keratosis": 8, "Cutaneous SCC": 4
    }
    qc, summary, contrasts = summarize_geo_counts(counts, metadata)
    assert len(qc) == 20
    assert {"Arid2", "Rora", "Mapk1"} == set(summary["Gene"])
    assert set(contrasts["contrast"]) == {"AK_vs_Normal", "cSCC_vs_AK", "cSCC_vs_Normal"}
    assert contrasts["interpretation"].str.contains("not differential-expression").all()


def test_geo_identifier_mismatch_does_not_invent_gene_mapping(tmp_path):
    data = {"feature_id": ["ENSMUSG00000000001", "ENSMUSG00000000002"]}
    for accession in SAMPLES:
        data[accession] = [100, 50]
    path = tmp_path / "counts.csv.gz"
    pd.DataFrame(data).to_csv(path, index=False, compression="gzip")
    counts, metadata = load_geo_counts(path)
    qc, summary, contrasts = summarize_geo_counts(counts, metadata)
    assert len(qc) == 20
    assert summary.empty
    assert contrasts.empty
