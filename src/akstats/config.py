"""Centralized project configuration."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
OUTPUT_DIR = ROOT / "outputs"

PANEL_URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-65687-y/MediaObjects/41467_2025_65687_MOESM8_ESM.xlsx"
VARIANTS_URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-65687-y/MediaObjects/41467_2025_65687_MOESM9_ESM.xlsx"

PANEL_SHEET = "v.3.0TargetGenes"
VARIANTS_SHEET = "TableS5"
PANEL_HEADER_ROW = 2
VARIANTS_HEADER_ROW = 1
CLONALITY_THRESHOLD = 0.02
MAF_COLUMNS = ("Normal_MAF", "AK_MAF", "SCC_MAF")

MAPK_MODULES = {
    "0. Chromatin (ARID2)": ["ARID2"],
    "1. RTK (upstream)": ["EGFR"],
    "2. RAS GTPases": ["HRAS", "KRAS", "NRAS", "MRAS", "RRAS", "RRAS2", "RASA1"],
    "3. MAP3K / RAF": ["ARAF", "BRAF", "RAF1", "MAP3K1", "MAP3K2", "MAP3K5", "MAP3K7", "MAP3K9"],
    "4. MAP2K / MEK": ["MAP2K1", "MAP2K2", "MAP2K4"],
    "5. MAPK effectors": ["MAPK1", "MAPK8IP1"],
}
GENE_TO_MODULE = {gene: module for module, genes in MAPK_MODULES.items() for gene in genes}
FUNCTIONAL_CLASSES = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Splice_Site", "In_Frame_Del", "Start_Codon_SNP", "De_novo_Start_InFrame",
}
