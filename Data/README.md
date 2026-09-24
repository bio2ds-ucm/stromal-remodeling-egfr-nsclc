# Data

Data files are **not versioned** in this repository. This README describes how
to obtain each dataset and where to place the files so the scripts find them.

## 1. Clinical cohorts (OSIRESP, OSIREAL, H12O TMA004)

Bruker GeoMx Digital Spatial Profiler data of EGFR-mutant NSCLC patient biopsies.

**Availability.** Data will be deposited in the Gene Expression Omnibus (GEO)
upon publication. The GEO accession number will be added here at that time.

**Expected layout after download:**

```
Data/
├── OSIRESP_cohort/
│   ├── Raw/
│   │   ├── GeoMx_Hs_CTA_v1.0.pkc
│   │   └── <DCC files ...>
│   ├── Clinical_annotations/OSIRESP_clinical_annotations.xlsx
│   └── Processed/                 # created by the processing scripts
├── OSIREAL_cohort/
│   ├── Raw/
│   ├── Clinical_annotations/OSIREAL_clinical_annotations.xlsx
│   └── Processed/
├── OSIRESP_RB_cohort/              # post-progression biopsies additionally collected and added to the OSIRESP cohort
│   ├── Raw/
│   └── Processed/
└── H12O_TMA004_cohort/
    ├── Raw/
    ├── Clinical_annotations/H12O_TMA004_clinical_annotations.xlsx
    └── Processed/
```

## 2. Cell line persistence data (GSE193258)

Publicly available bulk RNA-seq of EGFR-mutant NSCLC cell lines under
osimertinib exposure.

**Source.** GEO accession [GSE193258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193258).

Download the following file from GEO and place them in
`Data/Cell_line_models/Persistence/`:

- `GSE193258_raw_counts_GRCh38.p13_NCBI.tsv`
  
The per-cell-line `.rds` objects (`GSE193258_H1975.rds`, `GSE193258_HCC827.rds`,
`GSE193258_HCC2935.rds`, etc.) are produced by the scripts in
`Scripts/1_Data_processing/Cell_line_models/Persistence/`.

## 3. In-house cell line osimertinib-resistance data

RNA-seq of H1975 and HCC827 parental vs. osimertinib-resistant cell lines
generated for this study.

**Availability.** Data will be deposited in GEO upon publication. The
accession number will be added here at that time.

**Expected layout:**

```
Data/Cell_line_models/Osimertinib_resistance/
├── H1975_rawcounts.xlsx
├── HCC827_rawcounts.xlsx
├── H1975_normrlog_counts.rds       # created by the processing script
└── HCC827_normrlog_counts.rds      # created by the processing script
```

## 4. Gene set annotations

`Data/Bruker_CTA_assay_gene_sets.rds` — gene set definitions used for GSEA on
the CTA panel. Contact the corresponding author to obtain this file.
