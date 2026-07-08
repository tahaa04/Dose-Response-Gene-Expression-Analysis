# Dose-Response Gene Expression Analysis in Python

## Project Overview

This project analyzes dose-response gene expression data from the GEO dataset **GSE147379**, which examined transcriptional changes in **rat primary hepatocytes** after exposure to multiple chemicals. The dataset contains TempO-Seq RNA sequencing expression data along with sample-level metadata describing chemical exposure and dose.

The goal of this project is to:

* inspect raw GEO expression and metadata files,
* extract and link sample metadata with gene expression values,
* normalize raw expression counts using CPM,
* focus on a specific gene-chemical pair,
* and analyze whether gene expression shows a dose-response trend.

The final analysis focuses on the relationship between **S-bioallethrin dose** and **Acot1 expression**.

This project is best understood as a **data wrangling, exploratory data analysis, and statistical trend analysis project**. The final dose-response subset contains only six dose points, so the regression line is used only as a descriptive trend line.

## Biological Background

Chemical exposure can alter gene expression in a dose-dependent manner. In toxicogenomics, this idea is useful because changes in gene expression can provide early molecular evidence of how cells respond to chemical stress before visible disease outcomes appear.

The dataset used in this project, **GSE147379**, was designed around this idea. It measured gene expression changes in rat primary hepatocytes after exposure to 29 ToxCast chemicals across multiple doses. According to the GEO record, these chemicals were selected because in vitro transcriptional effects could be related to possible in vivo outcomes such as liver cancer and steatosis.

The broader biological motivation is therefore:

1. Chemical exposure can produce measurable transcriptional responses.
2. These responses may vary with dose.
3. Liver cells are especially relevant because the liver is a major site of chemical metabolism and toxicological response.
4. Genes involved in lipid metabolism can be useful markers for studying liver-related biological effects.

The gene selected for this project is **Acot1** (Acyl-CoA thioesterase 1), which is involved in fatty acyl-CoA metabolism and hepatic lipid regulation. Literature suggests that ACOT1 is connected to fatty acid metabolism.

Because the original dataset is based on rat primary hepatocytes and was motivated by chemical effects related to liver outcomes, Acot1 provides a biologically meaningful gene to examine in a small dose-response trend analysis.

## Dataset

**Dataset:** GSE147379
**Source:** NCBI Gene Expression Omnibus
**Organism:** Rattus norvegicus
**Experiment type:** Expression profiling by high-throughput sequencing
**Assay:** Rat TempO-Seq targeted RNA sequencing
**Cell type:** Rat primary hepatocytes
**Exposure format:** Chemical dose-response treatment

The raw data used in this project includes:

* a series matrix file containing sample metadata,
* an expression matrix file containing raw gene expression counts,
* processed files generated from the wrangling script.


## Project Workflow

### 1. Raw Data Inspection

The raw data EDA notebook inspects the original GEO files to understand their structure.

This step checks:

* whether the expression matrix loads correctly,
* whether the target gene `Acot1_7968` is present,
* what the sample columns look like,
* where chemical and dose metadata are stored,
* and whether metadata sample names match expression matrix sample columns.

This confirms that the expression data and sample metadata can be linked using sample identifiers.

### 2. Data Wrangling

The wrangling script extracts the gene and chemical subset needed for analysis.

The script:

* reads the GEO series matrix metadata,
* extracts sample name, chemical, and dose information,
* loads the raw expression count matrix,
* identifies the Acot1 expression row,
* links Acot1 counts with sample metadata,
* filters for S-bioallethrin samples,
* calculates total expression count per sample,
* normalizes Acot1 expression using CPM,
* and saves clean processed data files.

### 3. Trend Analysis

The trend analysis script uses the processed data to examine the relationship between S-bioallethrin dose and Acot1 expression.

Acot1 expression is analyzed using CPM (Counts Per Million), which normalizes the raw gene count relative to the total expression count in each sample. This makes expression values more comparable across samples.

The script calculates a small set of descriptive statistics, including Pearson correlation, Spearman correlation, and a simple linear trend line. The purpose is to summarize the observed dose-response pattern, not to build a predictive model.

The script also generates plots and saves the statistics, graphs, and written observations to the output/ folder.

## Repository Structure

```text
dose-response-gene-expression-analysis/
│
├── README.md
├── requirements.txt
│
├── data/
│   ├── raw_data/
│   │   ├── GSE147379_series_matrix.txt
│   │   └── GSE147379_Corton_29_Chem_data_matrix.txt
│   │
│   └── processed_data/
│       ├── acot1_sbioallethrin_wrangled.tsv
│       └── lr_data_sbio_acot.tsv
│
├── notebooks/
│   └── 01_raw_data_eda.ipynb
│
├── src/
│   ├── 01_data_wrangling.py
│   └── 02_trend_analysis.py
│
└── output/
    ├── summary_statistics.tsv
    ├── trend_observations.txt
    ├── dose_response_scatter.png
    ├── dose_response_line.png
    └── dose_response_regression_line.png
```

---

## Output Files

The analysis results are saved in the `output/` folder.

* `summary_statistics.tsv` — contains the calculated correlation and trend statistics.
* `trend_observations.txt` — contains a short written interpretation of the observed trend.
* `dose_response_scatter.png` — scatter plot of dose against normalized Acot1 expression.
* `dose_response_line.png` — line plot showing expression change across dose levels.
* `dose_response_regression_line.png` — scatter plot with a descriptive linear trend line.

The processed data files are saved in `data/processed_data/`.

* `acot1_sbioallethrin_wrangled.tsv` — full processed subset with sample metadata, raw counts, total expression count, normalized expression, and CPM.
* `lr_data_sbio_acot.tsv` — simplified file containing only `dose_uM` and `CPM` for trend analysis.

---

## Observed Trend

The processed data suggests that normalized **Acot1** expression generally increases as **S-bioallethrin** dose increases.

A fuller interpretation is provided in:

output/trend_observations.txt

The generated plots in the `output/` folder visually summarize this dose-response pattern.


## Limitations

This is a small exploratory analysis focused on one gene and one chemical. The final trend analysis uses only six dose points, so the regression line is descriptive and should not be interpreted as a predictive model or strong biological claim.

---

## How to Run

Install dependencies, then run the scripts in order from the project root.

```bash
pip install -r requirements.txt

python src/01_data_wrangling.py
python src/02_trend_analysis.py
```

The wrangling script reads the raw GEO files and saves the processed data to:

```text
data/processed_data/
```

The trend analysis script reads the processed data and saves the statistics, plots, and written observations to:

```text
output/
```

Optional raw data EDA notebook:

```bash
jupyter notebook notebooks/01_raw_data_eda.ipynb
```

## Skills Demonstrated

This project demonstrates:

* working with public GEO gene expression data,
* parsing metadata from a series matrix file,
* handling tab-separated expression matrices,
* linking sample metadata with expression data,
* normalizing sequencing count data,
* creating analysis-ready datasets,
* calculating correlation statistics,
* plotting dose-response trends,
* and writing interpretable biological observations from data.

---

## Literature and References

1. NCBI Gene Expression Omnibus. **GSE147379: Examination of gene expression changes after exposure to 29 chemicals in rat primary hepatocytes.**

2. Bushel, P. R. et al. **A comparison of the TempO-Seq S1500+ platform to RNA-Seq and microarray using rat liver mode of action samples.** *Frontiers in Genetics*, 2018.

3. Franklin, M. P. et al. **Acyl-CoA thioesterase 1 regulates PPARα to couple fatty acid flux with oxidative capacity during fasting.** *Diabetes*, 2017.
