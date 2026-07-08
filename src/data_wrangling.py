import csv
import pandas as pd
import numpy as np

#IMPORT FILES
series_file = "data/raw_data/GSE147379_series_matrix.txt"
expression_file = "data/raw_data/GSE147379_Corton_29_Chem_data_matrix.txt"

#INDEX TO ROWS WITH DESIRED DATA (SAMPLE DESC, CHEMICAL, DOSE) IN SERIES MATRIX FILE
with open(series_file, "r",) as f:
    reader = csv.reader(f, delimiter="\t")
    rows = list(reader)

for row in rows:
    if len(row)>0 and row[0] == "!Sample_description":
        sample_desc = row[1:]
        break

for row in rows:
    if len(row)>0 and row[0] == "!Sample_characteristics_ch1" and row[1]== "chemical: triclosan":
        sample_chem_row = row[1:]
        break

for row in rows:
    if len(row)>0 and row[0] == "!Sample_characteristics_ch1" and row[1]== "dose: 33.3 uM":
        sample_dose_row = row[1:]
        break

#CLEAN VALUES INTO LISTS

chem=[]
for value in sample_chem_row:
    chem.append(value.replace("chemical:", "").strip())

dose=[]
for value in sample_dose_row:
    dose.append(value.replace("dose:", "").replace("uM", "").strip())

    #CONVERT DOSE VALUE TO NUMBERS
dose = pd.to_numeric(dose, errors="coerce")  

#MAKE DATAFRAME WITH DESIRED METADATA VALUES (SAMPLE DESC, CHEMICAL, DOSE)
metadata = pd.DataFrame({
    "sample_column": sample_desc,
    "chemical": chem,
    "dose_uM": dose
})


#OPEN EXPRESSION FILE 
# (NOTE: COLUMN NAMES AFTER FIRST (ID REF) REFER TO EACH SAMPLE- SAMPLE DESC - WILL USE THIS LATER TO LINK METADATA AND EXPRESSION VALUES)
expression = pd.read_csv(expression_file, sep="\t")
sample_columns = expression.columns[1:]


            #    TO CHECK WHETHER METADATA AND CORRESPONDING GENE EXPRESSION VALUES ARE PRESENT FOR ALL SAMPLES:
            #    matching_samples = set(expression_sample_columns).intersection(set(metadata_sample_columns))
            #    print("Expression samples:", len(expression_sample_columns))
            #    print("Metadata samples:", len(metadata_sample_columns))
            #   print("Matching samples:", len(matching_samples))


#INDEX INTO TARGET GENE ROW IN EXPRESSION FILE AND ASSIGN VALUES TO A NEW SINGLE ROW DATAFRAME - gene_row
target_gene = "Acot1_7968"
gene_row = expression[expression["ID_REF"] == target_gene]
#ASSIGN VALUES OF ROW ZERO, AND COLUMNS (EXCEPT FIRST) OF gene_row-DATAFRAME TO SERIES-gene_counts
gene_counts = gene_row.iloc[0, 1:]
#ASSIGN SERIES-gene_counts values to DATAFRAME-gene_counts_df 
#zero index, add headings, convert vals to numbers
gene_counts_df = gene_counts.reset_index()
gene_counts_df.columns = ["sample_column", "raw_count"]
gene_counts_df["raw_count"] = pd.to_numeric(
    gene_counts_df["raw_count"],
    errors="coerce"
)
#MERGE DATAFRAMES-gene_counts_df AND metadata USING VALUES OF COLUMN sample_column - SAMPLE DESC OF EACH SAMPLE (UNIQUE) 
acot1_data = gene_counts_df.merge(
    metadata,
    on="sample_column",
    how="left"
)
#MAKE NEW DATAFRAME-total_counts_df WITH SAMPLE DESC-(sample_column) AND CORRESPONDING SUM OF THEIR ROWS (EXPRESSION VALUES)
sample_columns = expression.columns[1:]
sample_columns = expression.columns[1:]
expression[sample_columns] = expression[sample_columns].apply(
    pd.to_numeric,
    errors="coerce"
)
total_counts = expression[sample_columns].sum(axis=0)
total_counts_df = total_counts.reset_index()
total_counts_df.columns = ["sample_column", "total_expression_count"]

#MERGE DATAFRAMES- total_counts_df AND acot1_data USING COLUMN "sample_column"-SAMPLE DESC
acot1_data = acot1_data.merge(
    total_counts_df,
    on="sample_column",
    how="left"
)

#FILTER ROWS ONLY FOR TARGET CHEMICAL, SORT BY DOSE VALUES
target_chemical = "S-bioallethrin"
sbio_acot1 = acot1_data[
    acot1_data["chemical"].str.lower() == target_chemical.lower()
].copy()
sbio_acot1 = sbio_acot1.sort_values("dose_uM")

#CALCULATE AND ADD NORMALIZED AND CPM EXPRESSION COLUMNS
sbio_acot1["normalized_expression"] = (
    sbio_acot1["raw_count"] / sbio_acot1["total_expression_count"]
)
sbio_acot1["CPM"] = sbio_acot1["normalized_expression"] * 1_000_000

#ZERO INDEX
sbio_acot1 = sbio_acot1.reset_index(drop=True)

#MAKE NEW DATAFRAMAE WITH ONLY DOSE AND CPM
lr_data=sbio_acot1[["dose_uM","CPM"]].copy()



#EXPORT WRANGELD FILES

sbio_acot1.to_csv(
    "data/processed_data/acot1_sbioallethrin_wrangled.tsv",
    sep="\t",
    index=False
)

lr_data.to_csv(
    "data/processed_data/lr_data_sbio_acot.tsv",
    sep="\t",
    index=False
)