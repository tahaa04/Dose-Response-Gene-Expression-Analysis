from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr, linregress


# File paths

input_file = "data/processed_data/lr_data_sbio_acot.tsv"


# Load processed data

data = pd.read_csv(input_file, sep="\t")

x = data["dose_uM"]
y = data["CPM"]



# Correlation and trend statistics

# Pearson correlation:
# Measures strength of a linear relationship between dose and expression.
pearson_corr, _ = pearsonr(x, y)

# Spearman correlation:
# Measures whether expression generally increases or decreases as dose increases.
# Useful for monotonic trends that may not be perfectly linear.
spearman_corr, _ = spearmanr(x, y)

# Simple linear regression:
# Used here only as a descriptive trend line, not as a predictive model.
regression_result = linregress(x, y)

slope = regression_result.slope
intercept = regression_result.intercept
r_squared = regression_result.rvalue ** 2


# --------------------------------------------------
# Save summary statistics
# --------------------------------------------------

summary_stats = pd.DataFrame({
    "statistic": [
        "Pearson correlation",
        "Spearman correlation",
        "Linear regression slope",
        "Linear regression intercept",
        "R-squared"
    ],
    "value": [
        pearson_corr,
        spearman_corr,
        slope,
        intercept,
        r_squared
    ],
    "note": [
        "Measures strength of linear relationship between dose and CPM.",
        "Measures whether CPM generally increases as dose increases.",
        "Estimated change in CPM for each 1 uM increase in dose.",
        "Estimated CPM value when dose is 0 uM.",
        "Coefficient of determination; describes how much variation is explained by the linear trend."
    ]
})

summary_stats.to_csv(
    "output/summary_statistics.tsv",
    sep="\t",
    index=False
)



# Plot 1: Dose vs CPM scatter plot


plt.figure(figsize=(5, 3))

plt.scatter(data["dose_uM"], data["CPM"])

plt.xlabel("S-bioallethrin dose (uM)")
plt.ylabel("Acot1 expression (CPM)")
plt.title("Dose-response trend of Acot1 expression")

plt.tight_layout()
plt.savefig("output/dose_response_scatter.png", dpi=150)
plt.close()



# Plot 2: Dose-response line plot


plt.figure(figsize=(5, 3))

plt.plot(
    data["dose_uM"],
    data["CPM"],
    marker="o"
)

plt.xlabel("S-bioallethrin dose (uM)")
plt.ylabel("Acot1 expression (CPM)")
plt.title("Acot1 expression across S-bioallethrin doses")

plt.tight_layout()
plt.savefig("output/dose_response_line.png", dpi=150)
plt.close()



# Plot 3: Scatter plot with regression trend line


predicted_y = regression_result.intercept + regression_result.slope * x

plt.figure(figsize=(5, 3))

plt.scatter(data["dose_uM"], data["CPM"], label="Observed data")
plt.plot(data["dose_uM"], predicted_y, label="Linear trend line")

plt.xlabel("S-bioallethrin dose (uM)")
plt.ylabel("Acot1 expression (CPM)")
plt.title("Linear trend line for dose-response relationship")
plt.legend()

plt.tight_layout()
plt.savefig("output/dose_response_regression_line.png", dpi=150)
plt.close()


# --------------------------------------------------
# Save trend observations
# --------------------------------------------------

observations = f"""
Dose-Response Trend Analysis Summary

Gene of interest: Acot1
Chemical of interest: S-bioallethrin
Expression measure: CPM (Counts Per Million)

Statistics calculated:

Pearson correlation:
Measures the strength of a linear relationship between chemical dose and gene expression.

Spearman correlation:
Measures whether gene expression generally increases or decreases as dose increases.

R-squared:
Coefficient of determination. It describes how much variation in gene expression is explained by the linear trend line.

Results:

Pearson correlation: {pearson_corr:.3f}
Spearman correlation: {spearman_corr:.3f}
Linear regression slope: {slope:.3f}
Linear regression intercept: {intercept:.3f}
R-squared: {r_squared:.3f}

Observed trend:

Normalized Acot1 expression generally increases as S-bioallethrin dose increases.

The Spearman correlation is useful here because it checks whether the overall trend is increasing across dose levels. The regression line is used only as a descriptive trend line, not as a predictive machine learning model.

The observed pattern suggests a positive dose-response trend between S-bioallethrin concentration and Acot1 expression.
"""

with open("output/trend_observations.txt", "w", encoding="utf-8") as f:
    f.write(observations)


print("Trend analysis complete.")
print(f"Statistics saved to: {'output/summary_statistics.tsv'}")
print(f"Observations saved to: {'output/trend_observations.txt'}")
print("Plots saved to output folder.")