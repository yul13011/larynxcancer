# Oral Cancer Analysis: Propensity Score Weighting & Survival Analysis

**Author:** Eva Li, 2019

## Overview

This project analyzes the oral cancer (OC) outcomes between male and female patients using propensity score weighting and survival analysis techniques. The analysis examines patient characteristics, treatment patterns, and their association with patient survival using a comprehensive cancer registry dataset.

[📊 Click here to view the full report](https://yul13011.github.io/larynxcancer/Cancer_OC.html)

[🔬Click here to read the published research article](https://onlinelibrary.wiley.com/doi/abs/10.1002/hed.25897)

## Analysis Summary

### Key Methods
- **Propensity Score Weighting:** Used to balance baseline characteristics between treatment/comparison groups
- **Survival Analysis:** Kaplan-Meier (KM) curves to visualize survival outcomes
- **Data Cleaning:** Rigorous preprocessing to handle missing data and ensure data quality

### Analysis Steps

1. **Data Import & Cleaning**
   - Removed observations with missing outcome variables (death follow-up time & vital status)
   - Filtered tumor size categories for balance in propensity score analysis
   - Created subsets for oral cancer (OC) patients
   - Focused analysis on patients with primary tongue cancer

2. **Exploratory Survival Analysis**
   - Generated Kaplan-Meier curves pre-weighting

3. **Propensity Score Weighting**
   - Applied propensity score weighting separately for OC and tongue cancer subsets
   - Conducted balance checking to verify treatment group comparability

## Files

- **Cancer_OC.html** - Full interactive report with R code, analysis steps, visualizations, and results

## View the Analysis

Open `Cancer_OC.html` in a web browser to see the complete analysis including:
- R code for all data processing steps
- Data visualizations and survival curves
- Statistical summaries and methodology details

## Technologies

- **R** - Data analysis and visualization
- **RMarkdown** - Report generation
- **ggplot2 / base R** - Graphics

## Data

This analysis uses cancer registry data with variables including:
- Patient demographics and vital status
- Tumor characteristics (size, stage, primary site)
- Treatment information
- Follow-up time and outcome data

---

For questions or additional information about this analysis, please refer to the detailed code and documentation in `Cancer_OC.html`.
