# Time-Series Analysis of Meteorological Factors and Influenza Incidence in Kawasaki City, Japan

## Project Overview

This project evaluates and compares various time-series forecasting methods for the short-term prediction of daily influenza incidence in Kawasaki City, Japan. The study specifically investigates the impact of meteorological variables on influenza transmission and compares the performance of univariate models against models incorporating weather data.

## Team Members

- **Yongxi Feng** 
- **Jiayi Lyu** 
- **Keer Zheng** 
- **Gina Choi** 

## Research Questions

This project addresses the following key questions:

1. Which time-series forecasting methods provide the most accurate short-term forecasts of daily influenza incidence?
2. How do exponential smoothing methods (Simple Exponential Smoothing, Holt’s linear trend method, and Winters’ seasonal method) compare with ARIMA models in forecasting performance?
3. Does incorporating meteorological variables as explanatory predictors improve forecast accuracy compared with univariate time-series models?
4. Are there differences in forecastability between Total influenza cases, Influenza A, and Influenza B time series?

## Data Description

### Data Sources
- **Influenza Incidence Data**: Sourced from the Kawasaki City Infectious Disease Surveillance System. Includes mandatory daily reporting by all medical institutions using rapid diagnostic test kits.
- **Meteorological Data**: Obtained from a single monitoring station operated by the Japan Meteorological Agency in Kawasaki City.
- **Integrated Dataset**: The curated, analysis-ready dataset is available on [Zenodo](https://zenodo.org/records/15423347).

### Time Period and Frequency
- **Frequency**: Daily
- **Period**: March 2014 – April 2025 (approximately 11 years, ~4,079 observation days)

### Variables

#### Dependent (Outcome) Variables
The primary response variables are daily influenza incidence counts:
- Total daily influenza cases
- Daily Influenza A cases
- Daily Influenza B cases

Each series is analyzed as a separate univariate time series.

#### Independent (Explanatory) Variables
Daily meteorological variables considered as potential predictors:
- Mean temperature (°C)
- Minimum temperature (°C)
- Maximum temperature (°C)
- Relative humidity (%)
- Wind speed (m/s)
- Rainfall (mm)
- Sunshine duration (hours)
- Vapour pressure (hPa)

#### Control Variables
To account for temporal patterns:
- Lagged influenza counts (autoregressive terms)
- Seasonal components (e.g., annual seasonality)
- Trend components

## Methodology

The analytical approach is divided into five stages:

### Stage 1: Exploratory Analysis
- Visualization of time series to identify trends and seasonality.
- Interpretation of Sample Autocorrelation Function (ACF) and Partial Autocorrelation Function (PACF) to assess internal structure.

### Stage 2: Exponential Smoothing Models
Implementation of baseline forecasting rules:
- Simple Exponential Smoothing
- Holt’s Linear Exponential Smoothing
- Winters’ Seasonal Exponential Smoothing

### Stage 3: Regression with Explanatory Variables
- Development of regression models using lagged meteorological variables to forecast influenza cases.

### Stage 4: ARIMA and ARIMAX Modeling
- Fitting Autoregressive Integrated Moving Average (ARIMA) models using R (`arima` function).
- Extension to ARIMAX models to include meteorological covariates.
- Diagnostic checking via residual analysis (ACF/PACF of residuals, Ljung-Box tests).

### Stage 5: Model Comparison
- Comparative evaluation of all models based on forecast accuracy metrics (e.g., RMSE, MAE).
- Selection of the optimal model for each influenza category (Total, A, B).

## Expected Outcomes

1. **Methodological Application**: Direct application of Distributed Lag Non-linear Models (DLNM) and time-series techniques to quantify lagged and non-linear effects of weather on influenza.
2. **Public Health Insights**: Identification of key meteorological predictors for influenza in Kawasaki City, providing an evidence base for public health alerts.
3. **Reproducible Workflow**: Establishment of a robust, reproducible daily time-series analysis workflow that can serve as a template for future epidemiological studies.

## Project Structure

```text
.
├── data/
│   ├── raw/                # Raw data files (if not loading directly from Zenodo)
│   └── processed/          # Cleaned and merged dataset
├── src/
│   ├── exploratory_analysis.R  # Visualization and ACF/PACF plots
│   ├── exponential_smoothing.R # SES, Holt, Winters implementations
│   ├── regression_models.R     # Regression with weather variables
│   ├── arima_models.R          # ARIMA and ARIMAX fitting and diagnostics
│   └── comparison.R            # Model evaluation and comparison metrics
├── notebooks/              # Jupyter/RMarkdown notebooks for step-by-step analysis
├── results/                # Generated plots, tables, and forecast outputs
├── requirements.txt        # Python dependencies (if using Python)
├── packages.R              # R package dependencies
└── README.md
```

## Installation and Setup

### Prerequisites
- R (version 4.0 or higher) recommended, or Python 3.8+
- Access to the internet to download the dataset from Zenodo

### R Setup
Install required packages:
```R
install.packages(c("forecast", "tseries", "ggplot2", "dlnm", "tsibble", "feasts"))
```

### Python Setup (Optional)
If using Python:
```bash
pip install pandas numpy statsmodels scikit-learn matplotlib seaborn
```

## Usage

1. **Data Loading**:
   Download the dataset from the [Zenodo repository](https://zenodo.org/records/15423347) and place it in the `data/raw/` directory, or load it directly in the scripts using the provided URL.

2. **Run Analysis**:
   Execute the scripts in sequential order:
   ```bash
   Rscript src/exploratory_analysis.R
   Rscript src/exponential_smoothing.R
   Rscript src/regression_models.R
   Rscript src/arima_models.R
   Rscript src/comparison.R
   ```

3. **View Results**:
   Generated plots and performance tables will be saved in the `results/` directory.

## License

This project is for academic and research purposes. Data usage is subject to the terms of the original data providers (Kawasaki City and Japan Meteorological Agency).

## References

- Kawasaki City Infectious Disease Surveillance System.
- Japan Meteorological Agency.
- Zenodo Repository: [https://zenodo.org/records/15423347](https://zenodo.org/records/15423347)
