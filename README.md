README — Housing–Stock Volatility (Iran), 2016–2023

Reproducible EViews Program (DCC–GARCH + Cheung–Ng)
1) Overview

This repository reproduces monthly volatility linkages between Iran’s housing and stock markets (Apr 2016–Oct 2023) using:
ARMA (mean filtering) → univariate GARCH(1,1) → DCC(1,1) for dynamic correlation → Cheung–Ng variance causality (m=12).
All steps are scripted for EViews 13+.

2) Requirements

EViews 13 or later (MGARCH DCC support required)

Windows OS

Excel data: mgarch_gregorian_unmodified.xlsx (sheet gregorian)

Columns: date_gregorian (YYYY-MM), stock_index_close, tehran_price_per_m2_rial

3) Files

reproducible_results_eviews.prg — main EViews program

mgarch_gregorian_unmodified.xlsx — input data (place next to the PRG)

4) How to Run

Put both files in the same folder.

In EViews: File → Open → Program…, open reproducible_results_eviews.prg.

Run (▶ or F9).

Outputs will be generated:

Workfile: MGARCH_IR_EViews.wf1 (contains the SP_ALL spool and series)

Spool: SP_ALL (descriptives, ADF, ARMA, diagnostics, GARCH, DCC, Cheung–Ng)

Excel: results_tables_from_eviews.xlsx (supplementary tables: ARMA selection, DCC summary, Cheung–Ng summary & lags)

Figures:

Figure1_conditional_vols.png — Conditional volatilities (GARCH)

Figure2_dcc_correlation.png — Dynamic correlation (DCC)

If the data path differs, edit in the PRG:
%datapath = @runpath + "\mgarch_gregorian_unmodified.xlsx"

5) Mapping to the Paper

Main Results tables (Tables 1–8) are available in the spool outputs (ADF, ARMA estimates, GARCH(1,1) params, DCC params, Cheung–Ng summary).

Appendix tables (Q1 style) are exported to results_tables_from_eviews.xlsx:

Appendix Table A1 — ARMA Order Selection (AIC)

Appendix Table A2 — Residual Diagnostics (Q, Q², ARCH–LM)

Appendix Table A3 — DCC Correlation Summary

Appendix Table A4 — Cheung–Ng Lag-wise Correlations

Appendix Table A5 — Phillips–Perron Unit-Root Tests

Figures:

Figure 2 — Conditional Volatilities (GARCH(1,1))

Figure 3 — DCC Dynamic Correlation: Stock vs Housing

6) Reproducibility Notes

Estimation uses Gaussian QML with standard stationarity constraints.

Minor numerical differences may arise from solver tolerances/rounding; p-values shown as 0 or 1 reflect formatting/rounding.

Ensure EViews 13+ for mgarch(dcc, …) support.

7) Troubleshooting

Import errors: verify sheet/column names exactly as specified.

Locale/decimal: check Windows/Excel regional settings if decimals misread.

DCC convergence: re-run the program; close other heavy processes if needed.
