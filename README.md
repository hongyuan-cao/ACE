# Factor-Model-FAME：Factor Analysis based Methodology for Expression(FAME) data associates novel miRNAs with Ovarian Cancer

## Overview

This repository contains the MATLAB implementation of simulations and real data analysis in our manuscript: Factor Analysis based Methodology for Expression(FAME) data associates novel miRNAs with Ovarian Cancer.

All the simulations are in the folder "factor_one sample" for one-sample test, "factor_two sample(Pairwise Comparision)" and "factor_two sample(Transition to One-Sample Test)" for two-sample test. And we suggest use "factor_two sample(Transition to One-Sample Test)" rather than "factor_two sample(Pairwise Comparision)", it has the higher power.

The real data analysis are in the "real data analysis" folder containing the code, data and results visualization.

## Work flow

Here are some descriptions about the MATLAB .m files in the simulation:

- Example_one_sample_test.m and Example_two_sample_test.m: give examples for one-sample and two-sample test to show how to use the code.
- m0.m: main function file to get all important variables we need, such as FDR, rejection number and power.
- solveW.m: estimate factor model.
- solvet_hat00.m: function file to get FAME statistics and collect critical value which is computed by t_hat.m.
- pai.m: function file to get the estimation of the proportion of non-nulls. This method is based on Cao, H., Kosorok, M. (2011). Simultaneous critical values for t-tests in very high dimensions. Bernoulli, 17, 347–394.
- t_hat.m: function file to compute critical value.
- generate_Y.m: function file to generate data followed factor model.

For real data analysis:

- real_data_analysis.m: gives an example to run the analysis and print rejection number.
- plot_reject.m: visualization for the rejection number under varying FDR control level, and compare with CK, Fan, Storey method.
- poorly_expressed_plot.m: visualization for the miRNAs which are poorly expressed.
- vocano_plot.m: vocano plot compared with BH method.

And in the file /real data analysis/data/:

- data.benchmark.csv: gene expression data. The letter at the end is V and E in the first line indicate they are high-grade serous ovarian and endometrioid endometrial cancer samples, respectively. They both have n = 96 samples and p = 3523 markers for each of samples. And these 3523 markers in our data belong to 1347 miRNAs.
- Ebench.csv: endometrioid endometrial cancer samples
- Vbench.csv: high-grade serous ovarian cancer samples
