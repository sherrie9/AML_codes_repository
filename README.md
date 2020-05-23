
# Biomarker discovery in acute myeloid leukemia patients

*Disclaimer: All patients' identifiers in the codes are removed. Metadata csv files and intermediate RData files and some plots are not uploaded due to sensitivity reason. All information included in this readme is publishable under collaborator's consent and they are also published in my graduate thesis. This repository only contains codes for running analysis.*

## Problem Statement:
Discover bone marrow regeneration patterns in acute myeloid lymphoma patients who relapsed after chemotherapy

## Data:
Collected from multiple medical centers in Europe. Data consists of 5 tubes with 40 biomarkers in total + AML cohorts(n=20) + normal cohorts(n=20). 

|       |Markers|
|-------|------------------------------------------- |
|Tube1  | CD56,CD13,CD34,CD117,CD33,CD11B,HLADR,CD45|
|Tube2  | CD36,CD64,CD34,CD117,CD33,CD14,HLADR,CD45|
|Tube3  | 8 more markers|
|Tube4  | 8 more markers|
|Tube5  | 8 more markers|

The original data did not contain the same number of samples per group, roughly 4-6 samples per group. I did random sampling so that each group and subgroups contains 20 samples. This is to increase the statistical power for subsequent analysis. The random sampling function (**randomResampling**) can be found in [helperFunc_S.R](https://github.com/sherrie9/AML_codes_repository/blob/master/helperFunc_S.R)

|      |Normal| Day22 |Last before 2nd ind|Last before consolidation|
|------|------|-----------------------------|-----------------------------|----------------------------|
|Tube1 |20    |20 relasped + 20 non-relapsed|20 relasped + 20 non-relapsed|20 relasped + 20 non-relapsed|
|Tube2 |ditto |ditto                        |ditto                        |ditto|
|Tube3 |ditto |ditto                        |ditto                        |ditto|
|Tube4 |ditto |ditto                        |ditto                        |ditto|
|Tube5 |ditto |ditto                        |ditto                        |ditto|

## Study Design:
1. Identify significantly different cell populations at each time points for relapsed, non-relapsed and normal cohorts.
2. Identify cell populations with different regeneration trends over time between relapsed and non-relapsed patients.

## Methods:

### Data preprocessing:
Preprocessing involves data cleaning and data integration. Data was collected from different medical centers. Although each center tried to follow the same naming conventions for marker names and dye names and time points names, there were still minor differences that prevent them from being properly integrated. In [preprocess.R](https://github.com/sherrie9/AML_codes_repository/blob/master/preprocess.R), I changed names so that they are all consistent. At this point, each tube has one single integrated dataset containing 140 samples (sum of each row in the above table). Tubes cannot be integrated together because they contain different markers. Different combinations of markers require different analysis strategy as we will see later. 

Flow cytometry data requires compensation and transformation. Detailed steps are in [gating_preFaust_pretSne.R](https://github.com/sherrie9/AML_codes_repository/blob/master/gating_preFaust_pretSne.R). Inside this code, I also removed technical outlier data using [flowCut](https://github.com/jmeskas/flowCut), and removed dead cells and doublets.

### Cell Population Identification

Supervised automated gating is to customize codes for gating strategy provided for each tube. In [code_automatedGating.R](https://github.com/sherrie9/AML_codes_repository/blob/master/code_automatedGating.R), I wrote the code that produced the supervised gating. In code [gating_preFaust_pretSne.R], I wrote steps for doing the unsupervised analysis. The clustering of CD45+SSC- cells are using k means clustering implemented by flowPeak package. However, the resulting grouping requires human supervision and proper regroup and organization. 

![gating_image](https://github.com/sherrie9/AML_codes_repository/blob/master/Plots/gating.PNG)

### Biomarker Discovery
1. In code [flowType_pvalue.R](https://github.com/sherrie9/AML_codes_repository/blob/master/flowType_pvalue.R), I did pairwise t-test of 3^8 = 6561 (8 markers) cell populations from each tube for three cohorts; selecting populations with adjusted p value <0.05; In code [flowType-RchyOptimyx.R](https://github.com/sherrie9/AML_codes_repository/blob/master/flowType-Rchyoptimyx.R), I did marker optimization. This method corresponds to study design 1.

2. In code [regeneration_dynamics.R](https://github.com/sherrie9/AML_codes_repository/blob/master/regeneration_dynamics.R), I computed p value for difference in linear regression fitted to time series data for relapsed and non-relapsed patients. This corresponds to study design 2. 

## Results + Findings




