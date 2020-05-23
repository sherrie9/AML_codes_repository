
# Biomarker discovery in acute myeloid leukemia patients

*Disclaimer: All patients' identifiers in the codes are removed. Metadata csv files and intermediate RData files and some plots are not uploaded due to sensitivity reason. All information included in this readme is publishable under collaborator's consent and they are also published in my graduate thesis. This repository only contains codes for running analysis *

##Problem Statement:
Discover bone marrow regeneration patterns in acute myeloid lymphoma patients who relapsed after chemotherapy

##Data:
Collected from multiple medical centers in Europe. Data consists of 5 tubes with 40 biomarkers in total + AML cohorts(n=20) + normal cohorts(n=20). 

|       |Markers|
|-------|------------------------------------------- |
|Tube1  | CD56,CD13,CD34,CD117,CD33,CD11B,HLADR,CD45|
|Tube2  | CD36,CD64,CD34,CD117,CD33,CD14,HLADR,CD45|
|Tube3  | 8 more markers|
|Tube4  | 8 more markers|
|Tube5  | 8 more markers|

The original data did not contain the same number of samples per group, roughly 4-6 samples per group. I did random sampling so that each group and subgroups contains 20 samples. This is to increase the statistical power for subsequent analysis. 

|      |Normal| Day22 |Last before 2nd ind|Last before consolidation|
|------|------|-----------------------------|-----------------------------|----------------------------|
|Tube1 |20    |20 relasped + 20 non-relapsed|20 relasped + 20 non-relapsed|20 relasped + 20 non-relapsed|
|Tube2 |ditto |ditto                        |ditto                        |ditto|
|Tube3 |ditto |ditto                        |ditto                        |ditto|
|Tube4 |ditto |ditto                        |ditto                        |ditto|
|Tube5 |ditto |ditto                        |ditto                        |ditto|

