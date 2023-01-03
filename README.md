# Systematic assessment of prognostic molecular features across cancers
Balaji Santhanam, Panos Oikonomou, Saeed Tavazoie

(please send queries/comments to st2744 at columbia dot edu)

Go to section

* [Overview](#overview)
* [Data](#data)
* [Install R packages](#install-r-packages)
* [Calculate MPS of a module](#calculate-mps-of-a-module)
    *   [Existing module](#existing-module)
    *   [New user defined module](#new-user-defined-module)
* [Prognostic potential of single locus observations](#prognostic-potential-of-single-locus-observations)
* [Build random survival forests using PCMs](#build-random-survival-forests-using-pcms)
* [Run time](#run-time)

## Overview 

In this work, we present a computational framework for the unbiased discovery of perturbations in single-genes and in functionally coherent gene modules that are predictive of patient survival and to facilitate systematic comparisons with conventional clinically used features in cancer. We have quantified the prognostic potential of mutations, copy-number aberrations and expression changes at individual loci. Since cancer-drivers and other prominent genes including regulators function in coordination with their target gene modules, we used the dysregulation of these modules as a proxy for the perturbed activities of the associated drivers. A module is defined to be a collection of genes unified through an interpretable biological principle (e.g., VEGF signaling pathway or targets of let-7 microRNA). We devised a statistical measure to quantify the magnitude and direction of module perturbations in the transcriptomes of cancer patients. This measure, termed "Module Perturbation Score (MPS)" is indicative of the module's perturbed activity in an individual transcriptome within a cancer cohort. We have developed a computational framework that enables the discovery of prognostic cancer modules (PCMs) based on the ability of module perturbation scores to define significantly divergent patient survival trajectories.

## Data

A few notes on the data used:
- Transcriptome data for TCGA were obtained from Broad Institute's FIREBROWSE (doi:10.7908/C11G0KM9), from ICGC (https://dcc.icgc.org/releases/release_28) and ACRG (GSE62254).
- For all cancers, we utilize both overall survival and progression-free interval survival as outcome endpoints when available. For TCGA, these data were sourced from supplementary material in Liu et al., Cell 2018 (PMID: 29625055). For AML, PFS is not defined and thus ignored. For ICGC and ACRG datasets, clinical data were dowloaded from the same source as the transcriptome data.
- Unfortunately, all of the raw and/or preprocessed data are too large to be hosted online right now. We are working to address this. However, we have provided all of the module definitions used here so the analyses can be recreated. Additionally, we have also hosted the results of these analyses for the majority of the PCMs we identified at https://tavazoielab.c2b2.columbia.edu/PCMs/landing_page.html.


All analyses were performed using R (version 3.6.3). In order to explore our code, please install R from https://cran.r-project.org/. There are few libraries that may need to be installed, in case they are not already.

## Install R packages
Download the compressed file, uncompress and set that as your working directory.

```{r}
install.packages('vroom')
install.packages("foreach", repos="http://R-Forge.R-project.org")
install.packages("doMC", repos="http://R-Forge.R-project.org")
install.packages('ggplot2')
install.packages('grid')
install.packages('gridExtra')
install.packages('cowplot', dependencies=T)
install.packages('survminer', dependencies=T)

install.packages("caret", dependencies = T)
install.packages("randomForestSRC", dependencies = T)
```
Depending on what packages are already installed (and their dependencies), please set aside some time to install these R packages. Uncompress the compressed file and set that as the working directory.


### Calculate MPS of a module 

For this portion of the analysis, the file `pcm.R` holds all the function definitions, so please 'source' it. Alternatively, you could also copy paste the function definitions into your R session.

```{r}
source('pcm.R') ### ensure path is OK.
```


We have provided example SNV, CNA and gene expression data sets (pancreatic cancer from TCGA), an example module, corresponding clinical data, MPS computed for this module in pancreatic cancer patients and inferred cell-type compositions for M2 Macrophages and CD8+ T-cells from CIBERSORT (https://api.gdc.cancer.gov/data/b3df502e-3594-46ef-9f94-d041a20a0b9a). 

The code provides two functions for generating relevant plots. 
    Mode 1: plots for exisiting modules (in this case, 'MSigDBONC_52')
    Mode 2: plots for new module (user-defined) that needs to be a tab-delimited file with gene symbol in first column.


#### **Existing module**
```{r}
get_MPS_existingModule('pancreatic', 'MSigDBONC_52')
```

For a new user-defined module, we have provided an example tab-delimited file. For this new module, we compute the MPS score and then go on to plot the clinical characteristics associated with the module's perturbation scores. 


#### **New user defined module**
```{r}
get_MPS_newModule('sample_input_module', 'pancreatic', module_name='user_defined')
```

### Prognostic potential of single locus observations
For analyzing the prognostic abilities of single-gene measurements (SNVs, CNAs and expression) the function `getSurv_locus()` can be used. For example,
```{r}
getSurv_locus('p53', 'pancreatic')
```

All intermediate results and plots are contained in the `results` folder. The typical run time for a module in one cohort is less than 30 seconds.

### Build random survival forests using PCMs
For this portion of the analysis, the file `rsf.R` contains the code for running the analysis. Again, the data are too large to be hosted here. We have provided data to run this portion of the analysis and plot the results. We have included data for pancreatic cancer patients including module information for PCMs discovered in this cohort (for overall and progression-free interval survival), prominent mutations and copy-number aberration categories, inferred immune cell type frequencies (data source: https://api.gdc.cancer.gov/data/b3df502e-3594-46ef-9f94-d041a20a0b9a) and MPS computed for the PCMs discovered.

This code provides a way to build random survival forest models (https://cran.r-project.org/web/packages/randomForestSRC/randomForestSRC.pdf) trained on:-
- MPS of PCMs discovered in pancreatic cancer
- standard clinical factors (histological type, pathological stage, age groups)
- both MPS and standard clinical factors combined

These models are 10-fold cross-validated (variable name `folds`) and 10 such model instances (variable name `repeat_iterations`) are built. Survival comparisons are made by ordering patients based on predicted risk (model output `predicted`) and stratifying them into high and low risk categories for each of the three model categories listed above. Finally, comparisons of the hazard ratio and concordance index are plotted. All intermediate files and plots are located in the directory `results_rsf`.

#### Run time
Run times for the examples described here for the MPS quantification is about 30 seconds and for the random survival forest models is ~20 minutes. Times are obtained from tests on a 1.6 GHz Intel Core i5 16 GB 2133 MHz Macbook Air machine running macOS Mojave version 10.14.6.
