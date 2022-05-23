# BAMBI

## Overview
**B**enchmarking **A**nalysis of **M**icro**B**iome **I**nference methods

This repository contains the scripts to perform a large-scale benchmarking
effort of differential abundance testing methods for microbiome data that is
described in our [preprint](https://doi.org/10.1101/2022.05.09.491139).

The project relies heavily on [SIMBA](https://github.com/zellerlab/SIMBA), an
R package to simulate microbiome data with various methods, evaluate how well
they mimic reality, apply differential abundance testing methods on these simulations,
and finally evaluate the results from these differential abundance testing
methods. This repository contains the code that we used to perform different
benchmarks through `SIMBA`, relying on the
[`batchtools`](https://mllg.github.io/batchtools/) package for automation
and parallelization of computing jobs.
For reproduction purposes, we additionally make use of the
[`renv`](https://rstudio.github.io/renv/index.html) package, which provides
reproducible environments for R.

## Setup for Reproduction of Results

Since `BAMBI` is not an R package, the setup is a bit more complicated.
Reproducibility is ensured by using `renv` and `batchtools` is used for
communication with a high-performance computing infrastructure. We ran our
analyses on a SLURM cluster and on a Grid Engine architecture, but `batchtools`
is very versatile in which architectures it supports, including local machines 
and manually-defined computing clusters. We therefore encourage
you to check out the `batchtools` documentation if you are using a different
architecture.

* You will need R version 4.0.0 installed.
* Clone the [SIMBA](https://github.com/zellerlab/SIMBA) package which contains
functions for simulating data and applying differential abundance tests
* Install the `renv` v0.12.5 package from
[CRAN](https://cloud.r-project.org/web/packages/renv/index.html) with
`install.packages('renv')` to bootstrap the project-local library
contained in `renv.lock`
    * Restore the environment via `renv`, thereby creating the same environment
    (in terms of packages used and their versions) that we used in our
    analyses (you will be prompted by R the first time you open the terminal)
* Configure a `batchtools` template in the
[cluster_config](https://github.com/zellerlab/BAMBI/-/tree/master/cluster_config)
subdirectory
    * Default configuration is for SLURM HPC; details for other HPC
    environments available
    [here](https://mllg.github.io/batchtools/articles/batchtools.html#setup)
* Adjust the `.Rprofile` file within this repository (it is currently
    configured for our computing setting)
    * The `RENV_PATHS_ROOT` parameter indicates where `renv` stores packages
    * The `simba.loc` parameter indicates where the `SIMBA` package is located
    * The `temp.loc` parameter indicates where temporary files should be stored
* Download the requisite datasets to the data subdirectory by executing
`Rscript ./src/download_data.R`

## Instructions on Use

The scripts work on a per-simulation basis. That means that one script is
usually executed per simulation file. The scripts to create all the simulations
included in our analyses are stored in the `create_simulations` subdirectory.
In these example code chunks, we will use the simulations from McMurdie and
Holmes (MMH) as an example.

### 1. Create the simulation files
The simulations will be saved to `simulations/sim_MMH.h5` in
a [hierarchical data format](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)

```bash
Rscript ./create_simulations/sim_MMH.R
```
or (on a SLURM cluster)

```bash
cd create_simulations
sbatch -J MMH create_sim.sh sim_MMH.R
```
Expected output:
```
Loading SIMBA

+ Start checking data
++ BMI is interpreted as body mass index:
  will be converted to underweight/normal/overweight/obese
++ Age will be converted to a factor via quartiles
++ Library_Size will be converted to a factor via quartiles
++ meta-variables:
	Group
++ have only a single value across all samples and will be removed
++ meta-variables:
	Sample_ID
++ have too many values across samples and will be removed
++ Removing repeated measurements for individuals
+ Finished checking data

+ Start checking the simulation location
+ Finished checking the simulation location

+ Start filtering the data
++ Performing prevalence and abundance filtering
++ After filtering, the dataset contains 844 features and 880 samples
+ Finished filtering the data

+ Start checking the simulation parameters
+ Finished checking the simulation parameters

+ Save original data in the h5 file
+ Finished saving original data in the h5 file

+ Starting data generation using the method: McMurdie&Holmes
++ Create simulation template from data
++ Finished creating simulation template
+ Finished data generation   

+ Starting to check parameters for test idx creation
++ Remove test idx, if present
++ Test idx do not exist yet
+ Start test idx creation
+ Finished test idx creation  
```
_estimated time to create this specific simulation on a desktop computer: 3 mins_

_Please note: The time to create different simulations can vary dramatically
between methods_

_Please note: for some simulations (for example `sim_negbin_cor`), multiple
cores are needed, since the underlying code will call SPIEC-EASI to estimate
the underlying correlation structure, which is prohibitively slow when not
parallelized with multiple cores_

Also, please check out the [SIMBA](https://github.com/zellerlab/SIMBA) vignette
for more information about how the data is stored within the `.h5` file produced
by this script.

### 2. Check reality: Compare real and simulated data

In a second step, the simulated data are compared to the real data in various
measures, for example sparsity or feature variance. Additionally, PCoA are
computed to visualize the difference between real and simulated samples.
To do so, several jobs are created with `batchtools` when you call this script
with the simulation ID as parameters:

```bash
Rscript ./src/reality_checks_cluster.R sim_Weiss
```
The script also submits the jobs to the cluster (since it is usually only a
handful of jobs - as many as there are different effect sizes combinations).
After all the jobs have been completed, you can combine the results by calling:
```bash
Rscript ./src/reality_checks_combination.R sim_Weiss
```
Since this step demands quite a bit of memory (usually more than is available on
a submit node, we also have -again- a SLURM-specific script to send this compute
job to the cluster)
```bash
cd ./src
sbatch -J Weiss_reality reality_checks.sh
```
The results will be found in the subdirectory './reality_checks/sim_Weiss'

```
├── reality_checks
│   ├── sim_MMH      <- Reality check results from sim_MMH.
|   | ├── auc_all.tsv                       <- AUCs from machine learning and PERMANOVA
|   | ├── auc_plot.pdf                      <- AUC plotted as heatmap
|   | ├── effect_size.pdf                   <- Effect size measures plotted as scatter plot
|   | ├── pco_plots.pdf                     <- PCoA plots for different effect sizes
|   | ├── sparsity_plot.pdf                 <- Sparsity plotted as heatmap
|   | ├── sparsity.tsv                      <- Sparsity measures as table
|   | ├── variance.tsv                      <- Variance and effect size measures for all features
|   | ├── variance_by_group.pdf             <- Variance plotted as heatmap
|   | ├── variance_by_group_and_type.pdf    <- Variance plotted as heatmap
|   | └── variance_scatter_plots.pdf        <- Variance plotted as scatter plot
|   |
│   └── other_sim      <- Other simulations will get their own folder
```
### 3. Run differential abundance tests on the simulated data

In our benchmark, we applied various differential abundance testing methods to
the simulations, using the test indices created in step 1, so that each test
is applied to exactly the same data. Since there are many different simulated
or implanted effect sizes, several simulated sample sizes, and various methods
that all need a different amount of time, we end up with many different jobs to
send to the cluster. Here, `batchtools` comes to the rescue.

We first create a `batchtools` registry, in which we initialize all the jobs
that we want to have computed. Then, we send the jobs to the cluster (this part
requires a bit of manual oversight to check if some jobs have failed or if
all have finished).

For these tasks, we have two R scripts. First, we prepare the registry:
```{bash}
Rscript ./src/prepare_registry.R sim_MMH
```
Then, we can run the jobs with:
```{bash}
Rscript ./src/run_tests.R sim_MMH
```
Alternatively, to run the tests of a single method (for example ANCOM,
because there might be many for this specific method), we can also type:
```{bash}
Rscript ./src/run_tests.R sim_MMH ANCOM
```

Since the different tests have different time requirements, the jobs are chunked
differently across tests. Faster tests will have less jobs.
Therefore, for the `sim_MMH` simulation, we will have to run these tests:
```{R}
ALDEx2          980
ANCOM           980
corncob         980
distinct        980
ZIBSeq          980
ZIBSeq-sqrt     980
ANCOMBC         140
DESeq2          140
edgeR           140
metagenomeSeq   140
metagenomeSeq2  140
ZINQ            140
KS              70
limma           70
lm              70
wilcoxon        42
```    

For the implantation simulations, we have more repeats and also more effect
sizes, therefore we also have to run more tests.

For the confounded simulations, the setup is again slightly different, since
the resampling bias creates an additional parameter to be varied. Therefore,
there are dedicated scripts to prepare the `batchtools` registries for
those simulations, but the application of tests can be achieved with the same
script (`run_tests.sh`).

### 4. Evaluate test performance

Lastly, we can evaluate the performance of the results from the differential
abundance testing methods. For this, the P-values which are returned by the
various tools are compared with the information about which features were
actually used to carry a differential effect and the number of true positive
discoveries (TP), false positives (FP), true negative (TN), and false negatives
(FN) are computed with P=0.05 as cutoff. From those measures, the precision (PR),
the recall (R), and the false discovery rate (FDR) are inferred. Additionally,
the AUROC for the separation between true and background features is computed
using the P-values as predictor.

This can be achieved by:
```{bash}
Rscript ./src/eval_tests.R sim_MMH
```
Here, you can again evaluate by single testing methods (in this case ANCOM):
```{bash}
Rscript ./src/eval_tests.R sim_MMH ANCOM
```
The evaluation metrics are then stored in the subdirectory `test_results` and
can be read into R again:
```R
library(tidyverse)
test.results.MMH <- read_tsv('./test_results/sim_MMH.tsv')
head(test.results.MMH)
# A tibble: 6 x 17
    rep auroc    TP    FP    TN    FN    PR     R   FDR group subset job.id problem test  norm
  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <chr> <chr>   <dbl> <chr>   <chr> <chr>
1     1 0.5       0     0   761    84     0     0     0 ab1_… subse…      1 sim_MMH limma TSS  
2     2 0.501     0     0   761    84     0     0     0 ab1_… subse…      1 sim_MMH limma TSS  
3     3 0.5       0     0   761    84     0     0     0 ab1_… subse…      1 sim_MMH limma TSS  
4     4 0.499     0     0   761    84     0     0     0 ab1_… subse…      1 sim_MMH limma TSS  
5     5 0.5       0     0   761    84     0     0     0 ab1_… subse…      1 sim_MMH limma TSS  
6     6 0.5       0     0   761    84     0     0     0 ab1_… subse…      1 sim_MMH limma TSS  
# … with 2 more variables: time.running <dbl>, mem.used <dbl>
```

## Reproduction from Zenodo

We additionally provide the raw files from our analyses in a
[Zenodo repository](https://doi.org/0.5281/zenodo.6572226), so that you can
download those files (either the simulation files as `.h5` files or the
results from the differential abundance benchmarking) and run additional test
on them or compare to the results we obtained.

## Contact and Feedback

If you have any question or comment or if you run into any issue please:
- create an
[issue in this repository](https://github.com/zellerlab/BAMBI/issues/new) or
- email [Morgan Essex](mailto:Morgan.Essex@mdc-berlin.de) or
[Jakob Wirbel](mailto:jakob.wirbel@embl.de).

## License

`BAMBI` is distributed under the
[GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html) license.

## Citation

If you use `BAMBI`, please cite us by

> Wirbel J, Essex M, Foslund, SK Zeller G _Evaluation of microbiome
association models under realistic and confounded conditions_
**bioRxiv** (2022) https://doi.org/10.1101/2022.05.09.491139
