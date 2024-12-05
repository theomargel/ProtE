
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/ProtE.png" width="60%" />

One function to analyze them all! The Proteomics Eye (ProtE) establishes
an intuitive framework for the univariate analysis of label-free
proteomics data. By compiling all necessary data wrangling and
processing steps into the same function, ProtE automates all pairwise
statistical comparisons for a given categorical variable, returning to
the user performance quality metrics, measures to control for Type-I or
Type-II errors, and publication-ready visualizations.

ProtE is currently compatible with data generated by MaxQuant, DIA-NN
and Proteome Discoverer.

## How to install:

Step 1: Install R

To get started with R, first download and install the latest version of
R from the official CRAN website:

- Go to the [R download page](https://cran.r-project.org/).
- Select the appropriate version for your operating system:
  - **Windows**: Click on “Download R for Windows”.
  - **MacOS**: Click on “Download R for macOS”.
  - **Linux**: Follow the instructions based on your distribution.

Once downloaded, run the installer and follow the instructions to
complete the installation.

Step 2: Install RStudio

Next, you’ll need an Integrated Development Environment (IDE) to work
with R. The most popular IDE is RStudio.

- Go to the [RStudio download
  page](https://www.rstudio.com/products/rstudio/download/).
- Select “RStudio Desktop” and download the installer for your operating
  system.

Run the installer and follow the on-screen instructions to install
RStudio.

Step 3: Install RTools (For Windows Users)

RTools is necessary if you need to compile packages from source on
Windows, which is common when installing certain R packages.

- Go to the [RTools download
  page](https://cran.r-project.org/bin/windows/Rtools/).
- Download the version of RTools that matches your R version.
- Run the installer and follow the instructions.

Make sure to select the option that allows RTools to be added to your
system path during installation.

Step 4: Download the Package

Install the released ProtE version from CRAN:

    install.packages("ProtE")

Install the development version of ProtE from GitHub:

    install.packages("pak")
    pak::pak("theomargel/ProtE")

Then load its library with:


    library(ProtE)

## Function inputs

ProtE features 4 functions, each one tailored for a specific use case.

1.  `maximum_quantum()` accepts as input the MaxQuant generated file
    ProteinGroups.txt
2.  `dianno()` accepts as input either of the two DIA-NN (or the
    FragPipe - DIANN) output files pg_matrix.tsv or
    unique_genes_matrix.tsv
3.  `pd_single()` accepts as input the Proteome Discoverer output file
    that contains all sample protein intensities/abundances in one table
4.  `pd_multi()` accepts as input separate Proteome Discoverer protein
    intensity files

## How to use functions `maximum_quantum()`,`dianno()`,`pd_multi()`

All 3 functions expect the input file to be parsed in the parameter
`file`. To enable statistical analysis, in the input file, samples
(columns) belonging to the same group must be sorted next to each other.
For example, samples from an experiment with a 3-groups categorical
variable (control, treatment, compound) could be arranged such that:
first columns = Control samples, middle columns = Treatment samples,
last columns = Compound samples.

## Setting up the input file path:

Assuming a MaxQuant quantification has been performed, the file
ProteinGroups.txt can be fed to ProtE with the function
`maximum_quantum`.

Insert the file path of the ProteinGroups.txt in the `file` parameter.
To copy-paste the file path in Windows, firstly locate the desired file
inside your folders. Hold Shift and right-click the file, then select
“Copy as Path” from the context menu. Go to RStudio and click Ctrl+V or
right-click to paste the path.

Because usually the directories will be separated with a single
backlash, ensure to use forward slashes (/) for specifying paths or
adding a second backlash e.g:

    maximum_quantum(file = "C:\\Bioprojects\\BreastCancer\\Proteomics\\MaxQuant\\ProteinGroups.txt")

or

    maximum_quantum(file = "C:/Bioprojects/BreastCancer/Proteomics/MaxQuant/ProteinGroups.txt")

## Setting up `group_names` and number of `samples_per_group`

Group names are defined in the parameter `group_names` as a vector. The
order of the group names inside the vector must follow the order of the
groups by which the samples (columns) have been arranged in the input
proteomics file (from the left to the right). Same goes for the number
of samples of each group, which is defined again as a vector in the
parameter `samples_per_group`. In the following example there are 3
groups (Control,Treatment,Compound) with the Control group consisting of
10 samples the Treatment group of 12 samples and the Compound with 9:

``` r

maximum_quantum(
                    file = "C:\\Bioprojects\\BreastCancer\\Proteomics\\MaxQuant\\ProteinGroups.txt",
                    group_names = c("Control", "Treatment", "Compound"),
                    samples_per_group = c(10, 12, 9),
                    imputation = FALSE,
                    global_filtering = TRUE,
                    sample_relationship = "Independent",
                    filtering_value = 50,
                    normalization = FALSE,
                    parametric= FALSE,
                    significance = "p")
```

In the pairwise comparisons, nominators and denominators of the
FoldChange (and consequently the sign of Log2FoldChage) are defined
based on the order of the group names declared in the parameter
`group_names`. The general notion based on which FoldChange is
determined is: NextGroup/PreviousGroup. In our example the FoldChange
for every pairwise comparison will be set as: Treatment/Control,
Compound/Control and Compound/Treatment.

## How to use `pd_multi()`

`pd_multi` is tailored for the analysis of multiple **Proteome
Discoverer** (PD) exports, each one corresponding to a single sample. To
be able to use it, the user must save the PD exports to different
folders corresponding to the different groups of the variable that is
going to be analyzed. The paths to these folders are specified in the
parameter …:

    pd_multi(excel_file = "C:\\Bioprojects\\BreastCancer\\Proteomics\\PD\\Control",
                           "C:\\Bioprojects\\BreastCancer\\Proteomics\\PD\\Treatment",
                           "C:\\Bioprojects\\BreastCancer\\Proteomics\\PD\\Compound",
                        imputation = FALSE,
                        global_filtering = TRUE,
                        sample_relationship = "Independent",
                        filtering_value = 50,
                        normalization = FALSE,
                        parametric= FALSE,
                        significance = "p")

In the pairwise comparisons, nominators and denominators of the
FoldChange (and consequently the sign of Log2FoldChage) are defined
based on the order of the declared group folders in the `pd_multi`
function. Again, the general notion based on which FoldChange is
determined, is: NextGroup/PreviousGroup. In our imaginary example the
FoldChange for every pairwise comparison will be set as:
Treatment/Control, Compound/Control and Compound/Treatment.

## Summary of the ProtE pipeline

All 4 functions streamline the following process, which is reported in
more details in the `ProtE Guide` vignette.

1.  Normalization of proteomic intensity values (excluding DIA-NN files)

2.  Filtering based on the percentage of missing values.

3.  Imputation of missing data to ensure robust downstream analysis.

4.  Description fetching for DIA-NN input files: .pg_matrix.tsv

Once the data processing is complete, the package performs statistical
analysis for every pairwise comparison to identify significant protein
abundance differences between experimental groups. The results are
automatically exported as Excel files, and a range of visualizations is
generated to facilitate QC and interpretation. These include:

• Principal Component Analysis (PCA) plots for dimensionality reduction
and group comparison.

• Heatmap highlighting significant proteins.

• Protein rank-abundance and meanrank-sd scatterplots.

• Boxplots and violin plots to display data distribution and variability
across groups.

## Output Directory

The results from each function are saved in a folder named
**ProtE_Analysis**, which is created inside the last directory of the
provided file(s).

ProtE creates 3 sub-folders: • Data_processing, with the files of the
resulting data processing. • Statistical_Analysis, with the results of
the statistical tests. • Plots, with all the plots saved in pdf. format.

## How to cite ProtE

ProtE was compiled out of R scripts utilized for routine processing,
analysis and visualization of label-free MS data. These scripts were
built and maintained over the course of many years in the Proteomics lab
of Dr Antonia Vlahou, in the BRFAA institute (Greece). To continue
enriching and expanding the functionality of ProtE we kindly request
feedback from the users via..

In order to cite ProtE the following reference can be used:
