#' MaxQuant proteomics data analysis
#'
#' Processes the MaxQuant proteomics dataset and performs exploratory statistical analysis for a single categorical variable. Accepts as input the ProteinGroups.txt file.
#'
#' @param ... The specific path to the folder where the samples from each group are located. They are passed as unnamed arguments via "...".  Attention: Ensure paths use '/' as a directory separator.
#' @param file The whole path to the MaxQuant ProteinGroups.txt file. The folders in the file path must be separated either with the forward slashes (/), or with the double backslashes (\\). See the example for inserting correctly the file path.
#' @param group_names A character vector specifying group names. The order of the names should align with the order of the sample groups in the input tsv file.
#' @param samples_per_group A numerical vector giving the number of samples in each group. The order of the numbers should align with the order of the names in group_names.
#' @param normalization The specific method for normalizing the data.By default it is set to FALSE. Options are FALSE for no normalization of the data, "log2" for a simple log2 transformation, "Quantile" for a quantiles based normalization  and "Cyclic_Loess" for a Cyclic Loess normalization of the log2 data, "median" for a median one, "TIC" for Total Ion Current normalization, "VSN" for Variance Stabilizing Normalization and "PPM" for Parts per Million transformation of the data.
#' @param filtering_value The maximum allowable percentage of missing values for a protein. Proteins with missing values exceeding this percentage will be excluded from the analysis. By default it is set to 50.
#' @param global_filtering TRUE/FALSE. If TRUE, the per-protein percentage of missing values will be calculated across the entire dataset. If FALSE, it will be calculated separately for each group, allowing proteins to remain in the analysis if they meet the criteria within any group. By default it is set to TRUE.
#' @param imputation Imputes all remaining missing values. Available methods: "LOD" for assigning the dataset's Limit Of Detection (lowest protein intensity identified), "LOD/2", "Gaussian_LOD" for selecting random values from the normal distribution around LOD with sd= 0.2*LOD, "zeros" for simply assigning 0 to MVs, mean" for replacing missing values with the mean of each protein across the entire dataset, "kNN" for a k-nearest neighbors imputation using 5 neighbors (from the package VIM) and "missRanger" for a random forest based imputation using predictive mean matching (from the package missRanger). By default it is set to FALSE (skips imputation).
#' @param independent TRUE/FALSE If TRUE, the samples come from different populations, if FALSE they come from the same population (Dependent samples). By default, it is set to TRUE. If set to FALSE, the numbers given in the samples_per_group param must be equal to each other.
#' @param parametric TRUE/FALSE. Specifies the statistical tests that will be taken into account for creating the PCA plots and heatmap. By default it is set to FALSE (non-parametric).
#' @param significance "p" or "BH" Specifies which of the p-values (nominal vs BH adjusted for multiple hypothesis) will be taken into account for creating the PCA plots and the heatmap. By default it is set to "p" (nominal p-value).
#' @param bugs Either 0 to treat Proteome Discoverer bugs as Zeros (0) or "average" to convert them into the average of the protein between the samples. By default, it is set to 0. Bugs are referred to to the proteins with empty values inside a single-file analysis
#'
#' @return Returns the complete output of the exploratory analysis: i) The processed, or filtered/normalized data ii) Statistical output containing results for the parametric (limma+ANOVA) and non-parametric tests (Wilcoxon+Kruskal-Wallis+PERMANOVA), along with statistical tests for heteroscedasticity, iii) Quality metrics for the input samples iv) QC plots and exploratory visualizations.
#' @importFrom openxlsx write.xlsx  read.xlsx
#' @importFrom grDevices  dev.off bmp colorRampPalette
#' @importFrom stringr str_trunc
#' @importFrom dplyr select  group_by  group_modify everything  %>%
#' @importFrom tidyr gather pivot_longer
#' @importFrom broom tidy
#' @importFrom reshape2 melt
#' @importFrom ggpubr ggarrange
#' @importFrom ggplot2 ggplot ggsave geom_smooth geom_violin scale_color_gradient element_line theme_linedraw scale_fill_manual scale_color_manual aes geom_histogram element_rect geom_point xlab ylab ggtitle theme_bw theme_minimal theme element_text guides guide_legend geom_boxplot labs theme_classic element_blank geom_jitter position_jitter
#' @importFrom VIM kNN
#' @importFrom stats kruskal.test p.adjust prcomp complete.cases sd wilcox.test friedman.test rnorm bartlett.test model.matrix heatmap median na.omit
#' @importFrom forcats fct_inorder
#' @importFrom limma topTable eBayes contrasts.fit normalizeCyclicLoess normalizeCyclicLoess lmFit normalizeQuantiles duplicateCorrelation
#' @importFrom pheatmap pheatmap
#' @importFrom grid gpar
#' @importFrom car leveneTest
#' @importFrom missRanger missRanger
#' @importFrom utils read.delim
#' @importFrom vegan adonis2
#'
#' @examples
#' #Example of running the function with paths for two groups.
#' # The file path is a placeholder, replace it with an actual file.
#'
#'
#'
#' proteinGroups.txt <- system.file("extdata", "proteinGroups.txt", package = "ProtE")
#' # Copy the file to a temporary directory for CRAN checks
#' temp_file <- file.path(tempdir(), "proteinGroups.txt")
#' file.copy(proteinGroups.txt, temp_file, overwrite = TRUE)
#'
#' maximum_quantum(file = temp_file,
#'        group_names = c("Healthy","Control"),
#'        samples_per_group = c(4,4), filtering_value = 80)
#'
#'
#' @export

ProtE_analyse <-function(...,
                         file = NULL,
                         group_names = NULL,
                         samples_per_group = NULL,
                         normalization = FALSE,
                         imputation = FALSE,
                         global_filtering = TRUE,
                         independent = TRUE,
                         filtering_value = 50,
                         parametric = FALSE,
                         significance = "p",
                         description = FALSE,  bugs = 0)
{








}
