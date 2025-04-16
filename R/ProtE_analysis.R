#' Proteomics data processing, DE analysis, GSEA and visualization.
#'
#' Takes as input proteomics Data from DIA-NN, MaxQuant and ProteomeDiscoverer and it performs exploratory data analysis, providing different options for  data manipulation (normalization, filtering based on the missing values and imputation). It then proceeds to perform statistical analysis and gene set enrichment analysis, while creating exploratory plots such as relative log expression boxplots and violin plots, heatmaps, PCA and Volcano plots.
#'
#' @param pd_single_dir For Proteome Discoverer (PD) multiple-files' (format of a single file per sample) proteomic analysis. The specific path to the folder(s) where the Proteome Discoverer's samples from each group are located.  Attention: Ensure paths use '/' or '\\' as a directory separator.
#' @param file The whole path to the proteomics quantification data i.e. for software MaxQuant the ProteinGroups .txt or .xlsx file, for ProteomeDiscoverer the .pdResult derived .xlsx file and for DIA-NN the pg_matrix or the unique_genes_matrix .tsv or .xlsx files.  Attention: Ensure paths use '/' or '\\' as a directory separator.
#' @param metadata_file Not requisite if samples_per_group and group_names are defined. An excel or txt table that includes in column 1 the Samples' names in column 2 the experimental group of comparison and in the rest columns, covariates (e.g. age, sex, batch) that affect the limma statistical procedures.
#' @param group_names Not requisite if metadata_file is inserted. A character vector specifying group names. The order of the names should align with the order of the sample groups in the input tsv file.
#' @param samples_per_group Not requisite if metadata_file is inserted. A numerical vector giving the number of samples in each group. The order of the numbers should align with the order of the names in group_names.
#' @param normalization The specific method for normalizing the data.By default it is set to FALSE. Options are FALSE for no normalization of the data, "log2" for a simple log2 transformation, "Quantile" for a quantiles based normalization  and "Cyclic_Loess" for a Cyclic Loess normalization of the log2 data, "median" for a median one, "TIC" for Total Ion Current normalization, "VSN" for Variance Stabilizing Normalization and "PPM" for Parts per Million transformation of the data.
#' @param filtering_value The maximum allowable percentage of missing values for a protein. Proteins with missing values exceeding this percentage will be excluded from the analysis. By default it is set to 50.
#' @param global_filtering TRUE/FALSE. If TRUE, the per-protein percentage of missing values will be calculated across the entire dataset. If FALSE, it will be calculated separately for each group, allowing proteins to remain in the analysis if they meet the criteria within any group. By default it is set to TRUE.
#' @param imputation Imputes all remaining missing values. Available methods: "LOD" for assigning the dataset's Limit Of Detection (lowest protein intensity identified), "LOD/2", "Gaussian_LOD" for selecting random values from the normal distribution around LOD with sd= 0.2*LOD,  "Gaussian_mean_sd" for sampling values from the normal distribution of the mean value of each sample with its standard deviation, "zeros" for simply assigning 0 to MVs, mean" for replacing missing values with the mean of each protein across the entire dataset, "kNN" for a k-nearest neighbors imputation using 5 neighbors (from the package VIM) and "missRanger" for a random forest based imputation using predictive mean matching (from the package missRanger). By default it is set to FALSE (skips imputation).
#' @param independent TRUE/FALSE If TRUE, the samples come from different populations, if FALSE they come from the same population (Dependent samples). By default, it is set to TRUE. If set to FALSE, the numbers given in the samples_per_group param must be equal to each other.
#' @param parametric TRUE/FALSE. Specifies the statistical tests that will be taken into account for creating the PCA plots and heatmap. By default it is set to FALSE (non-parametric).
#' @param significance "p" or "p.adj" Specifies which of the p-values (nominal vs p.adjusted for multiple hypothesis) will be taken into account for creating the PCA plots, the heatmap and the Volcano plots. By default it is set to "p" (nominal p-value).
#' @param species Species name, such as "Homo sapiens" or "Mus musculus", for the genes symbols output of the msigdbr package that will be used for the GSEA. Use msigdbr_species() for available options.
#' @param subcollection Sub-collection abbreviation, (such as "CP:REACTOME" for Reactome pathways, "H" for Hallmark and "GO:BP" for GO Biological Processes) that will be used for the GSEA. Use msigdbr_collections() for the available options.
#' @param LFC The LogFC threshold (absolute value) that you want to be the cutoff for the significant up or downregulation of the volcano plots
#' @param p.adjust.method The method used to adjust the p-values for multiple testing. Options are "holm","hochberg","hommel","bonferroni","BH","BY","fdr". Default is "BH".
#'
#' @return Returns the complete output of the exploratory analysis: i) The processed, or filtered/normalized data ii) Statistical output containing results for the parametric (limma+ANOVA) and non-parametric tests (Wilcoxon+Kruskal-Wallis+PERMANOVA), along with statistical tests for heteroscedasticity, iii) Quality metrics for the input samples iv) QC plots and exploratory visualizations.
#' @importFrom openxlsx write.xlsx  read.xlsx createWorkbook addWorksheet writeData
#' @importFrom grDevices  dev.off bmp colorRampPalette hcl.colors
#' @importFrom dplyr select  group_by  group_modify everything  %>% arrange filter slice mutate
#' @importFrom tidyr gather pivot_longer
#' @importFrom broom tidy
#' @importFrom reshape2 melt
#' @importFrom ggpubr ggarrange
#' @importFrom queryup query_uniprot
#' @importFrom ggplot2 ggplot ggsave geom_smooth coord_cartesian ylim xlim geom_segment geom_vline scale_y_discrete scale_x_continuous scale_size_continuous scale_colour_manual geom_violin scale_color_gradient scale_x_discrete element_line theme_linedraw scale_fill_manual scale_color_manual aes geom_histogram element_rect geom_point xlab ylab ggtitle theme_bw theme_minimal theme element_text guides guide_legend geom_boxplot labs theme_classic element_blank geom_jitter position_jitter
#' @importFrom VIM kNN
#' @importFrom stats kruskal.test p.adjust prcomp complete.cases as.formula sd wilcox.test friedman.test rnorm bartlett.test model.matrix heatmap median na.omit
#' @importFrom forcats fct_inorder
#' @importFrom limma topTable eBayes contrasts.fit normalizeCyclicLoess normalizeCyclicLoess lmFit normalizeQuantiles duplicateCorrelation
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_block
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @importFrom car leveneTest
#' @importFrom missRanger missRanger
#' @importFrom utils read.delim install.packages update.packages
#' @importFrom vegan adonis2
#' @importFrom msigdbr msigdbr
#' @importFrom fgsea fgsea
#'
#' @examples
#' #Example of running the function with paths for two groups.
#' # The file path is a placeholder, replace it with an actual file.
#'
#'
#'
#' proteinGroups.txt <- system.file("extdata", "proteinGroups.txt", package = "ProtE")
#' # Copy the file to a temporary directory for CRAN checks
#' temp_file.txt <- file.path(tempdir(), "proteinGroups.txt")
#' file.copy(proteinGroups.txt, temp_file.txt, overwrite = TRUE)
#'
#' ProtE_analyse(file = temp_file.txt,
#'        group_names = c("Healthy","Control"),
#'        samples_per_group = c(4,4), filtering_value = 80)
#'
#'
#' @export

ProtE_analyse <-function(file = NULL,
                         pd_single_dir= NULL,
                         group_names = NULL,
                         samples_per_group = NULL,
                         normalization = FALSE,
                         imputation = FALSE,
                         global_filtering = TRUE,
                         independent = TRUE,
                         filtering_value = 50,
                         parametric = FALSE,
                         significance = "p",
                         metadata_file = NULL,
                         species = "Homo sapiens",
                         p.adjust.method = "BH",
                         subcollection = "CP:REACTOME",
                         LFC = 1)
{
 if (!requireNamespace("msigdbdf", quietly = TRUE)) {
  install.packages("msigdbdf", repos = "https://igordot.r-universe.dev")
 } else {update.packages("msigdbdf",repos = "https://igordot.r-universe.dev", ask = FALSE)}
 update.packages("msigdbr", ask = FALSE)

  Sample=group1=  Accession =Description =Symbol =X =p.value= Mean = SD=bartlett_result= size =Y =df4_wide= percentage=variable =.= g1.name =g2.name=key =value = Gene.Symbol = NES= Regulation = padj = pathway = NULL
  uqg = FALSE
  print("The ProtE process starts now!")
  if (!is.logical(global_filtering)) {
    stop("Argument 'global_filtering' must be a logical value (TRUE or FALSE).")
  }
  if (!is.logical(independent)) {
    stop("Argument 'independent' must be a logical value (TRUE or FALSE).")
  }
  if (!is.logical(parametric)) {
    stop("Argument 'parametric' must be a logical value (TRUE or FALSE).")
  }
  if (!significance %in% c("p","p.adj")) {
    stop("Argument 'significance' must be either 'p' or 'p.adj'.")
  }
  valid_species <- msigdbr::msigdbr_species()$species_name
  if (!species %in% valid_species) {
    stop("Argument 'species' must match a valid species from msigdbr::msigdbr_species(). Examples: 'Homo sapiens', 'Mus musculus'.")
  }
  if (!is.character(subcollection)) {
    stop("Argument 'subcollection' must be a character string (e.g., 'H', 'CP:REACTOME').")
  }
  if (!is.numeric(LFC)) {
    stop("Argument 'LFC' must be a numeric value.")
  }
  if (LFC < 0) {
    stop("Argument 'LFC' must be a non-negative number (absolute logFC threshold).")
  }
  if (!p.adjust.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {
    stop("Argument 'p.adjust.method' must be one of: 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'.")
  }
  if (!normalization %in% c(FALSE,"median","Quantile","log2","VSN", "Total_Ion_Current","Cyclic_Loess", "PPM") ){
    stop("An incompatible normalization method was provided, please check the available and run the function again.") }

  if (!imputation %in% c("kNN","missRanger","mean","zeros","LOD","Gaussian_LOD","Gaussian_mean_sd", FALSE))  {
    stop("incompatible imputation method was selected, please check the availables and run again.")}




  group_paths <- as.list(pd_single_dir)
  if (!is.null(metadata_file)){
    if(grepl("\\.txt$|\\.tsv$", metadata_file)) {
      metadata_df <- utils::read.delim2(metadata_file, header = TRUE, sep = "\t")
    }else if(grepl("\\.xlsx$", metadata_file)) {
      metadata_df <- openxlsx::read.xlsx(metadata_file)
    } else {stop("Error, the metadata file you provided is not on .txt , .tsv or .xlsx format")}
   if (ncol(metadata_df) < 2) stop("Incompatible metadatafile: less than 2 columns.")
   colnames(metadata_df)[1] <- "Samples"
   colnames(metadata_df)[2] <- "Group"
   if (sum(is.na(metadata_df$Group)) > 0) stop("empty values in group column of the metadata file.")
   if (any(colSums(is.na(metadata_df)) > 0)) {
     warning("Columns with at least one missing value in metadata will be omitted.")
   }
   metadata_df <- metadata_df[, !colSums(is.na(metadata_df)) > 0]

   metadata_df$Samples <- make.names(metadata_df$Samples, unique = TRUE)
     metadata_df <- metadata_df %>% arrange(.,Group)
    groups_number <- length(unique(metadata_df$Group))
    group_names <- unique(metadata_df$Group)
    for (i in 1:groups_number) {
      samples_per_group[i] <- sum(metadata_df$Group == group_names[i])
    }
    for (i in 1:groups_number) {
      assign(paste0("g",i,".name"),group_names[[i]])} }
  if(length(group_paths) > 0 & is.null(file)){
    dataspace <- data.frame()
    group_paths<- gsub( "\\\\", "/", group_paths)

    if (length(group_paths) == 1 ){
      if (is.null(metadata_file)) stop("With only one directory with single sample protemics files, a metadata_file must be present")
      file_names <- list.files(path = group_paths[[1]], pattern = "^[^~].*\\.xlsx$")
      for (j in seq_along(file_names)) {
        file_case <- openxlsx::read.xlsx(file.path(group_paths[[1]], file_names[j]), sheet = 1)
        dataspace <- rbind(dataspace, file_case[,grep("Accession|Description",colnames(file_case))])
        dataspace <- dataspace %>%
          filter(!is.na(Accession)) %>%
          filter(!duplicated(Accession)) %>%
          droplevels()}
      for (j in seq_along(file_names)){
        file_case <- openxlsx::read.xlsx(file.path(group_paths[[1]], file_names[j]), sheet = 1)
        dataspace <- merge(x = dataspace,y = file_case[,grep("Accession|Area|Abundance:",colnames(file_case))], by = "Accession" ,all.x = TRUE)
        colnames(dataspace)[length(colnames(dataspace))] <- file_names[j]
      }

      path <- dirname(group_paths[[1]])

    } else {
      groups_number <- length(group_paths)
      for (i in 1:groups_number) {
        assign(paste0("group",i),group_paths[[i]])}
      group_names <- basename(group_paths)
      for (i in 1:groups_number) {
        assign(paste0("g",i,".name"),group_names[[i]])}


      samples_per_group <- numeric(groups_number)

      for (i in seq_along(group_paths)) {
        file_names <- list.files(path = group_paths[[i]], pattern = "^[^~].*\\.xlsx$")
        samples_per_group[i] <- length(file_names)
        for (j in seq_along(file_names)) {
          file_case <- openxlsx::read.xlsx(file.path(group_paths[[i]], file_names[j]), sheet = 1)
          dataspace <- rbind(dataspace, file_case[,grep("Accession|Description",colnames(file_case))])
        }
      }
      dataspace <- dataspace %>%
        filter(!is.na(Accession)) %>%
        filter(!duplicated(Accession)) %>%
        droplevels()

      for (i in 1:groups_number) {
        file_names <- list.files(path = group_paths[[i]], pattern = "^[^~].*\\.xlsx$")
        for (j in seq_along(file_names)){
          file_case <- openxlsx::read.xlsx(file.path(group_paths[[i]], file_names[j]), sheet = 1)
          dataspace <- merge(x = dataspace,y = file_case[,grep("Accession|Area|Abundance:",colnames(file_case))], by = "Accession" ,all.x = TRUE)
          colnames(dataspace)[length(colnames(dataspace))] <- file_names[j]
        }
      }
      path <- dirname(group1)

    }
    if (!("Description" %in% colnames(dataspace))) {
      dataspace$Description <- "Not available"
    }

    colnames(dataspace) <- gsub(".xlsx", "", colnames(dataspace))
    path_res <- file.path(path, "ProtE_Analysis")

  }else {
    if (is.null(file)) stop("file argument must be provided if you are not using ProteomeDiscoverer multiple files.")
    if ((is.null(group_names) | is.null(samples_per_group)) & is.null(metadata_file))
      stop("For single file mode,either 'metadata_file' or 'group_names' and 'samples_per_group' must be provided.")
    if (is.null(metadata_file)){
      groups_number <- length(group_names)
      if (length(samples_per_group) != groups_number)
        stop("The length of 'samples_per_group' must match the length of 'group_names'.")
      for (i in 1:groups_number) {
        assign(paste0("g",i,".name"),group_names[[i]])}}

    if(grepl("\\.txt$|\\.tsv$", file)) {
      dataspace <- utils::read.delim2(file, header = TRUE, sep = "\t")
    }else if(grepl("\\.xlsx$", file)) {
      dataspace <- openxlsx::read.xlsx(file)
    } else {stop("Error, the file you provided is not on .txt , .tsv or .xlsx format")}

    path <- dirname(file)
    path_res <- file.path(path , "ProtE_Analysis")

  }
if (groups_number  == 1) stop("multiple groups should be inserted for the ProtE analysis.")

  if ("Accession" %in% colnames(dataspace) & any(grepl("^Abundance:", colnames(dataspace))) ) {
    print("ProteomeDiscoverer consesus file mode.")
    dataspace <- dataspace[,grep("Accession|Description|Abundance:",colnames(dataspace))]
    colnames(dataspace) <- gsub("Abundance:.|:.Sample", "", colnames(dataspace))
    dataspace <- dataspace[complete.cases(dataspace[,1]),]

    if (!("Description" %in% colnames(dataspace))) {
      dataspace$Description <- "Not available"
    }
  } else if ("Protein.IDs" %in% colnames(dataspace) & any(grepl("^Intensity", colnames(dataspace))) ) {
    print("MaxQuant file mode.")
    dataspace <- dataspace[!grepl("^;",dataspace$Protein.IDs),]
    dataspace <- dataspace[complete.cases(dataspace[,1]),]
    if("Reverse" %in% colnames(dataspace)) {
      dataspace <- dataspace[dataspace$Reverse == ""|is.na(dataspace$Reverse),]
      print("Removed REV_proteins: Reverse peptide Identifications. If you want to exclude the CON (Contaminants) proteins and proteins only identified by site, do so manually from the imported table.")}

    if (any(grepl("^LFQ intensity", colnames(dataspace)))){
      dataspace <- dataspace[,grep("Protein.IDs|Fasta.headers|LFQ.intensity",colnames(dataspace))]
      colnames(dataspace) <- gsub("LFQ.intensity.", "", colnames(dataspace))

       } else {
     dataspace <- dataspace[,grep("Protein.IDs|Fasta.headers|Intensity",colnames(dataspace))]
     dataspace$Intensity <- NULL
     colnames(dataspace) <- gsub("Intensity.", "", colnames(dataspace))
    }
    dataspace$Protein.IDs <- ifelse(
      grepl("sp\\|\\w+\\|", dataspace$Protein.IDs),
      sub(".*sp\\|(\\w+)\\|.*", "\\1", dataspace$Protein.IDs),
      dataspace$Protein.IDs
    )
    if ("Fasta.headers" %in% colnames(dataspace)){
      dataspace$Description <- dataspace$Fasta.headers
      dataspace$Fasta.headers <- NULL
    } else { dataspace$Description <- "Not available" }

    colnames(dataspace)[colnames(dataspace) == "Protein.IDs"] <- "Accession"
    dataspace[dataspace == 0] <- NA
  } else if(length(group_paths) > 0 & is.null(file)){
    print("ProteomeDiscoverer multiple files mode.")
  } else {
    if("Protein.Ids" %in% colnames(dataspace)) {
      print("DIA-NN mode.")
      dataspace <- dataspace[!grepl("^;",dataspace$Protein.Ids),]
      if("First.Protein.Description" %in% colnames(dataspace)){
        dataspace <- dataspace[,c(grep("Protein.Ids",colnames(dataspace)),6:ncol(dataspace))]
      }else{dataspace <- dataspace[,c(grep("Protein.Ids",colnames(dataspace)),5:ncol(dataspace))]}
      dataspace[, -1] <- lapply(dataspace[, -1], as.numeric)
      dataspace <- dataspace[complete.cases(dataspace[,1]),]
      col_names <- colnames(dataspace)
      col_names[-1] <- gsub("\\\\", "/", col_names[-1])
      col_names[-1] <- basename(col_names[-1])
      colnames(dataspace) <- col_names
      colnames(dataspace)[colnames(dataspace) == "Protein.Ids"] <- "Accession"
      dataspace$Description <- "Not available"

    } else if ("Genes" %in% colnames(dataspace)) {
      print("DIA-NN mode.")
      dataspace <- dataspace[!grepl("^;", dataspace[["Genes"]]), ]
      dataspace[, -1] <- lapply(dataspace[,-1], as.numeric)
      col_names <- colnames(dataspace)
      col_names[-1] <- gsub("\\\\", "/", col_names[-1])
      col_names[-1] <- basename(col_names[-1])
      colnames(dataspace) <- col_names
      dataspace <- dataspace[rowSums(!is.na(dataspace[,-1])) > 0, ]
      colnames(dataspace)[colnames(dataspace) == "Genes"] <- "Gene.Symbol"
      dataspace$Description <- rep("Unavailable for Genes matrices", nrow(dataspace))
      dataspace$Accession <- rep("Unavailable for Genes matrices", nrow(dataspace))
      uqg <- TRUE
    } else stop("Incompatible input file was inserted.")
  }

  colnames(dataspace) <- make.names(colnames(dataspace), unique = TRUE)
  rownames(dataspace) <- make.names(rownames(dataspace), unique = TRUE)


  if (!"Gene.Symbol" %in% colnames(dataspace)){
    if(all(dataspace$Description == "Not available")){
      print("The Description fetching from UniProt starts now, it might take some time depending on your Network speed.")
      dataspace$First_Accession <- sub("([^;]*).*", "\\1", dataspace$Accession)

      id_numbers <-   dataspace$First_Accession

      query <- list("accession_id" = id_numbers)
      conv_ID <- queryup::query_uniprot(query, columns = c("accession","id","gene_names","organism_name","organism_id","protein_name","protein_existence","sequence_version"), max_keys = 200, show_progress = TRUE )
      First_Accession <- conv_ID$Entry
      details<-paste0(conv_ID$`Protein names`," OS=",conv_ID$Organism, " OX=", conv_ID$`Organism (ID)`, " GN=",conv_ID$`Gene Names`, " PE=", conv_ID$`Protein existence`, " SV=", conv_ID$`Sequence version`)
      dfup <- cbind(First_Accession, details)
      dataspace <- merge(dataspace, dfup,
                         by = "First_Accession",
                         all.x = TRUE, sort = FALSE)
      dataspace$Description <- ifelse(is.na(dataspace$details),
                                      dataspace$Description,
                                      dataspace$details)
      dataspace$details <- NULL
      dataspace$First_Accession <- NULL
      dataspace <- dataspace[!duplicated(dataspace), ]
    }
  }else { print("Description fetching is not available for DIA-NN unique_genes matrices")}
  if (uqg == FALSE){
  dataspace$Gene.Symbol = sub(".*GN=(.*?) .*","\\1",dataspace$Description)}
  dataspace<-dataspace %>%  dplyr::select(Accession, Description,Gene.Symbol , everything())



  dataspace <- dataspace[rowSums(!is.na(dataspace[,-c(1:3)])) > 0, ]
  if (!is.null(metadata_file)){
    if (setequal(metadata_df$Samples, colnames(dataspace[,-c(1:3)])) == FALSE) stop("The samples provided in the metadata file do not match the samples provided in the proteomics file.")
    col_order <- match(metadata_df$Samples, colnames(dataspace[,-c(1:3)]))
    dataspace <- dataspace[, c(1:3, col_order + 3)]
  }

  print("Removing whichever proteins have only missing values in their abundances.")

  dir.create(path_res, showWarnings = FALSE)

  path_restat <- file.path(path_res, "Statistical_Analysis")
  path_resman <- file.path(path_res, "Data_processing")
  path_resplot <- file.path(path_res, "Plots")

  dir.create(path_restat, showWarnings = FALSE)
  dir.create(path_resman, showWarnings = FALSE)
  dir.create(path_resplot, showWarnings = FALSE)
  print("All files created will be stored in the folder ProtE_Analysis, that is located on the last directory of the first group that you input.")

  if(length(group_paths) > 0 & is.null(file)){
    mt_file_path <- file.path(path_resman, "Masterlist.xlsx")
    openxlsx::write.xlsx(dataspace, file = mt_file_path)
    print("Concatenated all data files to a single matrix, saved as Masterlist.xlsx")
  }

  if (sum(samples_per_group) != ncol(dataspace)-3) {stop("Error: The specified number of samples in the parameter samples_per_group does not align with the total samples in the input file.")}

  zero_per_sample <- colSums(is.na(dataspace[,-1:-3]))*100/nrow(dataspace)
  IDs <- colSums(!is.na(dataspace[,-1:-3]))
  dat.dataspace<-dataspace

  p <- function(x) {
    substr(x, 1, pmin(30, nchar(x)))
  }

  groups_list <- list("character")

  for (i in 1:groups_number) {
    groups_list[[i]] <- rep(group_names[i], times = samples_per_group[i])
  }

  groups_list_u <- unlist(groups_list)
  groups_list_f <- factor(groups_list_u, levels = unique(groups_list_u))

  Group <- groups_list_f


  Group2<-unique(groups_list_f)





  log.dataspace <- log(dataspace[,-c(1:3)]+1,2)
  melt.log.dataspace <- reshape2::melt(log.dataspace, id.vars = NULL)

  repvec <- as.data.frame(table(Group))$Freq * nrow(log.dataspace)
  storevec <- NULL
  storeres <- list()


  for (i in seq_along(Group2)){
    storevev <- rep(Group2[i], repvec[i])
    storeres[[i]] <- storevev
  }

  melt.log.dataspace$Group <- unlist(storeres)

  melt.log.dataspace$Group <- factor(melt.log.dataspace$Group, levels = Group2)

  qc.boxplots<-ggplot2::ggplot(melt.log.dataspace, aes(x=forcats::fct_inorder(variable), y=value, color=Group))+
    geom_boxplot(aes(color = Group),lwd=1, outlier.size=0.2, outlier.alpha = 0.2)+
    xlab("Sample")+
    ylab(expression(Log[2]~"Protein Abundance"))+
    theme_classic()+
    theme(text = element_text(size = 19),
          axis.text.x=element_text(size = 9, angle=90, vjust = 0.5, hjust = 0.5),
          axis.title.x=element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))+
    guides(color = guide_legend(override.aes = list(size = 1)))+
    scale_x_discrete(labels = p) +
    geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5, alpha = 0.5)


  suppressWarnings(ggplot2::ggsave("Boxplot_before_proccesing.bmp", plot = qc.boxplots,  path = path_resplot,
                  scale = 1, width = 12, height = 5, units = "in",
                  dpi = 300, limitsize = TRUE, bg = "white"))


  print("A boxplot showing the log2 Abundance of each protein, across the samples for the unprocessed data has been created as Boxplot_before_proccesing.bmp" )

  qc.violin<-ggplot2::ggplot(melt.log.dataspace , aes(x=forcats::fct_inorder(variable), y=value, color=Group))+
    geom_violin(aes(color = Group),lwd=1)+
    xlab("Sample")+
    ylab(expression(Log[2]~"Protein Abundance"))+
    theme_classic()+
    theme(text = element_text(size = 19),
          axis.text.x=element_text(size = 9, angle=90, vjust = 0.5, hjust = 0.5),
          axis.title.x=element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))+
    guides(color = guide_legend(override.aes = list(size = 1)))+
    scale_x_discrete(labels = p) +
    geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5, alpha = 0.5)

  suppressWarnings(ggplot2::ggsave("Violin_plot_before_processing.bmp", plot = qc.violin,  path = path_resplot,
                                   scale = 1, width = 12, height = 5, units = "in",
                                   dpi = 300, limitsize = TRUE, bg = "white"))
  print("A Violin Plot showing the log2 Protein Abundance of each protein, across the samples for the unprocessed data was created as Violin_plot_before_processing.bmp" )


  row_means <- rowMeans(log.dataspace, na.rm = TRUE)
  row_sds <- apply(log.dataspace, 1, sd, na.rm = TRUE)
  plot_data <- data.frame(Mean = row_means, SD = row_sds)
  meansd <- ggplot(plot_data, aes(x = Mean, y = SD)) +
    geom_point(alpha = 0.5, color = "blue") +
    geom_smooth(method = "loess", color = "red", se = FALSE) +
    theme_minimal() +
    labs(title = "Mean-SD Plot on the log2 unprocessed data", x = "Mean Expression", y = "Standard Deviation")
  suppressMessages(suppressWarnings(ggplot2::ggsave("Unprocessed_data_meanSdPlot.bmp", plot = meansd,  path = path_resplot,
                                                    scale = 1, width = 5, height = 4, units = "in",
                                                    dpi = 300, limitsize = TRUE)))
  print("Creating Mean-SD plot on the log2 unprocessed data, saved as Unprocessed_data_meanSdPlot.bmp")



  if (normalization == "PPM"){
    dataspace[, -1:-3] <- lapply(dataspace[, -1:-3], function(x) {
      sum_x <- sum(x, na.rm = TRUE)
      ifelse(is.na(x), NA, (x / sum_x) * 10^6)
    })

    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(dataspace, file = norm_file_path)
    print("Applying the selected normalization, saved as Normalized.xlsx")}

  if (normalization == "Quantile"){
    dataspace[, -1:-3]<- log(dataspace[, -1:-3]+1,2)
    dataspace[, -1:-3] <- limma::normalizeQuantiles(dataspace[, -1:-3])
    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(dataspace, file = norm_file_path)
    print("Applying Quantile normalization (the data was fristly log2 transformed) saved as Normalized.xlsx") }
  if (normalization == "log2"){
    dataspace[, -1:-3] <- log(dataspace[, -1:-3]+1,2)
    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(dataspace, file = norm_file_path)
    print("Applying the selected normalization, saved as Normalized.xlsx") }
  if (normalization == "Total_Ion_Current") {
    dataspace[, -1:-3] <- lapply(dataspace[, -1:-3], function(x) (x / sum(x, na.rm = TRUE)) * mean(colSums(dataspace[, -1:-2], na.rm = TRUE)))
    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(dataspace, file = norm_file_path)
    print("Applying the selected normalization, saved as Normalized.xlsx")}

  if (normalization == "Cyclic_Loess"){
    dataspace[, -1:-3] <- log(dataspace[, -1:-3]+1,2)
    dataspace[, -1:-3] <- limma::normalizeCyclicLoess(dataspace[, -1:-3])
    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(dataspace, file = norm_file_path)
    print("Applying Cyclic Loess normalization (the data was fristly log2 transformed), saved as Normalized.xlsx") }

  if ( normalization == "VSN") {
    dataspace[, -1:-3] <- suppressMessages(limma::normalizeVSN(dataspace[, -1:-3]))
    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(dataspace, file = norm_file_path)
    print("Applying the selected normalization, saved as Normalized.xlsx")
  }
  if (normalization == "median") {
    sample_medians <- apply(dataspace[, -1:-3], 2, median, na.rm = TRUE)
    dataspace[, -1:-3] <- sweep(dataspace[, -1:-3], 2, sample_medians, FUN = "/")
    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(dataspace, file = norm_file_path)
    print("Applying the selected normalization, saved as Normalized.xlsx")}

  if (normalization %in% c("median","VSN", "Total_Ion_Current","PPM") ){
    log.dataspace <- log(dataspace[,-c(1:3)]+1,2)
    row_means <- rowMeans(log.dataspace, na.rm = TRUE)
    row_sds <- apply(log.dataspace, 1, sd, na.rm = TRUE)
    plot_data <- data.frame(Mean = row_means, SD = row_sds)
    meansd <- ggplot(plot_data, aes(x = Mean, y = SD)) +
      geom_point(alpha = 0.5, color = "blue") +
      geom_smooth(method = "loess", color = "red", se = FALSE) +
      theme_minimal() +
      labs(title = "Mean-SD Plot on the log2 normalized data", x = "Mean Expression", y = "Standard Deviation")

    suppressMessages(suppressWarnings(ggplot2::ggsave("normalized_meanSdPlot.bmp", plot = meansd,  path = path_resplot,
                    scale = 1, width = 5, height = 4, units = "in",
                    dpi = 300, limitsize = TRUE)))

  } else if (normalization %in% c("log2","Quantile", "Cyclic_Loess")) {
    row_means <- rowMeans(dataspace[,-c(1:3)], na.rm = TRUE)
    row_sds <- apply(dataspace[,-c(1:3)], 1, sd, na.rm = TRUE)
    plot_data <- data.frame(Mean = row_means, SD = row_sds)
    meansd <- ggplot(plot_data, aes(x = Mean, y = SD)) +
      geom_point(alpha = 0.5, color = "blue") +
      geom_smooth(method = "loess", color = "red", se = FALSE) +
      theme_minimal() +
      labs(title = "Mean-SD Plot on the normalized data", x = "Mean Expression", y = "Standard Deviation")

    suppressMessages(suppressWarnings(ggplot2::ggsave("normalized_meanSdPlot.bmp", plot = meansd,  path = path_resplot,
                                                      scale = 1, width = 5, height = 4, units = "in",
                                                      dpi = 300, limitsize = TRUE)))
  }



  if (filtering_value < 0 && filtering_value > 100) {stop("Error: The filtering_value must be a number ranging from 0 to 100")}

  if (global_filtering == TRUE) {

    filtering_value <- 100- as.numeric(filtering_value)

    threshold <-  ceiling(sum(samples_per_group)-(sum(samples_per_group)*(as.numeric(filtering_value)/100))+0.00000000001)
  }

  if (global_filtering == FALSE) {
    threshold<-numeric(groups_number)
    filtering_value <- 100- as.numeric(filtering_value)
    for (i in 1:groups_number) {
      threshold[i] <-  ceiling(samples_per_group[i]-(samples_per_group[i])*(as.numeric(filtering_value)/100)+0.00000000001)}
  }



  coln <- list()
  case_last <- 3

  for (i in 1:groups_number) {
    case_last <- case_last + samples_per_group[i]
    coln[[i]] <- (case_last - samples_per_group[i] + 1):case_last
  }


  dataspace[is.na(dataspace)] <- 0

  for (j in 1:groups_number){

    for (i in c(1:nrow(dataspace))){
      a <- table(dataspace[i,coln[[j]]]==0)["TRUE"]
      if(is.na(a)){
        dataspace[i,paste0("Number_0_group",j)] <- 0
      } else{
        dataspace[i,paste0("Number_0_group",j)] <- table(dataspace[i,coln[[j]]]==0)["TRUE"]
      }
    }
  }
  dataspace$Number_0_all_groups <- rowSums(dataspace[,paste0("Number_0_group", 1:groups_number)])
  dataspace <- dataspace[dataspace$Number_0_all_groups < sum(samples_per_group),]

  if (global_filtering == TRUE) {
    dataspace <- dataspace[dataspace$Number_0_all_groups<threshold,]
    at_file_path <- file.path(path_resman, "Dataset_after_filtering.xlsx")
    openxlsx::write.xlsx(dataspace, file = at_file_path)
  }

  if (global_filtering == FALSE) {
    keep_rows <- rep(FALSE, nrow(dataspace))
    for (j in 1:groups_number) {
      keep_rows <- keep_rows | (dataspace[,paste0("Number_0_group", j)] < threshold[j])
    }
    dataspace <- dataspace[keep_rows, ]
    at_file_path <- file.path(path_resman, "Dataset_after_filtering.xlsx")
    openxlsx::write.xlsx(dataspace, file = at_file_path)}

  print("An excel file with the proteins remaining in the data after filtering for missing values, has been created as Dataset_after_filtering.xlsx")
  dataspace_0s<- dataspace
  dataspace[,paste0("Number_0_group", 1:groups_number)] <- NULL
  dataspace$Number_0_all_groups <- NULL


  zero_per_sample1 <- colSums(dataspace[,-1:-3] == 0)*100/nrow(dataspace)
  sample_names <- colnames(dataspace[,-1:-3])
  qc <- cbind(sample_names,IDs,zero_per_sample,zero_per_sample1)
  qc <- as.data.frame(qc)
  colnames(qc) <- c("Sample Name","Number of proteins detected in the sample","% of Missing values before filtering","% of Missing values after filtering")
  rownames(qc) <- NULL


  pre_dataspace <- dataspace

  if (imputation == "kNN") {
    print("kNN imputation starts now")
    dataspace[dataspace==0] <- NA
    dataspace[,-1:-3] <- VIM::kNN(dataspace[,-1:-3], imp_var = FALSE, k= 5)
    imp_file_path <- file.path(path_resman, "Dataset_Imputed.xlsx")
    openxlsx::write.xlsx(dataspace, file = imp_file_path)
  }
  if (imputation == "zeros"){
    imp_file_path <- file.path(path_resman, "Dataset_Imputed.xlsx")
    openxlsx::write.xlsx(dataspace, file = imp_file_path)    }
  if (imputation == "Gaussian_LOD") {
    dataspace[dataspace==0] <- NA
    LOD <- min(as.matrix(dataspace[,-1:-3]),na.rm = TRUE)
    replace_gaussian <- function(x) {
      ifelse(is.na(x), rnorm(length(x), mean = LOD, sd = LOD*0.1), x)
    }
    dataspace[,-1:-3] <- apply(dataspace[,-1:-3], 2, replace_gaussian)
    imp_file_path <- file.path(path_resman, "Dataset_Imputed.xlsx")
    openxlsx::write.xlsx(dataspace, file = imp_file_path)
  }
  if (imputation == "Gaussian_mean_sd") {
    dataspace[dataspace == 0] <- NA

    data_to_impute <- dataspace[, 4:ncol(dataspace)]

    for (i in seq_len(nrow(data_to_impute))) {
      row_data <- data_to_impute[i, ]

      if (any(is.na(row_data))) {
        row_data <- as.numeric(row_data)
        row_mean <- mean(row_data, na.rm = TRUE)
        sd_row <- sd(row_data, na.rm = TRUE)

        num_na <- sum(is.na(row_data))
        imputed_values <- rnorm(num_na, mean = row_mean, sd = sd_row * 0.3)
        imputed_values[imputed_values < 0] <- 0

        row_data[is.na(row_data)] <- imputed_values
        data_to_impute[i, ] <- row_data
      }
    }

    dataspace[, 4:ncol(dataspace)] <- data_to_impute

    imp_file_path <- file.path(path_resman, "Dataset_Imputed.xlsx")
    openxlsx::write.xlsx(dataspace, file = imp_file_path)
  }


  if (imputation == "mean"){
    dataspace[dataspace==0] <- NA

    data_to_impute <- dataspace[, 4:ncol(dataspace)]
    for (i in seq_len(nrow(data_to_impute))) {
      row_data <- data_to_impute[i, ]

      if (any(is.na(row_data))) {
        row_data  <- as.numeric(row_data)
        row_mean <- mean(row_data, na.rm = TRUE)
        row_data[is.na(row_data)] <- row_mean
        data_to_impute[i, ] <- row_data
      }
    }
    dataspace[, 4:ncol(dataspace)] <- data_to_impute
    imp_file_path <- file.path(path_resman, "Dataset_Imputed.xlsx")
    openxlsx::write.xlsx(dataspace, file = imp_file_path) }
  if (imputation == "LOD"){
    dataspace[dataspace==0] <- NA
    impute_value <- min(as.matrix(dataspace[,-1:-3]),na.rm = TRUE)
    dataspace[,-1:-3][is.na(dataspace[,-1:-3])]  <- impute_value
    imp_file_path <- file.path(path_resman, "Dataset_Imputed.xlsx")
    openxlsx::write.xlsx(dataspace, file = imp_file_path)  }
  if (imputation == "LOD/2"){
    dataspace[dataspace==0] <- NA
    impute_value <- min(as.matrix(dataspace[,-1:-3]),na.rm = TRUE)/2
    dataspace[,-1:-3][is.na(dataspace[,-1:-3])]  <- impute_value
    imp_file_path <- file.path(path_resman, "Dataset_Imputed.xlsx")
    openxlsx::write.xlsx(dataspace, file = imp_file_path)  }
  if(imputation == "missRanger"){
    print("missRanger imputation starts now")
    dataspace[dataspace==0] <- NA
    dataspace[,-1:-3] <- missRanger::missRanger(dataspace[,-1:-3])
    imp_file_path <- file.path(path_resman, "Dataset_Imputed.xlsx")
    openxlsx::write.xlsx(dataspace, file = imp_file_path)  }
  if(imputation == FALSE){    dataspace[dataspace==0] <- NA}


  if (imputation %in% c("kNN","missRanger","mean","Gaussian_mean_sd", "Gaussian_LOD"))    {
    pre_dataspace1<-pre_dataspace[,-1:-3]
    dataspace1<-dataspace[,-1:-3]
    imp.values<- dataspace1 - pre_dataspace1

    his_dataspace<-rbind(dataspace1,pre_dataspace1,imp.values)
    x = "Protein Abundance"
    if (normalization %in% c(FALSE,"median", "Total_Ion_Current", "PPM") ){
      his_dataspace<-log2(his_dataspace+1)
      x= expression(Log[2]~"Protein Abundance")}


    his_long <-tidyr::pivot_longer(his_dataspace, cols = everything())
    nrows<-nrow(his_long)
    his_long$Group <- rep(c("Final","Initial","Imputed"), each = (nrows/3))
    his_long_filtered <- his_long[his_long$value != 0,]
    his_long_filtered$Group <- factor(his_long_filtered$Group, levels = c("Final", "Initial", "Imputed"))


    imp_hist<- ggplot(his_long_filtered, aes(x = value, fill = Group, colour = Group)) +
      labs( x = "Protein Abundance", y = "Count") +
      scale_fill_manual(values = c("Final" = "#FF99FF", "Initial" = "#990000", "Imputed" = "#000033")) +
      scale_color_manual(values = c("Final" = "#FF99FF", "Initial" = "#990000", "Imputed" = "#000033")) +
      geom_histogram(alpha = 0.5, binwidth = 0.5,  position = "identity") +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white")) + theme(panel.background = element_rect(fill = "white"))

    imp_hist
    ggplot2::ggsave("Imputed_values_histogram.bmp", plot = imp_hist,  path = path_resplot,
                    scale = 1, width = 5, height = 4, units = "in",
                    dpi = 300, limitsize = TRUE)
    print("A frequency histogram of real and imputed values has been created as Imputed_values_histogram.bmp")

  }
  if (imputation %in% c("LOD/2","LOD","kNN","missRanger","mean","zeros","Gaussian_LOD")){
    print("An excel with the imputed missing values has been created as Dataset_Imputed.xlsx")

    dataspace_0s$percentage <- dataspace_0s$Number_0_all_groups*100/sum(samples_per_group)
    dataspace$percentage <- dataspace_0s$percentage
    dataspace$mean <- rowMeans(dataspace[,4:(3+sum(samples_per_group))])
    dataspace$log<-log2(dataspace$mean)
    dataspace$rank <- rank(-dataspace$mean)

    abund.plot <- ggplot(dataspace, aes(x = rank, y = log, colour = percentage)) +
      geom_point(size = 3, alpha = 0.8) +
      labs( x = "Proteins mean abundance rank", y= expression(Log[2]~"mean protein abundance")) +
      scale_color_gradient(low = "darkblue", high = "yellow",
                           name = "Imputations\nin each\nprotein\n(%)",limits = c(0,100-filtering_value)) +
      theme_linedraw()+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid = element_line(color = "grey80"),
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 9))
    abund.plot

    suppressWarnings(ggplot2::ggsave("Proteins_abundance_rank.bmp", plot = abund.plot ,  path = path_resplot,
                    scale = 1, width = 12, height = 5, units = "in",
                    dpi = 300, limitsize = TRUE, bg = "white"))
  }

  if (imputation == FALSE){
    dataspace_0s$percentage <- dataspace_0s$Number_0_all_groups*100/sum(samples_per_group)
    dataspace$percentage <- dataspace_0s$percentage
    dataspace[,4:(3+sum(samples_per_group))] <- lapply(dataspace[,4:(3+sum(samples_per_group))], as.numeric)
    dataspace$mean <- apply(dataspace[,4:(3+sum(samples_per_group))], 1, function(x) mean(x[!is.na(x)]))
    dataspace$log<-log2(dataspace$mean)
    dataspace$rank <- rank(-dataspace$mean)

    abund.plot <- ggplot(dataspace, aes(x = rank, y = log, colour = percentage)) +
      geom_point(size = 3, alpha = 0.8) +
      labs( x = "Proteins mean abundance rank", y= expression(Log[2]~"mean protein abundance")) +
      scale_color_gradient(low = "darkblue", high = "yellow",
                           name = "MVs\nin each\nprotein\n(%)",limits = c(0,100-filtering_value)) +
      theme_linedraw()+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid = element_line(color = "grey80"),
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 9))

    abund.plot

    suppressWarnings(ggplot2::ggsave("Proteins_abundance_rank.bmp", plot = abund.plot ,  path = path_resplot,
                    scale = 1, width = 12, height = 5, units = "in",
                    dpi = 300, limitsize = TRUE, bg = "white"))
    dataspace[dataspace==0] <- NA

  }
  print("A protein abundance rank order plot has been created as Proteins_abundance_rank.bmp")


  dataspace$percentage <- NULL
  dataspace$mean <- NULL
  dataspace$log<- NULL
  dataspace$rank <- NULL

  if (normalization %in% c(FALSE,"median","VSN", "Total_Ion_Current","PPM") ){
    log.dataspace <- log(dataspace[,-c(1:3)]+1,2)
  } else {
    log.dataspace <- dataspace[,-c(1:3)]
  }
  melt.log.dataspace <- reshape2::melt(log.dataspace, id.vars = NULL)

  repvec <- as.data.frame(table(Group))$Freq * nrow(log.dataspace)
  storevec <- NULL
  storeres <- list()


  for (i in seq_along(Group2)){
    storevev <- rep(Group2[i], repvec[i])
    storeres[[i]] <- storevev
  }

  melt.log.dataspace$Group <- unlist(storeres)

  melt.log.dataspace$Group <- factor(melt.log.dataspace$Group, levels = Group2)

  qc.boxplots<-ggplot2::ggplot(melt.log.dataspace, aes(x=forcats::fct_inorder(variable), y=value, color=Group))+
    geom_boxplot(aes(color = Group),lwd=1, outlier.size=0.2, outlier.alpha = 0.2)+
    xlab("Sample")+
    ylab(expression(Log[2]~"Protein Abundance"))+
    theme_classic()+
    theme(text = element_text(size = 19),
          axis.text.x=element_text(size = 9, angle=90, vjust = 0.5, hjust = 0.5),
          axis.title.x=element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))+
    guides(color = guide_legend(override.aes = list(size = 1)))+
    scale_x_discrete(labels = p) +
    geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5, alpha = 0.5)


  suppressWarnings(ggplot2::ggsave("After_Processing_Boxplot.bmp", plot = qc.boxplots,  path = path_resplot,
                                   scale = 1, width = 12, height = 5, units = "in",
                                   dpi = 300, limitsize = TRUE, bg = "white"))


  print("A boxplot showing the log2 Abundance of each protein, across the samples after processing has been created as After_Processing_Boxplot.bmp" )

  qc.violin<-ggplot2::ggplot(melt.log.dataspace, aes(x=forcats::fct_inorder(variable), y=value, color=Group))+
    geom_violin(aes(color = Group),lwd=1)+
    xlab("Sample")+
    ylab(expression(Log[2]~"Protein Abundance"))+
    theme_classic()+
    theme(text = element_text(size = 19),
          axis.text.x=element_text(size = 9, angle=90, vjust = 0.5, hjust = 0.5),
          axis.title.x=element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))+
    guides(color = guide_legend(override.aes = list(size = 1)))+
    scale_x_discrete(labels = p) +
    geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5, alpha = 0.5)

  suppressWarnings(ggplot2::ggsave("After_Processing_violin_plot.bmp", plot = qc.violin,  path = path_resplot,
                                   scale = 1, width = 12, height = 5, units = "in",
                                   dpi = 300, limitsize = TRUE, bg = "white"))
  print("A Violin Plot showing the log2 Protein Abundance of each protein, across the samples was created as After_Processing_violin_plot.bmp" )

  if(!is.null(metadata_file)){
    metadata_df <- metadata_df[, sapply(metadata_df, function(x) length(unique(x)) >= 2)]

     if (ncol(metadata_df) > 2){
      max_cols_to_process <- min(5, ncol(metadata_df) - 2)
      cov <- list()
      for (i in 1:max_cols_to_process) {
        if ( is.numeric(metadata_df[, i + 2]) == TRUE){
          cov[[i]] <- as.numeric(metadata_df[, i + 2])
        } else {
          cov[[i]] <- as.factor(metadata_df[, i + 2])
        }}
      cov_terms <- paste0("cov[[", 1:max_cols_to_process, "]]", collapse = " + ")
      formula <- as.formula(paste("~ 0 + groups_list_f +", cov_terms))
      mm <- model.matrix(formula)

      rank_mm <- qr(mm)$rank
      ncol_mm <- ncol(mm)
      if (rank_mm < ncol_mm) {
        qr_obj <- qr(mm)
        non_estimable <- colnames(mm)[-qr_obj$pivot[1:rank_mm]]
        non_estimable_short <- sub(".*\\]\\]", "", non_estimable)
        cov_remove <- unique(as.integer(sub("cov\\[\\[(\\d+)\\]\\].*", "\\1",
                                            non_estimable[grepl("cov\\[\\[", non_estimable)])))
        if (length(cov_remove) > 0) {
          cov_columns_to_remove <- unlist(lapply(cov_remove, function(i) {
            grep(paste0("cov\\[\\[", i, "\\]\\]"), colnames(mm), value = TRUE)
          }))
          mm <- mm[, !colnames(mm) %in% cov_columns_to_remove]
          excluded_metadata_cols <- colnames(metadata_df)[cov_remove + 2]
          warning("Non-estimable coefficients:", paste(non_estimable_short, collapse = ", "), "\n",
                  "Excluded entire covariates from the linear model based on:",
                  paste(excluded_metadata_cols, collapse = ", "), "\n")
        }
      }
      for (i in 1:max_cols_to_process) {
        names(cov)[[i]] <- colnames(metadata_df)[ i + 2]
        colnames(mm) <- gsub(paste0("cov\\[\\[", i, "\\]\\]") , names(cov)[[i]] , colnames(mm))
      }
      colnames(mm) <- gsub(paste0("cov\\[\\[", i, "\\]\\]") , paste0(names(cov)[[i]], "_") , colnames(mm))
     groups_coef <- grep("^groups_list_f", colnames(mm), value = TRUE)
     }else {
       mm <- model.matrix(~groups_list_f + 0)
       colnames(mm)<- group_names
     }
    } else {
      mm <- model.matrix(~groups_list_f + 0)
      colnames(mm)<- group_names
    }
  colnames(mm) <- make.names(colnames(mm))



  if (independent != FALSE && independent != TRUE){stop("Error: Sample relationship between groups should be defined in the parameter 'independent' either as TRUE or FALSE")}

  nndataspace<- dataspace[,-1:-3]
  if (imputation == FALSE ) {
    warning("To perform limma statistical analysis, without imputation NA coefficients are expected for certain proteins.")
  }
  nndataspace[nndataspace < 0] <- 0
  if (normalization %in% c(FALSE,"median", "VSN","Total_Ion_Current", "PPM") ){
    nndataspace <- log2(nndataspace+1)}
  if (independent == FALSE){
    if (length(unique(samples_per_group)) != 1){
      stop("Error: The paired group analysis requires an equal number of samples in each group. Please consider removing any samples that lack a corresponding matched pair.")
    }
    n = sum(samples_per_group)/groups_number
    pairing <- rep(1:n, each = groups_number)

    corfit <- limma::duplicateCorrelation(nndataspace, design = mm, block = pairing)
    fit <- limma::lmFit(nndataspace, mm, block = pairing, correlation = corfit$consensus.correlation)
  }
  if (independent == TRUE){
    fit <- limma::lmFit(nndataspace, mm)}
  fit<- limma::eBayes(fit)

  lima.res <- data.frame()

  if (!is.null(metadata_file)){
    if (ncol(metadata_df) > 2){
      if (groups_number == 2){
        coef_res <- limma::topTable(fit, adjust.method = p.adjust.method, number = Inf, sort.by = "none")
        colnames(coef_res)<-gsub("groups_list_f","",colnames(coef_res))
        colnames(coef_res)<- paste("coef",colnames(coef_res), sep = "_")
        coef_res <- coef_res[,1:ncol(mm)]
        lima.res <- coef_res
      }}}
  if (groups_number>2){
    anova_res <- data.frame()
    if (!is.null(metadata_file)){
      if (ncol(metadata_df) > 2){
        if (ncol(mm) > length(groups_coef)){
          coef_res <- data.frame()
          group_coef_indices <- match(grep("^groups_list_f", colnames(mm), value = TRUE), colnames(fit$coefficients))
          coef_res <- limma::topTable(fit, adjust.method = p.adjust.method, number = Inf, sort.by = "none")
          anova_res<- limma::topTable(fit, adjust.method = p.adjust.method, number = Inf, coef = group_coef_indices, sort.by = "none")
          colnames(coef_res)<-gsub("groups_list_f","",colnames(coef_res))
          colnames(coef_res)<- paste("coef",colnames(coef_res), sep = "_")
          coef_res <- coef_res[,1:ncol(mm)]
          anova_res <- anova_res[,-c(1:groups_number)]
          colnames(anova_res)<-paste("F-test",colnames(anova_res), sep = "_")
          anova_res <- cbind(coef_res,anova_res)
        } else {    anova_res<- limma::topTable(fit, adjust.method = p.adjust.method, number = Inf, sort.by = "none")
        colnames(anova_res)[-c(1:groups_number)]<-paste("F-test",colnames(anova_res)[-c(1:groups_number)], sep = "_")
        colnames(anova_res)[c(1:groups_number)]<- paste("coef",colnames(anova_res)[c(1:groups_number)], sep = "_")

        }
        }} else {    anova_res<- limma::topTable(fit, adjust.method = p.adjust.method, number = Inf, sort.by = "none")
        colnames(anova_res)[-c(1:groups_number)]<-paste("F-test",colnames(anova_res)[-c(1:groups_number)], sep = "_")
        colnames(anova_res)[c(1:groups_number)]<- paste("coef",colnames(anova_res)[c(1:groups_number)], sep = "_")

        }
  }

  for (i in 1:(groups_number-1)) {
    for (j in (i+1):groups_number) {
      comparison <- paste(colnames(mm)[i], "vs", colnames(mm)[j], sep = " ")
      comparison <- gsub("groups_list_f","",comparison)
      contrast_fref <- limma::makeContrasts(contrasts = paste0(colnames(mm)[i],"-",colnames(mm)[j]), levels = mm)
      fit2 <- limma::contrasts.fit(fit, contrast_fref)
      fit2 <- limma::eBayes(fit2)
      top_table<- limma::topTable(fit2, adjust.method = p.adjust.method, number = Inf, sort.by = "none")
      column_groups<- top_table[,c("logFC","AveExpr","t","P.Value","adj.P.Val","B")]
      colnames(column_groups)<-paste(colnames(column_groups), comparison, sep = "_")

      if (nrow(lima.res) == 0){
        lima.res <- column_groups
      } else {lima.res<- cbind(lima.res, column_groups)}
    }}

  if (groups_number>2){
    limma_dataspace <- cbind(anova_res,lima.res,dataspace)} else {limma_dataspace <- cbind(lima.res,dataspace)}
  ncollimma <- ncol(limma_dataspace) - ncol(dataspace) + 3
  limma_dataspace<-limma_dataspace %>%
    dplyr::select(Accession ,Description ,Gene.Symbol, everything())
  limma_dataspace <- limma_dataspace[,1:ncollimma]
  limma_file_path <- file.path(path_restat, "limma_statistics.xlsx")
  openxlsx::write.xlsx(limma_dataspace, file = limma_file_path)
  print("The limma output has been saved as limma_statistics.xlsx")

  data2 <- dataspace

  for (j in 1:groups_number) {
    data2[[paste0("Average_G", j)]] <- apply(data2[, coln[[j]]], 1, function(row) {
      if (all(is.na(row))) {
        0
      } else {
        mean(row, na.rm = TRUE)
      }
    })

    data2[[paste0("St_Dv_G", j)]] <- apply(data2[, coln[[j]]], 1, function(row) {
      if (sum(!is.na(row)) < 2) {
        0
      } else {
        sd(row, na.rm = TRUE)
      }
    })

    for (k in 1:j) {
      if (k < j) {
        for (i in 1:nrow(data2)) {
          values_k <- as.numeric(data2[i, coln[[k]]])
          values_j <- as.numeric(data2[i, coln[[j]]])

          if (independent == TRUE) {
            if (sum(!is.na(values_k)) > 0 & sum(!is.na(values_j)) > 0) {
              test_list <- stats::wilcox.test(
                values_k, values_j,
                exact = FALSE, paired = FALSE, na.rm = TRUE
              )
              data2[i, paste0("Wilcoxon_p_G", j, "vsG", k)] <- test_list$p.value
            } else {
              data2[i, paste0("Wilcoxon_p_G", j, "vsG", k)] <- NA
            }
          } else if (independent == FALSE) {
            paired_values <- na.omit(cbind(values_k, values_j))
            if (nrow(paired_values) > 0) {
              test_list <- stats::wilcox.test(
                paired_values[, 1], paired_values[, 2],
                exact = TRUE, paired = TRUE
              )
              data2[i, paste0("Wilcoxon_p_G", j, "vsG", k)] <- test_list$p.value
            } else {
              data2[i, paste0("Wilcoxon_p_G", j, "vsG", k)] <- NA
            }
          }
        }

        data2[[paste0("p.adj_p_G", j, "vsG", k)]] <- p.adjust(
          data2[[paste0("Wilcoxon_p_G", j, "vsG", k)]],
          method = p.adjust.method
        )

        avg_j <- data2[[paste0("Average_G", j)]]
        avg_k <- data2[[paste0("Average_G", k)]]
        data2[[paste0("Ratio_G", j, "vsG", k)]] <- ifelse(
          avg_k == 0, NA, avg_j / avg_k
        )
        data2[[paste0("Ratio_G", k, "vsG", j)]] <- ifelse(
          avg_j == 0, NA, avg_k / avg_j
        )
        if (normalization %in% c(FALSE,"median","VSN", "Total_Ion_Current","PPM") ){
          data2[[paste0("Log2_Ratio_G", j, "vsG", k)]] <- log2(
            data2[[paste0("Ratio_G", j, "vsG", k)]])
          data2[[paste0("Log2_Ratio_G", k, "vsG", j)]] <- log2(
            data2[[paste0("Ratio_G", k, "vsG", j)]])
        } else {
          data2[[paste0("Log2_Ratio_G", j, "vsG", k)]] <- avg_j - avg_k
          data2[[paste0("Log2_Ratio_G", k, "vsG", j)]] <- avg_k - avg_j
        }
      }
    }
  }

  for(i in 1:nrow(data2)){
    group_values <- list()
    for (j in 1:groups_number) {
      group_values[[j]]<- as.numeric(data2[i, coln[[j]]])
    }
    valid_groups <- group_values[sapply(group_values, function(x) sum(!is.na(x)) > 1)]
    if (length(valid_groups) >= 2) {
      bartlett_result <- bartlett.test(valid_groups)
     data2[i, "Bartlett_p"] <- bartlett_result$p.value
    leve_data <- data.frame(
      Value = unlist(valid_groups),
       Group = rep(seq_along(valid_groups), sapply(valid_groups, length))  )

   levene_result <- suppressWarnings(car::leveneTest(Value ~ as.factor(Group), data = leve_data, center = median))
   data2[i, "Levene_p"] <- levene_result$`Pr(>F)`[1]
    } else {
     data2[i, "Bartlett_p"] <- NA
     data2[i, "Levene_p"] <- NA
     }

  }

  only.data <- dataspace[,-c(1:3)]
  transposed_data <- t(only.data)
  metadata2 <- data.frame(group = groups_list_u)
  rownames(transposed_data) <- colnames(only.data)
  metadata2$samples <- colnames(only.data)
  adonis2_results <- vegan::adonis2(transposed_data ~ group, data = metadata2, method = "bray", na.rm = TRUE, permutations = 999)
  permanova_psF <- adonis2_results[1,4]
  permanova_pValue <- adonis2_results[1,5]
  permanova_modelR2 <- adonis2_results[1,3]
  data2[1, "PERMANOVA_PseudoF"] <- permanova_psF
  data2[1, "PERMANOVA_p"] <- permanova_pValue
  data2[1, "PERMANOVA_modelR2"] <- permanova_modelR2



  Ddataspace<-data2
  Fdataspace <- data2



  transposed_data<-data.frame(transposed_data)
  dataspace<-data.frame(dataspace)
  colnames(transposed_data)<-dataspace[,1]
  colnames(transposed_data) <- make.names(colnames(transposed_data), unique = TRUE)
  dataspace4<-cbind(Group,transposed_data)
  dataspace4$Group<-as.factor(dataspace4$Group)




  if (groups_number > 2) {
    if (independent == TRUE) {
      print("Mann-Whitney, Levene's and Bartlett's tests have been completed, calculating the Kruskal-Wallis p-values:")
      df3 <- dataspace4 %>% tidyr::gather(key, value, -Group)
      df4 <- df3 %>% dplyr::group_by(key)
      df4$value <- as.numeric(df4$value)
      df5 <- df4 %>%
        dplyr::group_by(key) %>%
        dplyr::group_modify(~ {
          valid_data <- .x[!is.na(.x$value), ]
          if (length(unique(valid_data$Group)) >= 2) {
            broom::tidy(kruskal.test(value ~ Group, data = valid_data))
          } else {
            data.frame(
              statistic = NA,
              p.value = NA,
              parameter = NA,
              method = "Kruskal-Wallis rank sum test"
            )
          }
        }, .keep = TRUE) %>%
        dplyr::select(key, p.value)
      data3 <- merge(Ddataspace, df5, by.x = colnames(Ddataspace)[1], by.y = "key", all.x = TRUE)
      data3 <- data3[match(Ddataspace[, 1], data3[, 1]), ]
      test_type <- "Kruskal_Wallis"


    }  else if (independent == FALSE) {
      print("Wilcoxon, Levene's and Bartlett's tests have been completed, performing Friedman test for paired samples:")
      df3 <- dataspace4 %>% tidyr::gather(key, value, -Group)
      df4 <- df3 %>% dplyr::group_by(key)
      df4$value <- as.numeric(df4$value)
      samples_fried <- rep(1:samples_per_group[1],nrow(dataspace)*groups_number)
      df4<- suppressMessages(cbind(df4, samples_fried ))
      colnames(df4)[ncol(df4)] <- "sample"
      keys <- unique(df4$key)
      p_values <- numeric(length(keys))
      for (i in seq_along(keys)) {
        gene_data <- df4[df4$key == keys[i], ]
        p_values[i] <- tryCatch({
          friedman.test(value ~ Group | sample, data = gene_data)$p.value
        }, error = function(e) NA)
      }
      Test.pvalue <- p_values
      test_type <- "Friedman"
      data3 <- cbind(Ddataspace, Test.pvalue)
    }

    colnames(data3)[ncol(data3)] <- paste0(test_type, ".pvalue")
    data3[[paste0(test_type, ".pvalue_adjusted")]] <- p.adjust(data3[[paste0(test_type, ".pvalue")]], method = p.adjust.method)
    Fdataspace <- data3
  }



  Fdataspace <- Fdataspace %>%
    dplyr::select(
      Accession, Description, Gene.Symbol,
      everything()
    )

  for (i in 1:groups_number){
    namesc<- colnames(Fdataspace)
    namesc<- gsub(paste0("G",i), get(paste0("g",i,".name")),namesc)
    colnames(Fdataspace)<-namesc
  }

  start_col <- 3 + as.numeric(sum(samples_per_group))
  Fdataspace <- Fdataspace[,-c(4:start_col)]

  stats_file_path <- file.path(path_restat, "traditional_statistics.xlsx")
  openxlsx::write.xlsx(Fdataspace, file = stats_file_path)
  print("The non-parametric statistical output along with tests for homoscedasticity have been saved as Statistical_analysis.xlsx")


  if (normalization %in% c(FALSE,"median","VSN", "Total_Ion_Current","PPM") ){
    log.dataspace <- log(dataspace[,-c(1:3)]+1,2)}
pca.log.dataspace <- log.dataspace

  if (imputation == FALSE ) {
    pca.log.dataspace[is.na(pca.log.dataspace)] <- 0
    hcluster <- FALSE
    warning("To perform principal component analysis, without imputation any missing values will be treated as 0s.")
  } else { hcluster <- TRUE }


  pca<-prcomp(t(pca.log.dataspace), scale=TRUE, center=TRUE)
  pca.data <- data.frame(Sample=rownames(pca$x),
                         X=pca$x[,1],
                         Y=pca$x[,2],
                         Group = Group)
  qc$PC1.score <- pca$x[,1]
  qc$PC2.score <-pca$x[,2]



  pca.var<-pca$sdev^2

  pca.var.per<-round(pca.var/sum(pca.var)*100,1)

  pca.data$Group<-factor(pca.data$Group, levels=Group2)

  pca.ent<-ggplot2::ggplot(data=pca.data, ggplot2::aes(x=X, y=Y, label=Sample))+
    geom_point(aes(color=Group), size = 2, alpha = 1)+
    xlab(paste("PC1 - ", pca.var.per[1], "%", sep=""))+
    ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+
    ggtitle("Complete set of proteins")+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5))+
    guides(color = guide_legend(override.aes = list(size=5)))+
    theme(legend.text = element_text(size = 12))+
    theme(legend.title = element_text(size = 12))+
    theme(legend.position="right")

  pca.ent

  ggplot2::ggsave("PCA_plot_alldata.bmp", plot = pca.ent,  path = path_resplot,
                  scale = 1, width = 5, height = 4, units = "in",
                  dpi = 300, limitsize = TRUE)
  print("PCA plot using all post-processing proteins has been created as PCA_plot_alldata.bmp")

  qc[,-1] <- lapply(qc[,-1], function(x) as.numeric(unlist(x)))
  qc[,-1]<-round(qc[,-1],3)
  qc_file_path <- file.path(path_restat, "Sample_QC.xlsx")
  openxlsx::write.xlsx(qc, file = qc_file_path)
  print("Sample quality metrics and association scores to the first two Principal Components have been saved as Sample_QC.xlsx")

  which.sig<-vector()
  groups_list_f <- factor(groups_list_f, levels=c(unique(groups_list_f)))

  if (parametric == TRUE) {
    tn <- "F-test_limma"
    if (significance == "p"){
      if (groups_number != 2){
        which.sig <- which(limma_dataspace$ANOVA_P.Value < 0.05)
      } else {(which.sig <- which(limma_dataspace[,grep("P.Value",colnames(limma_dataspace))] < 0.05))}
    }
    if (significance == "p.adj"){
      if (groups_number != 2){
        which.sig <- which(limma_dataspace$ANOVA_adj.P.Val < 0.05)
      } else {(which.sig <- which(limma_dataspace[,grep("adj.P.Val",colnames(limma_dataspace))] < 0.05))}
    }
  }

  if (parametric == FALSE) {
    if (independent == TRUE){
      tn <- "Kruskal"} else {tn <- "Friedman"}
    if (significance == "p"){
      if (groups_number != 2){
        which.sig <- which(Fdataspace$Kruskal_Wallis.pvalue < 0.05)
      }
    }
    if (significance == "p.adj"){
      if (groups_number != 2){
        which.sig <- which(Fdataspace$Kruskal_Wallis.pvalue_adjusted < 0.05)
      }
    }}
  group_colors <- hcl.colors(groups_number, palette = "Dark 3")
  names(group_colors) <- group_names
  if (length(which.sig) < 2){
    print("PCA and heatmap plots of the significant data cannot be generated since there are no significant proteins")
  }   else {
    log.dataspace.sig.m <- log.dataspace[which.sig,]
    zlog.dataspace.sig.m <- t(scale(t(log.dataspace.sig.m)))
    colnames(zlog.dataspace.sig.m) <- colnames(log.dataspace.sig.m)

    range_limit <- min(abs(min(zlog.dataspace.sig.m, na.rm = TRUE)), abs(max(zlog.dataspace.sig.m, na.rm = TRUE)))
    mycols <- circlize::colorRamp2(
      c(-range_limit, 0, range_limit),
      c("blue", "white", "red")
    )
    heatmap_data<- ComplexHeatmap::Heatmap(as.matrix(zlog.dataspace.sig.m),
                                           cluster_rows = hcluster,
                                           cluster_columns = FALSE,
                                           show_row_names = FALSE,
                                           show_column_names = FALSE,
                                           column_split = groups_list_f,
                                           top_annotation = ComplexHeatmap::HeatmapAnnotation(foo = anno_block(gp = gpar(fill = group_colors),
                                                                                                               labels = group_names, labels_gp = gpar(col = "white", fontface = "bold", fontsize = 14))),
                                           ,col = mycols, column_title = NULL,
                                           heatmap_legend_param = list(
                                             title = "Z-Score",
                                             color_bar = "continuous"
                                           ))

    bmp_file_path <- file.path(path_resplot, paste0(tn, "_significant_heatmap.bmp"))
    bmp(bmp_file_path,width = 1500, height = 1080, res = 150)
    ComplexHeatmap::draw(heatmap_data)
    dev.off()
    if (imputation == FALSE ) {
      log.dataspace.sig.m[is.na(log.dataspace.sig.m)] <- 0
    }

    pca<-prcomp(t(log.dataspace.sig.m), scale=TRUE, center=TRUE)

    pca.data <- data.frame(Sample=rownames(pca$x),
                           X=pca$x[,1],
                           Y=pca$x[,2],
                           Group = Group)

    pca.var<-pca$sdev^2
    pca.var.per<-round(pca.var/sum(pca.var)*100,1)

    pca.data$Group<-factor(pca.data$Group, levels=Group2)

    pca.sig<-ggplot2::ggplot(data=pca.data, aes(x=X, y=Y, label=Sample))+
      geom_point(aes(color=Group), size = 2, alpha = 1)+
      xlab(paste("PC1 - ", pca.var.per[1], "%", sep=""))+
      ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+
      ggtitle("Statistically significant proteins")+
      theme_bw()+
      theme(plot.title = element_text(hjust=0.5))+
      guides(color = guide_legend(override.aes = list(size=5)))+
      theme(legend.text = element_text(size = 12))+
      theme(legend.title = element_text(size = 12))+
      theme(legend.position="right")

    pca.sig

    ggplot2::ggsave(paste0(tn, "_PCA_plot_significant.bmp"), plot = pca.sig,  path = path_resplot,
                    scale = 1, width = 5, height = 4, units = "in",
                    dpi = 300, limitsize = TRUE)

  }

  which_sig <- list()
  volcano_list <- list()
  heatmap_cols <- list()
  if (parametric == TRUE) {
    if (significance  == "p"){
      for (i in 1:(ncol(mm)-1)) {
        for (j in (i+1):ncol(mm)) {
          comparison <- paste(colnames(mm)[i], "vs", colnames(mm)[j], sep = " ")
          which_sig[[comparison]] <- which(limma_dataspace[,grep(paste("P.Value", comparison, sep = "_"),colnames(limma_dataspace))] < 0.05)
          heatmap_cols[[comparison]] <- c(coln[[i]],coln[[j]])-3
          volcano.select <- c("Accession",paste0("logFC_", comparison) ,
                              paste0("P.Value_", comparison))
          volcano_list[[comparison]]<- limma_dataspace[,volcano.select ]
        }}
    }
    if (significance  == "p.adj"){
      significance_m = paste0(p.adjust.method," P")

      for (i in 1:(ncol(mm)-1)) {
        for (j in (i+1):ncol(mm)) {
          comparison <- paste(colnames(mm)[i], "vs", colnames(mm)[j], sep = " ")
          which_sig[[comparison]] <- which(limma_dataspace[,grep(paste("adj.P.Val", comparison, sep = "_"),colnames(limma_dataspace))] < 0.05)
          heatmap_cols[[comparison]] <- c(coln[[i]],coln[[j]])-3

          volcano.select <- c("Accession",paste0("logFC_", comparison) ,
                              paste0("adj.P.Val_", comparison))
          volcano_list[[comparison]]<- limma_dataspace[,volcano.select ]

        }}
    }
  }

  if (parametric == FALSE) {
    if (significance  == "p"){
      for (j in 2:groups_number) {
        for (k in 1:(j-1)) {
          comparison <- paste0("G", j, "vsG", k)
          which_sig[[comparison]] <- which(Ddataspace[,grep(paste0("Wilcoxon_p_", comparison),colnames(Ddataspace))] < 0.05)
          heatmap_cols[[comparison]] <- c(coln[[k]],coln[[j]])-3
          volcano.select <- c("Accession",paste0("Log2_Ratio_", comparison) ,
                              paste0("Wilcoxon_p_", comparison))
          volcano_list[[comparison]]<- Ddataspace[,volcano.select]

        }
      }
    }
    if (significance  == "p.adj"){
      significance_m = paste0(p.adjust.method," P")
      for (j in 2:groups_number) {
        for (k in 1:(j-1)) {
          comparison <- paste0("G", j, "vsG", k)
          which_sig[[comparison]] <- which(Ddataspace[,grep(paste0("p.adj_p_", comparison),colnames(Ddataspace))] < 0.05)
          heatmap_cols[[comparison]] <- c(coln[[k]],coln[[j]])-3

          volcano.select <- c("Accession",paste0("Log2_Ratio_", comparison) ,
                              paste0("p.adj_p_", comparison))
          volcano_list[[comparison]]<- Ddataspace[,volcano.select ]

        }
      }}
    for (i in 1:groups_number){
      names(which_sig)<- gsub(paste0("G",i), get(paste0("g",i,".name")),names(which_sig))
    }
  }

  for (i in 1:length(which_sig)) {
    if (length(which_sig[[i]]) <3 ){
      print(paste("No plots for the comparison ",names(which_sig[i]), ", will be created as there are not enough significant proteins"))
    } else {
      if (normalization %in% c(FALSE,"median","VSN", "Total_Ion_Current","PPM")){
        log.dataspace.sig <- log.dataspace[which_sig[[i]],]} else {
          log.dataspace.sig <-  dataspace[which_sig[[i]],-c(1:3)]
        }
      groups_list_h <- groups_list_f[heatmap_cols[[i]]]
      unique_groups_h <- unique(groups_list_h)
      block_colors_h <- group_colors[unique_groups_h]

      zlog.dataspace.sig <- t(scale(t(log.dataspace.sig)))
      colnames(zlog.dataspace.sig) <- colnames(log.dataspace.sig)
      zlog.dataspace.sig.h <- zlog.dataspace.sig[,heatmap_cols[[i]]]

      range_limit <- min(abs(min(zlog.dataspace.sig.h, na.rm = TRUE)), abs(max(zlog.dataspace.sig.h, na.rm = TRUE)))
      mycols <- circlize::colorRamp2(
        c(-range_limit, 0, range_limit),
        c("blue", "white", "red")
      )
      heatmap_data<- ComplexHeatmap::Heatmap(as.matrix(zlog.dataspace.sig.h),
                                             cluster_rows = hcluster,
                                             cluster_columns = FALSE,
                                             show_row_names = FALSE,
                                             show_column_names = FALSE,na_col = "grey",
                                             column_split = groups_list_h,
                                             top_annotation = ComplexHeatmap::HeatmapAnnotation(foo = anno_block(gp = gpar(fill =block_colors_h),
                                                                                                                 labels = unique_groups_h, labels_gp = gpar(col = "white", fontface = "bold", fontsize = 14))),
                                             ,col = mycols, column_title = NULL,
                                             heatmap_legend_param = list(
                                               title = "Z-Score",
                                               color_bar = "continuous"
                                             ))

      bmp_file_path <- file.path(path_resplot, paste0(names(which_sig[i]), "significant_heatmap.bmp"))
      bmp(bmp_file_path,width = 1500, height = 1080, res = 150)
      ComplexHeatmap::draw(heatmap_data)
      dev.off()


      volcano.P.vs.Con <- volcano_list[[i]]
      volcanocolors <- c("magenta", "navy", "#6D6D6D")
      max_y <- max(-log10(volcano.P.vs.Con[,3]), na.rm = TRUE) + 0.5
      names(volcanocolors) <- c("Significant upregulation","Significant downregulation", "Not significant")
      volcano.P.vs.Con$Regulation <- ifelse(volcano.P.vs.Con[,2]
                                            > LFC &
                                              volcano.P.vs.Con[,3]
                                            < 0.05, "Significant upregulation",
                                            ifelse(
                                              volcano.P.vs.Con[,2]
                                              < -LFC &
                                                volcano.P.vs.Con[,3]
                                              < 0.05,"Significant downregulation",
                                              "Not significant"))
      volcano.P.vs.Con$size <-
        ifelse(volcano.P.vs.Con$Regulation=="Not significant", 0.1,0.7)

      volcano.P.vs.Con <- volcano.P.vs.Con[!is.na(volcano.P.vs.Con$size), ]
      vol_plot <- ggplot(data = volcano.P.vs.Con, aes(x = volcano.P.vs.Con[,2],
                                                      y =
                                                        -log10(volcano.P.vs.Con[,3]), col = Regulation)) +
        geom_point(size= volcano.P.vs.Con$size) +
         scale_colour_manual(values = volcanocolors) +

        labs(title =  paste0(names(which_sig[i])),
             x =
               expression('Log'[2] * 'Fold Change')) +
        ylab(bquote('-Log'[10] * .(toupper(significance_m)) * ' value')) +
        guides(fill = guide_legend(title = NULL)) +
        theme_minimal() +
        theme(text = element_text(size = 14),
              plot.title = element_text(hjust = 0.5, face = "bold", size =
                                          16),
              axis.title.x = element_text(hjust = 0.5),
              legend.position = "none") +
        ylim(0, max_y) +
        xlim(-6, 6)

      suppressWarnings(ggplot2::ggsave(paste0(names(which_sig[i]),"_Volcano.bmp"), plot = vol_plot,  path = path_resplot,
                      scale = 1, width = 5, height = 4, units = "in",
                      dpi = 300, limitsize = TRUE))

    }}
  ##ENRICHMENT ANALYSIS
 # if (uqg == TRUE) stop("GSEA analysis unavailable for unique_genes matrices.")
  lfc_list <- list()
  dbgsea <- ifelse(species == "Mus musculus", "MM", "HS")
  if (subcollection == "H"){
    gene_sets <- msigdbr::msigdbr(species = species, db_species = dbgsea, collection = "H" )
  } else {
    gene_sets <- msigdbr::msigdbr(species = species, db_species = dbgsea, subcollection = subcollection)}
  gsea <- list()
  pathways <- split(gene_sets$gene_symbol, gene_sets$gs_name)

  for (j in 2:groups_number) {
    for (k in 1:(j-1)) {
      comparison <- paste0("G", j, "vsG", k)
      lfc_list[[comparison]] <- Ddataspace[,paste0("Log2_Ratio_",comparison)]
      names( lfc_list[[comparison]]) <- Ddataspace$Gene.Symbol
      valid_entries <- !is.na(names(lfc_list[[comparison]])) &
        names(lfc_list[[comparison]]) != "" &
        !is.na(lfc_list[[comparison]]) &
        is.finite(lfc_list[[comparison]])
      lfc_list[[comparison]] <- lfc_list[[comparison]][valid_entries]

      lfc_list[[comparison]] <- lfc_list[[comparison]][!duplicated(names(lfc_list[[comparison]]))]
      lfc_list[[comparison]] <- sort( lfc_list[[comparison]], decreasing = TRUE)
      gsea[[comparison]] <- suppressWarnings(fgsea(pathways = pathways, stats = lfc_list[[comparison]],
                                  eps = 0.0, minSize = 15, maxSize = 500))
    }}
  if ( nrow(gsea[[comparison]]) == 0 ) {
    warning("No pathways were found to be related with the features inserted. No GSEA will be conducted.")
  } else {

  for (i in 1:groups_number){
    names(gsea)<- gsub(paste0("G",i), get(paste0("g",i,".name")),names(gsea))
  }

  wb <- createWorkbook()

  for (name in names(gsea)) {
    addWorksheet(wb, name)
    writeData(wb, name, gsea[[name]])  # Write data to the sheet
  }
  gsea_file_path <- file.path(path_restat, "GSEA_results.xlsx")
  openxlsx::saveWorkbook(wb, file = gsea_file_path, overwrite =  TRUE)

  for (i in 1:length(gsea)){

    plot_data <- gsea[[i]] %>%
      filter(padj < 0.05) %>%
      arrange(NES)
    if (nrow(plot_data) >= 14) {
      plot_data <- plot_data %>% slice(c(1:7, (nrow(plot_data) - 6):nrow(plot_data)))
    }
    if (nrow(plot_data) == 0){
      warning(paste0("No gene sets came up significant (p.adjusted < 0.05) to the comparison ", names(gsea[i])))
    } else {
    if (subcollection == "CP:REACTOME") {
      plot_data <- plot_data %>% mutate(pathway = gsub("REACTOME_", "", pathway))
    }
      if (subcollection == "H") {
        plot_data <- plot_data %>% mutate(pathway = gsub("HALLMARK_", "", pathway))
      }
      if (subcollection == "GO:BP") {
        plot_data <- plot_data %>% mutate(pathway = gsub("GOBP_", "", pathway))
      }
    plot_data <-  plot_data %>% mutate( pathway = factor(pathway, levels = pathway))
    minNES <- floor(min(plot_data$NES))
    maxNES <- ceiling(max(plot_data$NES))

    comp_names <- strsplit(names(gsea[i]), "vs")[[1]]

    wrap_labels <- function(x) {
      sapply(x, function(label) {
        if (nchar(label) <= 28) return(label)
        lines <- character()
        remaining <- label
        while (nchar(remaining) > 28) {
          split_point <- 28 + regexpr("_", substr(remaining, 29, nchar(remaining)))[1]
          if (split_point < 28) split_point <- nchar(remaining)  # No _ found
          lines <- c(lines, substr(remaining, 1, split_point))
          remaining <- substr(remaining, split_point + 1, nchar(remaining))
        }
        if (nchar(remaining) > 0) lines <- c(lines, remaining)
        paste(lines, collapse = "\n")
      })
    }
    enr_plot <- ggplot(plot_data, aes(x = NES, y = pathway, size = size, color = padj)) +
      geom_point() + geom_segment(aes(xend = 0, yend = pathway), color = "black", size = 0.6, linetype = "dashed") +
      scale_color_gradient(low = "blue", high = "red", name = "P-adj",
                           limits = c(0, 0.05), breaks = seq(0, 0.05, by = 0.01)) +
      scale_size_continuous(name = "Size") +
      scale_x_continuous(limits = c(minNES, maxNES), breaks = seq(minNES, maxNES, by = 1)) +
      labs(x = "NES", y = "", title = paste0("GSEA Results for ", names(gsea[i]))) +
      theme_classic() +
      scale_y_discrete(labels = wrap_labels) + geom_vline(xintercept = 0, size = 0.2)+
      theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right"
      ) +
      labs(caption = paste("Enriched in", comp_names[2], "(Negative NES) and", comp_names[1], "(Positive NES)")) +
      theme(plot.caption = element_text(size = 8, face = "bold", hjust = 0.5))

  suppressWarnings(ggplot2::ggsave(paste0(names(gsea[i]),"_GSEA_plot.bmp"), plot = enr_plot,  path = path_resplot,
                    scale = 1, width = 8, height = 5, units = "in",
                    dpi = 300, limitsize = TRUE))

 } }
print("The ProtE analysis has been concluded.")
}
}

