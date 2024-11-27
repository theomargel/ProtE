#' MaxQuant proteomics data analysis
#'
#' It takes as input the Proteomics Data (output of Max Quant) in the format of an excel file that contains the information for each sample. Then it performs exploratory data analysis. The options for the data manipulation include different methods of normalization, and filtering based on the missing values per protein. Additionally, imputation of  the missing values can be performed, and a quality check with their percentage across every protein is provided. It then proceeds to perform statistical analysis using the Mann Whitney and the limma t-test for pairwise comparisons and  also Kruskal-Wallis and limma-ANOVA statistical tests,when there are more than 2 groups , while the pValues from the Levene and Bartlett statistical tests are also shown. The function also creates exploratory plots such as relative log espression boxplots and violin plots, heatmaps of the significant differentially expressed proteins and PCA plots.
#'
#' @param file The whole path to the ProteinGroups.txt file, that will be analyzed. Attention: Ensure to use forward slashes (/) for specifying paths.
#' @param group_names The names attributed to each different group. Insert in form of a vector. The order of the names should align with the order in the inserted ProteinGroups.txt file.
#' @param samples_per_group The number of samples attributed to each different group. Insert in form of a vector. The order of the number of groups should align with the order in the inserted ProteinGroups.txt file.
#' @param global_filtering TRUE/FALSE If TRUE the threshold for missing values filtering will be applied to the groups altogether, if FALSE it will be applied to each group separately.
#' @param imputation Imputation of the Missing Values. By default it is set to FALSE. Options are FALSE for no imputation implemented, "LOD" for assigning the lowest protein intensity identified to each MV and "LOD/2" to apply the half of it. Option "kNN" performs a default kNN imputation and "missRanger" a missRanger one. This 2 options are combined with a boxplot that visualizes the distribution of the log2 intensities of the imputed data compared to the initial ones.
#' @param sample_relationship Either "Independent" when the samples come from different populations or "Paired" when they come from the same. By default, it is set to "Independent". If "Paired" is selected the samples_per_group must be equal to each other
#' @param threshold_value The percentage of missing values per protein that will cause its omission. By default it is set to 50. (50 percent)
#' @param normalization The specific method for normalizing the data.By default it is set to FALSE. Options are FALSE for no normalization of the data, "log2" for a simple log2 transformation, "Quantile" for a quantiles based normalization, "median" for a median one, "TIC" for Total Ion Current normalization, "VSN" for Variance Stabilizing Normalization and "PPM" for Parts per Million transformation of the data.
#' @param parametric TRUE/FALSE Choose which statistical test will be taken into account when creating the optical statistical analysis (PCA plots, heatmap). By default it is set to FALSE (non Parametric)
#' @param significance pV or adj.pV Choose if the significant values for the PCA plots and the heatmap will derive from the unadjusted pValue or the adjusted pValue (Benjamini-Hochberg) of the comparison. By default it is set "pV" (pValue)
#'
#'
#' @return Excel files with the proteomic values that are optionally processed, via normalization, imputation and  filtering of proteins with a selected percentage of missing values. The result of the processing is visualized with an Protein Rank Abundance plot. PCA plots for all groups and for just their significant correlations are created. Furthermore violin and boxplots for the proteins of each sample is created and a heatmap for the significant proteins.
#' @importFrom openxlsx write.xlsx  read.xlsx
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom stringr str_trunc
#' @importFrom dplyr select  group_by  do everything  %>%
#' @importFrom tidyr gather pivot_longer
#' @importFrom broom tidy
#' @importFrom reshape2 melt
#' @importFrom ggpubr ggarrange
#' @importFrom ggplot2 ggplot ggsave geom_violin scale_color_gradient element_line theme_linedraw scale_fill_manual scale_color_manual aes geom_histogram element_rect geom_point xlab ylab ggtitle theme_bw theme_minimal theme element_text guides guide_legend geom_boxplot labs theme_classic element_blank geom_jitter position_jitter
#' @importFrom VIM kNN
#' @importFrom stats kruskal.test p.adjust prcomp sd wilcox.test model.matrix median na.omit
#' @importFrom forcats fct_inorder
#' @importFrom limma topTable eBayes contrasts.fit lmFit normalizeQuantiles duplicateCorrelation
#' @importFrom ComplexHeatmap HeatmapAnnotation anno_block draw Heatmap
#' @importFrom grid gpar
#' @importFrom car leveneTest
#' @importFrom vsn meanSdPlot
#' @importFrom missRanger missRanger
#' @importFrom utils read.delim
#'
#' @examples
#' #Example of running the function with paths for two groups.
#' # The file path is a placeholder, replace it with an actual file.
#' \dontrun{
#' maximum_quantum(
#'"C:/Users/User/Documents/ProteinGroups.txt",
#' groups_number = 2,
#' group_names = c("T0","T1"),
#' samples_per_group = c(3,3), threshold_value = 100)}
#'
#' @export

maximum_quantum <- function(file,
                            group_names,
                            samples_per_group,
                            imputation = FALSE,
                            global_filtering = TRUE,
                            sample_relationship = "Independent",
                            threshold_value = 50,
                            normalization = FALSE,
                      parametric= FALSE,
                      significance = "pV")
{
  group1 = group2 = Accession =Description =Symbol =X =Y = percentage=Sample= variable =.= g1.name =g2.name= g3.name =g4.name= g5.name =g6.name= g7.name= g8.name =g9.name =group3= group4= group5= group6 =group7= group8= group9 =key =value = NULL

message("The ProtE process starts now!")

  groups_number <- length(group_names)
  if (length(samples_per_group) != groups_number) {
    stop("The length of 'samples_per_group' must match 'groups_number'") }

  for (i in 1:groups_number) {
    assign(paste0("g",i,".name"),group_names[[i]])}

  dataspace <- read.delim(file, header = TRUE, sep = "\t")
  dataspace <- dataspace[!grepl("^;",dataspace$Protein.IDs),]
  dataspace <- dataspace[is.na(dataspace$Reverse), ]
  message("Removed REV_proteins: Reverse peptide Identifications")
  path <- dirname(file)
  path_res <- file.path(path , "MS_analysis")
  dir.create(path_res, showWarnings = FALSE)

  path_restat <- file.path(path_res, "Statistical_Analysis")
  path_resman <- file.path(path_res, "Data_processing")
  path_resplot <- file.path(path_res, "Plots")

  dir.create(path_restat, showWarnings = FALSE)
  dir.create(path_resman, showWarnings = FALSE)
  dir.create(path_resplot, showWarnings = FALSE)




  dataspace <- dataspace[,grep("Protein.IDs|Majority.protein.IDs|Intensity",colnames(dataspace))]
  dataspace$Protein.IDs <- ifelse(
    grepl("sp\\|\\w+\\|", dataspace$Protein.IDs),
    sub(".*sp\\|(\\w+)\\|.*", "\\1", dataspace$Protein.IDs),
    dataspace$Protein.IDs
  )
  dataspace$Majority.protein.IDs <- ifelse(
    grepl("sp\\|\\w+\\|\\w+", dataspace$Majority.protein.IDs),
    sub(".*sp\\|\\w+\\|([^;]+)","\\1", dataspace$Majority.protein.IDs),
    dataspace$Majority.protein.IDs)


  colnames(dataspace)[colnames(dataspace) == "Protein.IDs"] <- "Accession"
  colnames(dataspace)[colnames(dataspace) == "Majority.protein.IDs"] <- "Description"

  dataspace$Description[dataspace$Accession==dataspace$Description] = "Not available"

  colnames(dataspace) <- make.names(colnames(dataspace), unique = TRUE)
  rownames(dataspace) <- make.names(rownames(dataspace), unique = TRUE)


  dataspace <- dataspace[!dataspace$Intensity == 0,]
  dataspace$Intensity <- NULL
  colnames(dataspace) <- gsub("Intensity.", "", colnames(dataspace))

  if (sum(samples_per_group) != ncol(dataspace)-2) {stop("Error: Number of samples does not match the samples in the mastertable")}

  zero_per_sample <- colSums(is.na(dataspace[,-1:-2]))*100/nrow(dataspace)
  IDs <- colSums(!is.na(dataspace[,-1:-2]))


  if (normalization == "PPM"){
    dataspace[, -1:-2] <- lapply(dataspace[, -1:-2], function(x) {
      sum_x <- sum(x, na.rm = TRUE)
      ifelse(is.na(x), NA, (x / sum_x) * 10^6)
    })

    Gdataspace<-dataspace
    Gdataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Gdataspace$Description)
    Gdataspace$Symbol[Gdataspace$Symbol==Gdataspace$Description] = "Not available"
    Gdataspace<-Gdataspace %>%
      dplyr::select(Accession, Description, Symbol, everything())
    colnames(Gdataspace) <- gsub(".xlsx", "", colnames(Gdataspace))
    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(Gdataspace, file = norm_file_path)
    message("Applying the selected normalization, saved as Normalized.xlsx")}

  if (normalization == "Quantile"){
    dataspace[, -1:-2] <- log(dataspace[, -1:-2]+1,2)
    dataspace[, -1:-2] <- limma::normalizeQuantiles(dataspace[, -1:-2])
    Gdataspace<-dataspace
    Gdataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Gdataspace$Description)
    Gdataspace$Symbol[Gdataspace$Symbol==Gdataspace$Description] = "Not available"
    Gdataspace<-Gdataspace %>%
      dplyr::select(Accession, Description, Symbol, everything())
    colnames(Gdataspace) <- gsub(".xlsx", "", colnames(Gdataspace))
    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(Gdataspace, file = norm_file_path)
    message("Applying the selected normalization, saved as Normalized.xlsx") }
  if (normalization == "log2"){
    dataspace[, -1:-2] <- log(dataspace[, -1:-2]+1,2)
    Gdataspace<-dataspace
    Gdataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Gdataspace$Description)
    Gdataspace$Symbol[Gdataspace$Symbol==Gdataspace$Description] = "Not available"
    Gdataspace<-Gdataspace %>%
      dplyr::select(Accession, Description, Symbol, everything())
    colnames(Gdataspace) <- gsub(".xlsx", "", colnames(Gdataspace))
    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(Gdataspace, file = norm_file_path)
    message("Applying the selected normalization, saved as Normalized.xlsx") }
  if (normalization == "Total_Ion_Current") {
    dataspace[, -1:-2] <- lapply(dataspace[, -1:-2], function(x) (x / sum(x, na.rm = TRUE)) * mean(colSums(dataspace[, -1:-2], na.rm = TRUE)))
    Gdataspace<-dataspace
    Gdataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Gdataspace$Description)
    Gdataspace$Symbol[Gdataspace$Symbol==Gdataspace$Description] = "Not available"
    Gdataspace<-Gdataspace %>%
      dplyr::select(Accession, Description, Symbol, everything())
    colnames(Gdataspace) <- gsub(".xlsx", "", colnames(Gdataspace))
    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(Gdataspace, file = norm_file_path)
    message("Applying the selected normalization, saved as Normalized.xlsx")}
  if ( normalization == "VSN") {
    dataspace[, -1:-2] <- suppressMessages(limma::normalizeVSN(dataspace[, -1:-2]))
    Gdataspace<-dataspace
    Gdataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Gdataspace$Description)
    Gdataspace$Symbol[Gdataspace$Symbol==Gdataspace$Description] = "Not available"
    Gdataspace<-Gdataspace %>%
      dplyr::select(Accession, Description, Symbol, everything())
    colnames(Gdataspace) <- gsub(".xlsx", "", colnames(Gdataspace))
    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(Gdataspace, file = norm_file_path)
    message("Applying the selected normalization, saved as Normalized.xlsx")
  }
  if (normalization == "median") {
    sample_medians <- apply(dataspace[, -1:-2], 2, median, na.rm = TRUE)
    dataspace[, -1:-2] <- sweep(dataspace[, -1:-2], 2, sample_medians, FUN = "/")
    Gdataspace<-dataspace
    Gdataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Gdataspace$Description)
    Gdataspace$Symbol[Gdataspace$Symbol==Gdataspace$Description] = "Not available"
    Gdataspace<-Gdataspace %>%
      dplyr::select(Accession, Description, Symbol, everything())
    colnames(Gdataspace) <- gsub(".xlsx", "", colnames(Gdataspace))
    norm_file_path <- file.path(path_resman, "Normalized.xlsx")
    openxlsx::write.xlsx(Gdataspace, file = norm_file_path)
    message("Applying the selected normalization, saved as Normalized.xlsx")}

  if (normalization == FALSE ){
    dataspace <- dataspace
    sdrankplot_path <- file.path(path_resplot, "meanSdPlot.pdf")
    pdf(sdrankplot_path)
    suppressWarnings(vsn::meanSdPlot(as.matrix(dataspace[, -1:-2])))
    dev.off()
    message("Mean-SD plot of the data saved as meanSdPlot.pdf has been created.")
  } else {
    sdrankplot_path <- file.path(path_resplot, "meanSdPlot.pdf")
    pdf(sdrankplot_path)
    suppressWarnings(vsn::meanSdPlot(as.matrix(dataspace[, -1:-2])))
    dev.off()
    message("Mean-SD plot of the normalized data saved as meanSdPlot.pdf has been created.")
  }

  name_dataspace <-  dataspace[, -1:-2]


  if (threshold_value < 0 && threshold_value > 100) {stop("Error, you should add a threshold value number between 0 and 100")}

  if (global_filtering == TRUE) {

    threshold_value <- 100- as.numeric(threshold_value)

    threshold <-  ceiling(sum(samples_per_group)-(sum(samples_per_group)*(as.numeric(threshold_value)/100))+0.00000000001)
  }

  if (global_filtering == FALSE) {
    threshold<-numeric(groups_number)
    threshold_value <- 100- as.numeric(threshold_value)
    for (i in 1:groups_number) {
      threshold[i] <-  ceiling(samples_per_group[i]-(samples_per_group[i])*(as.numeric(threshold_value)/100)+0.00000000001)}
  }


  coln <- list()
  case_last <- 2

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
    at_file_path <- file.path(path_resman, "Dataset_filtering_applied.xlsx")
    openxlsx::write.xlsx(dataspace, file = at_file_path)
  }

  if (global_filtering == FALSE) {
    keep_rows <- rep(FALSE, nrow(dataspace))
    for (j in 1:groups_number) {
      keep_rows <- keep_rows | (dataspace[,paste0("Number_0_group", j)] < threshold[j])
    }
    dataspace <- dataspace[keep_rows, ]
    at_file_path <- file.path(path_resman, "Dataset_filtering_applied.xlsx")
    openxlsx::write.xlsx(dataspace, file = at_file_path)}
  message("An excel file with the proteins that have % of missing values at the selected threshold was created as Dataset_filtering_applied.xlsx")
  dataspace_0s<- dataspace
  dataspace[,paste0("Number_0_group", 1:groups_number)] <- NULL
  dataspace$Number_0_all_groups <- NULL



  zero_per_sample1 <- colSums(dataspace[,-1:-2] == 0)*100/nrow(dataspace)
  sample_names <- colnames(dataspace[,-1:-2])
  qc <- cbind(sample_names,IDs,zero_per_sample,zero_per_sample1)
  qc <- as.data.frame(qc)
  colnames(qc) <- c("Sample Name","Number of proteins detected in the sample","% of Missing values before filtering","% of Missing values after filtering")
  rownames(qc) <- NULL

  pre_dataspace <- dataspace

  if (imputation == "kNN") {
    message("kNN imputation starts now")
    dataspace[dataspace==0] <- NA
    dataspace[, -c(1, 2)] <- VIM::kNN(dataspace[, -c(1, 2)], imp_var = FALSE, k= 5)
    imp_file_path <- file.path(path_resman, "Dataset_Imputed.xlsx")
    openxlsx::write.xlsx(dataspace, file = imp_file_path)
  }
  if (imputation == "LOD"){
    dataspace[dataspace==0] <- NA
    impute_value <- min(as.matrix(dataspace[, -c(1, 2)]),na.rm = TRUE)
    dataspace[, -c(1, 2)][is.na(dataspace[, -c(1, 2)])]  <- impute_value
    imp_file_path <- file.path(path_resman, "Dataset_Imputed.xlsx")
    openxlsx::write.xlsx(dataspace, file = imp_file_path)  }
  if (imputation == "LOD/2"){
    dataspace[dataspace==0] <- NA
    impute_value <- min(as.matrix(dataspace[, -c(1, 2)]),na.rm = TRUE)/2
    dataspace[, -c(1, 2)][is.na(dataspace[, -c(1, 2)])]  <- impute_value
    imp_file_path <- file.path(path_resman, "Dataset_Imputed.xlsx")
    openxlsx::write.xlsx(dataspace, file = imp_file_path)  }
  if(imputation == "missRanger"){
    message("missRanger imputation starts now")
    dataspace[dataspace==0] <- NA
    dataspace[,-c(1,2)] <- missRanger::missRanger(dataspace[,-c(1,2)])
    imp_file_path <- file.path(path_resman, "Dataset_Imputed.xlsx")
    openxlsx::write.xlsx(dataspace, file = imp_file_path)  }
  if (imputation %in% c("kNN","missRanger"))    {
    pre_dataspace1<-pre_dataspace[,-1:-2]
    dataspace1<-dataspace[,-1:-2]
    imp.values<- dataspace1 - pre_dataspace1

    his_dataspace<-rbind(dataspace1,pre_dataspace1,imp.values)
    loghis_dataspace<-log2(his_dataspace+1)


    his_long <-tidyr::pivot_longer(loghis_dataspace, cols = everything())
    nrows<-nrow(his_long)
    his_long$Group <- rep(c("Final","Initial","Imputed"), each = (nrows/3))
    his_long_filtered <- his_long[his_long$value != 0,]
    his_long_filtered$Group <- factor(his_long_filtered$Group, levels = c("Final", "Initial", "Imputed"))

    imp_hist<- ggplot(his_long_filtered, aes(x = value, fill = Group, colour = Group)) +
      labs( x = expression(Log[2]~"Proteins Abundance"), y = "Count") +
      scale_fill_manual(values = c("Final" = "#FF99FF", "Initial" = "#990000", "Imputed" = "#000033")) +
      scale_color_manual(values = c("Final" = "#FF99FF", "Initial" = "#990000", "Imputed" = "#000033")) +
      geom_histogram(alpha = 0.5, binwidth = 0.3, position = "identity") +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white")) + theme(panel.background = element_rect(fill = "white"))

    imp_hist

    ggplot2::ggsave("Imputed_values_histogram.pdf", plot = imp_hist,  path = path_resplot,
                    scale = 1, width = 5, height = 4, units = "in",
                    dpi = 300, limitsize = TRUE)
    message("A plot named Imputed_values_histogram.pdf, showcasing the distribution of the imputed values was created.")
  }
  if (imputation %in% c("LOD/2","LOD","kNN","missRanger")){
    message("An excel with the imputed missing values was created as Dataset_Imputed.xlsx")

    dataspace_0s$percentage <- dataspace_0s$Number_0_all_groups*100/sum(samples_per_group)
    dataspace$percentage <- dataspace_0s$percentage
    dataspace$mean <- rowMeans(dataspace[,3:(2+sum(samples_per_group))])
    dataspace$log<-log2(dataspace$mean)
    dataspace$rank <- rank(-dataspace$mean)

    abund.plot <- ggplot(dataspace, aes(x = rank, y = log, colour = percentage)) +
      geom_point(size = 3, alpha = 0.8) +
      labs(title = "Protein Abundance Rank", x = "Rank", y = expression(Log[2] ~ "Proteins Abundance")) +
      scale_color_gradient(low = "darkblue", high = "yellow",
                           name = "Imputations\nin each\nprotein\n(%)") +
      theme_linedraw()+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid = element_line(color = "grey80"),
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 9))
    abund.plot

    ggplot2::ggsave("Proteins_abundance_rank.pdf", plot = abund.plot ,  path = path_resplot,
                    scale = 1, width = 12, height = 5, units = "in",
                    dpi = 300, limitsize = TRUE, bg = "white")
  }

  if (imputation == FALSE){dataspace <- dataspace

  dataspace_0s$percentage <- dataspace_0s$Number_0_all_groups*100/sum(samples_per_group)
  dataspace$percentage <- dataspace_0s$percentage
  dataspace$mean <- apply(dataspace[, 3:(2+sum(samples_per_group))], 1, function(x) mean(x[x != 0]))
  dataspace$log<-log2(dataspace$mean)
  dataspace$rank <- rank(-dataspace$mean)

  abund.plot <- ggplot(dataspace, aes(x = rank, y = log, colour = percentage)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = "Protein Abundance Rank", x = "Rank", y = expression(Log[2] ~ "Proteins Abundance")) +
    scale_color_gradient(low = "darkblue", high = "yellow",
                         name = "MVs\nin each\nprotein\n(%)") +
    theme_linedraw()+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid = element_line(color = "grey80"),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9))

  abund.plot

  ggplot2::ggsave("Proteins_abundance_rank.pdf", plot = abund.plot ,  path = path_resplot,
                  scale = 1, width = 12, height = 5, units = "in",
                  dpi = 300, limitsize = TRUE, bg = "white")
  }
  message("The plot named Proteins_abundance_rank.pdf, which depicts the proteins abundance rank and their percentage of missing values was created.")
  dataspace$percentage <- NULL
  dataspace$mean <- NULL
  dataspace$log<- NULL
  dataspace$rank <- NULL

  groups_list <- list("character")
  for (i in 1:groups_number) {
    groups_list[[i]] <- rep(group_names[i], times = samples_per_group[i])
  }

  groups_list_u <- unlist(groups_list)
  groups_list_f <- factor(groups_list_u, levels = unique(groups_list_u))

  mm <- model.matrix(~groups_list_f + 0)
  colnames(mm)<- group_names

  nndataspace<- dataspace[,-1:-2]
  nndataspace[is.na(nndataspace)] <- 0
  nndataspace <- log2(nndataspace+1)

  if (sample_relationship == "Paired"){
    if (length(unique(samples_per_group)) != 1){
      message("Error: Paired comparison did not happen correctly because you have input different number of samples in some groups")
    }
    n = sum(samples_per_group)/groups_number
    pairing <- rep(1:n, each = groups_number)

    corfit <- limma::duplicateCorrelation(nndataspace, design = mm, block = pairing)
    fit <- limma::lmFit(nndataspace, mm, block = pairing, correlation = corfit$consensus.correlation)
  }
  if (sample_relationship == "Independent"){
    fit <- limma::lmFit(nndataspace, mm)}
  fit<- limma::eBayes(fit)
  if (groups_number>2){
    anova_res<- limma::topTable(fit, adjust.method = "BH", number = Inf)
    colnames(anova_res)<-paste("ANOVA",colnames(anova_res), sep = "_")
    anova_res<- anova_res[,-c(1:groups_number)]}
  lima.res <- data.frame()
  for (i in 1:(ncol(mm)-1)) {
    for (j in (i+1):ncol(mm)) {
      comparison <- paste(colnames(mm)[i], "vs", colnames(mm)[j], sep = " ")
      contrast_fref <- limma::makeContrasts(contrasts = paste0(colnames(mm)[i],"-",colnames(mm)[j]), levels = mm)
      fit2 <- limma::contrasts.fit(fit, contrast_fref)
      fit2 <- limma::eBayes(fit2)
      top_table<- limma::topTable(fit2, adjust.method = "BH", number = Inf)
      column_groups<- top_table[,c("logFC","AveExpr","t","P.Value","adj.P.Val","B")]
      colnames(column_groups)<-paste(colnames(column_groups), comparison, sep = "_")

      if (nrow(lima.res) == 0){
        lima.res <- column_groups
      } else {lima.res<- cbind(lima.res, column_groups)}
    }}

  if (groups_number>2){
    limma_dataspace <- cbind(anova_res,lima.res,dataspace)} else {limma_dataspace <- cbind(lima.res,dataspace)}

  ncollimma <- ncol(limma_dataspace) - ncol(dataspace) + 2
  limma_dataspace<-limma_dataspace %>%
    dplyr::select(Accession, any_of(c("Protein.Names","Description")), dplyr::everything())
  limma_dataspace <- limma_dataspace[,1:ncollimma]
  limma_file_path <- file.path(path_restat, "Dataset_limma.test.xlsx")
  openxlsx::write.xlsx(limma_dataspace, file = limma_file_path)
  message("The statistics from the limma parametric tests are showcased in the created Dataset_limma.test.xlsx file.")

  if (sample_relationship != "Paired" && sample_relationship != "Independent"){stop("Error. You need to assign sample_relationship = 'Paired' or 'Indepedent'")}
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

          if (sample_relationship == "Independent") {
            if (sum(!is.na(values_k)) > 0 & sum(!is.na(values_j)) > 0) {
              test_list <- stats::wilcox.test(
                values_k, values_j,
                exact = FALSE, paired = FALSE, na.rm = TRUE
              )
              data2[i, paste0("MW_G", j, "vsG", k)] <- test_list$p.value
            } else {
              data2[i, paste0("MW_G", j, "vsG", k)] <- NA
            }
          } else if (sample_relationship == "Paired") {
            paired_values <- na.omit(cbind(values_k, values_j))
            if (nrow(paired_values) > 0) {
              test_list <- stats::wilcox.test(
                paired_values[, 1], paired_values[, 2],
                exact = FALSE, paired = TRUE
              )
              data2[i, paste0("MW_G", j, "vsG", k)] <- test_list$p.value
            } else {
              data2[i, paste0("MW_G", j, "vsG", k)] <- NA
            }
          }
        }

        data2[[paste0("BH_p_G", j, "vsG", k)]] <- p.adjust(
          data2[[paste0("MW_G", j, "vsG", k)]],
          method = "BH"
        )

        avg_j <- data2[[paste0("Average_G", j)]]
        avg_k <- data2[[paste0("Average_G", k)]]
        data2[[paste0("Ratio_G", j, "vsG", k)]] <- ifelse(
          avg_k == 0, NA, avg_j / avg_k
        )
        data2[[paste0("Log2_Ratio_G", j, "vsG", k)]] <- log2(
          data2[[paste0("Ratio_G", j, "vsG", k)]]+1
        )
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

      levene_result <- car::leveneTest(Value ~ as.factor(Group), data = leve_data, center = median)
      data2[i, "Levene_p"] <- levene_result$`Pr(>F)`[1]      } else {
        data2[i, "Bartlett_p"] <- NA
        data2[i, "Levene_p"] <- NA
      }

  }



  Ddataspace<-data2
  Ddataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Ddataspace$Description)
  Ddataspace$Symbol[Ddataspace$Symbol==Ddataspace$Description] = "Not available"
  Ddataspace<-Ddataspace %>%
    dplyr::select(Accession, Description, Symbol, everything())
  Fdataspace<-Ddataspace


  Group<-list()
  times<-vector()
  for (i in (1:groups_number)){
    times<-samples_per_group[i]
    Group[[i]] <- rep(paste0("G",i), times)
    times<-NULL}

  Group<-unlist(Group)
  dataspace3<-t(dataspace[,-c(1:2)])
  dataspace3<-data.frame(dataspace3)
  dataspace<-data.frame(dataspace)
  colnames(dataspace3)<-dataspace[,1]
  dataspace4<-cbind(Group,dataspace3)
  dataspace4$Group<-as.factor(dataspace4$Group)

  if (groups_number>2){
    message("Mann-Whitney, Levene, Bartlett tests done, now calculating the Kruskal-Wallis test's results:")

    df3 <- dataspace4 %>% tidyr::gather(key, value, -Group)
    df4 <- df3 %>% dplyr::group_by(key)
    df4$value<-as.numeric(df4$value)
    df5 <- df4 %>% dplyr::do(broom::tidy(kruskal.test(x= .$value, g = .$Group)))
    Kruskal_Wallis.pvalue <- df5$p.value
    data3<-cbind(Ddataspace,Kruskal_Wallis.pvalue)
    data3$Kruskal_Wallis.pvalue_BH.adjusted<- p.adjust(data3$Kruskal_Wallis.pvalue, method = "BH")

    Fdataspace<-data3
    Fdataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Fdataspace$Description)
    Fdataspace$Symbol[Fdataspace$Symbol==Fdataspace$Description] = "Not available"
    Fdataspace<-Fdataspace %>%
      dplyr::select(Accession, Description, Symbol, everything())

  }
  for (i in 1:groups_number){
    namesc<- colnames(Fdataspace)
    namesc<- gsub(paste0("G",i), get(paste0("g",i,".name")),namesc)
    colnames(Fdataspace)<-namesc
  }

  colnames(Fdataspace) <- gsub(".xlsx", "", colnames(Fdataspace))
  start_col <- 3 + as.numeric(sum(samples_per_group))
  Fdataspace <- Fdataspace[,-c(4:start_col)]


  stats_file_path <- file.path(path_restat, "Statistics.xlsx")
  openxlsx::write.xlsx(Fdataspace, file = stats_file_path)
  message("An excel with the statistical tests for the normalized data was created as Statistical_analysis.xlsx")



  Group <- groups_list_f

  Group2<-unique(groups_list_f)

  log.dataspace <- log(dataspace[,-c(1:2)]+1,2)


  pca<-prcomp(t(log.dataspace), scale=TRUE, center=TRUE)

  pca.data <- data.frame(Sample=rownames(pca$x),
                         X=pca$x[,1],
                         Y=pca$x[,2],
                         Group = Group)
  qc$PC1.score <- pca$x[,1]
  qc$PC2.score <-pca$x[,2]
  if (groups_number == 2){
    Group<-list()
    times<-vector()
    for (i in (1:length(samples_per_group))){
      times<-samples_per_group[i]
      Group[[i]] <- rep(paste0("G",i), times)
      times<-NULL
    }
    Group<-unlist(Group)
    Group<-gsub("G1", g1.name, Group)
    Group<-gsub("G2", g2.name, Group)
  }
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

  ggplot2::ggsave("PCA_plot_alldata.pdf", plot = pca.ent,  path = path_resplot,
                  scale = 1, width = 5, height = 4, units = "in",
                  dpi = 300, limitsize = TRUE)
  message("PCA plot using all data was created as PCA_plot_alldata.pdf")


  which.sig<-vector()
  if (parametric == TRUE) {
    if (significance == "pV"){
      if (groups_number != 2){
        which.sig <- which(limma_dataspace$ANOVA_P.Value < 0.05)
      } else {(which.sig <- which(limma_dataspace[,grep("P.Value",colnames(limma_dataspace))] < 0.05))}
    }
    if (significance == "adj.pV"){
      if (groups_number != 2){
        which.sig <- which(limma_dataspace$ANOVA_adj.P.Val < 0.05)
      } else {(which.sig <- which(limma_dataspace[,grep("adj.P.Val",colnames(limma_dataspace))] < 0.05))}
    }
  }

  if (parametric == FALSE) {
    if (significance == "pV"){
      if (groups_number != 2){
        which.sig <- which(Fdataspace$Kruskal_Wallis.pvalue < 0.05)
      } else {(which.sig <- which(Ddataspace$MW_G2vsG1 < 0.05))}
    }
    if (significance == "adj.pV"){
      if (groups_number != 2){
        which.sig <- which(Fdataspace$Kruskal_Wallis.pvalue_BH.adjusted < 0.05)
      } else {(which.sig <- which(Ddataspace$BH_p_G2vsG1 < 0.05))}
    }}

  if (length(which.sig) == 0){
    message("There are no significant proteins, to create a PCA plot with them and a heatmap")
    qc[,-1] <- lapply(qc[,-1], function(x) as.numeric(unlist(x)))
    qc[,-1]<-round(qc[,-1],3)

    qc_file_path <- file.path(path_restat, "Quality_check.xlsx")
    openxlsx::write.xlsx(qc, file = qc_file_path)    }   else {
      log.dataspace.sig <- log.dataspace[which.sig,]


      zlog.dataspace.sig <- t(scale(t(log.dataspace.sig)))
      colnames(zlog.dataspace.sig) <- colnames(log.dataspace.sig)

      mycols <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
      heatmap_data<- ComplexHeatmap::Heatmap(as.matrix(zlog.dataspace.sig),
                                             cluster_rows = TRUE,
                                             cluster_columns = FALSE,
                                             show_row_names = FALSE,
                                             show_column_names = FALSE,
                                             column_split = groups_list_f,
                                             top_annotation = ComplexHeatmap::HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:(groups_number+1)),
                                                                                                                 labels = group_names, labels_gp = gpar(col = "black", fontsize = 10))),
                                             col = mycols, column_title = NULL,
                                             heatmap_legend_param = list(
                                               title = "Z-Score",
                                               color_bar = "continuous"
                                             ))
      pdf_file_path <- file.path(path_resplot, "heatmap.pdf")
      pdf(pdf_file_path, width = 7.37, height = 6.09)
      ComplexHeatmap::draw(heatmap_data)
      dev.off()

      message("A heatmap with the significantly differentially expressed proteins was created as heatmap.pdf")


      pca<-prcomp(t(log.dataspace.sig), scale=TRUE, center=TRUE)

      pca.data <- data.frame(Sample=rownames(pca$x),
                             X=pca$x[,1],
                             Y=pca$x[,2],
                             Group = Group)
      qc$PC1.score.Significant <- pca$x[,1]
      qc$PC2.score.Significant <-pca$x[,2]
      qc[,-1] <- lapply(qc[,-1], function(x) as.numeric(unlist(x)))
      qc[,-1]<-round(qc[,-1],3)
      qc_file_path <- file.path(path_restat, "Quality_check.xlsx")
      openxlsx::write.xlsx(qc, file = qc_file_path)

      message("An excel file named Quality_check.xlsx, that provides information on the missing values and the Principal Component score for each sample was created")
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

      ggplot2::ggsave("PCA_plot_significant.pdf", plot = pca.sig,  path = path_resplot,
                      scale = 1, width = 5, height = 4, units = "in",
                      dpi = 300, limitsize = TRUE)
      message("PCA plot with the significant data was created as PCA_plot_significant.pdf" )
      a<-ggpubr::ggarrange(pca.ent, pca.sig, nrow = 1, ncol=2,
                           common.legend = TRUE, legend = "bottom")

      ggplot2::ggsave("PCA_plots_combined.pdf", plot = a,  path = path_resplot,
                      scale = 1, width = 8, height = 4.5, units = "in",
                      dpi = 300, limitsize = TRUE)
      message ("The 2 PCA plots are combined in PCA_plots_combined.pdf")
    }
  p<- function(x) {
    sapply(x, function(label) {
      truncated_label <- substr(label, nchar(label) - 24, nchar(label))
      truncated_label
    })
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

  if (imputation == FALSE) {
    qc.boxplots<-ggplot2::ggplot(melt.log.dataspace, aes(x=forcats::fct_inorder(variable), y=value, color=Group))+
      geom_boxplot(aes(color = Group),lwd=1, outlier.size=0.2, outlier.alpha = 0.2)+
      xlab("Sample")+
      ylab(expression(Log[2]~"Protein Abundance"))+
      theme_classic()+
      theme(text = element_text(size = 19),
            axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 0.5),
            axis.title.x=element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"))+
      guides(color = guide_legend(override.aes = list(size = 1)))+
      scale_x_discrete(labels = p) +
      geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5, alpha = 0.5)

    qc.boxplots

    ggplot2::ggsave("Boxplot_withZeros.pdf", plot = qc.boxplots,  path = path_resplot,
                    scale = 1, width = 12, height = 5, units = "in",
                    dpi = 300, limitsize = TRUE, bg = "white")

  }
  melt.log.dataspace.na <- melt.log.dataspace
  melt.log.dataspace.na$value[melt.log.dataspace.na$value == 0] <- NA
  melt.log.dataspace.na$Group <- factor(melt.log.dataspace.na$Group, levels = Group2)
  is.factor(melt.log.dataspace.na$variable)

  qc.boxplots.na<-ggplot2::ggplot(melt.log.dataspace.na, aes(x=forcats::fct_inorder(variable), y=value, color=Group))+
    geom_boxplot(aes(color = Group),lwd=1, outlier.size=0.2, outlier.alpha = 0.2)+
    xlab("Sample")+
    ylab(expression(Log[2]~"Protein Abundance"))+
    theme_classic()+
    theme(text = element_text(size = 19),
          axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 0.5),
          axis.title.x=element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))+
    guides(color = guide_legend(override.aes = list(size = 1)))+
    scale_x_discrete(labels = p) +
    geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5, alpha = 0.5)

  qc.boxplots.na
  if (imputation == FALSE) {
    ggplot2::ggsave("Boxplot_withoutZeros.pdf", plot = qc.boxplots.na, path = path_resplot,
                    scale = 1, width = 12, height = 5, units = "in",
                    dpi = 300, limitsize = TRUE, bg = "white")
  } else
  {
    ggplot2::ggsave("Boxplot.pdf", plot = qc.boxplots.na, path = path_resplot,
                    scale = 1, width = 12, height = 5, units = "in",
                    dpi = 300, limitsize = TRUE, bg = "white")
  }
  message("A boxplot showing the ", expression(Log[2]~"Abundance")," of each protein, across the samples was created as Boxplot.pdf" )

  qc.violin<-ggplot2::ggplot(melt.log.dataspace.na, aes(x=forcats::fct_inorder(variable), y=value, color=Group))+
    geom_violin(aes(color = Group),lwd=1)+
    xlab("Sample")+
    ylab(expression(Log[2]~"Protein Abundance"))+
    theme_classic()+
    theme(text = element_text(size = 19),
          axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 0.5),
          axis.title.x=element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))+
    guides(color = guide_legend(override.aes = list(size = 1)))+
    scale_x_discrete(labels = p) +
    geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5, alpha = 0.5)

  ggplot2::ggsave("Violin_plot.pdf", plot = qc.violin,  path = path_resplot,
                  scale = 1, width = 12, height = 5, units = "in",
                  dpi = 300, limitsize = TRUE, bg = "white")
  message("A Violin Plot showing the ", expression(Log[2]~"Abundance")," of each protein, across the samples was created as Violin_plot.pdf" )
  message("The analysis was created. The results are saved inside the MS_analysis folder. Thanks for your patience!")

}
