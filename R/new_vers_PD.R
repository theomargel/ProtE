#' Proteome Discoverer analysis
#'
#' It takes Proteomics Data from samples in different groups, in the format they are created by Proteome Discoverer (PD). It concatenates the Accession IDs, Descriptions and Areas from different PD export files into a master table and performs exploratory data analysis. The function outputs a normalized Parts Per Million protein dataset along with descriptive statistics and results of significance testing. The script also creates exploratory plots such as relative log espression boxplots and PCA plots.
#'
#' @param excel_file The whole path to the excel .xlsx file, that will be analysed. Attention: Add '/' between the directories.
#' @param group_names The names attributed to each different group. Insert in form of a vector. The order of the names should align with the order in the inserted excel file.
#' @param case_number The number of samples attributed to each different group. Insert in form of a vector. The order of the number of groups should align with the order in the inserted excel file.
#' @param global_threshold TRUE/FALSE If TRUE threshold for missing values will be applied to the groups altogether, if FALSE to each group seperately
#' @param imputation TRUE/FALSE Data imputation using kNN classification or assigning missing values as 0.
#' @param MWtest Either "Paired" for a Wilcoxon Signed-rank test or "Independent" for a Mann-Whitney U test.
#' @param threshold_value The percentage of missing values per protein that will cause its deletion
#' @param normalization PPM, TIC etc.

#'
#' @return Excel files with the proteomic values from all samples, processed with normalization and imputation and substraction of samples with high number of missing values. PCA plots for all or for just the significant correlations, and boxplots for the proteins of each sample.
#' @importFrom openxlsx write.xlsx  read.xlsx
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom dplyr select  group_by  do everything  %>%
#' @importFrom tidyr gather pivot_longer
#' @importFrom broom tidy
#' @importFrom reshape2 melt
#' @importFrom ggpubr ggarrange
#' @importFrom ggplot2 ggplot ggsave geom_violin scale_color_gradient element_line theme_linedraw scale_fill_manual scale_color_manual aes geom_histogram element_rect geom_point xlab ylab ggtitle theme_bw theme_minimal theme element_text guides guide_legend geom_boxplot labs theme_classic element_blank geom_jitter position_jitter
#' @importFrom VIM kNN
#' @importFrom stringr str_trunc
#' @importFrom stats kruskal.test p.adjust prcomp sd wilcox.test model.matrix heatmap median
#' @importFrom forcats fct_inorder
#' @importFrom limma topTable eBayes contrasts.fit lmFit normalizeQuantiles
#' @importFrom ComplexHeatmap HeatmapAnnotation anno_block draw Heatmap
#' @importFrom grid gpar
#' @importFrom missRanger missRanger
#'
#' @examples #' # Example of running the function with paths for two groups.
#' #Do not add if (interactive()){} condition in your code
#' if (interactive()){
#' new_user_inputs(excel_file = "C:/Users/User/Documents/itern/new_version_PD.xlsx",
#' groups_number = 2,group_names = c("exp_2020","exp_2024"),
#'  case_number = c(18,8), imputation = FALSE, threshold_value = 100)}
#'
#' @export

new_user_inputs <- function(excel_file,
                            group_names,
                            case_number,
                        imputation = TRUE,
                        global_threshold = TRUE,
                        MWtest = "Independent",
                        threshold_value = 50
                        ,normalization = FALSE)
{
  group1 = group2 = Accession =Description =Symbol =X =Y = percentage=Sample= variable =.= g1.name =g2.name= g3.name =g4.name= g5.name =g6.name= g7.name= g8.name =g9.name =group3= group4= group5= group6 =group7= group8= group9 =key =value = NULL

  groups_number <- length(group_names)
  if (length(case_number) != groups_number) {
    stop("The length of 'case_number' must match 'groups_number'") }
#  if (groups_number>9){stop("You can add up to 9 groups")}

  for (i in 1:groups_number) {
    assign(paste0("g",i,".name"),group_names[[i]])}

 dataspace <- openxlsx::read.xlsx(excel_file)

   path <- dirname(excel_file)
  path_res <- file.path(path , "MS_analysis_new")
  dir.create(path_res, showWarnings = FALSE)

  #modify PD masterlist to exclude unwanted stats
dataspace <- dataspace[,grep("Accession|Description|Abundance:",colnames(dataspace))]

colnames(dataspace) <- gsub("Abundance:.|:.Sample", "", colnames(dataspace))
dataspace <- dataspace[rowSums(!is.na(dataspace[,-c(1,2)])) > 0, ]

zero_per_sample <- colSums(is.na(dataspace[,-1:-2]))*100/nrow(dataspace)
IDs <- colSums(!is.na(dataspace[,-1:-2]))

if (normalization == FALSE ){
  dataspace <- dataspace
}
if (normalization == "PPM"){
  #normalize PPm
  dataspace[, -1:-2] <- lapply(dataspace[, -1:-2], function(x) {
    sum_x <- sum(x, na.rm = TRUE)  # Sum of the column, ignoring NAs
    ifelse(is.na(x), NA, (x / sum_x) * 10^6)  # NA=0 , normalize the rest
  })}

if (normalization == "Quantile"){
  dataspace[, -1:-2] <- log(dataspace[, -1:-2]+1,2)
  dataspace[, -1:-2] <- limma::normalizeQuantiles(dataspace[, -1:-2])
}

if (normalization == "log2"){
  dataspace[, -1:-2] <- log(dataspace[, -1:-2]+1,2)
}
if (normalization == "Total_Ion_Current") {
  dataspace[, -1:-2] <- lapply(dataspace[, -1:-2], function(x) (x / sum(x, na.rm = TRUE)) * mean(colSums(dataspace[, -1:-2], na.rm = TRUE)))
}
if (normalization == "median") {
  sample_medians <- apply(dataspace[, -1:-2], 2, median, na.rm = TRUE)
  dataspace[, -1:-2] <- sweep(dataspace[, -1:-2], 2, sample_medians, FUN = "/")
}
name_dataspace <-  dataspace[, -1:-2]
dat.dataspace<-dataspace


if (threshold_value < 0 && threshold_value > 100) {stop("Error, you should add a threshold value number between 0 and 100")}

if (global_threshold == TRUE) {

  threshold_value <- 100- as.numeric(threshold_value)
  # Iterate over each group and calculate the new threshold

  threshold <-  ceiling(sum(case_number)-(sum(case_number)*(as.numeric(threshold_value)/100))+0.00000000001)
}

if (global_threshold == FALSE) {
  threshold<-numeric(groups_number)
  threshold_value <- 100- as.numeric(threshold_value)
  # Iterate over each group and calculate the new threshold
  for (i in 1:groups_number) {
    threshold[i] <-  ceiling(case_number[i]-(case_number[i])*(as.numeric(threshold_value)/100)+0.00000000001)}
}



# Create identifier variables for the thhreshold and statistics
coln <- list()
case_last <- 2

for (i in 1:groups_number) {
  case_last <- case_last + case_number[i]
  coln[[i]] <- (case_last - case_number[i] + 1):case_last
}


setwd(path_res)
Gdataspace<-dataspace
Gdataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Gdataspace$Description)
Gdataspace$Symbol[Gdataspace$Symbol==Gdataspace$Description] = "Not available"
Gdataspace<-Gdataspace %>%
  dplyr::select(Accession, Description, Symbol, everything())
colnames(Gdataspace) <- gsub(".xlsx", "", colnames(Gdataspace))
openxlsx::write.xlsx(Gdataspace, file = "Normalized.xlsx")
message("Applying Parts Per Million normalization, saved as Normalized.xlsx")

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
#omiting empty rows
dataspace <- dataspace[dataspace$Number_0_all_groups < sum(case_number),]

if (global_threshold == TRUE) {
  openxlsx::write.xlsx(dataspace,file = "Dataset_before_threshold.xlsx")
  dataspace <- dataspace[dataspace$Number_0_all_groups<threshold,]
  openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")
}

if (global_threshold == FALSE) {
  openxlsx::write.xlsx(dataspace,file = "Dataset_before_threshold.xlsx")
  keep_rows <- rep(FALSE, nrow(dataspace))
  for (j in 1:groups_number) {
    keep_rows <- keep_rows | (dataspace[,paste0("Number_0_group", j)] < threshold[j])
  }
  dataspace <- dataspace[keep_rows, ]
  openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")}
message("An excel file with the proteins that have % of missing values below the threshold was created as Dataset_threshold_applied.xlsx")
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


##imputation KNN
if (imputation == "kNN") {
  dataspace[dataspace==0] <- NA
  dataspace[, -c(1, 2)] <- VIM::kNN(dataspace[, -c(1, 2)], imp_var = FALSE, k= 5)
  openxlsx::write.xlsx(dataspace,file = "Dataset_Imputed.xlsx")
}
if (imputation == "LOD"){
  dataspace[dataspace==0] <- NA
  impute_value <- min(as.matrix(dataspace[, -c(1, 2)]),na.rm = TRUE)
  dataspace[, -c(1, 2)][is.na(dataspace[, -c(1, 2)])]  <- impute_value
  openxlsx::write.xlsx(dataspace,file = "Dataset_Imputed.xlsx")
}
if (imputation == "LOD/2"){
  dataspace[dataspace==0] <- NA
  impute_value <- min(as.matrix(dataspace[, -c(1, 2)]),na.rm = TRUE)/2
  dataspace[, -c(1, 2)][is.na(dataspace[, -c(1, 2)])]  <- impute_value
  openxlsx::write.xlsx(dataspace,file = "Dataset_Imputed.xlsx")
}
if(imputation == "missRanger"){
  dataspace[dataspace==0] <- NA
  dataspace[,-c(1,2)] <- missRanger::missRanger(dataspace[,-c(1,2)])
  openxlsx::write.xlsx(dataspace,file = "Dataset_Imputed.xlsx")
}
if (imputation %in% c("kNN","missRanger"))    {
    pre_dataspace1<-pre_dataspace[,-1:-2]
    dataspace1<-dataspace[,-1:-2]
    imp.values<- dataspace1 - pre_dataspace1

    his_dataspace<-rbind(dataspace1,pre_dataspace1,imp.values)
    loghis_dataspace<-log2(his_dataspace+1)

    #his_long <- reshape2::melt(loghis_dataspace)

    his_long <-tidyr::pivot_longer(loghis_dataspace, cols = everything())
    nrows<-nrow(his_long)
    his_long$Group <- rep(c("Final","Initial","Imputed"), each = (nrows/3))
    his_long_filtered <- his_long[his_long$value != 0,]
    his_long_filtered$Group <- factor(his_long_filtered$Group, levels = c("Final", "Initial", "Imputed"))

    imp_hist<- ggplot(his_long_filtered, aes(x = value, fill = Group, colour = Group)) +
      labs( x = expression(Log[2]~"Parts per Million"), y = "Count") +
      scale_fill_manual(values = c("Final" = "#FF99FF", "Initial" = "#990000", "Imputed" = "#000033")) +
      scale_color_manual(values = c("Final" = "#FF99FF", "Initial" = "#990000", "Imputed" = "#000033")) +
      geom_histogram(alpha = 0.5, binwidth = 0.3, position = "identity") +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white")) + theme(panel.background = element_rect(fill = "white"))

    imp_hist

    ggplot2::ggsave("Imputed_values_histogram.pdf", plot = imp_hist,  path = path_res,
                    scale = 1, width = 5, height = 4, units = "in",
                    dpi = 300, limitsize = TRUE)

    message("An excel with the imputed missing values was created as Dataset_Imputed.xlsx and a histogram documentating these values")
    if (imputation %in% c("LOD/2","LOD","kNN")){    #create histogramm for imputed values

      dataspace_0s$percentage <- dataspace_0s$Number_0_all_groups*100/sum(case_number)
      dataspace$percentage <- dataspace_0s$percentage
      dataspace$mean <- rowMeans(dataspace[,3:(2+sum(case_number))])
      dataspace$log<-log2(dataspace$mean)
      dataspace$rank <- rank(-dataspace$mean)

      abund.plot <- ggplot(dataspace, aes(x = rank, y = log, colour = percentage)) +
        geom_point(size = 3, alpha = 0.8) +
        labs(title = "Protein Abundance Rank", x = "Rank", y = expression(Log[2] ~ "Parts per Million")) +
        scale_color_gradient(low = "darkblue", high = "yellow",
                             name = "Imputations\nin each\nprotein\n(%)") +
        theme_linedraw()+
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              panel.grid = element_line(color = "grey80"),
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 9))
      abund.plot

      ggplot2::ggsave("Proteins_abundance_rank.pdf", plot = abund.plot ,  path = path_res,
                      scale = 1, width = 12, height = 5, units = "in",
                      dpi = 300, limitsize = TRUE, bg = "white")
    }}

if (imputation == FALSE){dataspace <- dataspace

dataspace_0s$percentage <- dataspace_0s$Number_0_all_groups*100/sum(case_number)
dataspace$percentage <- dataspace_0s$percentage
dataspace$mean <- apply(dataspace[, 3:(2+sum(case_number))], 1, function(x) mean(x[x != 0]))
dataspace$log<-log2(dataspace$mean)
dataspace$rank <- rank(-dataspace$mean)

abund.plot <- ggplot(dataspace, aes(x = rank, y = log, colour = percentage)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "Protein Abundance Rank", x = "Rank", y = expression(Log[2] ~ "Parts per Million")) +
  scale_color_gradient(low = "darkblue", high = "yellow",
                       name = "MVs\nin each\nprotein\n(%)") +
  #theme_minimal(base_size = 15)  +
  theme_linedraw()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_line(color = "grey80"),  # Make grids more visible
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9))

abund.plot

ggplot2::ggsave("Proteins_abundance_rank.pdf", plot = abund.plot ,  path = path_res,
                scale = 1, width = 12, height = 5, units = "in",
                dpi = 300, limitsize = TRUE, bg = "white")
}

dataspace$percentage <- NULL
dataspace$mean <- NULL
dataspace$log<- NULL
dataspace$rank <- NULL

groups_for_test <- NULL
for (i in 1:groups_number) {
  groups_for_test <- factor(c(as.character(groups_for_test), rep(group_names[i], times = case_number[i])))
}
mm <- model.matrix(~groups_for_test + 0)
colnames(mm)<- group_names

nndataspace<- dataspace[,-1:-2]
nndataspace <- log2(nndataspace+1)

fit <- limma::lmFit(nndataspace, mm)
fit<- limma::eBayes(fit)
if (groups_number>2){
  anova_res<- limma::topTable(fit, adjust.method = "BH", number = Inf)
  colnames(anova_res)<-paste("ANOVA",colnames(anova_res), sep = "_")
anova_res<- anova_res[,-c(1:groups_number)]}
lima.res <- data.frame()
message("ebayes.")
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
  limma_dataspace <- cbind(anova_res,lima.res,dataspace)}
else {limma_dataspace <- cbind(lima.res,dataspace)}

limma_dataspace<-limma_dataspace %>%
  dplyr::select(Accession, Description, dplyr::everything())
openxlsx::write.xlsx(limma_dataspace, file = "Dataset_limma.test.xlsx")
message("limma test was created")
#####
###- Mann-Whitney and Kruskal-Wallis starts here! - ###
if (MWtest != "Paired" && MWtest != "Independent"){stop("Error. You need to assign MWtest = 'Paired' or 'Indepedent'")}
### 1 ### Specify file for statistical analysis
data2 <- dataspace

# Initialize Average, Standard Deviation, and p-values for all groups
for (j in 1:groups_number) {
  data2[[paste0("Average_G", j)]] <- rowMeans(data2[,coln[[j]]])
  data2[[paste0("St_Dv_G", j)]] <- apply(data2[,coln[[j]]], 1, sd)

  # Perform Mann-Whitney U test comparisons for pairs
  for (k in 1:j) {
    if (k < j) {
      if (MWtest == "Independent") {
        for (i in 1:nrow(data2)) {
          test_list <- stats::wilcox.test(as.numeric(data2[i,coln[[k]]]),
                                          as.numeric(data2[i,coln[[j]]]),
                                          exact = FALSE, paired = FALSE)
          data2[i, paste0("MW_G", j, "vsG", k)] <- test_list$p.value
        }
      } else if (MWtest == "Paired") {
        for (i in 1:nrow(data2)) {
          test_list <- stats::wilcox.test(as.numeric(data2[i,coln[[k]]]),
                                          as.numeric(data2[i,coln[[j]]]),
                                          exact = FALSE, paired = TRUE)
          data2[i, paste0("MW_G", j, "vsG", k)] <- test_list$p.value
        }
      }

      # Adjust the p-values
      data2[[paste0("BH_p_G", j, "vsG", k)]] <- p.adjust(data2[[paste0("MW_G", j, "vsG", k)]], method = "BH")

      # Calculate the ratio and log2 ratio
      data2[[paste0("Ratio_G", j, "vsG", k)]] <- data2[[paste0("Average_G", j)]] / data2[[paste0("Average_G", k)]]
      data2[[paste0("Log2_Ratio.G", j, "vsG", k)]] <- log2(data2[[paste0("Ratio_G", j, "vsG", k)]])
    }
  }
}


Ddataspace<-data2
Ddataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Ddataspace$Description)
Ddataspace$Symbol[Ddataspace$Symbol==Ddataspace$Description] = "Not available"
Ddataspace<-Ddataspace %>%
  dplyr::select(Accession, Description, Symbol, everything())
Fdataspace<-Ddataspace


#if (groups_number != 2){
# Create grouping variable for the Kruskal test


Group<-list()
times<-vector()
for (i in (1:groups_number)){
  times<-case_number[i]
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
  #Do the Kruskal-Wallis test
  df3 <- dataspace4 %>% tidyr::gather(key, value, -Group)
  df4 <- df3 %>% dplyr::group_by(key)
  df4$value<-as.numeric(df4$value)
  df5 <- df4 %>% dplyr::do(broom::tidy(kruskal.test(x= .$value, g = .$Group)))
  Kruskal_Wallis.pvalue <- df5$p.value
  data3<-cbind(Ddataspace,Kruskal_Wallis.pvalue)
  data3$Kruskal_Wallis.pvalue_BH.adjusted<- p.adjust(data3$Kruskal_Wallis.pvalue, method = "BH")

  #Create gene symbols and write the data
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
start_col <- 4 + as.numeric(sum(case_number))
Fdataspace <- Fdataspace %>%
  dplyr::select(1:3, start_col:ncol(Fdataspace), 4:(start_col - 1))

openxlsx::write.xlsx(Fdataspace, file = "Normalized_stats.xlsx")
message("An excel with the statistical tests for the normalized data was created as Normalized_stats.xlsx")



Group <- groups_for_test

Group2<-unique(groups_for_test)

log.dataspace <- log(dataspace[,-c(1:2)]+1,2)

# PCA of the entire data

pca<-prcomp(t(log.dataspace), scale=TRUE, center=FALSE)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group = Group)
qc$PC1.score <- pca$x[,1]
qc$PC2.score <-pca$x[,2]
if (groups_number == 2){
  Group<-list()
  times<-vector()
  for (i in (1:length(case_number))){
    times<-case_number[i]
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
  #scale_colour_manual(values=cbbPalette)+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+
  #ylim(-60,60)+
  ggtitle("Complete set of proteins")+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))+
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_text(size = 12))+
  theme(legend.position="right")

pca.ent

ggplot2::ggsave("PCA_plot_alldata.pdf", plot = pca.ent,  path = path_res,
                scale = 1, width = 5, height = 4, units = "in",
                dpi = 300, limitsize = TRUE)
message("PCA plot using all data was created as PCA_plot_alldata.pdf")
# PCA of the significant data. If number of groups = 2, the script uses the
# unadjusted Mann-Whitney test; else, it uses the unadjusted Kruskal-Wallis test.

which.sig<-vector()
if (groups_number != 2){
  which.sig <- which(Fdataspace$Kruskal_Wallis.pvalue < 0.05)
} else {(which.sig <- which(Ddataspace$MW_G2vsG1 < 0.05))}

which(Ddataspace$MW_G2vsG1< 0.05)
if (length(which.sig) == 0){
  message("There are no significant proteins, to create a PCA plot with them and a heatmap")
  qc[,-1] <- lapply(qc[,-1], function(x) as.numeric(unlist(x)))
  qc[,-1]<-round(qc[,-1],3)

   openxlsx::write.xlsx(qc,file = "Quality_check.xlsx")
}
else {
log.dataspace.sig <- log.dataspace[which.sig,]


zlog.dataspace.sig <- t(scale(t(log.dataspace.sig)))
colnames(zlog.dataspace.sig) <- colnames(log.dataspace.sig)

mycols <- grDevices::colorRampPalette(c("blue", "white", "red"))(100)
heatmap_data<- ComplexHeatmap::Heatmap(as.matrix(zlog.dataspace.sig),
                                       cluster_rows = TRUE,
                                       cluster_columns = TRUE ,
                                       show_row_names = FALSE,
                                       show_column_names = FALSE,
                                       column_split = groups_for_test,
                                       top_annotation = ComplexHeatmap::HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:(groups_number+1)),
                                                                                                           labels = group_names, labels_gp = gpar(col = "white", fontsize = 10))),
                                       col = mycols, column_title = NULL,
                                       heatmap_legend_param = list(
                                         title = "Z-Score",
                                         color_bar = "continuous"
                                       ))
pdf("heatmap.pdf", width = 7.37, height = 6.09)
ComplexHeatmap::draw(heatmap_data)
dev.off()



pca<-prcomp(t(log.dataspace.sig), scale=TRUE, center=FALSE)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group = Group)
qc$PC1.score.Significant <- pca$x[,1]
qc$PC2.score.Significant <-pca$x[,2]
qc[,-1] <- lapply(qc[,-1], function(x) as.numeric(unlist(x)))
qc[,-1]<-round(qc[,-1],3)
openxlsx::write.xlsx(qc,file = "Quality_check.xlsx")

pca.var<-pca$sdev^2
pca.var.per<-round(pca.var/sum(pca.var)*100,1)

pca.data$Group<-factor(pca.data$Group, levels=Group2)

pca.sig<-ggplot2::ggplot(data=pca.data, aes(x=X, y=Y, label=Sample))+
  geom_point(aes(color=Group), size = 2, alpha = 1)+
  #scale_colour_manual(values=cbbPalette)+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep=""))+
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep=""))+
  #ylim(-60,60)+
  ggtitle("Statistically significant proteins")+
  theme_bw()+
  theme(plot.title = element_text(hjust=0.5))+
  guides(color = guide_legend(override.aes = list(size=5)))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_text(size = 12))+
  theme(legend.position="right")

pca.sig

ggplot2::ggsave("PCA_plot_significant.pdf", plot = pca.sig,  path = path_res,
                scale = 1, width = 5, height = 4, units = "in",
                dpi = 300, limitsize = TRUE)
message("PCA plot with the significant data was created as PCA_plot_significant.pdf" )

a<-ggpubr::ggarrange(pca.ent, pca.sig, nrow = 1, ncol=2,
                     common.legend = TRUE, legend = "bottom")

ggplot2::ggsave("PCA_plots_combined.pdf", plot = a,  path = path_res,
                scale = 1, width = 8, height = 4.5, units = "in",
                dpi = 300, limitsize = TRUE)
message ("The 2 PCA plots are combined in PCA_plots_combined.pdf")
}
# Quality check - boxplots of data distribution
p<- function(x) {
  sapply(x, function(label) {
    # Get the last 25 characters from each label
    truncated_label <- substr(label, nchar(label) - 24, nchar(label))
    truncated_label
  })
}
melt.log.dataspace <- reshape2::melt(log.dataspace)
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
    #scale_colour_manual(values=colors)+
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

  ggplot2::ggsave("QC_dataDistribution_withZeros.pdf", plot = qc.boxplots,  path = path_res,
                  scale = 1, width = 12, height = 5, units = "in",
                  dpi = 300, limitsize = TRUE, bg = "white")

}
melt.log.dataspace.na <- melt.log.dataspace
melt.log.dataspace.na$value[melt.log.dataspace.na$value == 0] <- NA
melt.log.dataspace.na$Group <- factor(melt.log.dataspace.na$Group, levels = Group2)
is.factor(melt.log.dataspace.na$variable)

qc.boxplots.na<-ggplot2::ggplot(melt.log.dataspace.na, aes(x=forcats::fct_inorder(variable), y=value, color=Group))+
  geom_boxplot(aes(color = Group),lwd=1, outlier.size=0.2, outlier.alpha = 0.2)+
  #scale_colour_manual(values=colors)+
  xlab("Sample")+
  ylab(expression(Log[2]~"Protein Abundance"))+
  theme_classic()+
  theme(text = element_text(size = 19),
        axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 0.5),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 1)))+
  scale_x_discrete(labels = p) +
  #geom_dotplot(aes(color = Group), binaxis='y', stackdir='center', dotsize=0.1, stackgroups = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5, alpha = 0.5)

qc.boxplots.na
if (imputation == FALSE) {
  ggplot2::ggsave("QC_dataDistribution_NoZeros.pdf", plot = qc.boxplots.na, path = path_res,
                  scale = 1, width = 12, height = 5, units = "in",
                  dpi = 300, limitsize = TRUE, bg = "white")
}
else
{
  ggplot2::ggsave("QC_dataDistribution.pdf", plot = qc.boxplots.na, path = path_res,
                  scale = 1, width = 12, height = 5, units = "in",
                  dpi = 300, limitsize = TRUE, bg = "white")
}
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

ggplot2::ggsave("Violin_plot.pdf", plot = qc.violin,  path = path_res,
                scale = 1, width = 12, height = 5, units = "in",
                dpi = 300, limitsize = TRUE, bg = "white")
message("The Boxplots for each sample have been created!!")

}
