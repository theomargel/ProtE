#' Get data from 2 different sample groups from protein discoverer and create a dataspace wwith all the protein IDs found.
#'
#' Creates an a excel file
#'
#' @param group1  path to the folder where the samples from group 1 are
#' @param group2 path to the folder where the samples from group 2 are
#'
#'
#' @return excel file
#'
#' @examples user_inputs(group1, group2)
#' @export
user_inputs <- function(group1, group2, imputation = TRUE)
  {
group1<- gsub( "\\\\", "/", group1)
group2<-  gsub( "\\\\", "/", group2)

g1.name <- basename(group1)
g2.name <- basename(group2)


#group number is now 2
groups_number <- 2
#create the dataspace for all the data
dataspace <- data.frame()
#assign the excel files to a list
file_names_g1<-list.files(path=group1,pattern="*.xlsx")

setwd(group1)
for (i in 1:length(file_names_g1)) {
  file_case <- readxl::read_excel(paste(group1,file_names_g1[i],sep = "/"), sheet = 1)
  dataspace <- rbind(dataspace,file_case[,1:2])
}
if (groups_number == 2 | groups_number == 3 | groups_number == 4 | groups_number==5 | groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
  file_names_g2<-list.files(path=group2,pattern="*.xlsx")
  setwd(group2)
  for (i in 1:length(file_names_g2)) {
    file_case <- readxl::read_excel(paste(group2,file_names_g2[i],sep = "/"), sheet = 1)
    dataspace <- rbind(dataspace,file_case[,1:2])
  }
}
dataspace <- dataspace[!is.na(dataspace$Accession),]
dupl <- duplicated(dataspace[,1])
dataspace_no_dupl <- dataspace[!dupl,]
dataspace_no_dupl <- droplevels(dataspace_no_dupl)
ML <- data.frame(dataspace_no_dupl, stringsAsFactors = FALSE)
### Merge all areas from the files
dataspace <- dataspace_no_dupl

setwd(group1)
for (i in 1:length(file_names_g1)) {
  file_case <- readxl::read_excel(paste(group1,file_names_g1[i],sep = "/"), sheet = 1)
  file_case[file_case==0] <-1
  file_case[is.na(file_case)]<-0
  dataspace <- merge(x = dataspace,y = file_case[,c(1,9)], by = "Accession" ,all.x = TRUE)
  colnames(dataspace)[length(colnames(dataspace))] <- file_names_g1[i]
}


##Add area from group 2
if (groups_number == 2 | groups_number == 3 | groups_number == 4 | groups_number == 5 | groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
  setwd(group2)
  for (i in 1:length(file_names_g2)) {
    file_case <- readxl::read_excel(paste(group2,file_names_g2[i],sep = "/"), sheet = 1)
    file_case[file_case==0] <-1
    file_case[is.na(file_case)]<-0
    dataspace <- merge(x = dataspace,y = file_case[,c(1,9)], by = "Accession" ,all.x = TRUE)
    colnames(dataspace)[length(colnames(dataspace))] <- file_names_g2[i]
  }
}

colnames(dataspace) <- gsub(".xlsx", "", colnames(dataspace))


    path_res <- readline(prompt = "Specify folder where you want to save the results: e.g.: C:/Users/User/Documents/")
  if (!dir.exists(path_res)) {
    stop("The specified folder does not exist.")
  }
    setwd(path_res)
    openxlsx::write.xlsx(dataspace, file = "Masterlist.xlsx")
    message("the excel of the list with all proteomics data was created as Masterlist.xlsx")
    #normalize PPm
    dataspace[, -1:-2] <- lapply(dataspace[, -1:-2], function(x) {
      sum_x <- sum(x, na.rm = TRUE)  # Sum of the column, ignoring NAs
      ifelse(is.na(x), NA, (x / sum_x) * 10^6)  # NA=0 , normalize the rest
    })

    dat.dataspace<-dataspace

    #assign values to case number
    case_number <- numeric(groups_number)

    for (i in 1:groups_number) {
      case_number[i] <- length(get(paste0("file_names_g",i)))
    }



    threshold_value <- readline ("Set a percentage threshold for missing values.
                         Proteins with missing values greater than this threshold, will be deleted.
                           (Enter a number, e.g., 50. Write 'D' for default, which is 40%):")
      if (threshold_value < 0 && threshold_value > 100 && threshold_value != "D") {stop("Error, you should add D for default or a number between 0 and 100")}

      threshold<-numeric(groups_number)

      if(threshold_value == "D"){
        for (i in 1:groups_number){threshold[i] <- ceiling((case_number[i]-(case_number[i])*60/100)+0.00000000001)}
      } else {
        threshold_value <- 100- as.numeric(threshold_value)
        # Iterate over each group and calculate the new threshold
        for (i in 1:groups_number) {
          threshold[i] <-  ceiling(case_number[i]-(case_number[i])*(as.numeric(threshold_value)/100)+0.00000000001)}
      }


    method_number <- readline ("How to treat Proteome Discoverer's bugs (blank values)?
  0 = as zeros
  1 = as group averages
  ")
    if (method_number != 0 && method_number != 1) {stop("Error, you should add 0 or 1")}

    # Create identifier variables for the thhreshold and statistics

    if (groups_number==2){
      control_last <-(3+case_number[1]-1)
      coln <- c(3:control_last)
      case_last <- (control_last+case_number[2])
      coln2 <- c((control_last+1):case_last)
    }

    # assign average of group to discoverer bugs!
    if (method_number==1){
      dat.dataspace[dat.dataspace==0] <- 1
      dat.dataspace[is.na(dat.dataspace)] <- 0

      if (groups_number==2){
        dat.data.1 <- dat.dataspace[,coln]
        dat.data.2 <- dat.dataspace[,coln2]

        rep.data.1 <- dataspace[,coln]
        rep.data.2 <- dataspace[,coln2]

        m1<-rowMeans(rep.data.1)
        m2<-rowMeans(rep.data.2)

        idx1 <- dat.data.1 == 1
        idx2 <- dat.data.2 == 1

        tmp1 <- idx1 * m1
        tmp2 <- idx2 * m2

        rep.data.1[idx1] <- tmp1[idx1]
        rep.data.2[idx2] <- tmp2[idx2]

        dataspace <- data.frame(dat.dataspace[,c(1:2)], rep.data.1, rep.data.2)
      }}

     Gdataspace<-dataspace
    Gdataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Gdataspace$Description)
    Gdataspace$Symbol[Gdataspace$Symbol==Gdataspace$Description] = "Not available"
    Gdataspace<-Gdataspace %>%
      dplyr::select(Accession, Description, Symbol, everything())



    colnames(Gdataspace) <- gsub(".xlsx", "", colnames(Gdataspace))
    openxlsx::write.xlsx(Gdataspace, file = "Normalized.xlsx")





    if (groups_number==2){


      #Count the numberof 0 for group1
      dataspace[is.na(dataspace)] <- 0
      for (i in c(1:length(dataspace[,1]))){
        a <- table(dataspace[i,coln]==0)["TRUE"]
        if(is.na(a)){
          dataspace[i,"Number_0_group1"] <- 0
        } else{
          dataspace[i,"Number_0_group1"] <- table(dataspace[i,coln]==0)["TRUE"]
        }
      }

      #Count the number of 0 for group 2
      for (i in c(1:length(dataspace[,1]))){
        a <- table(dataspace[i,coln2]==0)["TRUE"]
        if(is.na(a)){
          dataspace[i,"Number_0_group2"] <- 0
        } else{
          dataspace[i,"Number_0_group2"] <- table(dataspace[i,coln2]==0)["TRUE"]
        }
        }




       setwd<-path_res
      openxlsx::write.xlsx(dataspace, file = "Dataset_before_threshold.xlsx")

      #Delete row with x or more zeros
      dataspace <- dataspace[dataspace$Number_0_group1<threshold[1] | dataspace$Number_0_group2<threshold[2],]


      #write.table(dataspace, file=paste(path_res,"Dataset_threshold_applied.xlsx",sep=""), dec=".",sep="\t", row.names=FALSE)
      openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")
      dataspace$Number_0_group1 <- NULL
      dataspace$Number_0_group2 <- NULL}

##imputation KNN
if (imputation == TRUE) {
dataspace[dataspace==0] <- NA
dataspace <- VIM::kNN(dataspace, imp_var = FALSE)
openxlsx::write.xlsx(dataspace,file = "Dataset_IMPUTED.xlsx")}
else {dataspace <- dataspace }

### - Mann-Whitney and Kruskal-Wallis starts here! - ###

### 1 ### Specify file for statistical analysis
data2 <- dataspace

if (groups_number==2){
  data2$Average_G1 <- rowMeans(data2[,coln])
  data2$Average_G2 <- rowMeans(data2[,coln2])
  data2$St_Dv_G1 <- apply(data2[,coln], 1, sd)
  data2$St_Dv_G2 <- apply(data2[,coln2], 1, sd)
  ### Calculate unadjusted p-value
  for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=TRUE)
    data2[i,"MW_G2vsG1"]<-test_list[[3]]
  }
  #### adjust the p-values
  data2$BH_p_G2vsG1 <- stats::p.adjust(data2$MW_G2vsG1, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G2vsG1 <- (data2$Average_G2)/(data2$Average_G1)
  data2$Log2_Ratio.G2vsG1 <- log2(data2$Ratio_G2vsG1)
  #data$subtraction <- (data$average_case)-(data$average_control)


  Ddataspace<-data2
  Ddataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Ddataspace$Description)
  Ddataspace$Symbol[Ddataspace$Symbol==Ddataspace$Description] = "Not available"
  Ddataspace<-Ddataspace %>%
    dplyr::select(Accession, Description, Symbol, everything())
  Fdataspace<-Ddataspace
}

if (groups_number != 2){


  # Create grouping variable for the Kruskal test


  Group<-list()
  times<-vector()
  for (i in (1:length(case_number))){
    times<-case_number[i]
    Group[[i]] <- rep(paste0("G",i), times)
    times<-NULL
  }
  Group<-unlist(Group)



  dataspace3<-t(dataspace[,-c(1:2)])
  dataspace3<-as.data.frame(dataspace3)

  colnames(dataspace3)<-dataspace[,1]

  #dataspace3 <- dataspace3 %>% mutate_if(is.factor, as.numeric)

  dataspace4<-cbind(Group,dataspace3)
  dataspace4$Group<-as.factor(dataspace4$Group)


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

if (groups_number == 2){
  namesc<-colnames(Fdataspace)
  namesc<-gsub("G1", g1.name, namesc)
  namesc<-gsub("G2", g2.name, namesc)
  colnames(Fdataspace)<-namesc
}

colnames(Fdataspace) <- gsub(".xlsx", "", colnames(Fdataspace))

write.xlsx(Fdataspace, file = "Normalized_stats.xlsx")
message("normalized stats ok!")

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



Group2<-unique(Group)

log.dataspace <- log(dataspace[,-c(1:2)]+1)

# PCA of the entire data

pca<-prcomp(t(log.dataspace), scale=FALSE, center=FALSE)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group = Group)


pca.var<-pca$sdev^2
pca.var.per<-round(pca.var/sum(pca.var)*100,1)

pca.data$Group<-factor(pca.data$Group, levels=Group2)

pca.ent<-ggplot2::ggplot(data=pca.data, aes(x=X, y=Y, label=Sample))+
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

ggplot2::ggsave("PCA_plot_alldata.tiff", plot = pca.ent, device = "tiff", path = path_res,
       scale = 1, width = 5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

# PCA of the significant data. If number of groups = 2, the script uses the
# unadjusted Mann-Whitney test; else, it uses the unadjusted Kruskal-Wallis test.

which.sig<-vector()
if (groups_number != 2){
  which.sig <- which(Fdataspace$Kruskal_Wallis.pvalue < 0.05)
} else {(which.sig <- which(Ddataspace$MW_G2vsG1 < 0.05))}

which(Ddataspace$MW_G2vsG1< 0.05)
log.dataspace.sig <- log.dataspace[which.sig,]


pca<-prcomp(t(log.dataspace.sig), scale=FALSE, center=FALSE)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group = Group)


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


ggplot2::ggsave("PCA_plot_significant.tiff", plot = pca.sig, device = "tiff", path = path_res,
       scale = 1, width = 5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)


a<-ggpubr::ggarrange(pca.ent, pca.sig, nrow = 1, ncol=2,
             common.legend = TRUE, legend = "bottom")



ggplot2::ggsave("PCA_plots_combined.tiff", plot = a, device = "tiff", path = path_res,
       scale = 1, width = 8, height = 4.5, units = "in",
       dpi = 300, limitsize = TRUE)

# Quality check - boxplots of data distribution

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


qc.boxplots<-ggplot2::ggplot(melt.log.dataspace, aes(x=fct_inorder(variable), y=value, color=Group))+
  geom_boxplot(aes(color = Group),lwd=1, outlier.size=0.2, outlier.alpha = 0.2)+
  #scale_colour_manual(values=colors)+
  xlab("Sample")+
  ylab("Log parts per million")+
  theme_classic()+
  theme(text = element_text(size = 19),
        axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 0.5),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 1)))+
  #geom_dotplot(aes(color = Group), binaxis='y', stackdir='center', dotsize=0.1, stackgroups = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5, alpha = 0.5)


qc.boxplots

ggplot2::ggsave("QC_dataDistribution_withZeros.tiff", plot = qc.boxplots, device = "tiff", path = path_res,
       scale = 1, width = 12, height = 5, units = "in",
       dpi = 300, limitsize = TRUE, bg = "white")


melt.log.dataspace.na <- melt.log.dataspace

melt.log.dataspace.na$value[melt.log.dataspace.na$value == 0] <- NA

melt.log.dataspace.na$Group <- factor(melt.log.dataspace.na$Group, levels = Group2)

is.factor(melt.log.dataspace.na$variable)

qc.boxplots.na<-ggplot2::ggplot(melt.log.dataspace.na, aes(x=fct_inorder(variable), y=value, color=Group))+
  geom_boxplot(aes(color = Group),lwd=1, outlier.size=0.2, outlier.alpha = 0.2)+
  #scale_colour_manual(values=colors)+
  xlab("Sample")+
  ylab("Log parts per million")+
  theme_classic()+
  theme(text = element_text(size = 19),
        axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 0.5),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 1)))+
  #geom_dotplot(aes(color = Group), binaxis='y', stackdir='center', dotsize=0.1, stackgroups = FALSE)+
  geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5, alpha = 0.5)


qc.boxplots.na

ggplot2::ggsave("QC_dataDistribution_NoZeros.tiff", plot = qc.boxplots.na, device = "tiff", path = path_res,
       scale = 1, width = 12, height = 5, units = "in",
       dpi = 300, limitsize = TRUE, bg = "white")




}
