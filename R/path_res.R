#' Proteome Discoverer analysis
#'
#' It takes Proteomics Data from samples in different groups, in the format they are created by Proteome Discoverer (PD). It concatenates the Accession IDs, Descriptions and Areas from different PD export files into a master table and performs exploratory data analysis. The function outputs a normalized Parts Per Million protein dataset along with descriptive statistics and results of significance testing. The script also creates exploratory plots such as relative log espression boxplots and PCA plots.
#'
#' @param ... The specific path to the folder where the samples from each group are located. They are passed as unnamed arguments via "...".  Attention: Add '/' between the directories.
#' @param global_threshold TRUE/FALSE If TRUE threshold for missing values will be applied to the groups altogether, if FALSE to each group seperately
#' @param imputation TRUE/FALSE Data imputation using kNN classification or assigning missing values as 0.
#' @param MWtest Either "Paired" for a Wilcoxon Signed-rank test or "Independent" for a Mann-Whitney U test.
#' @param threshold_value The percentage of missing values per protein that will cause its omittion.
#' @param bugs Either 0 to treat Proteome Discoverer bugs as Zeros (0) or "average" to convert them into the average of the protein between the samples.
#'
#' @return Excel files with the proteomic values from all samples, processed with normalization and imputation and substraction of samples with high number of missing values. PCA plots for all or for just the significant correlations, and boxplots for the proteins of each sample.
#' @importFrom openxlsx write.xlsx  read.xlsx
#' @importFrom dplyr select  group_by  do everything  %>%
#' @importFrom tidyr gather pivot_longer
#' @importFrom broom tidy
#' @importFrom reshape2 melt
#' @importFrom ggpubr ggarrange
#' @importFrom ggplot2 ggplot ggsave scale_color_gradient element_line theme_linedraw scale_fill_manual scale_color_manual aes geom_histogram element_rect geom_point xlab ylab ggtitle theme_bw theme_minimal theme element_text guides guide_legend geom_boxplot labs theme_classic element_blank geom_jitter position_jitter
#' @importFrom VIM kNN
#' @importFrom stats kruskal.test p.adjust prcomp sd wilcox.test model.matrix
#' @importFrom forcats fct_inorder
#' @importFrom limma topTable eBayes contrasts.fit lmFit
#'
#' @examples #' # Example of running the function with paths for two groups.
#' #Do not add if (interactive()){} condition in your code
#' if (interactive()){
#' user_inputs(
#'   "C:/Users/User/Documents/T0_samples",
#'   "C:/Users/User/Documents/T0_samples",
#'   MWtest = "Paired",
#'   imputation = TRUE,
#'   global_threshold = TRUE
#' )}
#'
#' @export

user_inputs <- function(...,
                        imputation = TRUE,
                        global_threshold = TRUE,
                        MWtest = "Independent",
                        threshold_value = 50,
                        bugs = 0)
  {
  group1 = group2 = Accession =Description =Symbol =X =Y = percentage=Sample= variable =.= g1.name =g2.name= g3.name =g4.name= g5.name =g6.name= g7.name= g8.name =g9.name =group3= group4= group5= group6 =group7= group8= group9 =key =value = NULL

group_paths <- list(...)
groups_number <- length(group_paths)
if (groups_number>9){stop("You can add up to 9 groups")}

group_paths<- gsub( "\\\\", "/", group_paths)
for (i in 1:groups_number) {
  assign(paste0("group",i),group_paths[[i]])
}
group_names <- basename(group_paths)
for (i in 1:groups_number) {
assign(paste0("g",i,".name"),group_names[[i]])}

#group number is now 2

#create the dataspace for all the data
dataspace <- data.frame()

print("this has succesfully been installed")

#assign An excel files to a list
file_names_g1<-list.files(path=group1,pattern="*.xlsx")

setwd(group1)
for (i in 1:length(file_names_g1)) {
  file_case <- openxlsx::read.xlsx(paste(group1,file_names_g1[i],sep = "/"), sheet = 1)
  dataspace <- rbind(dataspace,file_case[,1:2])
}
if (groups_number == 2 | groups_number == 3 | groups_number == 4 | groups_number==5 | groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
  file_names_g2<-list.files(path=group2,pattern="*.xlsx")
  setwd(group2)
  for (i in 1:length(file_names_g2)) {
    file_case <- openxlsx::read.xlsx(paste(group2,file_names_g2[i],sep = "/"), sheet = 1)
    dataspace <- rbind(dataspace,file_case[,1:2])
  }
  # group 3 protein IDs
  if (groups_number==3 | groups_number == 4 | groups_number==5 | groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
    file_names_g3<-list.files(path=group3,pattern="*.xlsx")
    setwd(group3)
    for (i in 1:length(file_names_g3)) {
      file_case <- openxlsx::read.xlsx(paste(group3,file_names_g3[i],sep = "/"), sheet = 1)
      dataspace <- rbind(dataspace,file_case[,1:2])
    }
  }

  # group 4 protein IDs
  if (groups_number==4 | groups_number==5 | groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
    file_names_g4<-list.files(path=group4,pattern="*.xlsx")
    setwd(group4)
    for (i in 1:length(file_names_g4)) {
      file_case <- openxlsx::read.xlsx(paste(group4,file_names_g4[i],sep = "/"), sheet = 1)
      dataspace <- rbind(dataspace,file_case[,1:2])
    }
  }



  # group 5 protein IDs
  if (groups_number==5 | groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
    file_names_g5<-list.files(path=group5,pattern="*.xlsx")
    setwd(group5)
    for (i in 1:length(file_names_g5)) {
      file_case <- openxlsx::read.xlsx(paste(group5,file_names_g5[i],sep = "/"), sheet = 1)
      dataspace <- rbind(dataspace,file_case[,1:2])
    }
  }

  # group 6 protein IDs
  if (groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
    file_names_g6<-list.files(path=group6,pattern="*.xlsx")
    setwd(group6)
    for (i in 1:length(file_names_g6)) {
      file_case <- openxlsx::read.xlsx(paste(group6,file_names_g6[i],sep = "/"), sheet = 1)
      dataspace <- rbind(dataspace,file_case[,1:2])
    }
  }

  # group 7 protein IDs
  if (groups_number==7 | groups_number==8 | groups_number==9){
    file_names_g7<-list.files(path=group7,pattern="*.xlsx")
    setwd(group7)
    for (i in 1:length(file_names_g7)) {
      file_case <- openxlsx::read.xlsx(paste(group7,file_names_g7[i],sep = "/"), sheet = 1)
      dataspace <- rbind(dataspace,file_case[,1:2])
    }
  }

  # group 8 protein IDs
  if (groups_number==8 | groups_number==9){
    file_names_g8<-list.files(path=group8,pattern="*.xlsx")
    setwd(group8)
    for (i in 1:length(file_names_g8)) {
      file_case <- openxlsx::read.xlsx(paste(group8,file_names_g8[i],sep = "/"), sheet = 1)
      dataspace <- rbind(dataspace,file_case[,1:2])
    }
  }

  # group 9 protein IDs
  if (groups_number==9){
    file_names_g9<-list.files(path=group9,pattern="*.xlsx")
    setwd(group9)
    for (i in 1:length(file_names_g9)) {
      file_case <- openxlsx::read.xlsx(paste(group9,file_names_g9[i],sep = "/"), sheet = 1)
      dataspace <- rbind(dataspace,file_case[,1:2])
    }
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
  file_case <- openxlsx::read.xlsx(paste(group1,file_names_g1[i],sep = "/"), sheet = 1)
  file_case[file_case==0] <-1
  file_case[is.na(file_case)]<-0
  dataspace <- merge(x = dataspace,y = file_case[,c(1,9)], by = "Accession" ,all.x = TRUE)
  colnames(dataspace)[length(colnames(dataspace))] <- file_names_g1[i]
}


##Add area from group 2
if (groups_number == 2 | groups_number == 3 | groups_number == 4 | groups_number == 5 | groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
  setwd(group2)
  for (i in 1:length(file_names_g2)) {
    file_case <- openxlsx::read.xlsx(paste(group2,file_names_g2[i],sep = "/"), sheet = 1)
    file_case[file_case==0] <-1
    file_case[is.na(file_case)]<-0
    dataspace <- merge(x = dataspace,y = file_case[,c(1,9)], by = "Accession" ,all.x = TRUE)
    colnames(dataspace)[length(colnames(dataspace))] <- file_names_g2[i]
  }
}
##Add area from group 3
if (groups_number==3 | groups_number == 4 | groups_number == 5 | groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
  setwd(group3)
  for (i in 1:length(file_names_g3)) {
    file_case <- openxlsx::read.xlsx(paste(group3,file_names_g3[i],sep = "/"), sheet = 1)
    file_case[file_case==0] <-1
    file_case[is.na(file_case)]<-0
    dataspace <- merge(x = dataspace,y = file_case[,c(1,9)], by = "Accession" ,all.x = TRUE)
    colnames(dataspace)[length(colnames(dataspace))] <- file_names_g3[i]
  }
}

##Add area from group 4
if (groups_number==4 | groups_number==5 | groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
  setwd(group4)
  for (i in 1:length(file_names_g4)) {
    file_case <- openxlsx::read.xlsx(paste(group4,file_names_g4[i],sep = "/"), sheet = 1)
    file_case[file_case==0] <-1
    file_case[is.na(file_case)]<-0
    dataspace <- merge(x = dataspace,y = file_case[,c(1,9)], by = "Accession" ,all.x = TRUE)
    colnames(dataspace)[length(colnames(dataspace))] <- file_names_g4[i]
  }
}

##Add area from group 5
if (groups_number==5 | groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
  setwd(group5)
  for (i in 1:length(file_names_g5)) {
    file_case <- openxlsx::read.xlsx(paste(group5,file_names_g5[i],sep = "/"), sheet = 1)
    file_case[file_case==0] <-1
    file_case[is.na(file_case)]<-0
    dataspace <- merge(x = dataspace,y = file_case[,c(1,9)], by = "Accession" ,all.x = TRUE)
    colnames(dataspace)[length(colnames(dataspace))] <- file_names_g5[i]
  }
}

##Add area from group 6
if (groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
  setwd(group6)
  for (i in 1:length(file_names_g6)) {
    file_case <- openxlsx::read.xlsx(paste(group6,file_names_g6[i],sep = "/"), sheet = 1)
    file_case[file_case==0] <-1
    file_case[is.na(file_case)]<-0
    dataspace <- merge(x = dataspace,y = file_case[,c(1,9)], by = "Accession" ,all.x = TRUE)
    colnames(dataspace)[length(colnames(dataspace))] <- file_names_g6[i]
  }
}

##Add area from group 7
if (groups_number==7 | groups_number==8 | groups_number==9){
  setwd(group7)
  for (i in 1:length(file_names_g7)) {
    file_case <- openxlsx::read.xlsx(paste(group7,file_names_g7[i],sep = "/"), sheet = 1)
    file_case[file_case==0] <-1
    file_case[is.na(file_case)]<-0
    dataspace <- merge(x = dataspace,y = file_case[,c(1,9)], by = "Accession" ,all.x = TRUE)
    colnames(dataspace)[length(colnames(dataspace))] <- file_names_g7[i]
  }
}

##Add area from group 8
if (groups_number==8 | groups_number==9){
  setwd(group8)
  for (i in 1:length(file_names_g8)) {
    file_case <- openxlsx::read.xlsx(paste(group8,file_names_g8[i],sep = "/"), sheet = 1)
    file_case[file_case==0] <-1
    file_case[is.na(file_case)]<-0
    dataspace <- merge(x = dataspace,y = file_case[,c(1,9)], by = "Accession" ,all.x = TRUE)
    colnames(dataspace)[length(colnames(dataspace))] <- file_names_g8[i]
  }
}

##Add area from group 9
if (groups_number==9){
  setwd(group9)
  for (i in 1:length(file_names_g9)) {
    file_case <- openxlsx::read.xlsx(paste(group9,file_names_g9[i],sep = "/"), sheet = 1)
    file_case[file_case==0] <-1
    file_case[is.na(file_case)]<-0
    dataspace <- merge(x = dataspace,y = file_case[,c(1,9)], by = "Accession" ,all.x = TRUE)
    colnames(dataspace)[length(colnames(dataspace))] <- file_names_g9[i]
  }
}

colnames(dataspace) <- gsub(".xlsx", "", colnames(dataspace))
path_g1 <- dirname(group1)
path_res <- file.path(path_g1, "MS_analysis")
dir.create(path_res, showWarnings = FALSE)

#assign values to case number
case_number <- numeric(groups_number)

for (i in 1:groups_number) {
  case_number[i] <- length(get(paste0("file_names_g",i)))
}

    setwd(path_res)
    openxlsx::write.xlsx(dataspace, file = "Masterlist.xlsx")
    message("An excel of the list with all proteomics data was created as Masterlist.xlsx")

    zero_per_sample <- colSums(is.na(dataspace[,-1:-2]))*100/nrow(dataspace)
    IDs <- colSums(!is.na(dataspace[,-1:-2]))

    #normalize PPm
    dataspace[, -1:-2] <- lapply(dataspace[, -1:-2], function(x) {
      sum_x <- sum(x, na.rm = TRUE)  # Sum of the column, ignoring NAs
      ifelse(is.na(x), NA, (x / sum_x) * 10^6)  # NA=0 , normalize the rest
    })
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


    if (bugs != 0 && bugs != "average") {stop("Error, you should assign bugs as 0 or average")}

    # Create identifier variables for the thhreshold and statistics

    if (groups_number >= 2){
      control_last <-(3+case_number[1]-1)
      coln <- c(3:control_last)
      case_last <- (control_last+case_number[2])
      coln2 <- c((control_last+1):case_last)
    }

    if (groups_number >=3){
      case2_last <- (case_last+case_number[3])
      coln3 <- c((case_last+1):case2_last)
    }


    if (groups_number>=4){
      case3_last <- (case2_last+case_number[4])
      coln4 <- c((case2_last+1):case3_last)
    }

    if (groups_number>=5){
      case4_last <- (case3_last+case_number[5])
      coln5 <- c((case3_last+1):case4_last)
    }


    if (groups_number>=6){
      case5_last <- (case4_last+case_number[6])
      coln6 <- c((case4_last+1):case5_last)
    }
    if (groups_number>=7){
      case6_last <- (case5_last+case_number[7])
      coln7 <- c((case5_last+1):case6_last)
    }

    if (groups_number>=8){
      case7_last <- (case6_last+case_number[8])
      coln8 <- c((case6_last+1):case7_last)
    }

    if (groups_number>=9){
      case8_last <- (case7_last+case_number[9])
      coln9 <- c((case7_last+1):case8_last)
    }


    # assign average of group to discoverer bugs!
    if (bugs== "average"){
      dat.dataspace[dat.dataspace==0] <- 1
      dat.dataspace[is.na(dat.dataspace)] <- 0

      if (groups_number>=2){
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
      }
      if (groups_number>=3){

        dat.data.3 <- dat.dataspace[,coln3]

        rep.data.3 <- dataspace[,coln3]

        m3<-rowMeans(rep.data.3)

        idx3 <- dat.data.3 == 1

        tmp3 <- idx3 * m3

        rep.data.3[idx3] <- tmp3[idx3]

        dataspace <- data.frame(dat.dataspace[,c(1:2)], rep.data.1, rep.data.2, rep.data.3)
      }

      if (groups_number>=4){

        dat.data.4 <- dat.dataspace[,coln4]

        rep.data.4 <- dataspace[,coln4]

        m4<-rowMeans(rep.data.4)


        idx4 <- rep.data.4 == 1

        tmp4 <- idx4 * m4

        rep.data.4[idx4] <- tmp4[idx4]

        dataspace <- data.frame(dat.dataspace[,c(1:2)], rep.data.1, rep.data.2, rep.data.3, rep.data.4)
      }

      if (groups_number>=5){

        dat.data.5 <- dat.dataspace[,coln5]

        rep.data.5 <- dataspace[,coln5]

        m5<-rowMeans(rep.data.5)

        idx5 <- dat.data.5 == 1

        tmp5 <- idx5 * m5

        rep.data.5[idx5] <- tmp5[idx5]

        dataspace <- data.frame(dat.dataspace[,c(1:2)], rep.data.1, rep.data.2, rep.data.3, rep.data.4, rep.data.5)
      }

      if (groups_number>=6){

        dat.data.6 <- dat.dataspace[,coln6]

        rep.data.6 <- dataspace[,coln6]

        m6<-rowMeans(rep.data.6)

        idx6 <- dat.data.6 == 1

        tmp6 <- idx6 * m6

        rep.data.6[idx6] <- tmp6[idx6]

        dataspace <- data.frame(dat.dataspace[,c(1:2)], rep.data.1, rep.data.2, rep.data.3,
                                rep.data.4, rep.data.5, rep.data.6)
      }

      if (groups_number>=7){

        dat.data.7 <- dat.dataspace[,coln7]

        rep.data.7 <- dataspace[,coln7]

        m7<-rowMeans(rep.data.7)

        idx7 <- dat.data.7 == 1

        tmp7 <- idx7 * m7

        rep.data.7[idx7] <- tmp7[idx7]

        dataspace <- data.frame(dat.dataspace[,c(1:2)], rep.data.1, rep.data.2, rep.data.3,
                                rep.data.4, rep.data.5, rep.data.6, rep.data.7)
      }

      if (groups_number>=8){

        dat.data.8 <- dat.dataspace[,coln8]


        rep.data.8 <- dataspace[,coln8]

        m8<-rowMeans(rep.data.8)

        idx8 <- dat.data.8 == 1

        tmp8 <- idx8 * m8

        rep.data.8[idx8] <- tmp8[idx8]

        dataspace <- data.frame(dat.dataspace[,c(1:2)], rep.data.1, rep.data.2, rep.data.3,
                                rep.data.4, rep.data.5, rep.data.6, rep.data.7, rep.data.8)
      }


      if (groups_number==9){

        dat.data.9 <- dat.dataspace[,coln9]

        rep.data.9 <- dataspace[,coln9]

        m9<-rowMeans(rep.data.9)

        idx9 <- dat.data.9 == 1

        tmp9 <- idx9 * m9

        rep.data.9[idx9] <- tmp9[idx9]

        dataspace <- data.frame(dat.dataspace[,c(1:2)], rep.data.1, rep.data.2, rep.data.3,
                                rep.data.4, rep.data.5, rep.data.6, rep.data.7, rep.data.8,
                                rep.data.9)
      }


    }


     Gdataspace<-dataspace
    Gdataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Gdataspace$Description)
    Gdataspace$Symbol[Gdataspace$Symbol==Gdataspace$Description] = "Not available"
    Gdataspace<-Gdataspace %>%
      dplyr::select(Accession, Description, Symbol, everything())



    colnames(Gdataspace) <- gsub(".xlsx", "", colnames(Gdataspace))
    openxlsx::write.xlsx(Gdataspace, file = "Normalized.xlsx")
    message("An excel file with the normalized values for the proteomics data was created as Normalized.xlsx")
    if (groups_number>=1){
      control_last <-(3+case_number[1]-1)
      coln <- c(3:control_last)

      #Count the numberof 0 for group1
      dataspace[is.na(dataspace)] <- 0
      for (i in c(1:length(dataspace[,1]))){
        a <- table(dataspace[i,coln]==0)["TRUE"]
        if(is.na(a)){
          dataspace[i,"Number_0_group1"] <- 0
        } else{
          dataspace[i,"Number_0_group1"] <- table(dataspace[i,coln]==0)["TRUE"]
        }
      }}
    if (groups_number==1){
      setwd<-path_res
      write.xlsx(dataspace, file = "Dataset_before_threshold.xlsx")
      #write.xlsx(dataspace, "C:/Users/Raf/Desktop/Cell lines GFP/high_con/test/Dataset_before_threshold.xlsx")
      #Delete row with x or more zeros
      dataspace <- dataspace[dataspace$Number_0_group1<threshold[1],]

      write.xlsx(dataspace, file = "Dataset_threshold_applied.xlsx")
      message("An excel file with the proteins that have % of missing values below the threshold was created as Dataset_threshold_applied.xlsx")
      #write.xlsx(dataspace, "G:/LC-MS Analysis Normalized Area/weeks_6/Outputs/Dataset_threshold_applied.xlsx")
      dataspace_0s<- dataspace

      dataspace$Number_0_group1 <- NULL

    }
    if (groups_number>=2){

      #Count the number of 0 for group 2
      for (i in c(1:length(dataspace[,1]))){
        a <- table(dataspace[i,coln2]==0)["TRUE"]
        if(is.na(a)){
          dataspace[i,"Number_0_group2"] <- 0
        } else{
          dataspace[i,"Number_0_group2"] <- table(dataspace[i,coln2]==0)["TRUE"]
        }
      }
 dataspace$Number_0_all_groups <- dataspace$Number_0_group1 + dataspace$Number_0_group2}

    if (groups_number==2){
           openxlsx::write.xlsx(dataspace, file = "Dataset_before_threshold.xlsx")

       setwd<-path_res
      if (global_threshold == TRUE) {
        dataspace <- dataspace[dataspace$Number_0_all_groups<threshold,]
        openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")}

      if (global_threshold == FALSE) {
      #Delete row with x or more zeros
      dataspace <- dataspace[dataspace$Number_0_group1<threshold[1] | dataspace$Number_0_group2<threshold[2],]

      #write.table(dataspace, file=paste(path_res,"Dataset_threshold_applied.xlsx",sep=""), dec=".",sep="\t", row.names=FALSE)
      openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")}
      message("An excel file with the proteins that have % of missing values below the threshold was created as Dataset_threshold_applied.xlsx")
      dataspace_0s<- dataspace

      dataspace$Number_0_group1 <- NULL
      dataspace$Number_0_group2 <- NULL
      dataspace$Number_0_all_groups <- NULL }

    if (groups_number>=3){

      #Count the number of 0 for group 3
      for (i in c(1:length(dataspace[,1]))){
        a <- table(dataspace[i,coln3]==0)["TRUE"]
        if(is.na(a)){
          dataspace[i,"Number_0_group3"] <- 0
        } else{
          dataspace[i,"Number_0_group3"] <- table(dataspace[i,coln3]==0)["TRUE"]
        }
      }


     dataspace$Number_0_all_groups <- dataspace$Number_0_group1 + dataspace$Number_0_group2 +dataspace$Number_0_group3
    }
    if (groups_number ==3){

       openxlsx::write.xlsx(dataspace, file = "Dataset_before_threshold.xlsx")

      setwd<-path_res
      if (global_threshold == TRUE) {
        dataspace <- dataspace[dataspace$Number_0_all_groups<threshold,]
        openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")}

      if (global_threshold == FALSE) {
  #Delete row with x or more zeros
      dataspace <- dataspace[dataspace$Number_0_group1<threshold[1] | dataspace$Number_0_group2<threshold[2]| dataspace$Number_0_group3<threshold[3],]

      #write.table(dataspace, file=paste(path_res,"Dataset_threshold_applied.xlsx",sep=""), dec=".",sep="\t", row.names=FALSE)
      openxlsx::write.xlsx(dataspace, file = "Dataset_threshold_applied.xlsx")
      }
      message("An excel file with the proteins that have % of missing values below the threshold was created as Dataset_threshold_applied.xlsx")
      dataspace_0s<- dataspace

      dataspace$Number_0_group1 <- NULL
      dataspace$Number_0_group2 <- NULL
      dataspace$Number_0_group3 <- NULL
       dataspace$Number_0_all_groups <- NULL
    }

    if (groups_number>=4){

      #Count the number of 0 for group 4
      for (i in c(1:length(dataspace[,1]))){
        a <- table(dataspace[i,coln4]==0)["TRUE"]
        if(is.na(a)){
          dataspace[i,"Number_0_group4"] <- 0
        } else{
          dataspace[i,"Number_0_group4"] <- table(dataspace[i,coln4]==0)["TRUE"]
        }
      }
      dataspace$Number_0_all_groups <- dataspace$Number_0_group1 + dataspace$Number_0_group2 +dataspace$Number_0_group3 +dataspace$Number_0_group4
    }
    if (groups_number==4){
       setwd<-path_res
       openxlsx::write.xlsx(dataspace, file = "Dataset_before_threshold.xlsx")

      if (global_threshold == TRUE) {
        dataspace <- dataspace[dataspace$Number_0_all_groups<threshold,]
        openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")}

      if (global_threshold == FALSE) {

      #Delete row with x or more zeros
      dataspace <- dataspace[dataspace$Number_0_group1<threshold[1] | dataspace$Number_0_group2<threshold[2]| dataspace$Number_0_group3<threshold[3] | dataspace$Number_0_group4<threshold[4],]

      #write.table(dataspace, file=paste(path_res,"Dataset_threshold_applied.xlsx",sep=""), dec=".",sep="\t", row.names=FALSE)
      openxlsx::write.xlsx(dataspace, file = "Dataset_threshold_applied.xlsx")
      }
      message("An excel file with the proteins that have % of missing values below the threshold was created as Dataset_threshold_applied.xlsx")
      dataspace_0s<- dataspace

      dataspace$Number_0_group1 <- NULL
      dataspace$Number_0_group2 <- NULL
      dataspace$Number_0_group3 <- NULL
      dataspace$Number_0_group4 <- NULL
       dataspace$Number_0_all_groups <- NULL
    }
    if (groups_number >=5){

      #Count the number of 0 for group 5
      for (i in c(1:length(dataspace[,1]))){
        a <- table(dataspace[i,coln5]==0)["TRUE"]
        if(is.na(a)){
          dataspace[i,"Number_0_group5"] <- 0
        } else{
          dataspace[i,"Number_0_group5"] <- table(dataspace[i,coln5]==0)["TRUE"]
        }
      }

      dataspace$Number_0_all_groups <- dataspace$Number_0_group1 + dataspace$Number_0_group2 +dataspace$Number_0_group3 +dataspace$Number_0_group4 +dataspace$Number_0_group5
    }
    if (groups_number == 5){

       openxlsx::write.xlsx(dataspace, file = "Dataset_before_threshold.xlsx")

      setwd<-path_res
      if (global_threshold == TRUE) {
        dataspace <- dataspace[dataspace$Number_0_all_groups<threshold,]
        openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")}

      if (global_threshold == FALSE) {

      #Delete row with x or more zeros
      dataspace <- dataspace[dataspace$Number_0_group1<threshold[1] | dataspace$Number_0_group2<threshold[2]| dataspace$Number_0_group3<threshold[3] | dataspace$Number_0_group4<threshold[4] | dataspace$Number_0_group5<threshold[5],]

      #write.table(dataspace, file=paste(path_res,"Dataset_threshold_applied.xlsx",sep=""), dec=".",sep="\t", row.names=FALSE)
      openxlsx::write.xlsx(dataspace, file = "Dataset_threshold_applied.xlsx")}
      message("An excel file with the proteins that have % of missing values below the threshold was created as Dataset_threshold_applied.xlsx")
      dataspace_0s<- dataspace

      dataspace$Number_0_group1 <- NULL
      dataspace$Number_0_group2 <- NULL
      dataspace$Number_0_group3 <- NULL
      dataspace$Number_0_group4 <- NULL
      dataspace$Number_0_group5 <- NULL
       dataspace$Number_0_all_groups <- NULL
    }

    if (groups_number>=6){

          #Count the number of 0 for group 6
      for (i in c(1:length(dataspace[,1]))){
        a <- table(dataspace[i,coln6]==0)["TRUE"]
        if(is.na(a)){
          dataspace[i,"Number_0_group6"] <- 0
        } else{
          dataspace[i,"Number_0_group6"] <- table(dataspace[i,coln6]==0)["TRUE"]
        }
      }
      dataspace$Number_0_all_groups <- dataspace$Number_0_group1 + dataspace$Number_0_group2 +dataspace$Number_0_group3 +dataspace$Number_0_group4 +dataspace$Number_0_group5 +dataspace$Number_0_group6
    }
    if (groups_number == 6){

    openxlsx::write.xlsx(dataspace, file = "Dataset_before_threshold.xlsx")

      setwd<-path_res
      if (global_threshold == TRUE) {
        dataspace <- dataspace[dataspace$Number_0_all_groups<threshold,]
        openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")}

      if (global_threshold == FALSE) { #Delete row with x or more zeros
      dataspace <- dataspace[dataspace$Number_0_group1<threshold[1] | dataspace$Number_0_group2<threshold[2]| dataspace$Number_0_group3<threshold[3] | dataspace$Number_0_group4<threshold[4] | dataspace$Number_0_group5<threshold[5] | dataspace$Number_0_group6<threshold[6],]

      #write.table(dataspace, file=paste(path_res,"Dataset_threshold_applied.xlsx",sep=""), dec=".",sep="\t", row.names=FALSE)
      openxlsx::write.xlsx(dataspace, file = "Dataset_threshold_applied.xlsx")}
      message("An excel file with the proteins that have % of missing values below the threshold was created as Dataset_threshold_applied.xlsx")
      dataspace_0s<- dataspace

      dataspace$Number_0_group1 <- NULL
      dataspace$Number_0_group2 <- NULL
      dataspace$Number_0_group3 <- NULL
      dataspace$Number_0_group4 <- NULL
      dataspace$Number_0_group5 <- NULL
      dataspace$Number_0_group6 <- NULL
       dataspace$Number_0_all_groups <- NULL
    }

    if (groups_number>=7){

      #Count the number of 0 for group 7
      for (i in c(1:length(dataspace[,1]))){
        a <- table(dataspace[i,coln7]==0)["TRUE"]
        if(is.na(a)){
          dataspace[i,"Number_0_group7"] <- 0
        } else{
          dataspace[i,"Number_0_group7"] <- table(dataspace[i,coln7]==0)["TRUE"]
        }
      }
      dataspace$Number_0_all_groups <- dataspace$Number_0_group1 + dataspace$Number_0_group2 +dataspace$Number_0_group3 +dataspace$Number_0_group4 +dataspace$Number_0_group5 +dataspace$Number_0_group6 +dataspace$Number_0_group7
    }
    if (groups_number == 7){
      openxlsx::write.xlsx(dataspace, file = "Dataset_before_threshold.xlsx")

      setwd<-path_res
      if (global_threshold == TRUE) {
        dataspace <- dataspace[dataspace$Number_0_all_groups<threshold,]
        openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")}

      if (global_threshold == FALSE) { #Delete row with x or more zeros
      dataspace <- dataspace[dataspace$Number_0_group1<threshold[1] | dataspace$Number_0_group2<threshold[2]| dataspace$Number_0_group3<threshold[3] | dataspace$Number_0_group4<threshold[4] | dataspace$Number_0_group5<threshold[5] | dataspace$Number_0_group6<threshold[6] | dataspace$Number_0_group7<threshold[7],]

      #write.table(dataspace, file=paste(path_res,"Dataset_threshold_applied.xlsx",sep=""), dec=".",sep="\t", row.names=FALSE)
      openxlsx::write.xlsx(dataspace, file = "Dataset_threshold_applied.xlsx")}
      message("An excel file with the proteins that have % of missing values below the threshold was created as Dataset_threshold_applied.xlsx")
      dataspace_0s<- dataspace

      dataspace$Number_0_group1 <- NULL
      dataspace$Number_0_group2 <- NULL
      dataspace$Number_0_group3 <- NULL
      dataspace$Number_0_group4 <- NULL
      dataspace$Number_0_group5 <- NULL
      dataspace$Number_0_group6 <- NULL
      dataspace$Number_0_group7 <- NULL
       dataspace$Number_0_all_groups <- NULL
    }

    if (groups_number>=8){

      #Count the number of 0 for group 8
      for (i in c(1:length(dataspace[,1]))){
        a <- table(dataspace[i,coln8]==0)["TRUE"]
        if(is.na(a)){
          dataspace[i,"Number_0_group8"] <- 0
        } else{
          dataspace[i,"Number_0_group8"] <- table(dataspace[i,coln8]==0)["TRUE"]
        }
      }
      dataspace$Number_0_all_groups <- dataspace$Number_0_group1 + dataspace$Number_0_group2 +dataspace$Number_0_group3 +dataspace$Number_0_group4 +dataspace$Number_0_group5 +dataspace$Number_0_group6 +dataspace$Number_0_group7 +dataspace$Number_0_group8
    }
    if (groups_number == 8){

      openxlsx::write.xlsx(dataspace, file = "Dataset_before_threshold.xlsx")

      setwd<-path_res
      if (global_threshold == TRUE) {
        dataspace <- dataspace[dataspace$Number_0_all_groups<threshold,]
        openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")}

      if (global_threshold == FALSE) {  #Delete row with x or more zeros
      dataspace <- dataspace[dataspace$Number_0_group1<threshold[1] | dataspace$Number_0_group2<threshold[2]| dataspace$Number_0_group3<threshold[3] | dataspace$Number_0_group4<threshold[4] | dataspace$Number_0_group5<threshold[5] | dataspace$Number_0_group6<threshold[6] | dataspace$Number_0_group7<threshold[7] | dataspace$Number_0_group8<threshold[8],]


      #write.table(dataspace, file=paste(path_res,"Dataset_threshold_applied.xlsx",sep=""), dec=".",sep="\t", row.names=FALSE)
      openxlsx::write.xlsx(dataspace, file = "Dataset_threshold_applied.xlsx")}
      message("An excel file with the proteins that have % of missing values below the threshold was created as Dataset_threshold_applied.xlsx")
      dataspace_0s<- dataspace

      dataspace$Number_0_group1 <- NULL
      dataspace$Number_0_group2 <- NULL
      dataspace$Number_0_group3 <- NULL
      dataspace$Number_0_group4 <- NULL
      dataspace$Number_0_group5 <- NULL
      dataspace$Number_0_group6 <- NULL
      dataspace$Number_0_group7 <- NULL
      dataspace$Number_0_group8 <- NULL
       dataspace$Number_0_all_groups <- NULL
    }

    if (groups_number==9){
  #Count the number of 0 for group 9
      for (i in c(1:length(dataspace[,1]))){
        a <- table(dataspace[i,coln9]==0)["TRUE"]
        if(is.na(a)){
          dataspace[i,"Number_0_group9"] <- 0
        } else{
          dataspace[i,"Number_0_group9"] <- table(dataspace[i,coln9]==0)["TRUE"]
        }
      }

    dataspace$Number_0_all_groups <- dataspace$Number_0_group1 + dataspace$Number_0_group2 +dataspace$Number_0_group3 +dataspace$Number_0_group4 +dataspace$Number_0_group5 +dataspace$Number_0_group6 +dataspace$Number_0_group7 +dataspace$Number_0_group8 +dataspace$Number_0_group9
    openxlsx::write.xlsx(dataspace, file = "Dataset_before_threshold.xlsx")

    setwd<-path_res
    if (global_threshold == TRUE) {
      dataspace <- dataspace[dataspace$Number_0_all_groups<threshold,]
      dataspace_0s<- dataspace
      openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")}

    if (global_threshold == FALSE) {  #Delete row with x or more zeros
      dataspace <- dataspace[dataspace$Number_0_group1<threshold[1] | dataspace$Number_0_group2<threshold[2]| dataspace$Number_0_group3<threshold[3] | dataspace$Number_0_group4<threshold[4] | dataspace$Number_0_group5<threshold[5] | dataspace$Number_0_group6<threshold[6] | dataspace$Number_0_group7<threshold[7] | dataspace$Number_0_group8<threshold[8] | dataspace$Number_0_group9<threshold[9],]
     dataspace_0s<- dataspace
      #write.table(dataspace, file=paste(path_res,"Dataset_threshold_applied.xlsx",sep=""), dec=".",sep="\t", row.names=FALSE)
      openxlsx::write.xlsx(dataspace, file = "Dataset_threshold_applied.xlsx")}
      message("An excel file with the proteins that have % of missing values below the threshold was created as Dataset_threshold_applied.xlsx")
      dataspace_0s<- dataspace

      dataspace$Number_0_group1 <- NULL
      dataspace$Number_0_group2 <- NULL
      dataspace$Number_0_group3 <- NULL
      dataspace$Number_0_group4 <- NULL
      dataspace$Number_0_group5 <- NULL
      dataspace$Number_0_group6 <- NULL
      dataspace$Number_0_group7 <- NULL
      dataspace$Number_0_group8 <- NULL
      dataspace$Number_0_group9 <- NULL
       dataspace$Number_0_all_groups <- NULL
    }
    zero_per_sample1 <- colSums(dataspace[,-1:-2] == 0)*100/nrow(dataspace)
    sample_names <- colnames(dataspace[,-1:-2])
    qc <- cbind(sample_names,zero_per_sample,zero_per_sample1,IDs)
    colnames(qc) <- c("Sample Name","% of Missing values before filtering","% of Missing values after filtering","Number of proteins detected in the sample")
    rownames(qc) <- NULL
    openxlsx::write.xlsx(qc,file = "Quality_check.xlsx")

pre_dataspace <- dataspace

##imputation KNN
if (imputation == TRUE) {
dataspace[dataspace==0] <- NA
dataspace <- VIM::kNN(dataspace, imp_var = FALSE, k= 5)
openxlsx::write.xlsx(dataspace,file = "Dataset_Imputed.xlsx")

#create histogramm for imputed values
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

ggplot2::ggsave("Imputed_values_histogram.tiff", plot = imp_hist, device = "tiff", path = path_res,
                scale = 1, width = 5, height = 4, units = "in",
                dpi = 300, limitsize = TRUE)

message("An excel with the imputed missing values was created as Dataset_Imputed.xlsx and a histogram documentating these values")

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
  #theme_minimal(base_size = 15)  +
  theme_linedraw()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_line(color = "grey80"),  # Make grids more visible
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9))
abund.plot

ggplot2::ggsave("Proteins_abundance_rank.tiff", plot = abund.plot , device = "tiff", path = path_res,
                scale = 1, width = 12, height = 5, units = "in",
                dpi = 300, limitsize = TRUE, bg = "white")

}
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

ggplot2::ggsave("Proteins_abundance_rank.tiff", plot = abund.plot , device = "tiff", path = path_res,
                scale = 1, width = 12, height = 5, units = "in",
                dpi = 300, limitsize = TRUE, bg = "white")
}



dataspace$percentage <- NULL
dataspace$mean <- NULL
dataspace$log<- NULL
dataspace$rank <- NULL
groups_for_test<-NULL

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
colnames(anova_res)<-paste("ANOVA",colnames(anova_res), sep = "_")}

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
openxlsx::write.xlsx(limma_dataspace, file = "Dataset_limma.t-test.xlsx")
message("limma test was created")
#####
###- Mann-Whitney and Kruskal-Wallis starts here! - ###
if (MWtest != "Paired" && MWtest != "Independent"){stop("Error. You need to assign MWtest = 'Paired' or 'Indepedent'")}
### 1 ### Specify file for statistical analysis
data2 <- dataspace

if (groups_number>=2){
  data2$Average_G1 <- rowMeans(data2[,coln])
  data2$Average_G2 <- rowMeans(data2[,coln2])
  data2$St_Dv_G1 <- apply(data2[,coln], 1, sd)
  data2$St_Dv_G2 <- apply(data2[,coln2], 1, sd)
  ### Calculate unadjusted p-value
 if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G2vsG1"]<-test_list[[3]]
  }  }
   if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=TRUE)
    data2[i,"MW_G2vsG1"]<-test_list[[3]]
  }}


  #### adjust the p-values
  data2$BH_p_G2vsG1 <- stats::p.adjust(data2$MW_G2vsG1, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G2vsG1 <- (data2$Average_G2)/(data2$Average_G1)
  data2$Log2_Ratio.G2vsG1 <- log2(data2$Ratio_G2vsG1)
  #data$subtraction <- (data$average_case)-(data$average_control)

}

if (groups_number>=3){

  data2$Average_G3 <- rowMeans(data2[,coln3])
  data2$St_Dv_G3 <- apply(data2[,coln3], 1, sd)
  ### Calculate unadjusted p-value

      if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
        test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
        data2[i,"MW_G3vsG1"]<-test_list[[3]]
      }  }
        if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
      test_list<-wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln3]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G3vsG1"]<-test_list[[3]]}
    }

  #### adjust the p-values
  data2$BH_p_G3vsG1 <- p.adjust(data2$MW_G3vsG1, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G3vsG1 <- (data2$Average_G3)/(data2$Average_G1)
  data2$Log2_Ratio.G3vsG1 <- log2(data2$Ratio_G3vsG1)
  #data$subtraction <- (data$average_case)-(data$average_control)

    if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
        test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
        data2[i,"MW_G3vsG2"]<-test_list[[3]]
      }  }
        if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
      test_list<-wilcox.test(as.numeric(data2[i,coln2]),as.numeric(data2[i,coln3]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G3vsG2"]<-test_list[[3]]}
    }

  #### adjust the p-values
  data2$BH_p_G3vsG2 <- p.adjust(data2$MW_G3vsG2, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G3vsG2 <- (data2$Average_G3)/(data2$Average_G2)
  data2$Log2_Ratio.G3vsG2 <- log2(data2$Ratio_G3vsG2)
  #data$subtraction <- (data$average_case)-(data$average_control)

}

if (groups_number>=4){
  data2$Average_G4 <- rowMeans(data2[,coln4])
  data2$St_Dv_G4 <- apply(data2[,coln4], 1, sd)

    if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
        test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
        data2[i,"MW_G4vsG1"]<-test_list[[3]]
      }  }
        if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
      test_list<-wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln4]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G4vsG1"]<-test_list[[3]]}
    }

  #### adjust the p-values
  data2$BH_p_G4vsG1 <- p.adjust(data2$MW_G4vsG1, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G4vsG1 <- (data2$Average_G4)/(data2$Average_G1)
  data2$Log2_Ratio.G4vsG1 <- log2(data2$Ratio_G4vsG1)
  #data$subtraction <- (data$average_case)-(data$average_control)


    if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
        test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
        data2[i,"MW_G4vsG2"]<-test_list[[3]]
      }  }
        if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
      test_list<-wilcox.test(as.numeric(data2[i,coln2]),as.numeric(data2[i,coln4]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G4vsG2"]<-test_list[[3]]}
    }

  #### adjust the p-values
  data2$BH_p_G4vsG2 <- p.adjust(data2$MW_G4vsG2, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G4vsG2 <- (data2$Average_G4)/(data2$Average_G2)
  data2$Log2_Ratio.G4vsG2 <- log2(data2$Ratio_G4vsG2)
  #data$subtraction <- (data$average_case)-(data$average_control)

      if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
        test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
        data2[i,"MW_G4vsG3"]<-test_list[[3]]
      }  }
        if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
      test_list<-wilcox.test(as.numeric(data2[i,coln3]),as.numeric(data2[i,coln4]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G4vsG3"]<-test_list[[3]]}
    }

  #### adjust the p-values
  data2$BH_p_G4vsG3 <- p.adjust(data2$MW_G4vsG3, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G4vsG3 <- (data2$Average_G4)/(data2$Average_G3)
  data2$Log2_Ratio.G4vsG3 <- log2(data2$Ratio_G4vsG3)
  #data$subtraction <- (data$average_case)-(data$average_control)

}

if (groups_number>=5){
  data2$Average_G5 <- rowMeans(data2[,coln5])
  data2$St_Dv_G5 <- apply(data2[,coln5], 1, sd)

   if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
        test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
        data2[i,"MW_G5vsG1"]<-test_list[[3]]
      }  }
        if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
      test_list<-wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln5]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G5vsG1"]<-test_list[[3]]}
    }

  #### adjust the p-values
  data2$BH_p_G5vsG1 <- p.adjust(data2$MW_G5vsG1, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G5vsG1 <- (data2$Average_G5)/(data2$Average_G1)
  data2$Log2_Ratio.G5vsG1 <- log2(data2$Ratio_G5vsG1)
  #data$subtraction <- (data$average_case)-(data$average_control)

   if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
        test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
        data2[i,"MW_G5vsG2"]<-test_list[[3]]
      }  }
        if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
      test_list<-wilcox.test(as.numeric(data2[i,coln2]),as.numeric(data2[i,coln5]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G5vsG2"]<-test_list[[3]]}
    }

  #### adjust the p-values
  data2$BH_p_G5vsG2 <- p.adjust(data2$MW_G5vsG2, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G5vsG2 <- (data2$Average_G5)/(data2$Average_G2)
  data2$Log2_Ratio.G5vsG2 <- log2(data2$Ratio_G5vsG2)
  #data$subtraction <- (data$average_case)-(data$average_control)


   if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
        test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
        data2[i,"MW_G5vsG3"]<-test_list[[3]]
      }  }
        if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
      test_list<-wilcox.test(as.numeric(data2[i,coln3]),as.numeric(data2[i,coln5]), exact=FALSE, paired=TRUE,)
        data2[i,"MW_G5vsG3"]<-test_list[[3]]}
    }

  #### adjust the p-values
  data2$BH_p_G5vsG3 <- p.adjust(data2$MW_G5vsG3, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G5vsG3 <- (data2$Average_G5)/(data2$Average_G3)
  data2$Log2_Ratio.G5vsG3 <- log2(data2$Ratio_G5vsG3)
  #data$subtraction <- (data$average_case)-(data$average_control)

   if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
        test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
        data2[i,"MW_G5vsG4"]<-test_list[[3]]
      } }
        if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
      test_list<-wilcox.test(as.numeric(data2[i,coln4]),as.numeric(data2[i,coln5]), exact=FALSE, paired=TRUE,)
        data2[i,"MW_G5vsG4"]<-test_list[[3]]}
    }

  #### adjust the p-values
  data2$BH_p_G5vsG4 <- p.adjust(data2$MW_G5vsG4, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G5vsG4 <- (data2$Average_G5)/(data2$Average_G4)
  data2$Log2_Ratio.G5vsG4 <- log2(data2$Ratio_G5vsG4)
  #data$subtraction <- (data$average_case)-(data$average_control)

}

if (groups_number>=6){
  data2$Average_G6 <- rowMeans(data2[,coln6])
  data2$St_Dv_G6 <- apply(data2[,coln6], 1, sd)

    if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
        test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
        data2[i,"MW_G6vsG1"]<-test_list[[3]]
      }  }
        if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
      test_list<-wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln6]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G6vsG1"]<-test_list[[3]]}
    }

  #### adjust the p-values
  data2$BH_p_G6vsG1 <- p.adjust(data2$MW_G6vsG1, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G6vsG1 <- (data2$Average_G6)/(data2$Average_G1)
  data2$Log2_Ratio.G6vsG1 <- log2(data2$Ratio_G6vsG1)
  #data$subtraction <- (data$average_case)-(data$average_control)

  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G6vsG2"]<-test_list[[3]]
  }  }
    if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
      test_list<-wilcox.test(as.numeric(data2[i,coln2]),as.numeric(data2[i,coln6]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G6vsG2"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G6vsG2 <- p.adjust(data2$MW_G6vsG2, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G6vsG2 <- (data2$Average_G6)/(data2$Average_G2)
  data2$Log2_Ratio.G6vsG2 <- log2(data2$Ratio_G6vsG2)
  #data$subtraction <- (data$average_case)-(data$average_control)


  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G6vsG3"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){  test_list<-wilcox.test(as.numeric(data2[i,coln3]),as.numeric(data2[i,coln6]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G6vsG3"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G6vsG3 <- p.adjust(data2$MW_G6vsG3, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G6vsG3 <- (data2$Average_G6)/(data2$Average_G3)
  data2$Log2_Ratio.G6vsG3 <- log2(data2$Ratio_G6vsG3)
  #data$subtraction <- (data$average_case)-(data$average_control)

  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G6vsG4"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){  test_list<-wilcox.test(as.numeric(data2[i,coln4]),as.numeric(data2[i,coln6]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G6vsG4"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G6vsG4 <- p.adjust(data2$MW_G6vsG4, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G6vsG4 <- (data2$Average_G6)/(data2$Average_G4)
  data2$Log2_Ratio.G6vsG4 <- log2(data2$Ratio_G6vsG4)
  #data$subtraction <- (data$average_case)-(data$average_control)

  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G6vsG5"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){ test_list<-wilcox.test(as.numeric(data2[i,coln5]),as.numeric(data2[i,coln6]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G6vsG5"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G6vsG5 <- p.adjust(data2$MW_G6vsG5, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G6vsG5 <- (data2$Average_G6)/(data2$Average_G5)
  data2$Log2_Ratio.G6vsG5 <- log2(data2$Ratio_G6vsG5)
  #data$subtraction <- (data$average_case)-(data$average_control)

}

if (groups_number>=7){
  data2$Average_G7 <- rowMeans(data2[,coln7])
  data2$St_Dv_G7 <- apply(data2[,coln7], 1, sd)

  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G7vsG1"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){ test_list<-wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln7]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G7vsG1"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G7vsG1 <- p.adjust(data2$MW_G7vsG1, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G7vsG1 <- (data2$Average_G7)/(data2$Average_G1)
  data2$Log2_Ratio.G7vsG1 <- log2(data2$Ratio_G7vsG1)
  #data$subtraction <- (data$average_case)-(data$average_control)


  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G7vsG2"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){ test_list<-wilcox.test(as.numeric(data2[i,coln2]),as.numeric(data2[i,coln7]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G7vsG2"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G7vsG2 <- p.adjust(data2$MW_G7vsG2, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G7vsG2 <- (data2$Average_G7)/(data2$Average_G2)
  data2$Log2_Ratio.G7vsG2 <- log2(data2$Ratio_G7vsG2)
  #data$subtraction <- (data$average_case)-(data$average_control)


  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G7vsG3"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){  test_list<-wilcox.test(as.numeric(data2[i,coln3]),as.numeric(data2[i,coln7]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G7vsG3"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G7vsG3 <- p.adjust(data2$MW_G7vsG3, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G7vsG3 <- (data2$Average_G7)/(data2$Average_G3)
  data2$Log2_Ratio.G7vsG3 <- log2(data2$Ratio_G7vsG3)
  #data$subtraction <- (data$average_case)-(data$average_control)

  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G7vsG4"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){  test_list<-wilcox.test(as.numeric(data2[i,coln4]),as.numeric(data2[i,coln7]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G7vsG4"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G7vsG4 <- p.adjust(data2$MW_G7vsG4, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G7vsG4 <- (data2$Average_G7)/(data2$Average_G4)
  data2$Log2_Ratio.G7vsG4 <- log2(data2$Ratio_G7vsG4)
  #data$subtraction <- (data$average_case)-(data$average_control)
   if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G7vsG5"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){  test_list<-wilcox.test(as.numeric(data2[i,coln5]),as.numeric(data2[i,coln7]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G7vsG5"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G7vsG5 <- p.adjust(data2$MW_G7vsG5, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G7vsG5 <- (data2$Average_G7)/(data2$Average_G5)
  data2$Log2_Ratio.G7vsG5 <- log2(data2$Ratio_G7vsG5)
  #data$subtraction <- (data$average_case)-(data$average_control)
  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G7vsG6"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){ test_list<-wilcox.test(as.numeric(data2[i,coln6]),as.numeric(data2[i,coln7]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G7vsG6"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G7vsG6 <- p.adjust(data2$MW_G7vsG6, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G7vsG6 <- (data2$Average_G7)/(data2$Average_G6)
  data2$Log2_Ratio.G7vsG6 <- log2(data2$Ratio_G7vsG6)
  #data$subtraction <- (data$average_case)-(data$average_control)

 }

if (groups_number>=8){
  data2$Average_G8 <- rowMeans(data2[,coln8])
  data2$St_Dv_G8 <- apply(data2[,coln8], 1, sd)
  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G8vsG1"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){  test_list<-wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln8]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G8vsG1"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G8vsG1 <- p.adjust(data2$MW_G8vsG1, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G8vsG1 <- (data2$Average_G8)/(data2$Average_G1)
  data2$Log2_Ratio.G8vsG1 <- log2(data2$Ratio_G8vsG1)
  #data$subtraction <- (data$average_case)-(data$average_control)

 if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G8vsG2"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){
      test_list<-wilcox.test(as.numeric(data2[i,coln2]),as.numeric(data2[i,coln8]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G8vsG2"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G8vsG2 <- p.adjust(data2$MW_G8vsG2, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G8vsG2 <- (data2$Average_G8)/(data2$Average_G2)
  data2$Log2_Ratio.G8vsG2 <- log2(data2$Ratio_G8vsG2)
  #data$subtraction <- (data$average_case)-(data$average_control)


  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G8vsG3"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){ test_list<-wilcox.test(as.numeric(data2[i,coln3]),as.numeric(data2[i,coln8]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G8vsG3"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G8vsG3 <- p.adjust(data2$MW_G8vsG3, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G8vsG3 <- (data2$Average_G8)/(data2$Average_G3)
  data2$Log2_Ratio.G8vsG3 <- log2(data2$Ratio_G8vsG3)
  #data$subtraction <- (data$average_case)-(data$average_control)

 if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G8vsG4"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){    test_list<-wilcox.test(as.numeric(data2[i,coln4]),as.numeric(data2[i,coln8]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G8vsG4"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G8vsG4 <- p.adjust(data2$MW_G8vsG4, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G8vsG4 <- (data2$Average_G8)/(data2$Average_G4)
  data2$Log2_Ratio.G8vsG4 <- log2(data2$Ratio_G8vsG4)
  #data$subtraction <- (data$average_case)-(data$average_control)
  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G8vsG5"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){   test_list<-wilcox.test(as.numeric(data2[i,coln5]),as.numeric(data2[i,coln8]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G8vsG5"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G8vsG5 <- p.adjust(data2$MW_G8vsG5, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G8vsG5 <- (data2$Average_G8)/(data2$Average_G5)
  data2$Log2_Ratio.G8vsG5 <- log2(data2$Ratio_G8vsG5)
  #data$subtraction <- (data$average_case)-(data$average_control)
 if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G8vsG6"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){   test_list<-wilcox.test(as.numeric(data2[i,coln6]),as.numeric(data2[i,coln8]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G8vsG6"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G8vsG6 <- p.adjust(data2$MW_G8vsG6, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G8vsG6 <- (data2$Average_G8)/(data2$Average_G6)
  data2$Log2_Ratio.G8vsG6 <- log2(data2$Ratio_G8vsG6)
  #data$subtraction <- (data$average_case)-(data$average_control)

  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G8vsG7"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){   test_list<-wilcox.test(as.numeric(data2[i,coln7]),as.numeric(data2[i,coln8]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G8vsG7"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G8vsG7 <- p.adjust(data2$MW_G8vsG7, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G8vsG7 <- (data2$Average_G8)/(data2$Average_G7)
  data2$Log2_Ratio.G8vsG7 <- log2(data2$Ratio_G8vsG7)
  #data$subtraction <- (data$average_case)-(data$average_control)

 }


if (groups_number==9){
  data2$Average_G9 <- rowMeans(data2[,coln9])
  data2$St_Dv_G9 <- apply(data2[,coln9], 1, sd)
   if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G9vsG1"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){   test_list<-wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln9]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G9vsG1"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G9vsG1 <- p.adjust(data2$MW_G9vsG1, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G9vsG1 <- (data2$Average_G9)/(data2$Average_G1)
  data2$Log2_Ratio.G9vsG1 <- log2(data2$Ratio_G9vsG1)
  #data$subtraction <- (data$average_case)-(data$average_control)
 if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G9vsG2"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){test_list<-wilcox.test(as.numeric(data2[i,coln2]),as.numeric(data2[i,coln9]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G9vsG2"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G9vsG2 <- p.adjust(data2$MW_G9vsG2, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G9vsG2 <- (data2$Average_G9)/(data2$Average_G2)
  data2$Log2_Ratio.G9vsG2 <- log2(data2$Ratio_G9vsG2)
  #data$subtraction <- (data$average_case)-(data$average_control)
  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G9vsG3"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){ test_list<-wilcox.test(as.numeric(data2[i,coln3]),as.numeric(data2[i,coln9]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G9vsG3"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G9vsG3 <- p.adjust(data2$MW_G9vsG3, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G9vsG3 <- (data2$Average_G9)/(data2$Average_G3)
  data2$Log2_Ratio.G9vsG3 <- log2(data2$Ratio_G9vsG3)
  #data$subtraction <- (data$average_case)-(data$average_control)

  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G9vsG4"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){  test_list<-wilcox.test(as.numeric(data2[i,coln4]),as.numeric(data2[i,coln9]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G9vsG4"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G9vsG4 <- p.adjust(data2$MW_G9vsG4, method = "BH")
  #### calculate the ratio, use the subtraction (insted of ratio) only when the values of the data are zero centered
  data2$Ratio_G9vsG4 <- (data2$Average_G9)/(data2$Average_G4)
  data2$Log2_Ratio.G9vsG4 <- log2(data2$Ratio_G9vsG4)
  #data$subtraction <- (data$average_case)-(data$average_control)
  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G9vsG5"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){  test_list<-wilcox.test(as.numeric(data2[i,coln5]),as.numeric(data2[i,coln9]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G9vsG5"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G9vsG5 <- p.adjust(data2$MW_G9vsG5, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G9vsG5 <- (data2$Average_G9)/(data2$Average_G5)
  data2$Log2_Ratio.G9vsG5 <- log2(data2$Ratio_G9vsG5)
  #data$subtraction <- (data$average_case)-(data$average_control)


  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G9vsG6"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){  test_list<-wilcox.test(as.numeric(data2[i,coln6]),as.numeric(data2[i,coln9]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G9vsG6"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G9vsG6 <- p.adjust(data2$MW_G9vsG6, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G9vsG6 <- (data2$Average_G9)/(data2$Average_G6)
  data2$Log2_Ratio.G9vsG6 <- log2(data2$Ratio_G9vsG6)
  #data$subtraction <- (data$average_case)-(data$average_control)

 if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G9vsG7"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){    test_list<-wilcox.test(as.numeric(data2[i,coln7]),as.numeric(data2[i,coln9]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G9vsG7"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G9vsG7 <- p.adjust(data2$MW_G9vsG7, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G9vsG7 <- (data2$Average_G9)/(data2$Average_G7)
  data2$Log2_Ratio.G9vsG7 <- log2(data2$Ratio_G9vsG7)
  #data$subtraction <- (data$average_case)-(data$average_control)

  if (MWtest == "Independent"){ for (i in c(1:length(data2[,1]))){
    test_list<-stats::wilcox.test(as.numeric(data2[i,coln]),as.numeric(data2[i,coln2]), exact=FALSE, paired=FALSE)
    data2[i,"MW_G9vsG8"]<-test_list[[3]]
  }  }
  if (MWtest == "Paired"){ for (i in c(1:length(data2[,1]))){   test_list<-wilcox.test(as.numeric(data2[i,coln8]),as.numeric(data2[i,coln9]), exact=FALSE, paired=TRUE,)
      data2[i,"MW_G9vsG8"]<-test_list[[3]]
    }
  }
  #### adjust the p-values
  data2$BH_p_G9vsG8 <- p.adjust(data2$MW_G9vsG8, method = "BH")
  #### calculate the ratio, use the subtraction (instead of ratio) only when the values of the data are zero centered
  data2$Ratio_G9vsG8 <- (data2$Average_G9)/(data2$Average_G8)
  data2$Log2_Ratio.G9vsG8 <- log2(data2$Ratio_G9vsG8)
  #data$subtraction <- (data$average_case)-(data$average_control)

}

  Ddataspace<-data2
  Ddataspace$Symbol = sub(".*GN=(.*?) .*","\\1",Ddataspace$Description)
  Ddataspace$Symbol[Ddataspace$Symbol==Ddataspace$Description] = "Not available"
  Ddataspace<-Ddataspace %>%
    dplyr::select(Accession, Description, Symbol, everything())


  #Ddataspace <- replace(Ddataspace, is.nan(Ddataspace), "Not Applicable")


  Fdataspace<-Ddataspace


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

if (groups_number >= 2){
  namesc<-colnames(Fdataspace)
  namesc<-gsub("G1", g1.name, namesc)
  namesc<-gsub("G2", g2.name, namesc)
  colnames(Fdataspace)<-namesc
}

if (groups_number >= 3){
  namesc<-colnames(Fdataspace)
  namesc<-gsub("G3", g3.name, namesc)
  colnames(Fdataspace)<-namesc
}

if (groups_number >= 4){
  namesc<-colnames(Fdataspace)
  namesc<-gsub("G4", g4.name, namesc)
  colnames(Fdataspace)<-namesc
}

if (groups_number >= 5){
  namesc<-colnames(Fdataspace)
  namesc<-gsub("G5", g5.name, namesc)
  colnames(Fdataspace)<-namesc
}

if (groups_number >= 6){
  namesc<-colnames(Fdataspace)
  namesc<-gsub("G6", g6.name, namesc)
  colnames(Fdataspace)<-namesc
}

if (groups_number >= 7){
  namesc<-colnames(Fdataspace)
  namesc<-gsub("G7", g7.name, namesc)
  colnames(Fdataspace)<-namesc
}


if (groups_number >= 8){
  namesc<-colnames(Fdataspace)
  namesc<-gsub("G8", g8.name, namesc)
  colnames(Fdataspace)<-namesc
}

if (groups_number >= 9){
  namesc<-colnames(Fdataspace)
  namesc<-gsub("G9", g9.name, namesc)
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

log.dataspace <- log(dataspace[,-c(1:2)]+1)

# PCA of the entire data

pca<-prcomp(t(log.dataspace), scale=TRUE, center=FALSE)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group = Group)

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

ggplot2::ggsave("PCA_plot_alldata.tiff", plot = pca.ent, device = "tiff", path = path_res,
       scale = 1, width = 5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)
message("PCA plot using all data was created as PCA_plot_alldata.tiff")
# PCA of the significant data. If number of groups = 2, the script uses the
# unadjusted Mann-Whitney test; else, it uses the unadjusted Kruskal-Wallis test.

which.sig<-vector()
if (groups_number != 2){
  which.sig <- which(Fdataspace$Kruskal_Wallis.pvalue < 0.05)
} else {(which.sig <- which(Ddataspace$MW_G2vsG1 < 0.05))}

which(Ddataspace$MW_G2vsG1< 0.05)
log.dataspace.sig <- log.dataspace[which.sig,]


pca<-prcomp(t(log.dataspace.sig), scale=TRUE, center=FALSE)

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
message("PCA plot with the significant data was created as PCA_plot_significant.tiff" )

#Prec.var<- barplot(pca.var.per, xlab= "Principal Component", ylab= "Principal Variation")

a<-ggpubr::ggarrange(pca.ent, pca.sig, nrow = 1, ncol=2,
             common.legend = TRUE, legend = "bottom")

#hhedhewd

ggplot2::ggsave("PCA_plots_combined.tiff", plot = a, device = "tiff", path = path_res,
       scale = 1, width = 8, height = 4.5, units = "in",
       dpi = 300, limitsize = TRUE)
message ("The 2 PCA plots are combined in PCA_plots_combined.tiff")
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

if (imputation == FALSE) {
qc.boxplots<-ggplot2::ggplot(melt.log.dataspace, aes(x=forcats::fct_inorder(variable), y=value, color=Group))+
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

}
melt.log.dataspace.na <- melt.log.dataspace

melt.log.dataspace.na$value[melt.log.dataspace.na$value == 0] <- NA

melt.log.dataspace.na$Group <- factor(melt.log.dataspace.na$Group, levels = Group2)

is.factor(melt.log.dataspace.na$variable)

qc.boxplots.na<-ggplot2::ggplot(melt.log.dataspace.na, aes(x=forcats::fct_inorder(variable), y=value, color=Group))+
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
if (imputation == FALSE) {
ggplot2::ggsave("QC_dataDistribution_NoZeros.tiff", plot = qc.boxplots.na, device = "tiff", path = path_res,
       scale = 1, width = 12, height = 5, units = "in",
       dpi = 300, limitsize = TRUE, bg = "white")
}
else
{
  ggplot2::ggsave("QC_dataDistribution.tiff", plot = qc.boxplots.na, device = "tiff", path = path_res,
                  scale = 1, width = 12, height = 5, units = "in",
                  dpi = 300, limitsize = TRUE, bg = "white")
}
message("The Boxplots for each sample have been created!!")



}
