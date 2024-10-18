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
user_inputs <- function(group1, group2)
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


       setwd<-path_res
      openxlsx::write.xlsx(dataspace, file = "Dataset_before_threshold.xlsx")

      #Delete row with x or more zeros
      dataspace <- dataspace[dataspace$Number_0_group1<threshold[1] | dataspace$Number_0_group2<threshold[2],]


      #write.table(dataspace, file=paste(path_res,"Dataset_threshold_applied.xlsx",sep=""), dec=".",sep="\t", row.names=FALSE)
      openxlsx::write.xlsx(dataspace,file = "Dataset_threshold_applied.xlsx")

}


