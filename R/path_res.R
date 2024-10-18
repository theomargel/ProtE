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
      ifelse(is.na(x), 0, (x / sum_x) * 10^6)  # NA=0 , normalize the rest
    })
    openxlsx::write.xlsx(dataspace, file = "normalized_list.xlsx")
    message("the excel of the normalized list was created")
}

