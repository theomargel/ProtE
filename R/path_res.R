#' Get data from 2 different sample groups from protein discoverer and create a dataspace wwith all the protein IDs found.
#'
#' Creates an a excel file
#'
#' @param group1  path to the folder where the samples from group 2 are
#' @param group2 path to the folder where the samples from group 2 are
#'
#'
#' @return the group name for group1
#'
#' @examples user_inputs(group1, group2)
#' @export
user_inputs <- function(group1, group2)
  {
g1.name <- basename(group1)
return(g1.name)
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
    path_res <- readline(prompt = "Specify folder where you want to save the results: e.g.: C:/Users/User/Documents/")
  if (!dir.exists(path_res)) {
    stop("The specified folder does not exist.")
  }
    setwd(path_res)
    openxlsx::write.xlsx(dataspace, file = "dataspace.xlsx")
    message("the excel was created")
}

