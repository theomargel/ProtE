user_inputs <- function(path_for_group1, path_for_group2) {
g1.name <- basename(path_for_group1)
g2.name <- basename(path_for_group2)



   dataspace <- data.frame()
file_names_g1<-list.files(path=group1,pattern="*.xlsx")
setwd(group1)
for (i in 1:length(file_names_g1)) {
  file_case <- read_excel(paste(group1,file_names_g1[i],sep = "/"), sheet = 1)
  dataspace <- rbind(dataspace,file_case[,1:2])
}

# group 2 protein IDs
if (groups_number == 2 | groups_number == 3 | groups_number == 4 | groups_number==5 | groups_number==6 | groups_number==7 | groups_number==8 | groups_number==9){
  file_names_g2<-list.files(path=group2,pattern="*.xlsx")
  setwd(group2)
  for (i in 1:length(file_names_g2)) {
    file_case <- read_excel(paste(group2,file_names_g2[i],sep = "/"), sheet = 1)
    dataspace <- rbind(dataspace,file_case[,1:2])
  }
}
    path_res <- readline(prompt = "Specify folder where you want to save the results: e.g.: C:/Users/User/Documents/")
  if (!dir.exists(path_res)) {
    stop("The specified folder does not exist.")
  }
}

