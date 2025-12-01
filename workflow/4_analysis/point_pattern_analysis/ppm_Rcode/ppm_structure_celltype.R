library(spatstat)
library(readxl)
library(openxlsx)
savepath = "path/to/struct_celltype"

main_datapath = "path/to/metadata/struct_celltype"
distancemap_path = 'path/to/distance_map'

struct_grp <- list("Bone", "Fat", "Arteriole", "Sinusoid")

for (g in 1:length(struct_grp)){
  theDatalist <- list.files(sprintf("%s/%s",distancemap_path,struct_grp[g]))
  subsavepath = sprintf('%s/%s',savepath,struct_grp[g])
  if (!file.exists(subsavepath)) {
    # Directory or file does not exist, so you can create it
    dir.create(subsavepath, recursive = TRUE)  # For directories
  }
  for (f in 1:length(theDatalist)){
  #for (f in 28:28){
    print(sprintf('%s: %d out of %d : %s started',struct_grp[g],f,length(theDatalist),theDatalist[f]))
    Datapath = sprintf("%s/%s.csv", main_datapath,strsplit(theDatalist[f],".csv"))
    theData <- read.csv(Datapath, check.names = FALSE)
    u_data  <- theData[, c('Cell Type',"x","y")]
    cell_type = unique(theData['Cell Type'])
    colname <- paste("Distance to", struct_grp[g])
    df <- data.frame(
      "Cell Type" = character(),
      stringsAsFactors = FALSE,  # Avoids converting character vectors to factors
      check.names = FALSE
    )
    df[[colname]] <- numeric()
    
    dist_map = read.csv(sprintf('%s/%s/%s.csv',distancemap_path,struct_grp[g],strsplit(theDatalist[f],".csv")[[1]]))
    covariate_im <- as.im(dist_map, c(as.numeric(max(dist_map$x)),as.numeric(max(dist_map$y)),1,1))
    
    
    for (i in 1:nrow(cell_type)){
      if (cell_type[[1]][i] == 'CD69'){
        next
      }
      indices <- which(theData$"Cell Type"==cell_type[[1]][i], arr.ind = TRUE)
      sub_u_data <- u_data[indices,]
      sub_u_data$x = (sub_u_data$x+1)/20
      sub_u_data$y = (sub_u_data$y+1)/20
      win <- owin(xrange = range(dist_map$x), yrange = range(dist_map$y))
      points_ppp <- ppp(x = sub_u_data$x, y = sub_u_data$y, window = win)
      fit_ <- ppm(points_ppp ~ covariate_im)
      
      if (!colname %in% names(df)) df[[colname]] <- numeric()
      
      new_row <- list()
      new_row[["Cell Type"]] <- cell_type[[1]][i]
      new_row[[colname]]     <- fit_$coef["covariate_im"][[1]]
      new_row_df <- as.data.frame(new_row, stringsAsFactors = FALSE, check.names = FALSE)
      
      df <- rbind(df, new_row_df)
      
      
    }
    write.csv(df, file = sprintf("%s/%s",subsavepath, theDatalist[f]),row.names = FALSE)
    print(sprintf('%d out of %d : %s done',f,length(theDatalist),theDatalist[f]))
  }
}
