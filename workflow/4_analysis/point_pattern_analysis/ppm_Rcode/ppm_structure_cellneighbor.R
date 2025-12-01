library(spatstat)
library(readxl)
library(openxlsx)

savepath = "path/to/struct_cellneighbor"
data_folder = "path/to/metadata/struct_cellneighbor"
distancemap_path = 'path/to/distance_map'
struct_grp <- list("Bone", "Fat", "Arteriole", "Sinusoid")

theDatalist <- list.files(data_folder)

for (f in 1:length(struct_grp)){
  strct_name = struct_grp[f]
  subsavepath = sprintf('%s/%s',savepath,strct_name)
  if (!file.exists(subsavepath)) {
    # Directory or file does not exist, so you can create it
    dir.create(subsavepath, recursive = TRUE)  # For directories
  }
  for (ff in 1:length(theDatalist)){
    if (file.exists(sprintf('%s/%s',subsavepath,theDatalist[ff][1]))) {
      next
    }
    theData <- read.csv(sprintf('%s/%s',data_folder,theDatalist[ff][1]), check.names = FALSE)
    u_data  <- theData[, c('Cluster Type',"x","y")]
    cluster_type = unique(theData['Cluster Type'])
    #cov_data <- theData[, c(sprintf('Distance to %s',paste0(toupper(substr(strct_name, 1, 1)), substr(strct_name, 2, nchar(strct_name)))))]
    
    colname <- paste("Distance to", struct_grp[f])
    
    df <- data.frame(
      "Cluster Type" = character(),
      stringsAsFactors = FALSE,  # Avoids converting character vectors to factors
      check.names = FALSE
    )
    
    df[[colname]] <- numeric()
    
    distmap_file = sprintf('%s/%s/%s.csv',distancemap_path,strct_name,strsplit(theDatalist[ff],".csv")[[1]])
    
    if (!file.exists(distmap_file)) {
      next
    }
    
    strct_dist_map = read.csv(distmap_file)
    covariate_im_strct <- as.im(strct_dist_map, c(as.numeric(max(strct_dist_map$x)),as.numeric(max(strct_dist_map$y)),1,1))
    
    for (i in 1:nrow(cluster_type)){
      indices <- which(theData$"Cluster Type"==cluster_type[[1]][i], arr.ind = TRUE)
      sub_u_data <- u_data[indices,]
      sub_u_data$x = (sub_u_data$x+1)/20
      sub_u_data$y = (sub_u_data$y+1)/20
      if (max(sub_u_data$x) > max(strct_dist_map$x)){
        browser()
      }
      if (max(sub_u_data$y) > max(strct_dist_map$y)){
        browser()
      }
      win <- owin(xrange = range(strct_dist_map$x), yrange = range(strct_dist_map$y))
      points_ppp <- ppp(x = sub_u_data$x, y = sub_u_data$y, window = win)
      fit_strct <- ppm(points_ppp ~ covariate_im_strct)
      
      new_row <- list()
      new_row[["Cluster Type"]] <- cluster_type[[1]][i]
      new_row[[colname]]     <- fit_strct$coef["covariate_im_strct"][[1]]
      new_row_df <- as.data.frame(new_row, stringsAsFactors = FALSE, check.names = FALSE)
      
      df <- rbind(df, new_row_df)
      
    write.csv(df, file = sprintf("%s/%s",subsavepath, theDatalist[ff]),row.names = FALSE)
    print(sprintf('%s: %d out of %d : %s done',strct_name,ff,length(theDatalist),theDatalist[ff]))
    
  }
  }
}
