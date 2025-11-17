library(spatstat)
library(readxl)
library(openxlsx)

savepath = "path/to/output"
data_folder = "path/to/input_data"
strct_dist_map = read.csv(sprintf('path/to/distance_map/%s/%s.csv',strct_name,strsplit(theDatalist[ff],".xlsx")[[1]]))
the_strct_list <- list.files(data_folder)
for (f in 1:length(the_strct_list)){
  strct_name = the_strct_list[f][1]
  subsavepath = sprintf('%s/%s',savepath,strct_name)
  if (!file.exists(subsavepath)) {
    # Directory or file does not exist, so you can create it
    dir.create(subsavepath, recursive = TRUE)  # For directories
  }
  theDatalist <- list.files(sprintf("%s/%s", data_folder,strct_name))
  for (ff in 1:length(theDatalist)){
    if (file.exists(sprintf('%s/%s',subsavepath,theDatalist[ff][1]))) {
      next
    }
    theData <- read_excel(sprintf('%s/%s/%s',data_folder,strct_name,theDatalist[ff][1]))
    u_data  <- theData[, c('Cluster Type',"x","y")]
    cluster_type = unique(theData['Cluster Type'])
    cov_data <- theData[, c(sprintf('Distance to %s',paste0(toupper(substr(strct_name, 1, 1)), substr(strct_name, 2, nchar(strct_name)))))]
    
    df <- data.frame(
      cluster_type = character(),
      distance_to = numeric(),
      stringsAsFactors = FALSE  # Avoids converting character vectors to factors
    )
    names(df)[2] <- sprintf('distance_to_%s', strct_name)
    
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
      if (strct_name == 'arteriole'){
        df <- rbind(df, data.frame(name = cluster_type[[1]][i], distance_to_arteriole = fit_strct$'coef'['covariate_im_strct'][[1]]))
      }
      if (strct_name == 'bone'){
        df <- rbind(df, data.frame(name = cluster_type[[1]][i], distance_to_bone = fit_strct$'coef'['covariate_im_strct'][[1]]))
      }
      if (strct_name == 'fat'){
        df <- rbind(df, data.frame(name = cluster_type[[1]][i], distance_to_fat = fit_strct$'coef'['covariate_im_strct'][[1]]))
      }
      if (strct_name == 'sinusoid'){
        df <- rbind(df, data.frame(name = cluster_type[[1]][i], distance_to_sinusoid = fit_strct$'coef'['covariate_im_strct'][[1]]))
      }
    }
    write.xlsx(df, file = sprintf("%s/%s",subsavepath, theDatalist[ff]))
    print(sprintf('%s: %d out of %d : %s done',strct_name,ff,length(theDatalist),theDatalist[ff]))
    
  }
}
