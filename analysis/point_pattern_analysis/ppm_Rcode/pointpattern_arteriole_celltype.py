library(spatstat)
library(readxl)
library(openxlsx)

data_path = "path/to/data" # From get_metadata_for_cell_to_artsinu.py
distmap_path = "path/to/distancemap"
savepath = "path/to/output"

if (!file.exists(savepath)) {
  # Directory or file does not exist, so you can create it
  dir.create(savepath, recursive = TRUE)  # For directories
}
theDatalist <- list.files(distmap_path)
for (f in 1:length(theDatalist)){
  print(sprintf('%d out of %d : %s started',f,length(theDatalist),theDatalist[f]))
  Datapath = sprintf("%s/%s.xlsx", data_path, strsplit(theDatalist[f],".csv"))
  theData <- read_excel(Datapath)
  u_data  <- theData[, c('Cell Type',"x","y")]
  cell_type = unique(theData['Cell Type'])
  df <- data.frame(
    cell_type = character(),
    # distance_to_arteriole = numeric(),
    distance_to_sinusoid = numeric(),
    stringsAsFactors = FALSE  # Avoids converting character vectors to factors
  )
  art_dist_map = read.csv(sprintf('%s/%s.csv',data_path,strsplit(theDatalist[f],".csv")[[1]]))
  covariate_im_art <- as.im(art_dist_map, c(as.numeric(max(art_dist_map$x)),as.numeric(max(art_dist_map$y)),1,1))
  
  for (i in 1:nrow(cell_type)){
    if (cell_type[[1]][i] == 'CD69'){
      next
    }
    indices <- which(theData$"Cell Type"==cell_type[[1]][i], arr.ind = TRUE)
    sub_u_data <- u_data[indices,]
    sub_u_data$x = (sub_u_data$x+1)/20
    sub_u_data$y = (sub_u_data$y+1)/20
    # win <- owin(xrange = range(art_dist_map$x), yrange = range(art_dist_map$y))
    win <- owin(xrange = range(sinu_dist_map$x), yrange = range(sinu_dist_map$y))
    points_ppp <- ppp(x = sub_u_data$x, y = sub_u_data$y, window = win)
    # fit_art <- ppm(points_ppp ~ covariate_im_art)
    fit_sinu <- ppm(points_ppp ~ covariate_im_sinu)
    # df <- rbind(df, data.frame(name = cell_type[[1]][i], distance_to_arteriole = fit_art$'coef'['covariate_im_art'][[1]]))
    df <- rbind(df, data.frame(name = cell_type[[1]][i], distance_to_sinusoid = fit_sinu$'coef'['covariate_im_sinu'][[1]]))
  }
  write.xlsx(df, file = sprintf("%s/%s",savepath, theDatalist[f]))
  print(sprintf('%d out of %d : %s done',f,length(theDatalist),theDatalist[f]))
}
