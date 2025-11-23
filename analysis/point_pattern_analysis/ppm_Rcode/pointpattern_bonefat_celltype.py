library(spatstat)
library(readxl)
library(openxlsx)

savepath = "path/to/output"
input_data_path = "path/to/input_data"
theDatalist <- list.files(input_data_path)
distmap_path = 'path/to/distance_map'
for (f in 1:length(theDatalist)){
  print(sprintf('%d out of %d : %s started',f,length(theDatalist),theDatalist[f]))
  Datapath = sprintf("%s/%s.", input_data_path, theDatalist[f])
  theData <- read_excel(Datapath)
  u_data  <- theData[, c('Cell Type',"x","y")]
  cell_type = unique(theData['Cell Type'])
  cov_data <- theData[, c('Distance to Bone',"Distance to Fat")]
  df <- data.frame(
    cell_type = character(),
    distance_to_bone = numeric(),
    stringsAsFactors = FALSE  # Avoids converting character vectors to factors
  )
  
  for (ct in 1:nrow(cell_type)){
    
    if (cell_type[[1]][ct] == 'CD69'){
      next
    }
    
    if (!file.exists(sprintf('savepath/%s',cell_type[[1]][ct]))) {
      # Directory or file does not exist, so you can create it
      dir.create(savepath, recursive = TRUE)  # For directories
    }
    
    cell_dist_map = read.csv(sprintf('%s/%s/%s.csv',distmap_path,cell_type[[1]][ct],strsplit(theDatalist[f],".xlsx")[[1]]))
    covariate_im_cell <- as.im(cell_dist_map, c(as.numeric(max(cell_dist_map$x)),as.numeric(max(cell_dist_map$y)),1,1))
    
    for (i in 1:nrow(cell_type)){
      indices <- which(theData$"Cell Type"==cell_type[[1]][i], arr.ind = TRUE)
      sub_u_data <- u_data[indices,]
      sub_u_data$x = (sub_u_data$x+1)/20
      sub_u_data$y = (sub_u_data$y+1)/20
      win <- owin(xrange = range(cell_dist_map$x), yrange = range(cell_dist_map$y))
      points_ppp <- ppp(x = sub_u_data$x, y = sub_u_data$y, window = win)
      fit_cell <- ppm(points_ppp ~ covariate_im_cell)
      df <- rbind(df, data.frame(name = cell_type[[1]][i], distance_to_bone = fit_cell$'coef'['covariate_im_cell'][[1]]))
    }
  write.xlsx(df, file = sprintf("%s/%s",savepath, theDatalist[f]))
  print(sprintf('%d out of %d : %s done',f,length(theDatalist),theDatalist[f]))
  }
}
