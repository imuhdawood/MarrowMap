library(spatstat)
library(readxl)
library(openxlsx)

datapath = "path/to/data" #From get_metadata
distmap_path = "path/to/distancemap"
savepath = "path/to/output"
theDatalist <- list.files(datapath)

for (f in 1:length(theDatalist)){
  print(sprintf('%d out of %d : %s started',f,length(theDatalist),theDatalist[f]))
  Datapath = sprintf("%s/%s", datapath,theDatalist[f])
  theData <- read_excel(Datapath)
  u_data  <- theData[, c('Cell Type',"x","y")]
  cell_type = unique(theData['Cell Type'])
  
  for (ct in 1:nrow(cell_type)){
    
    if (cell_type[[1]][ct] == 'CD69'){
      next
    }
    
    if (cell_type[[1]][ct] == 'MNP'){
      next
    }
    
    cell_name = cell_type[[1]][ct]
    
    if (cell_name == 'Granulocyte/mast'){
      cell_name = 'Granulocyte_mast'
    }
    
    subsavepath = sprintf('%s/%s',savepath,cell_name)
    
    if (file.exists(sprintf('%s/%s.xlsx',subsavepath,strsplit(theDatalist[f],".xlsx")[[1]]))) {
      print(sprintf('%s already exists',cell_type[[1]][ct]))
      next
    }
    
    if (!file.exists(subsavepath)) {
      # Directory or file does not exist, so you can create it
      dir.create(subsavepath, recursive = TRUE)  # For directories
    }
    
    df <- data.frame(
      cell_type = character(),
      distance_to_cell = numeric(),
      stringsAsFactors = FALSE  # Avoids converting character vectors to factors
    )
    
    cell_dist_map = read.csv(sprintf('%s/%s/%s.csv',distmap_path,cell_name,strsplit(theDatalist[f],".xlsx")[[1]]))
    covariate_im_cell <- as.im(cell_dist_map, c(as.numeric(max(cell_dist_map$x)),as.numeric(max(cell_dist_map$y)),1,1))
    
    for (i in 1:nrow(cell_type)){
      if (cell_type[[1]][i] == 'CD69'){
        next
      }
      
      if (cell_type[[1]][i] == 'MNP'){
        next
      }
      
      indices <- which(theData$"Cell Type"==cell_type[[1]][i], arr.ind = TRUE)
      sub_u_data <- u_data[indices,]
      sub_u_data$x = (sub_u_data$x+1)/20
      sub_u_data$y = (sub_u_data$y+1)/20
      win <- owin(xrange = range(cell_dist_map$x), yrange = range(cell_dist_map$y))
      points_ppp <- ppp(x = sub_u_data$x, y = sub_u_data$y, window = win)
      fit_cell <- ppm(points_ppp ~ covariate_im_cell)
      df <- rbind(df, data.frame(name = cell_type[[1]][i], distance_to_cell = fit_cell$'coef'['covariate_im_cell'][[1]]))
    }
    write.xlsx(df, file = sprintf("%s/%s",subsavepath, theDatalist[f]))
    print(sprintf('%s done',cell_type[[1]][ct]))
  }
  print(sprintf('%d out of %d : %s done',f,length(theDatalist),theDatalist[f]))
}
