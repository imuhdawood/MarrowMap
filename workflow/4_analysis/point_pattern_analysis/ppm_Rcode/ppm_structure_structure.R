library(spatstat)
library(readxl)
library(openxlsx)
savepath = "path/to/struct_struct"
input_path = "path/to/struct_boundary_points"
distance_map_path = "path/to/distance_map_structure"
if (!file.exists(savepath)) {
  # Directory or file does not exist, so you can create it
  dir.create(savepath, recursive = TRUE)  # For directories
}
theDatalist <- list.files(input_path)
for (f in 1:length(theDatalist)){
  print(sprintf('%d out of %d : %s started',f,length(theDatalist),theDatalist[f]))
  Datapath = sprintf("%s/%s.", input_path, theDatalist[f])
  theData_bone <- read_excel(Datapath, sheet = "Bone")
  theData_fat <- read_excel(Datapath, sheet = "Fat")
  theData_art <- read_excel(Datapath, sheet = "Arteriole")
  theData_sinu <- read_excel(Datapath, sheet = "Sinusoid")

  thepatchdim <- read_excel(sprintf("path/to/sample_dim/%s.xlsx", strsplit(theDatalist[f],".xlsx")[[1]]))
  
  filtered_bone <- theData_bone[theData_bone$x >= thepatchdim$min_x & theData_bone$x <= thepatchdim$max_x & theData_bone$y >= thepatchdim$min_y & theData_bone$y <= thepatchdim$max_y, ]
  filtered_fat <- theData_fat[theData_fat$x >= thepatchdim$min_x & theData_fat$x <= thepatchdim$max_x & theData_fat$y >= thepatchdim$min_y & theData_fat$y <= thepatchdim$max_y, ]
  filtered_art <- theData_art[theData_art$x >= thepatchdim$min_x & theData_art$x <= thepatchdim$max_x & theData_art$y >= thepatchdim$min_y & theData_art$y <= thepatchdim$max_y, ]
  filtered_sinu <- theData_sinu[theData_sinu$x >= thepatchdim$min_x & theData_sinu$x <= thepatchdim$max_x & theData_sinu$y >= thepatchdim$min_y & theData_sinu$y <= thepatchdim$max_y, ]
  
  u_data_bone  <- filtered_bone[, c("x","y")]
  u_data_fat  <- filtered_fat[, c("x","y")]
  u_data_art  <- filtered_art[, c("x","y")]
  u_data_sinu  <- filtered_sinu[, c("x","y")]
  
  df_tobone <- data.frame(
    structure_type = character(),
    beta_1 = numeric(),
    stringsAsFactors = FALSE  # Avoids converting character vectors to factors
  )
  
  df_tofat <- data.frame(
    structure_type = character(),
    beta_1 = numeric(),
    stringsAsFactors = FALSE  # Avoids converting character vectors to factors
  )
  
  
  bone_dist_map = read.csv(sprintf('%s/bone/%s.csv',distance_map_path,strsplit(theDatalist[f],".xlsx")[[1]]))
  
  thepatchdim$min_x = as.integer((thepatchdim$min_x+1)/20)
  thepatchdim$max_x = as.integer((thepatchdim$max_x+1)/20)
  thepatchdim$min_y = as.integer((thepatchdim$min_y+1)/20)
  thepatchdim$max_y = as.integer((thepatchdim$max_y+1)/20)

  filtered_bone_map <- bone_dist_map[bone_dist_map$x >= thepatchdim$min_x & bone_dist_map$x <= thepatchdim$max_x & bone_dist_map$y >= thepatchdim$min_y & bone_dist_map$y <= thepatchdim$max_y, ]
  
  fat_dist_map = read.csv(sprintf('%s/fat/%s.csv',distance_map_path,strsplit(theDatalist[f],".xlsx")[[1]]))
  
  filtered_fat_map <- fat_dist_map[fat_dist_map$x >= thepatchdim$min_x & fat_dist_map$x <= thepatchdim$max_x & fat_dist_map$y >= thepatchdim$min_y & fat_dist_map$y <= thepatchdim$max_y, ]
  
  art_path = sprintf('%s/arteriole/%s.csv',distance_map_path,strsplit(theDatalist[f],".xlsx")[[1]])
  if (file.exists(art_path)) {
    df_toart <- data.frame(
      structure_type = character(),
      beta_1 = numeric(),
      stringsAsFactors = FALSE  # Avoids converting character vectors to factors
    )
    art_dist_map = read.csv(art_path)
    
    filtered_art_map <- art_dist_map[art_dist_map$x >= thepatchdim$min_x & art_dist_map$x <= thepatchdim$max_x & art_dist_map$y >= thepatchdim$min_y & art_dist_map$y <= thepatchdim$max_y, ]
    
    covariate_im_art <- as.im(filtered_art_map, c(as.numeric(max(filtered_art_map$x)),as.numeric(max(filtered_art_map$y)),1,1))
  }
  sinu_path = sprintf('%s/%s.csv',distance_map_path,strsplit(theDatalist[f],".xlsx")[[1]])
  if (file.exists(sinu_path)) {
    df_tosinu <- data.frame(
      structure_type = character(),
      beta_1 = numeric(),
      stringsAsFactors = FALSE  # Avoids converting character vectors to factors
    )
    sinu_dist_map = read.csv(sinu_path)
    
    filtered_sinu_map <- sinu_dist_map[sinu_dist_map$x >= thepatchdim$min_x & sinu_dist_map$x <= thepatchdim$max_x & sinu_dist_map$y >= thepatchdim$min_y & sinu_dist_map$y <= thepatchdim$max_y, ]
    
    covariate_im_sinu <- as.im(filtered_sinu_map, c(as.numeric(max(filtered_sinu_map$x)),as.numeric(max(filtered_sinu_map$y)),1,1))
  }
  covariate_im_bone <- as.im(filtered_bone_map, c(as.numeric(max(filtered_bone_map$x)),as.numeric(max(filtered_bone_map$y)),1,1))
  covariate_im_fat <- as.im(filtered_fat_map, c(as.numeric(max(filtered_fat_map$x)),as.numeric(max(filtered_fat_map$y)),1,1))
  
  u_data_bone$x = as.integer((u_data_bone$x+1)/20)
  u_data_bone$x[u_data_bone$x == 0] <- 1
  u_data_bone$y = as.integer((u_data_bone$y+1)/20)
  u_data_bone$y[u_data_bone$y == 0] <- 1
  u_data_bone = unique(u_data_bone)
  u_data_fat$x = as.integer((u_data_fat$x+1)/20)
  u_data_fat$x[u_data_fat$x == 0] <- 1
  u_data_fat$y = as.integer((u_data_fat$y+1)/20)
  u_data_fat$y[u_data_fat$y == 0] <- 1
  u_data_fat = unique(u_data_fat)
  if (nrow(u_data_art)!=0){
    u_data_art$x = as.integer((u_data_art$x+1)/20)
    u_data_art$x[u_data_art$x == 0] <- 1
    u_data_art$y = as.integer((u_data_art$y+1)/20)
    u_data_art$y[u_data_art$y == 0] <- 1
    u_data_art = unique(u_data_art)
  }
  if (nrow(u_data_sinu)!=0){
    u_data_sinu$x = as.integer((u_data_sinu$x+1)/20)
    u_data_sinu$x[u_data_sinu$x == 0] <- 1
    u_data_sinu$y = as.integer((u_data_sinu$y+1)/20)
    u_data_sinu$y[u_data_sinu$y == 0] <- 1
    u_data_sinu = unique(u_data_sinu)
  }
  
  win_bone <- owin(xrange = range(filtered_bone_map$x), yrange = range(filtered_bone_map$y))
  points_ppp_bone_to_bone <- ppp(x = u_data_bone$x, y = u_data_bone$y, window = win_bone)
  fit_bone_to_bone <- ppm(points_ppp_bone_to_bone ~ covariate_im_bone)
  points_ppp_fat_to_bone <- ppp(x = u_data_fat$x, y = u_data_fat$y, window = win_bone)
  fit_fat_to_bone <- ppm(points_ppp_fat_to_bone ~ covariate_im_bone)
  win_fat <- owin(xrange = range(filtered_fat_map$x), yrange = range(filtered_fat_map$y))
  points_ppp_fat_to_fat <- ppp(x = u_data_fat$x, y = u_data_fat$y, window = win_fat)
  fit_fat_to_fat <- ppm(points_ppp_fat_to_fat ~ covariate_im_fat)
  points_ppp_bone_to_fat <- ppp(x = u_data_bone$x, y = u_data_bone$y, window = win_fat)
  fit_bone_to_fat <- ppm(points_ppp_bone_to_fat ~ covariate_im_fat)
  if (file.exists(art_path)) {
    win_art <- owin(xrange = range(filtered_art_map$x), yrange = range(filtered_art_map$y))
    points_ppp_art_to_art <- ppp(x = u_data_art$x, y = u_data_art$y, window = win_art)
    fit_art_to_art <- ppm(points_ppp_art_to_art ~ covariate_im_art)
    points_ppp_bone_to_art <- ppp(x = u_data_bone$x, y = u_data_bone$y, window = win_art)
    fit_bone_to_art <- ppm(points_ppp_bone_to_art ~ covariate_im_art)
    points_ppp_fat_to_art <- ppp(x = u_data_fat$x, y = u_data_fat$y, window = win_art)
    fit_fat_to_art <- ppm(points_ppp_fat_to_art ~ covariate_im_art)
  }
  
  if (file.exists(sinu_path)) {
    win_sinu <- owin(xrange = range(filtered_sinu_map$x), yrange = range(filtered_sinu_map$y))
    points_ppp_sinu_to_sinu <- ppp(x = u_data_sinu$x, y = u_data_sinu$y, window = win_sinu)
    fit_sinu_to_sinu <- ppm(points_ppp_sinu_to_sinu ~ covariate_im_sinu)
    points_ppp_bone_to_sinu <- ppp(x = u_data_bone$x, y = u_data_bone$y, window = win_sinu)
    fit_bone_to_sinu <- ppm(points_ppp_bone_to_sinu ~ covariate_im_sinu)
    points_ppp_fat_to_sinu <- ppp(x = u_data_fat$x, y = u_data_fat$y, window = win_sinu)
    fit_fat_to_sinu <- ppm(points_ppp_fat_to_sinu ~ covariate_im_sinu)
  }
    
  if (nrow(u_data_art)!=0){
    points_ppp_art_to_bone <- ppp(x = u_data_art$x, y = u_data_art$y, window = win_bone)
    fit_art_to_bone <- ppm(points_ppp_art_to_bone ~ covariate_im_bone)
    points_ppp_art_to_fat <- ppp(x = u_data_art$x, y = u_data_art$y, window = win_fat)
    fit_art_to_fat <- ppm(points_ppp_art_to_fat ~ covariate_im_fat)
    if (file.exists(sinu_path)) {
      points_ppp_art_to_sinu <- ppp(x = u_data_art$x, y = u_data_art$y, window = win_sinu)
      fit_art_to_sinu <- ppm(points_ppp_art_to_sinu ~ covariate_im_sinu)
    }
  }
  if (nrow(u_data_sinu)!=0){
    points_ppp_sinu_to_bone <- ppp(x = u_data_sinu$x, y = u_data_sinu$y, window = win_bone)
    fit_sinu_to_bone <- ppm(points_ppp_sinu_to_bone ~ covariate_im_bone)
    points_ppp_sinu_to_fat <- ppp(x = u_data_sinu$x, y = u_data_sinu$y, window = win_fat)
    fit_sinu_to_fat <- ppm(points_ppp_sinu_to_fat ~ covariate_im_fat)
    if (file.exists(art_path)) {
      points_ppp_sinu_to_art <- ppp(x = u_data_sinu$x, y = u_data_sinu$y, window = win_art)
      fit_sinu_to_art <- ppm(points_ppp_sinu_to_art ~ covariate_im_art)
    }
  }
  
  df_tobone <- rbind(df_tobone, data.frame(name = 'Bone', beta_1 = fit_bone_to_bone$'coef'['covariate_im_bone'][[1]]))
  df_tobone <- rbind(df_tobone, data.frame(name = 'Fat', beta_1 = fit_fat_to_bone$'coef'['covariate_im_bone'][[1]]))
  if (nrow(u_data_art)!=0){
    df_tobone <- rbind(df_tobone, data.frame(name = 'Arteriole', beta_1 = fit_art_to_bone$'coef'['covariate_im_bone'][[1]]))
  }
  if (nrow(u_data_sinu)!=0){
    df_tobone <- rbind(df_tobone, data.frame(name = 'Sinusoid', beta_1 = fit_sinu_to_bone$'coef'['covariate_im_bone'][[1]]))
  }
  
  df_tofat <- rbind(df_tofat, data.frame(name = 'Bone', beta_1 = fit_bone_to_fat$'coef'['covariate_im_fat'][[1]]))
  df_tofat <- rbind(df_tofat, data.frame(name = 'Fat', beta_1 = fit_fat_to_fat$'coef'['covariate_im_fat'][[1]]))
  if (nrow(u_data_art)!=0){
    df_tofat <- rbind(df_tofat, data.frame(name = 'Arteriole', beta_1 = fit_art_to_fat$'coef'['covariate_im_fat'][[1]]))
  }
  if (nrow(u_data_sinu)!=0){
    df_tofat <- rbind(df_tofat, data.frame(name = 'Sinusoid', beta_1 = fit_sinu_to_fat$'coef'['covariate_im_fat'][[1]]))
  }
  
  if (file.exists(art_path)){
    df_toart <- rbind(df_toart, data.frame(name = 'Bone', beta_1 = fit_bone_to_art$'coef'['covariate_im_art'][[1]]))
    df_toart <- rbind(df_toart, data.frame(name = 'Fat', beta_1 = fit_fat_to_art$'coef'['covariate_im_art'][[1]]))
    df_toart <- rbind(df_toart, data.frame(name = 'Arteriole', beta_1 = fit_art_to_art$'coef'['covariate_im_art'][[1]]))
    if (nrow(u_data_sinu)!=0){
      df_toart <- rbind(df_toart, data.frame(name = 'Sinusoid', beta_1 = fit_sinu_to_art$'coef'['covariate_im_art'][[1]]))
    }
  }
  
  if (file.exists(sinu_path)){
    df_tosinu <- rbind(df_tosinu, data.frame(name = 'Bone', beta_1 = fit_bone_to_sinu$'coef'['covariate_im_sinu'][[1]]))
    df_tosinu <- rbind(df_tosinu, data.frame(name = 'Fat', beta_1 = fit_fat_to_sinu$'coef'['covariate_im_sinu'][[1]]))
    if (nrow(u_data_art)!=0){
      df_tosinu <- rbind(df_tosinu, data.frame(name = 'Arteriole', beta_1 = fit_art_to_sinu$'coef'['covariate_im_sinu'][[1]]))
    }
    df_tosinu <- rbind(df_tosinu, data.frame(name = 'Sinusoid', beta_1 = fit_sinu_to_sinu$'coef'['covariate_im_sinu'][[1]]))
  }
  
  wb <- createWorkbook()
  
  addWorksheet(wb, "To Bone")
  writeData(wb, "To Bone", df_tobone, rowNames = FALSE)
  
  addWorksheet(wb, "To Fat")
  writeData(wb, "To Fat", df_tofat, rowNames = FALSE)
  
  if (file.exists(art_path)){
    addWorksheet(wb, "To Arteriole")
    writeData(wb, "To Arteriole", df_toart, rowNames = FALSE)
  }
  
  if (file.exists(sinu_path)){
    addWorksheet(wb, "To Sinusoid")
    writeData(wb, "To Sinusoid", df_tosinu, rowNames = FALSE)
  }
  
  # Save the workbook
  saveWorkbook(wb, sprintf("%s/%s",savepath, theDatalist[f]), overwrite = TRUE)
  
  print(sprintf('%d out of %d : %s done',f,length(theDatalist),theDatalist[f]))
}
