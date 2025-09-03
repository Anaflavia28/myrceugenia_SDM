########################################################################
###### Modeling rare species (flexsdm - Velazco et al. 2022) ###########
########################################################################
### by: Ana Flávia Augustin in 2025 - Myrceugenia joinvillensis 
## current and future scenarios 

######################## Pre-modeling #################################
#######################################################################

### 1 Require packages 
library(flexsdm)
library(terra)
library(dplyr)
library(readr)
library(raster)
library(openxlsx)

# Directories
coord_path <- "C:/Users/LENOVO/Desktop/Review FLORA/Testes/flexdm"
env_path <- "C:/Users/LENOVO/Desktop/Review FLORA/Testes/flexdm/env_present"
env_project <- "C:/Users/LENOVO/Desktop/Review FLORA/Testes/flexdm/env_projections"

### Environmental variables
env_files <- list.files(env_path, pattern = ".asc$", full.names = TRUE)
env_stack <- rast(env_files)

# Upload one env layer to get the CRS
env_files <- list.files(env_path, pattern = ".asc$", full.names = TRUE)
env_rast <- rast(env_files[1])
env_rast

### Occurrences data
occ <- read.delim(file.path(coord_path, "occ_joinvillensis_2.txt"))
colnames(occ) <- tolower(gsub("[^a-zA-Z0-9_]", "", colnames(occ)))
print(colnames(occ))

occ

### Calibration area ########################
## (buffer 300 km)
ca <- calib_area(
  data = occ,
  x = "x",
  y = "y",
  method = c("buffer", width = 300000),
  crs = crs(env_stack)
)


# Visualize calibration area and occurrences data
plot(env_rast[[1]], main = "Área de Calibração sobre variável ambiental")   # background
plot(ca, add = TRUE, col = NA, border = "blue", lwd = 2)                    # buffer
points(occ$x, occ$y, pch = 20, col = "red")                       # occurrences

# Collinearity reduction on predictors
# Upload env of current scenario 
env_files <- list.files(env_path, pattern = ".asc$", full.names = TRUE)
env_stack <- rast(env_files)

# Prepare paths for each scenario
cenarios <- list.dirs(env_project, full.names = TRUE, recursive = FALSE)
cenarios

env_project <- lapply(cenarios, function(caminho) {
  list.files(caminho, pattern = ".asc$", full.names = TRUE)
})

names(env_project) <- basename(cenarios)

env_project

# Extract the main dir. of sub-folders 
dir_w_proj <- dirname(cenarios[1])

save_pca_path <- "C:/Users/LENOVO/Desktop/Review FLORA/Testes/flexdm/env_pca_proj"
dir.create(save_pca_path, showWarnings = FALSE)

# PCA - current and projected for different scenarios 
env_pca <- correct_colinvar(
  env_layer = env_stack,
  method = "pca",
  proj = dir_w_proj,        
  save_proj = save_pca_path,          
  restric_to_region = NULL,
  restric_pca_proj = FALSE
)

env_pca 

# Saving PCs in .tif format ## optional
pca_dir <- file.path(env_path, "pca_resultado")
dir.create(pca_dir, showWarnings = FALSE)

for (i in 1:nlyr(env_pca$env_layer)) {
  writeRaster(
    env_pca$env_layer[[i]],
    filename = file.path(pca_dir, paste0("PC", i, ".tif")),
    overwrite = TRUE
  )
}

write_csv(as.data.frame(env_pca$coefficients), file.path(pca_dir, "pca_coefficients.csv"))
write_csv(env_pca$cumulative_variance, file.path(pca_dir, "pca_variance.csv"))

### Create pseudo-absences and/or background points

# checking how many pixels - to use 10-20% as background points (Guillera-Arroita et al 2014)
masked_rast <- mask(env_pca$env_layer[[1]], ca)

# Count number of cells not NA (i.e. valid cells)
num_valid_cells <- global(!is.na(masked_rast), "sum", na.rm = TRUE)
print(num_valid_cells)

# create background points 
rlayer_base <- env_stack[[1]]  # or any specific layer
plot(rlayer_base)

occ$pr_ab <- 1
occ_p <- occ %>% dplyr::filter(pr_ab == 1)
nrow(occ_p) 

bg <- sample_background(
  data = occ_p,
  x = "x",
  y = "y",
  n = 12000, # 10-20% of the valid cells (PC1 = 64.000)
  method = "random",
  rlayer = rlayer_base,
  calibarea = ca,
  sp_name = NULL,
)

bg

plot(rlayer_base[[1]], main = "Background Points")
points(bg[, c("x", "y")], pch = 20, col = "blue")
points(occ[, c("x", "y")], pch = 20, col = "red")


### Bind a presences and pseudo-absences or bg
occ_bg <- bind_rows(occ, bg)

occ_bg

### Partitioning
occ_partition <- part_random(
  data = occ_bg,
  pr_ab = "pr_ab",
  method = c(method = "loocv") 
)

occ_partition

#### Extracting env values (ex: PCs)
names(env_pca$env_layer) # it must be PC1, PC2...

occ_partition_vars <- sdm_extract(
  data = occ_partition,
  x = "x",
  y = "y",
  env_layer = env_pca$env_layer,
  variables = names(env_pca$env_layer)
)

occ_partition_vars

## save 
write.table(occ_partition_vars,
            file = "occ_partition_vars.txt",
            sep = "\t",         
            row.names = FALSE,    
            quote = FALSE)  

####################### Modeling #####################################
######################################################################

##### Maxent with hyper-parameters test ########
gridtest <-
  expand.grid(
    regmult = seq(1, 2.5, 3),
    classes = c("l", "lq")
  )

max_tun <- tune_max(
  data = occ_partition_vars,
  response = "pr_ab",
  predictors = c("PC1", "PC2", "PC3", "PC4"),
  predictors_f = NULL,
  partition = ".part",
  background = occ_partition_vars, ### background points included
  grid = gridtest,
  thr = "equal_sens_spec", 
  metric = "TSS",
  clamp = TRUE,
  pred_type = "cloglog",
  n_cores = 1 
)

max_tun$performance

## predicting models
max_tun_pred <- sdm_predict(
  models = max_tun,
  pred = env_pca$env_layer,
  thr = "equal_sens_spec",
  con_thr = FALSE,
  predict_area = NULL,
)

max_tun_pred

## plotting models
list_models <- terra::rast(max_tun_pred)
plot(list_models)

## predicting models for future scenarios #############

################ Scenario 1 ###########################
pca_scenario1_path <- "C:/Users/LENOVO/Desktop/Review FLORA/Testes/flexdm/env_pca_proj/scenario_1/pcs.tif"

pca_stack_1 <- rast(pca_scenario1_path)

names(pca_stack_1) <- c("PC1", "PC2", "PC3", "PC4")

max_pred_scenario1 <- sdm_predict(
  models = max_tun,
  pred = pca_stack_1,
  thr = "equal_sens_spec",
  con_thr = FALSE,
  predict_area = NULL
)

max_pred_scenario1

## plotting models
list_models <- terra::rast(max_pred_scenario1)
plot(list_models)

################ Scenario 2 ########################
pca_scenario2_path <- "C:/Users/LENOVO/Desktop/Review FLORA/Testes/flexdm/env_pca_proj/scenario_2/pcs.tif"

pca_stack_2 <- rast(pca_scenario2_path)

names(pca_stack_2) <- c("PC1", "PC2", "PC3", "PC4")

max_pred_scenario2 <- sdm_predict(
  models = max_tun,
  pred = pca_stack_2,
  thr = "equal_sens_spec",
  con_thr = FALSE,
  predict_area = NULL
)

max_pred_scenario2

## plotting models
list_models <- terra::rast(max_pred_scenario2)
plot(list_models)

################ Scenario 3 ###########################
pca_scenario3_path <- "C:/Users/LENOVO/Desktop/Review FLORA/Testes/flexdm/env_pca_proj/scenario_3/pcs.tif"

pca_stack_3 <- rast(pca_scenario3_path)

names(pca_stack_3) <- c("PC1", "PC2", "PC3", "PC4")

max_pred_scenario3 <- sdm_predict(
  models = max_tun,
  pred = pca_stack_3,
  thr = "equal_sens_spec",
  con_thr = FALSE,
  predict_area = NULL
)

max_pred_scenario3

## plotting models
list_models <- terra::rast(max_pred_scenario3)
plot(list_models)

################ Scenario 4 ###########################
pca_scenario4_path <- "C:/Users/LENOVO/Desktop/Review FLORA/Testes/flexdm/env_pca_proj/scenario_4/pcs.tif"

pca_stack_4 <- rast(pca_scenario4_path)

names(pca_stack_4) <- c("PC1", "PC2", "PC3", "PC4")

max_pred_scenario4 <- sdm_predict(
  models = max_tun,
  pred = pca_stack_4,
  thr = "equal_sens_spec",
  con_thr = FALSE,
  predict_area = NULL
)

## plotting models
list_models <- terra::rast(max_pred_scenario4)
plot(list_models)

####################### post modeling ################################
######################################################################

#### save a spreadsheet with performances ######
model_eval <- sdm_summarize(models = list(max_tun))
model_eval

write.xlsx(model_eval, file = file.path(coord_path, "Maxent_perf_join.xlsx"))

#### Saving rasters in Geotiff ################

### Maxent

unique(values(max_tun_pred$max[["equal_sens_spec"]]))

raster_continuo <- max_tun_pred$max[["max"]]
raster_binario <- max_tun_pred$max[["equal_sens_spec"]]

writeRaster(raster_continuo, "pres_cont_joi.tif", overwrite = TRUE)
writeRaster(raster_binario, "pres_bin_joi.tif", overwrite = TRUE)

### scenario 1
max_pred_scenario1

unique(values(max_pred_scenario1$max[["equal_sens_spec"]]))

raster_continuo <- max_pred_scenario1$max[["max"]]
raster_binario <- max_pred_scenario1$max[["equal_sens_spec"]]

writeRaster(raster_continuo, "scenario1_cont_joi.tif", overwrite = TRUE)
writeRaster(raster_binario, "scenario1_bin_joi.tif", overwrite = TRUE)

## scenario 2
max_pred_scenario2

unique(values(max_pred_scenario2$max[["equal_sens_spec"]]))

raster_continuo <- max_pred_scenario2$max[["max"]]
raster_binario <- max_pred_scenario2$max[["equal_sens_spec"]]

writeRaster(raster_continuo, "scenario2_cont_joi.tif", overwrite = TRUE)
writeRaster(raster_binario, "scenario2_bin_joi.tif", overwrite = TRUE)

## scenario 3
max_pred_scenario3

unique(values(max_pred_scenario3$max[["equal_sens_spec"]]))

raster_continuo <- max_pred_scenario3$max[["max"]]
raster_binario <- max_pred_scenario3$max[["equal_sens_spec"]]

writeRaster(raster_continuo, "scenario3_cont_joi.tif", overwrite = TRUE)
writeRaster(raster_binario, "scenario3_bin_joi.tif", overwrite = TRUE)

## scenario 4
max_pred_scenario4

unique(values(max_pred_scenario4$max[["equal_sens_spec"]]))

raster_continuo <- max_pred_scenario4$max[["max"]]
raster_binario <- max_pred_scenario4$max[["equal_sens_spec"]]

writeRaster(raster_continuo, "scenario4_cont_joi.tif", overwrite = TRUE)
writeRaster(raster_binario, "scenario4_bin_joi.tif", overwrite = TRUE)

#####################################################################################