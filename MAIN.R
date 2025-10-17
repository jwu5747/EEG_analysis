# ============================================================================
# EEG CORRELATION ANALYSIS - ORGANIZED CODE
# ============================================================================

# ----------------------------------------------------------------------------
# LIBRARIES
# ----------------------------------------------------------------------------
library('R.matlab')
library(HellCor)
library(doParallel)
library(ggplot2)
library(reshape2)

# ----------------------------------------------------------------------------
# GLOBAL VARIABLES
# ----------------------------------------------------------------------------

# Electrode names (63 electrodes)
electrodes <- c(
  'Fp1', 'Fz', 'F3', 'F7', 'FT9', 'FC5', 'FC1', 'C3', 'T7', 'TP9',
  'CP5', 'CP1', 'Pz', 'P3', 'P7', 'O1', 'Oz', 'O2', 'P4', 'P8',
  'TP10', 'CP6', 'CP2', 'Cz', 'C4', 'T8', 'FT10', 'FC6', 'FC2', 'F4',
  'F8', 'Fp2', 'AF7', 'AF3', 'AFz', 'F1', 'F5', 'FT7', 'FC3', 'FCz',
  'C1', 'C5', 'TP7', 'CP3', 'P1', 'P5', 'PO7', 'PO3', 'POz', 'PO4',
  'PO8', 'P6', 'P2', 'CP4', 'TP8', 'C6', 'C2', 'FC4', 'FT8', 'F6',
  'F2', 'AF4', 'AF8'
)



# File paths
base_dir <- "D:/ds/test/" #This is where all our files will be saved
matlab_dir <- "D:/ds/participant/" #This is the directory containing all participant files

participant.names <- list.files(path = matlab_dir)
participant <- list.files(path = matlab_dir, full.names = TRUE)

# ----------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ----------------------------------------------------------------------------

# Flatten nested list structure
flatten <- function(x) {
  x[[1]] <- x[[1]][[1]]
  return(x)
}

# Reshape matrix by combining rows
reshape <- function(mat, factor = 4) {
  if (1000 %% factor == 0) {
    n <- dim(mat)[1]
    m <- dim(mat)[2]
    new.mat <- matrix(NA, nrow = n / factor, ncol = m * factor)
    
    for (i in 1:(n / factor)) {
      start_row <- (i - 1) * factor + 1
      new.mat[i, ] <- c(mat[start_row, ], mat[start_row + 1, ], 
                        mat[start_row + 2, ], mat[start_row + 3, ])
    }
    return(new.mat)  
  } else {
    print("Error: 1000 must be divisible by factor")
  }
}

# ----------------------------------------------------------------------------
# DATA EXTRACTION FUNCTIONS
# ----------------------------------------------------------------------------

# Extract amplitude data from specific channel and latency
pull_data <- function(dat, channel, lat) {
  amplitudes <- dat[channel][[1]][lat, ]
  return(amplitudes)
}

# ----------------------------------------------------------------------------
# CORRELATION MEASUREMENT FUNCTIONS
# ----------------------------------------------------------------------------

# Measure correlation between two electrodes at specified latencies
measure_cor <- function(dat, elec1, elec2, lat1, lat2, type = 'hcor') {
  X <- pull_data(dat, elec1, lat1)
  Y <- pull_data(dat, elec2, lat2)
  
  if (type == 'hcor') {
    measure <- HellCor(X, Y)$Hcor  
  }
  if (type == 'pearson') {
    measure <- cor(X, Y)
  }
  if (type == 'spearman') {
    measure <- unname(cor.test(X, Y, method = "spearman")$estimate)
  }
  return(measure)
}

# ----------------------------------------------------------------------------
# CORRELATION MATRIX SAVING FUNCTIONS
# ----------------------------------------------------------------------------

# Save full correlation matrix for all electrode pairs across time
save_cormatrix <- function(dat, save_dir, cor_type) {
  n <- 250 #this has been hard-coded, note that we reshape our 1000 x n into 250 x (4n) where n is number of trials (or in this case tones)

  for (l in 1:62) {
    channel1 <- electrodes[l]
    
    for (k in (l + 1):63) {
      channel2 <- electrodes[k]
      filename <- paste(channel1, ';', channel2, '_tgt', '.csv', sep = "")
      file_path <- paste0(save_dir, filename)
      
      if (file.exists(file_path)) {
        cat("File", filename, "already exists. Skipping...\n")
        next
      }
      
      x <- foreach(i = 1:n, .combine = 'rbind') %:% 
        foreach(j = 1:n, .combine = 'c', .packages = 'HellCor', 
                .export = c('measure_cor', 'pull_data')) %dopar% {
                  measure_cor(dat, elec1 = l, elec2 = k, lat1 = i, 
                              lat2 = j, type = cor_type)
                }
      write.table(x, file = file_path, col.names = FALSE)
    }
  }
}

# ----------------------------------------------------------------------------
# PARALLEL PROCESSING SETUP
# ----------------------------------------------------------------------------

cores_used <- detectCores() - 1
cl <- makeCluster(cores_used)
registerDoParallel(cl)

# ----------------------------------------------------------------------------
# MAIN ANALYSIS - Target
# ----------------------------------------------------------------------------


for (i in length(participant)) {
  tmp <- readMat(participant[i])
  datatgt <- tmp$datatgt[2][[1]]
  
  # Flatten and reshape data
  for (j in 1:63) {
    datatgt[j] <- flatten(datatgt[j])
  }
  for (k in 1:63) {
    datatgt[[k]] <- reshape(datatgt[[k]])
  }
  
  # Create output directory and save correlation matrix
  save_dir <- paste0(base_dir, participant.names[i], "/")
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  save_cormatrix(datatgt, save_dir, "hcor")
}


# ----------------------------------------------------------------------------
# CLEANUP (remember to run this when done)
# ----------------------------------------------------------------------------

stopCluster(cl)