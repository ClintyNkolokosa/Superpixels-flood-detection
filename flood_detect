# ============================================================================
# SUPERPIXEL-BASED ANOMALY DETECTION - ALL BANDS, MULTIPLE CI LEVELS
# Processing VV, VH, and RVI bands with 90%, 95%, and 99% confidence intervals
# Methods: IQR, GESD, and Isolation Forest
# Demo: Bangula wetland, lower Shire valley, Malawi
# Clinty 2025
# ============================================================================

library(terra)
library(supercells)    # Superpixels
library(dplyr)
library(tibble)
library(anomalize)     # for anomaly detection using IQR and GESD
library(tibbletime)
library(future)
library(future.apply)
library(isotree)       # for anomaly detection using iForest

# ============================================================================
# SETUP
# ============================================================================

data_dir <- "/home/cln3/SAR/Anomaly detection/S1_Wetland_Project"
#setwd(data_dir)

# Create output directory
output_base <- "superpixel_anomaly_results_multiband"
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

# Get all TIF files
tif_files <- list.files(pattern = "\\.tif$", full.names = TRUE)
dates <- as.Date(gsub(".tif", "", basename(tif_files)))
tif_files <- tif_files[order(dates)]
dates <- sort(dates)

cat(sprintf("Found %d images\n", length(tif_files)))
cat(sprintf("Date range: %s to %s\n", as.character(min(dates)), as.character(max(dates))))

# Load reference raster for dimensions
ref_raster <- rast(tif_files[1])
n_pixels <- ncell(ref_raster)
cat(sprintf("Original pixels: %s\n", format(n_pixels, big.mark=",")))
cat(sprintf("Bands per image: %d (1=VV, 2=VH, 3=RVI)\n", nlyr(ref_raster)))

# Configuration
bands_to_process <- list(
  VV = 1,
  VH = 2,
  RVI = 3
)

confidence_levels <- list(
  CI90 = 0.10,
  CI95 = 0.05,
  CI99 = 0.01
)

methods <- c("iqr", "gesd", "iforest")
# ============================================================================
# ANOMALY DETECTION FUNCTIONS
# ============================================================================

# Function: IQR/GESD anomaly detection
detect_superpixel_anomalies <- function(sp_ts, dates, method, alpha, max_anoms) {
  valid_idx <- !is.na(sp_ts)
  n_valid <- sum(valid_idx)
  
  if (n_valid < 24 || length(unique(sp_ts[valid_idx])) < 5) {
    return(rep(0, length(sp_ts)))
  }
  
  tryCatch({
    df <- tibble(date = dates[valid_idx], value = sp_ts[valid_idx]) %>%
      as_tbl_time(index = date) %>%
      time_decompose(value, method = "stl", frequency = "auto",
                     trend = "auto", message = FALSE) %>%
      anomalize(remainder, method = method, alpha = alpha, max_anoms = max_anoms)
    
    result <- rep(0, length(sp_ts))
    for (i in 1:nrow(df)) {
      date_idx <- which(dates == df$date[i])
      if (length(date_idx) > 0 && df$anomaly[i] == "Yes") {
        if (df$remainder[i] < df$remainder_l1[i]) {
          result[date_idx] <- 1  # Low anomaly
        } else {
          result[date_idx] <- 2  # High anomaly
        }
      }
    }
    return(result)
  }, error = function(e) {
    return(rep(0, length(sp_ts)))
  })
}

# Function: Isolation Forest anomaly detection
detect_superpixel_iforest <- function(sp_ts, dates, alpha) {
  valid_idx <- !is.na(sp_ts)
  n_valid <- sum(valid_idx)
  
  if (n_valid < 24 || length(unique(sp_ts[valid_idx])) < 5) {
    return(rep(0, length(sp_ts)))
  }
  
  tryCatch({
    df <- tibble(date = dates[valid_idx], value = sp_ts[valid_idx]) %>%
      as_tbl_time(index = date) %>%
      time_decompose(value, method = "stl", frequency = "auto",
                     trend = "auto", message = FALSE)
    
    iso_data <- data.frame(remainder = df$remainder)
    contamination <- alpha
    
    iso_model <- isolation.forest(
      iso_data,
      ntrees = 100,
      sample_size = min(256, nrow(iso_data)),
      prob_pick_pooled_gain = 1.0
    )
    
    iso_scores <- predict(iso_model, iso_data, type = "score")
    threshold <- quantile(iso_scores, probs = 1 - contamination)
    
    result <- rep(0, length(sp_ts))
    for (i in 1:nrow(df)) {
      date_idx <- which(dates == df$date[i])
      if (length(date_idx) > 0 && iso_scores[i] > threshold) {
        result[date_idx] <- 1
      }
    }
    return(result)
  }, error = function(e) {
    return(rep(0, length(sp_ts)))
  })
}

# ============================================================================
# PARALLEL PROCESSING SETUP
# ============================================================================

n_cores <- 12 # edit according to your machine
plan(multisession, workers = n_cores)
options(future.globals.maxSize = 2000 * 1024^2)
cat(sprintf("\nParallel processing: %d cores\n", n_cores))

# ============================================================================
# MAIN PROCESSING LOOP - BY BAND
# ============================================================================

processing_times <- list()
overall_start <- Sys.time()

for (band_name in names(bands_to_process)) {
  band_idx <- bands_to_process[[band_name]]
  
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat(sprintf("PROCESSING BAND: %s (Band %d)\n", band_name, band_idx))
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  band_start <- Sys.time()
  
  # ============================================================================
  # STEP 1: CREATE SUPERPIXELS FROM MEAN IMAGE
  # ============================================================================
  
  cat("\n[1/5] Loading band and computing temporal mean...\n")
  band_stack <- rast()
  for (i in seq_along(tif_files)) {
    r <- rast(tif_files[i])
    band_layer <- r[[band_idx]]
    band_stack <- c(band_stack, band_layer)
    if (i %% 10 == 0) cat(sprintf("  Loaded %d/%d\n", i, length(tif_files)))
  }
  
  cat("Computing temporal mean...\n")
  band_mean <- app(band_stack, fun = mean, na.rm = TRUE)
  
  # Generate superpixels
  cat("\n[2/5] Generating superpixels...\n")
  superpixel_start <- Sys.time()
  
  k_value <- 10000
  superpixels <- supercells(
    x = band_mean,
    k = k_value,
    compactness = 1,
    iter = 10
  )
  
  superpixel_time <- as.numeric(difftime(Sys.time(), superpixel_start, units = "secs"))
  n_superpixels <- nrow(superpixels)
  reduction_factor <- n_pixels / n_superpixels
  
  cat(sprintf("  Generated %d superpixels in %.1f seconds\n", n_superpixels, superpixel_time))
  cat(sprintf("  Reduction factor: %.1fx\n", reduction_factor))
  
  # Save superpixel map
  superpixel_dir <- file.path(output_base, band_name)
  dir.create(superpixel_dir, showWarnings = FALSE, recursive = TRUE)
  
  superpixel_raster <- rasterize(vect(superpixels), band_mean, field = "supercells")
  writeRaster(superpixel_raster, 
              file.path(superpixel_dir, "superpixel_map.tif"), 
              overwrite = TRUE)
  
  # ============================================================================
  # STEP 2: EXTRACT TIME SERIES
  # ============================================================================
  
  cat("\n[3/5] Extracting time series for each superpixel...\n")
  extract_start <- Sys.time()
  
  superpixel_vect <- vect(superpixels)
  superpixel_ts <- matrix(NA, nrow = length(dates), ncol = n_superpixels)
  
  for (i in seq_along(dates)) {
    r <- rast(tif_files[i])
    band_layer <- r[[band_idx]]
    extracted <- extract(band_layer, superpixel_vect, fun = mean, na.rm = TRUE)
    superpixel_ts[i, ] <- extracted[, 2]
    
    if (i %% 10 == 0) cat(sprintf("  Processed %d/%d dates\n", i, length(dates)))
  }
  
  extract_time <- as.numeric(difftime(Sys.time(), extract_start, units = "secs"))
  cat(sprintf("Extraction complete in %.1f seconds\n", extract_time))
  
  # ============================================================================
  # STEP 3: ANOMALY DETECTION FOR ALL CI LEVELS AND METHODS
  # ============================================================================
  
  cat("\n[4/5] Running anomaly detection...\n")
  
  for (ci_name in names(confidence_levels)) {
    alpha <- confidence_levels[[ci_name]]
    max_anoms <- 0.20
    
    cat(sprintf("\n  Processing %s (alpha=%.2f):\n", ci_name, alpha))
    
    for (method in methods) {
      cat(sprintf("    [%s] Processing %d superpixels...", 
                  toupper(method), n_superpixels))
      method_start <- Sys.time()
      
      if (method == "iforest") {
        results <- future_lapply(1:n_superpixels, function(i) {
          detect_superpixel_iforest(superpixel_ts[, i], dates, alpha)
        }, future.seed = TRUE)
      } else {
        results <- future_lapply(1:n_superpixels, function(i) {
          detect_superpixel_anomalies(superpixel_ts[, i], dates, method, alpha, max_anoms)
        }, future.seed = TRUE)
      }
      
      method_time <- as.numeric(difftime(Sys.time(), method_start, units = "secs"))
      cat(sprintf(" %.1fs\n", method_time))
      
      # Convert to matrix
      results_matrix <- do.call(cbind, results)
      
      # ============================================================================
      # STEP 4: CONVERT TO RASTERS AND SAVE
      # ============================================================================
      
      output_dir <- file.path(superpixel_dir, ci_name, method)
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      
      cat(sprintf("      Saving %d rasters...", length(dates)))
      save_start <- Sys.time()
      
      for (i in seq_along(dates)) {
        date_str <- format(dates[i], "%Y%m%d")
        anomaly_values <- results_matrix[i, ]
        superpixels$anomaly <- anomaly_values
        
        anomaly_raster <- rasterize(
          vect(superpixels), 
          band_mean, 
          field = "anomaly"
        )
        
        output_file <- file.path(output_dir, 
                                 sprintf("%s_%s_%s_%s.tif", 
                                         method, band_name, ci_name, date_str))
        writeRaster(anomaly_raster, output_file, 
                    overwrite = TRUE, datatype = "INT1U")
      }
      
      save_time <- as.numeric(difftime(Sys.time(), save_start, units = "secs"))
      cat(sprintf(" %.1fs\n", save_time))
      
      # Store timing
      timing_key <- sprintf("%s_%s_%s", band_name, ci_name, method)
      processing_times[[timing_key]] <- list(
        detection = method_time,
        saving = save_time
      )
    }
  }
  
  band_time <- as.numeric(difftime(Sys.time(), band_start, units = "mins"))
  cat(sprintf("\n%s processing complete in %.2f minutes\n", band_name, band_time))
}

# Clean up parallel processing
plan(sequential)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

overall_time <- as.numeric(difftime(Sys.time(), overall_start, units = "mins"))

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("PROCESSING COMPLETE - FINAL SUMMARY\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat(sprintf("\nTotal processing time: %.2f minutes (%.2f hours)\n", 
            overall_time, overall_time/60))

cat(sprintf("\nOutput structure:\n"))
cat(sprintf("  Base directory: %s/\n", output_base))
cat(sprintf("  ├── VV/\n"))
cat(sprintf("  │   ├── superpixel_map.tif\n"))
cat(sprintf("  │   ├── CI90/\n"))
cat(sprintf("  │   │   ├── iqr/     (%d rasters)\n", length(dates)))
cat(sprintf("  │   │   ├── gesd/    (%d rasters)\n", length(dates)))
cat(sprintf("  │   │   └── iforest/ (%d rasters)\n", length(dates)))
cat(sprintf("  │   ├── CI95/ (same structure)\n"))
cat(sprintf("  │   └── CI99/ (same structure)\n"))
cat(sprintf("  ├── VH/ (same structure as VV)\n"))
cat(sprintf("  └── RVI/ (same structure as VV)\n"))

cat(sprintf("\nTotal files created:\n"))
total_rasters <- length(bands_to_process) * length(confidence_levels) * 
  length(methods) * length(dates)
total_maps <- length(bands_to_process)
cat(sprintf("  Anomaly rasters: %d\n", total_rasters))
cat(sprintf("  Superpixel maps: %d\n", total_maps))
cat(sprintf("  Total files: %d\n", total_rasters + total_maps))

cat(sprintf("\nProcessing breakdown:\n"))
for (band_name in names(bands_to_process)) {
  cat(sprintf("  %s:\n", band_name))
  for (ci_name in names(confidence_levels)) {
    for (method in methods) {
      timing_key <- sprintf("%s_%s_%s", band_name, ci_name, method)
      if (!is.null(processing_times[[timing_key]])) {
        times <- processing_times[[timing_key]]
        cat(sprintf("    %s-%s: detection=%.1fs, saving=%.1fs\n", 
                    ci_name, method, times$detection, times$saving))
      }
    }
  }
}

cat(sprintf("\nComparison with pixel-based approach:\n"))
cat(sprintf("  Original method: ~108 hours per band\n"))
cat(sprintf("  Superpixel method: %.2f hours for all 3 bands\n", overall_time/60))
cat(sprintf("  Speedup: ~%.0fx faster\n", (108*3*60) / overall_time))

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Script complete! Check output directory for results.\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# ============================================================================
# OPTIONAL: CREATE SAMPLE VISUALIZATIONS
# ============================================================================

cat("\nCreating sample visualization for first anomaly date...\n")

# Find a date with anomalies in VV band
sample_dir <- file.path(output_base, "VV", "CI90", "iqr")
sample_files <- list.files(sample_dir, pattern = "\\.tif$", full.names = TRUE)

if (length(sample_files) > 0) {
  # Load first few rasters and check for anomalies
  for (f in sample_files[1:min(10, length(sample_files))]) {
    r <- rast(f)
    if (max(values(r), na.rm = TRUE) > 0) {
      # Found anomaly
      date_str <- gsub(".*_(\\d{8})\\.tif", "\\1", basename(f))
      
      cat(sprintf("Example date with anomalies: %s\n", date_str))
      
      # Create comparison plot
      png(file.path(output_base, "sample_comparison.png"), 
          width = 3000, height = 2000, res = 150)
      par(mfrow = c(3, 3), mar = c(2,2,3,2))
      
      for (band_name in names(bands_to_process)) {
        for (method in methods) {
          file_path <- file.path(output_base, band_name, "CI90", method,
                                 sprintf("%s_%s_CI90_%s.tif", 
                                         method, band_name, date_str))
          if (file.exists(file_path)) {
            r <- rast(file_path)
            plot(r, main = sprintf("%s - %s - CI90", band_name, toupper(method)),
                 col = c("grey90", "blue", "red"))
          }
        }
      }
      
      dev.off()
      cat(sprintf("Sample visualization saved: %s/sample_comparison.png\n", output_base))
      break
    }
  }
}

cat("\nAll processing complete!\n")
