# Superpixels-flood-detection
Superpixel-Based Flood Anomaly Detection for Sentinel-1 SAR Time Series

**Authors:** Clinton Nkolokosa  
**Contact:** [clinton.nkolokosa@stir.ac.uk]  

---

## Overview

This repository contains an **unsupervised, superpixel-based framework** for detecting flood-related anomalies in **Sentinel-1 Ground Range Detected (GRD)** data. The workflow processes multi-band SAR imagery (VV, VH, and RVI) and identifies temporal anomalies using multiple statistical approaches.  

Key features:  
- **Bands processed:** VV, VH, RVI  
- **Anomaly detection methods:** IQR, GESD, Isolation Forest  
- **Confidence intervals:** 90%, 95%, 99%  
- **Output:** Superpixel maps and time series anomaly rasters  

The method significantly reduces computation time compared to pixel-based approaches, enabling rapid processing of **hundreds of SAR images**.

---

## Workflow

1. **Load SAR images:** Raster data in `.tif` format.  
2. **Compute temporal mean:** Generates a base image for superpixel segmentation.  
3. **Superpixel segmentation:** Using `supercells` to reduce computational load while preserving spatial structure.  
4. **Extract superpixel time series:** Average backscatter per superpixel over all dates.  
5. **Anomaly detection:**  
   - **IQR/GESD:** STL-decomposed residuals  
   - **Isolation Forest:** Robust detection on residuals  
   - Multiple confidence levels (90%, 95%, 99%) applied  
6. **Export results:**  
   - Superpixel map (`superpixel_map.tif`)  
   - Anomaly rasters per date, method, CI, and band  

Optional: Sample visualizations can be created for quick inspection of anomalies.

---

## Installation

This script requires R (≥4.3) and the following packages:

```r
install.packages(c(
  "terra", "supercells", "dplyr", "tibble", 
  "anomalize", "tibbletime", "future", "future.apply", "isotree"
))
