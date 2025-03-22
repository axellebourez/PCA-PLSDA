# PCA / PLSDA Shiny Application

## Overview
This repository contains a Shiny application for performing **PCA (Principal Component Analysis)** and **PLSDA (Partial Least Squares Discriminant Analysis)** on CSV datasets. The application provides an intuitive interface for uploading data, selecting variables, and visualizing results in both 2D and 3D plots.

## Features
- **Data Upload:** Easily import CSV files with customizable sample and annotation columns.
- **Data Preprocessing:** Remove unwanted samples and automatically handle missing values by imputing with column means.
- **Analysis Options:** 
  - **PCA:** Visualize principal components with interactive scatter plots.
  - **PLSDA:** Perform discriminant analysis with group differentiation and permutation tests for model validation.
- **Visualizations:** 
  - Interactive 2D plots with ellipses for group confidence regions.
  - 3D scatter plots with ellipsoids for advanced visualization.
  - Scree plots for explained variance and loadings plots to display variable contributions.
- **Interactivity:** Built-in tooltips and dynamic controls using Plotly.

## Requirements
- R (version 4.x or higher recommended)
- R Packages:
  - shiny
  - ggplot2
  - mixOmics
  - DT
  - colourpicker
  - plotly
  - ellipse
  - ggrepel

You can install the required packages using the following commands in R:

```R
install.packages(c("shiny", "ggplot2", "DT", "colourpicker", "plotly", "ellipse", "ggrepel"))
install.packages("mixOmics")  # mixOmics may require installation from Bioconductor
