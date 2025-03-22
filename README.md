# PCA-PLSDA
This R script is an interactive Shiny application designed for performing PCA (Principal Component Analysis) and PLSDA (Partial Least Squares Discriminant Analysis) on user-provided CSV data. It provides a user-friendly interface where users can:

Upload Data: Easily load CSV files containing sample data.

Data Preprocessing: Choose the sample and annotation columns, and optionally remove specific samples.

Analysis Options: Select between PCA and PLSDA, along with scaling methods (Standard or Pareto).

Visualization: Generate interactive 2D and 3D plots using ggplot2 and Plotly. For 2D plots, the app can overlay ellipses representing the confidence regions for each group, while for 3D plots, it computes ellipsoids.

Results & Interpretation: View contributions (scree plots) and loadings plots to understand how variables contribute to the analyses.

Statistical Testing: For PLSDA, the app runs a permutation test to assess model significance.

The script makes use of several powerful R libraries including shiny, ggplot2, mixOmics, DT, colourpicker, plotly, ellipse, and ggrepel to ensure comprehensive data analysis and interactive visualizations.

