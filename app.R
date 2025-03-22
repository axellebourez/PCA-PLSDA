# Auteur : Axelle Bourez
# Date : 2025-02-22
# Licence : MIT (voir fichier LICENSE pour les détails)
# Description : PCA et PLSDA interactive

library(shiny)
library(ggplot2)
library(mixOmics)       
library(DT)
library(colourpicker)
library(plotly)
library(ellipse)
library(ggrepel)

## Fonction pour calculer l'ellipse 2D
computeEllipse <- function(data, level = 0.95) {
  if(nrow(data) < 2) return(NULL)
  cov_mat <- cov(data)
  if(abs(det(cov_mat)) < 1e-8) {
    cov_mat <- cov_mat + diag(1e-6, ncol(data))
  }
  ell <- ellipse(cov_mat, centre = colMeans(data), level = level)
  ell_df <- as.data.frame(ell)
  colnames(ell_df) <- c("x", "y")
  ell_df
}

computeEllipsoid3D <- function(data, level = 0.95, n = 30) {
  if(ncol(data) != 3) stop("Les données doivent comporter 3 colonnes pour l'ellipsoïde 3D.")
  m <- colMeans(data)
  cov_mat <- cov(data)
  if(abs(det(cov_mat)) < 1e-8) {
    cov_mat <- cov_mat + diag(1e-6, ncol(data))
  }
  eig <- eigen(cov_mat)
  V <- eig$vectors
  lambda <- eig$values
  scale_factor <- sqrt(qchisq(level, df = 3))
  u <- seq(0, pi, length.out = n)
  v <- seq(0, 2*pi, length.out = n)
  xs <- outer(sin(u), cos(v))
  ys <- outer(sin(u), sin(v))
  zs <- outer(cos(u), rep(1, n))
  xs <- sqrt(lambda[1]) * scale_factor * xs
  ys <- sqrt(lambda[2]) * scale_factor * ys
  zs <- sqrt(lambda[3]) * scale_factor * zs
  pts <- cbind(as.vector(xs), as.vector(ys), as.vector(zs))
  pts_rot <- t(V %*% t(pts))
  pts_rot <- sweep(pts_rot, 2, m, "+")
  x_mat <- matrix(pts_rot[,1], n, n)
  y_mat <- matrix(pts_rot[,2], n, n)
  z_mat <- matrix(pts_rot[,3], n, n)
  list(x = x_mat, y = y_mat, z = z_mat)
}

ui <- fluidPage(
  titlePanel("PCA / PLSDA"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Télécharger le fichier CSV", accept = c(".csv")),
      uiOutput("selectColsUI"),
      uiOutput("removeSamplesUI"),
      radioButtons("analysis", "Type d'analyse",
                   choices = c("PCA", "PLSDA"), selected = "PCA"),
      radioButtons("scaling", "Méthode de scaling",
                   choices = c("Standard", "Pareto"), selected = "Standard"),
      helpText("La colonne d'annotation est optionnelle. Si 'None' est sélectionnée, tous les points seront affichés en bleu."),
      selectInput("plot_type", "Type de graphique", choices = c("2D", "3D"), selected = "2D"),
      conditionalPanel(
        condition = "input.plot_type == '3D'",
        selectInput("comp_z", "Composante Z", choices = 1:4, selected = 3)
      ),
      selectInput("comp_x", "Composante X", choices = 1:4, selected = 1),
      selectInput("comp_y", "Composante Y", choices = 1:4, selected = 2),
      checkboxInput("add_ellipse", "Ajouter ellipses", value = TRUE),
      conditionalPanel(
        condition = "input.analysis == 'PLSDA'",
        numericInput("nperm", "Nombre de permutations", value = 100, min = 1)
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Graphique", plotlyOutput("analysisPlot")),
        tabPanel("Contributions", plotlyOutput("screePlot")),
        tabPanel("Loadings", plotlyOutput("loadingsPlot")
        )
      ),
      verbatimTextOutput("perm_result")
    )
  )
)

server <- function(input, output, session) {
  
  dataRaw <- reactive({
    req(input$file)
    read.csv2(input$file$datapath, header = TRUE, stringsAsFactors = FALSE)
  })
  
  output$selectColsUI <- renderUI({
    req(dataRaw())
    cols <- names(dataRaw())
    tagList(
      selectInput("sample_col", "Colonne des échantillons", choices = cols, selected = cols[1]),
      selectInput("annotation_col", "Colonne d'annotation", 
                  choices = c("None", setdiff(cols, cols[1])), selected = "None")
    )
  })
  
  output$removeSamplesUI <- renderUI({
    req(dataRaw(), input$sample_col)
    df <- dataRaw()
    sampleNames <- make.unique(as.character(df[[input$sample_col]]))
    selectizeInput("remove_samples", "Supprimer ces échantillons (optionnel):", 
                   choices = sampleNames, selected = NULL, multiple = TRUE)
  })
  
  parsedData <- reactive({
    req(dataRaw(), input$sample_col)
    df <- dataRaw()
    sampleCol <- as.character(input$sample_col)
    annotationCol <- as.character(input$annotation_col)
    
    if (!(sampleCol %in% names(df))) {
      showNotification("La colonne d'échantillons sélectionnée n'existe pas.", type = "error")
      return(NULL)
    }
    sampleNames <- make.unique(as.character(df[[sampleCol]]))
    
    if(!is.null(input$remove_samples)) {
      keepRows <- !(sampleNames %in% input$remove_samples)
      sampleNames <- sampleNames[keepRows]
      df <- df[keepRows, , drop = FALSE]
    }
    
    if (!is.null(annotationCol) && annotationCol != "None") {
      if (!(annotationCol %in% names(df))) {
        showNotification("La colonne d'annotation sélectionnée n'existe pas.", type = "error")
        return(NULL)
      }
      annotation <- as.factor(df[[annotationCol]])
      featCols <- setdiff(names(df), c(sampleCol, annotationCol))
    } else {
      annotation <- NULL
      featCols <- setdiff(names(df), sampleCol)
    }
    
    if (length(featCols) < 1) {
      showNotification("Aucune colonne de features n'est disponible pour l'analyse.", type = "error")
      return(NULL)
    }
    
    features <- df[, featCols, drop = FALSE]
    
    features <- data.frame(lapply(features, function(x) {
      x_num <- suppressWarnings(as.numeric(as.character(x)))
      if (any(is.na(x_num)) && !all(is.na(x_num))) {
        x_num[is.na(x_num)] <- mean(x_num, na.rm = TRUE)
      }
      x_num
    }))
    
    keep <- sapply(features, function(x) {
      s <- sd(x, na.rm = TRUE)
      if (is.na(s)) return(FALSE) else return(s > 0)
    })
    features <- features[, keep, drop = FALSE]
    
    list(samples = sampleNames, features = features, annotation = annotation)
  })
  
  analysisData <- reactive({
    req(parsedData())
    feats <- parsedData()$features
    if(ncol(feats) < 1) return(NULL)
    if(input$scaling == "Standard") {
      X <- scale(feats, center = TRUE, scale = TRUE)
    } else {
      X_center <- scale(feats, center = TRUE, scale = FALSE)
      sds <- apply(feats, 2, sd, na.rm = TRUE)
      X <- sweep(X_center, 2, sqrt(sds), FUN = "/")
    }
    list(X = X, samples = parsedData()$samples, annotation = parsedData()$annotation)
  })
  
  output$analysisPlot <- renderPlotly({
    req(analysisData())
    ad <- analysisData()
    X <- ad$X
    if(is.null(X) || ncol(X)==0) return(plotly_empty())
    samples <- ad$samples
    annot <- ad$annotation
    comp_x <- as.numeric(input$comp_x)
    comp_y <- as.numeric(input$comp_y)
    
    if(input$plot_type == "2D") {
      if(input$analysis == "PCA") {
        pcaRes <- prcomp(X, center = FALSE, scale. = FALSE)
        scores <- pcaRes$x
        explained <- (pcaRes$sdev^2) / sum(pcaRes$sdev^2) * 100
        dfPlot <- data.frame(Sample = samples,
                             PCx = scores[, comp_x],
                             PCy = scores[, comp_y])
        if(!is.null(annot)) dfPlot$Annotation <- annot
        x_lab <- paste0("PC", comp_x, " (", round(explained[comp_x], 1), "%)")
        y_lab <- paste0("PC", comp_y, " (", round(explained[comp_y], 1), "%)")
        if(!is.null(annot)) {
          p <- ggplot(dfPlot, aes(x = PCx, y = PCy, color = Annotation, text = Sample)) +
            geom_point(size = 3) +
            labs(x = x_lab, y = y_lab, title = "PCA") +
            theme_minimal() +
            scale_color_discrete(name = input$annotation_col)
        } else {
          p <- ggplot(dfPlot, aes(x = PCx, y = PCy, text = Sample)) +
            geom_point(size = 3, color = "blue") +
            labs(x = x_lab, y = y_lab, title = "PCA") +
            theme_minimal()
        }
        if(input$add_ellipse && !is.null(annot)) {
          groups <- levels(annot)
          ellipseData <- do.call(rbind, lapply(groups, function(g) {
            subData <- dfPlot[dfPlot$Annotation == g, c("PCx", "PCy")]
            if(nrow(subData) < 2) return(NULL)
            ell_df <- computeEllipse(subData, level = 0.95)
            if(is.null(ell_df)) return(NULL)
            ell_df$Annotation <- g
            ell_df
          }))
          if(!is.null(ellipseData)) {
            p <- p + geom_polygon(data = ellipseData, aes(x = x, y = y, fill = Annotation, group = Annotation),
                                  inherit.aes = FALSE, alpha = 0.2, color = NA)
          }
        }
        ggplotly(p, tooltip = "text")
      } else { 
        req(annot)
        validate(
          need(length(unique(annot)) > 1, "La PLSDA nécessite au moins deux groupes dans l'annotation.")
        )
        ncomp <- max(comp_x, comp_y, 4)
        plsdaRes <- plsda(X, annot, ncomp = ncomp)
        scores <- plsdaRes$variates$X
        explained <- plsdaRes$prop$X * 100
        dfPlot <- data.frame(Sample = samples,
                             T_x = scores[, comp_x],
                             T_y = scores[, comp_y],
                             Annotation = annot)
        x_lab <- paste0("T", comp_x, " (", round(explained[comp_x], 1), "%)")
        y_lab <- paste0("T", comp_y, " (", round(explained[comp_y], 1), "%)")
        p <- ggplot(dfPlot, aes(x = T_x, y = T_y, color = Annotation, text = Sample)) +
          geom_point(size = 3) +
          labs(x = x_lab, y = y_lab, title = "PLSDA") +
          theme_minimal() +
          scale_color_discrete(name = input$annotation_col)
        if(input$add_ellipse) {
          groups <- levels(annot)
          ellipseData <- do.call(rbind, lapply(groups, function(g) {
            subData <- dfPlot[dfPlot$Annotation == g, c("T_x", "T_y")]
            if(nrow(subData) < 2) return(NULL)
            ell_df <- computeEllipse(subData, level = 0.95)
            if(is.null(ell_df)) return(NULL)
            ell_df$Annotation <- g
            ell_df
          }))
          if(!is.null(ellipseData)) {
            p <- p + geom_polygon(data = ellipseData, aes(x = x, y = y, fill = Annotation, group = Annotation),
                                  inherit.aes = FALSE, alpha = 0.2, color = NA)
          }
        }
        ggplotly(p, tooltip = "text")
      }
    } else { 
      req(input$comp_z)
      comp_z <- as.numeric(input$comp_z)
      if(input$analysis == "PCA") {
        pcaRes <- prcomp(X, center = FALSE, scale. = FALSE)
        scores <- pcaRes$x
        explained <- (pcaRes$sdev^2) / sum(pcaRes$sdev^2) * 100
        dfPlot <- data.frame(Sample = samples,
                             PCx = scores[, comp_x],
                             PCy = scores[, comp_y],
                             PCz = scores[, comp_z])
        if(!is.null(annot)) dfPlot$Annotation <- annot
        p <- plot_ly(dfPlot, x = ~PCx, y = ~PCy, z = ~PCz, type = 'scatter3d', mode = 'markers',
                     marker = list(size = 5),
                     text = ~Sample,
                     color = if(!is.null(annot)) ~Annotation else I("blue"))
        if(input$add_ellipse && !is.null(annot)) {
          groups <- levels(annot)
          for(g in groups) {
            subData <- dfPlot[dfPlot$Annotation == g, c("PCx", "PCy", "PCz")]
            if(nrow(subData) < 2) next
            ell <- computeEllipsoid3D(as.matrix(subData), level = 0.95, n = 30)
            p <- add_surface(p, x = ell$x, y = ell$y, z = ell$z, opacity = 0.2,
                             showscale = FALSE, inherit = FALSE)
          }
        }
        p
      } else {
        req(annot, input$comp_z)
        comp_z <- as.numeric(input$comp_z)
        validate(
          need(length(unique(annot)) > 1, "La PLSDA nécessite au moins deux groupes dans l'annotation.")
        )
        ncomp <- max(comp_x, comp_y, comp_z, 4)
        plsdaRes <- plsda(X, annot, ncomp = ncomp)
        scores <- plsdaRes$variates$X
        explained <- plsdaRes$prop$X * 100
        dfPlot <- data.frame(Sample = samples,
                             T_x = scores[, comp_x],
                             T_y = scores[, comp_y],
                             T_z = scores[, comp_z],
                             Annotation = annot)
        p <- plot_ly(dfPlot, x = ~T_x, y = ~T_y, z = ~T_z, type = 'scatter3d', mode = 'markers',
                     marker = list(size = 5),
                     text = ~Sample,
                     color = ~Annotation)
        p
      }
    }
  })
  
  output$screePlot <- renderPlotly({
    req(analysisData())
    ad <- analysisData()
    if(input$analysis == "PCA") {
      pcaRes <- prcomp(ad$X, center = FALSE, scale. = FALSE)
      explained <- (pcaRes$sdev^2) / sum(pcaRes$sdev^2) * 100
      n <- min(length(explained), 4)
      dfScree <- data.frame(Composante = paste0("PC", 1:n), Contribution = explained[1:n])
      p <- ggplot(dfScree, aes(x = Composante, y = Contribution, fill = Composante)) +
        geom_bar(stat = "identity") +
        labs(title = "Contribution des composantes (PCA)", y = "Pourcentage d'explication", x = "Composantes") +
        theme_minimal() +
        theme(legend.position = "none")
      ggplotly(p)
    } else {
      plsdaRes <- plsda(ad$X, ad$annotation, ncomp = 4)
      explained <- plsdaRes$prop$X * 100
      n <- min(length(explained), 4)
      dfScree <- data.frame(Composante = paste0("T", 1:n), Contribution = explained[1:n])
      p <- ggplot(dfScree, aes(x = Composante, y = Contribution, fill = Composante)) +
        geom_bar(stat = "identity") +
        labs(title = "Contribution des composantes (PLSDA)", y = "Pourcentage d'explication", x = "Composantes") +
        theme_minimal() +
        theme(legend.position = "none")
      ggplotly(p)
    }
  })
  
  output$loadingsPlot <- renderPlotly({
    req(analysisData())
    ad <- analysisData()
    X <- ad$X
    
    if(input$analysis == "PCA") {
      pcaRes <- prcomp(X, center = FALSE, scale. = FALSE)
      loadings <- as.data.frame(pcaRes$rotation)
      loadings$variable <- rownames(loadings)
      explained <- (pcaRes$sdev^2) / sum(pcaRes$sdev^2) * 100
      p <- ggplot(loadings, aes(x = PC1, y = PC2, label = variable)) +
        geom_point(size = 3, color = "steelblue") +
        geom_text_repel() +
        labs(x = paste0("PC1 (", round(explained[1], 1), "%)"),
             y = paste0("PC2 (", round(explained[2], 1), "%)"),
             title = "Loadings Plot (PC1 vs PC2) - PCA") +
        theme_light()
      ggplotly(p)
    } else {
      req(ad$annotation)
      ncomp <- max(2, 4)
      plsdaRes <- plsda(X, ad$annotation, ncomp = ncomp)
      loadings <- as.data.frame(plsdaRes$loadings$X)
      loadings$variable <- rownames(loadings)
      explained <- plsdaRes$prop$X * 100
      p <- ggplot(loadings, aes(x = comp1, y = comp2, label = variable)) +
        geom_point(size = 3, color = "steelblue") +
        geom_text_repel() +
        labs(x = paste0("comp1 (", round(explained[1], 1), "%)"),
             y = paste0("comp2 (", round(explained[2], 1), "%)"),
             title = "Loadings Plot (Comp1 vs Comp2) - PLSDA") +
        theme_light()
      ggplotly(p)
    }
  })
  
  output$perm_result <- renderPrint({
    if(input$analysis == "PLSDA") {
      req(analysisData())
      ad <- analysisData()
      X <- ad$X
      annot <- ad$annotation
      if(is.null(annot)) {
        cat("La PLSDA nécessite une colonne d'annotation avec au moins deux niveaux.")
        return()
      }
      validate(
        need(length(unique(annot)) > 1, "La PLSDA nécessite au moins deux groupes dans l'annotation.")
      )
      if(input$scaling == "Standard") {
        X_scaled <- scale(X, center = TRUE, scale = TRUE)
      } else {
        X_centered <- scale(X, center = TRUE, scale = FALSE)
        sds <- apply(X, 2, sd, na.rm = TRUE)
        X_scaled <- sweep(X_centered, 2, sqrt(sds), FUN = "/")
      }
      comp_x <- as.numeric(input$comp_x)
      comp_y <- as.numeric(input$comp_y)
      ncomp <- max(comp_x, comp_y, 4)
      p_val <- permutation_test_plsda(X_scaled, annot, ncomp, input$nperm)
      cat("P-value du test de permutation : ", p_val)
    }
  })
  
  permutation_test_plsda <- function(X, Y, ncomp, nperm) {
    original_model <- plsda(X, Y, ncomp = ncomp)
    original_pred <- predict(original_model, X)
    original_class <- original_pred$class$max.dist[, ncomp]
    original_error <- mean(original_class != Y)
    count <- 0
    for(i in 1:nperm) {
      Y_perm <- sample(Y)
      perm_model <- plsda(X, Y_perm, ncomp = ncomp)
      perm_pred <- predict(perm_model, X)
      perm_class <- perm_pred$class$max.dist[, ncomp]
      perm_error <- mean(perm_class != Y_perm)
      if(perm_error <= original_error) count <- count + 1
    }
    (count + 1) / (nperm + 1)
  }
}

shinyApp(ui, server)
