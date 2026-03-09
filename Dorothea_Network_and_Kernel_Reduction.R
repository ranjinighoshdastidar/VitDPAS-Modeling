### DOROTHEA MATRIX GENERATION ###
# Load required packages
library(igraph)
library(dplyr)
library(ggplot2)
library(svglite)

# Set target genes
genes_of_interest <- c("JUNB", "ETS2", "FOS", "NR4A2", "KLF10", "CSRNP1", "HIF1A")

# Load Dorothea regulons
library(dorothea)
data(dorothea_hs, package = "dorothea")

# Filter for A and B confidence levels
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B"))

# Extract interactions for our genes of interest
tf_interactions <- regulons %>%
  dplyr::filter(tf %in% genes_of_interest)

target_interactions <- regulons %>%
  dplyr::filter(target %in% genes_of_interest)

# Combine interactions
direct_interactions <- rbind(tf_interactions, target_interactions) %>%
  distinct()

# Get secondary TFs and limit interactions
secondary_tfs <- setdiff(unique(c(direct_interactions$tf, direct_interactions$target)), genes_of_interest)

secondary_interactions <- regulons %>%
  dplyr::filter((tf %in% genes_of_interest & target %in% secondary_tfs) |
                  (tf %in% secondary_tfs & target %in% genes_of_interest)) %>%
  arrange(desc(abs(mor))) %>%
  head(50)

# Combine all interactions
all_interactions <- rbind(direct_interactions, secondary_interactions) %>%
  distinct()

# Create igraph network
network <- graph_from_data_frame(
  d = all_interactions %>% dplyr::dplyr::select(tf, target, mor),
  directed = TRUE
)

# Add node attributes
V(network)$type <- ifelse(V(network)$name %in% genes_of_interest, "primary", "secondary")

# Add edge attributes
E(network)$type <- ifelse(E(network)$mor > 0, "activation", "inhibition")
E(network)$weight <- abs(E(network)$mor)

# Simplify network for visualization
simplify_for_vis <- function(g, max_nodes = 40) {
  if(vcount(g) <= max_nodes) return(g)
  
  primary_nodes <- V(g)[V(g)$type == "primary"]$name
  remaining_slots <- max_nodes - length(primary_nodes)
  
  if(remaining_slots <= 0) {
    return(induced_subgraph(g, V(g)[V(g)$name %in% primary_nodes]))
  }
  
  secondary_nodes <- V(g)[V(g)$type == "secondary"]
  if(length(secondary_nodes) == 0) return(g)
  
  secondary_data <- data.frame(
    name = secondary_nodes$name,
    degree = degree(g, v = secondary_nodes)
  )
  
  secondary_data <- secondary_data[order(secondary_data$degree, decreasing = TRUE), ]
  top_secondaries <- head(secondary_data$name, remaining_slots)
  
  nodes_to_keep <- c(primary_nodes, top_secondaries)
  return(induced_subgraph(g, V(g)[V(g)$name %in% nodes_to_keep]))
}

# Function to plot network with circular layout
plot_network <- function(g, layout_type, title) {
  g <- simplify_for_vis(g)
  
  node_colors <- c("primary" = "#E41A1C", "secondary" = "#377EB8")
  edge_colors <- c("activation" = "#4DAF4A", "inhibition" = "#984EA3")
  
  layout <- layout_in_circle(g)
  
  vertices_df <- data.frame(
    name = V(g)$name,
    type = V(g)$type,
    x = layout[,1],
    y = layout[,2],
    stringsAsFactors = FALSE
  )
  
  edges_df <- data.frame(
    from = get.edgelist(g)[,1],
    to = get.edgelist(g)[,2],
    type = E(g)$type,
    weight = E(g)$weight,
    stringsAsFactors = FALSE
  )
  
  edges_with_coords <- edges_df %>%
    left_join(vertices_df, by = c("from" = "name")) %>%
    rename(x_from = x, y_from = y) %>%
    left_join(vertices_df, by = c("to" = "name")) %>%
    rename(x_to = x, y_to = y)
  
  p <- ggplot() +
    geom_segment(
      data = edges_with_coords,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to, color = type, size = weight),
      arrow = arrow(length = unit(2, 'mm')),
      alpha = 0.6
    ) +
    geom_point(
      data = vertices_df,
      aes(x = x, y = y, color = type, size = ifelse(type == "primary", 5, 3))
    ) +
    geom_text(
      data = vertices_df,
      aes(x = x, y = y, label = name),
      size = 2.5,
      vjust = -0.8,
      fontface = "bold.italic"  # Bold italic gene names
    ) +
    scale_color_manual(values = c(node_colors, edge_colors)) +
    scale_size_continuous(range = c(0.3, 2)) +
    labs(title = title) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right"
    )
  
  return(p)
}

# Create circular network plot
circular_plot <- plot_network(network, "circle", "Circular Layout")

# Save as PNG
ggsave("network_circle_layout.png", circular_plot, width = 8, height = 6, dpi = 100)
print("Saved PNG: network_circle_layout.png")

# Save as SVG
ggsave("network_circle_layout.svg", circular_plot, width = 8, height = 6, device = svglite)
print("Saved SVG: network_circle_layout.svg")


# Script to extract subnetwork interactions for specific genes
# And create an adjacency matrix with only +1 or -1 values

# Load required packages
library(igraph)
library(dplyr)
library(tidyr)
library(svglite)

# List of genes for the subnetwork
genes_of_interest <- c(
  "ETS2", "FOS", "HIF1A", "ATF2", "CREB1", "EGR1", "ELK1", "ESR1", "ETS1", 
  "FLI1", "FOXM1", "JUN", "MITF", "MYC", "NFKB1", "SP1", "STAT1", "STAT3", 
  "TP53", "TWIST1", "ANGPT2", "CSF2", "FOSL1", "IL12B", "MMP3", "MMP9", 
  "TERT", "BCL2L1", "CCL2", "CCND1", "CDKN1A", "FGF2", "IL1B", "NOS2", 
  "NR3C1", "PLA2G4A", "PLAUR", "NR4A2", "KLF10", "JUNB")

# Load Dorothea regulons
if (!requireNamespace("dorothea", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("dorothea")
}

library(dorothea)
data(dorothea_hs, package = "dorothea")

# Filter for A and B confidence levels to reduce noise
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B"))

# Extract only the interactions between genes of interest
subnetwork_interactions <- regulons %>%
  dplyr::filter(tf %in% genes_of_interest & target %in% genes_of_interest) %>%
  # Keep only essential columns
  dplyr::dplyr::select(tf, target, mor)

# Create a simplified version with just +1 or -1
simplified_interactions <- subnetwork_interactions %>%
  dplyr::mutate(simplified_mor = ifelse(mor > 0, 1, -1)) %>%
  dplyr::dplyr::select(tf, target, simplified_mor)

# Save the simplified interactions
print(paste("Found", nrow(simplified_interactions), "direct interactions between the genes of interest"))
write.csv(simplified_interactions, "subnetwork_gene_interactions1.csv", row.names = FALSE)
print("Saved subnetwork_gene_interactions.csv")

# Create an adjacency matrix with values +1/-1
# Initialize matrix with zeros
adjacency_matrix <- matrix(
  0, 
  nrow = length(genes_of_interest), 
  ncol = length(genes_of_interest),
  dimnames = list(genes_of_interest, genes_of_interest)
)

# Fill matrix with +1/-1 values based on interactions
for (i in 1:nrow(simplified_interactions)) {
  row <- simplified_interactions[i, ]
  adjacency_matrix[row$tf, row$target] <- row$simplified_mor
}

# Save adjacency matrix
write.csv(adjacency_matrix, "subnetwork_adjacency_matrix1.csv")
print("Saved subnetwork_adjacency_matrix.csv")

# Create a visualization of the matrix
print("Creating heatmap visualization...")

# Create bold and italicized gene names for labels
bold_italic_genes <- paste0("bolditalic(", genes_of_interest, ")")

# Create a heatmap of the adjacency matrix
png("subnetwork_heatmap2.png", width = 1000, height = 800)
# Use custom color palette for -1, 0, 1
colors <- c("#984EA3", "white", "#4DAF4A")  # purple for negative, green for positive
breaks <- c(-1, -0.01, 0.01, 1)
image(
  x = 1:length(genes_of_interest),
  y = 1:length(genes_of_interest),
  z = t(adjacency_matrix)[,length(genes_of_interest):1],  # Transpose and reverse for correct orientation
  col = colorRampPalette(colors)(3),
  breaks = breaks,
  xlab = "", ylab = "",
  axes = FALSE,
  main = "Gene Regulatory Interactions (+1 activation, -1 inhibition)"
)

# Add gene names as axis labels with bold italics
axis(1, at = 1:length(genes_of_interest), labels = parse(text = bold_italic_genes), las = 2, cex.axis = 0.7)
axis(2, at = 1:length(genes_of_interest), labels = parse(text = rev(bold_italic_genes)), las = 2, cex.axis = 0.7)
dev.off()
print("Saved subnetwork_heatmap.png")

# Create SVG version
svglite("subnetwork_heatmap2.svg", width = 10, height = 8)
image(
  x = 1:length(genes_of_interest),
  y = 1:length(genes_of_interest),
  z = t(adjacency_matrix)[,length(genes_of_interest):1],
  col = colorRampPalette(colors)(3),
  breaks = breaks,
  xlab = "", ylab = "",
  axes = FALSE,
  main = "Gene Regulatory Interactions (+1 activation, -1 inhibition)"
)
axis(1, at = 1:length(genes_of_interest), labels = parse(text = bold_italic_genes), las = 2, cex.axis = 0.7)
axis(2, at = 1:length(genes_of_interest), labels = parse(text = rev(bold_italic_genes)), las = 2, cex.axis = 0.7)
dev.off()
print("Saved subnetwork_heatmap.svg")

print("Processing complete!")

# Additional output: Network visualization
library(igraph)

# Create graph from the simplified interactions
g <- graph_from_data_frame(
  d = simplified_interactions,
  vertices = genes_of_interest,
  directed = TRUE
)

# Set edge colors based on +1/-1
E(g)$color <- ifelse(E(g)$simplified_mor > 0, "#4DAF4A", "#984EA3")

# Create SVG visualization with bold and italicized gene names
svglite("subnetwork_graph2.svg", width = 10, height = 8)
plot(g, 
     layout = layout_with_fr(g),
     vertex.color = "#E41A1C",
     vertex.size = 7,
     vertex.label.cex = 0.7,
     vertex.label.font = 4,  # This sets bold italic font (4 = bold italic)
     edge.arrow.size = 0.5,
     main = "Gene Regulatory Subnetwork")
dev.off()
print("Saved subnetwork_graph.svg")

# Create circular layout visualization with bold and italicized gene names
svglite("subnetwork_graph_circular.svg", width = 10, height = 8)
plot(g, 
     layout = layout_in_circle(g),
     vertex.color = "#E41A1C",
     vertex.size = 7,
     vertex.label.cex = 0.7,
     vertex.label.font = 4,  # This sets bold italic font (4 = bold italic)
     edge.arrow.size = 0.5,
     main = "Gene Regulatory Subnetwork (Circular Layout)")
dev.off()
print("Saved subnetwork_graph_circular.svg")






### DOROTHEA SUBNETWORK VISUALIZATION ##
# Cross-Platform Gene Regulatory Network Analysis with Dorothea
# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "dorothea", "viper", "igraph", "ggplot2", "ggraph",
  "dplyr", "RColorBrewer", "pheatmap", "svglite")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("dorothea", "viper")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

options(expressions = 500000)
options(java.parameters = "-Xmx4g")

gc_timing <- 0
force_gc <- function() {
  gc_timing <<- gc_timing + 1
  if (gc_timing %% 5 == 0) {
    gc()
    cat("Running garbage collection...\n")
  }
}

safe_execute <- function(expr, error_message) {
  tryCatch(
    expr,
    error = function(e) {
      cat(paste0("\nERROR: ", error_message, "\n"))
      cat("Original error: ", e$message, "\n")
      return(NULL)
    }
  )
}

dir.create("png_output", showWarnings = FALSE)
dir.create("svg_output", showWarnings = FALSE)

# ── Gene sets ─────────────────────────────────────────────────────────────────
# Core 7 starting genes = INPUT genes (blue fill)
genes_of_interest <- c("JUNB", "ETS2", "FOS", "NR4A2", "KLF10", "CSRNP1", "HIF1A")

# INPUT genes: the 7 starting genes — blue fill
input_genes <- c("JUNB", "ETS2", "FOS", "NR4A2", "KLF10", "CSRNP1", "HIF1A")

# INTERMEDIATE genes: TFs added in Stage 2/3 (yellow fill)
# These are the TF-row genes visible in Stage 2/3 of the network figure
# that are NOT in the input set and NOT pure output targets
intermediate_genes <- c(
  "ATF2", "CREB1", "EGR1", "ELK1", "ESR1", "ETS1",
  "FLI1", "FOSL1", "FOXM1", "JUN", "MITF", "MYC",
  "NFKB1", "NR3C1", "SP1", "STAT1", "STAT3", "TP53", "TWIST1"
)

# OUTPUT genes: pure targets (pink fill) — bottom circle row in Stage 2/3
# These have no outgoing regulatory edges in the network
output_genes <- c(
  "ANGPT2", "BCL2L1", "CCL2", "CCND1", "CDKN1A",
  "CSF2", "FGF2", "IL12B", "IL1B", "MMP3", "MMP9",
  "NOS2", "PLA2G4A", "PLAUR", "TERT"
)

# Direct Vitamin D targets — red asterisk on label
# (genes with a red star marker in the network figure)
vitd_direct <- c(
  "JUNB", "ETS2", "FOS", "HIF1A",
  "EGR1", "ELK1", "MYC", "NFKB1", "SP1",
  "STAT3", "TP53", "FOXM1", "TWIST1",
  "CCL2", "CDKN1A", "IL1B", "MMP9", "TERT"
)

# Genes that act as TFs — triangle shape; all others get circle shape
tfs_in_network <- c(
  "JUNB", "ETS2", "FOS", "NR4A2", "KLF10", "HIF1A",
  "ATF2", "CREB1", "EGR1", "ELK1", "ESR1", "ETS1",
  "FLI1", "FOSL1", "FOXM1", "JUN", "MITF", "MYC",
  "NFKB1", "NR3C1", "SP1", "STAT1", "STAT3", "TP53",
  "TWIST1", "CSRNP1"
)

# ── Colours ───────────────────────────────────────────────────────────────────
COL_ACT    <- "#4DAF4A"   # green  — activation
COL_INH    <- "#E8837A"   # salmon — inhibition
COL_INPUT  <- "#377EB8"   # blue   — input genes (7 starting genes)
COL_INTERM <- "#F5C518"   # yellow — intermediate genes
COL_OUTPUT <- "#FF69B4"   # pink   — output/target genes

# ── Helper: per-node display attributes ──────────────────────────────────────
node_display <- function(gene_name) {
  # Shape: triangle for TFs, circle for non-TFs
  shape <- if (gene_name %in% tfs_in_network) 24L else 21L
  
  # Fill colour by gene role
  fill <- if      (gene_name %in% input_genes)        COL_INPUT
  else if (gene_name %in% intermediate_genes)  COL_INTERM
  else if (gene_name %in% output_genes)        COL_OUTPUT
  else                                          COL_INTERM  # fallback = intermediate
  
  # Label: add * for direct VitD targets
  label <- if (gene_name %in% vitd_direct) paste0(gene_name, "*") else gene_name
  
  list(shape = shape, fill = fill, label = label)
}

# ── Load Dorothea ─────────────────────────────────────────────────────────────
cat("Loading Dorothea human regulons database...\n")
data(dorothea_hs, package = "dorothea")

regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B"))
cat("Starting with", nrow(regulons), "high-confidence regulatory interactions from Dorothea.\n")

tf_interactions     <- regulons %>% dplyr::filter(tf     %in% genes_of_interest)
target_interactions <- regulons %>% dplyr::filter(target %in% genes_of_interest)
direct_interactions <- rbind(tf_interactions, target_interactions) %>% distinct()
cat("Found", nrow(direct_interactions), "direct interactions involving genes of interest.\n")
force_gc()

secondary_tfs <- setdiff(unique(c(direct_interactions$tf, direct_interactions$target)),
                         genes_of_interest)
cat("Found", length(secondary_tfs), "secondary transcription factors.\n")

secondary_interactions <- regulons %>%
  dplyr::filter((tf %in% genes_of_interest & target %in% secondary_tfs) |
                  (tf %in% secondary_tfs & target %in% genes_of_interest)) %>%
  arrange(desc(abs(mor))) %>%
  head(50)

all_interactions <- rbind(direct_interactions, secondary_interactions) %>% distinct()
cat("Network size:", length(unique(c(all_interactions$tf, all_interactions$target))),
    "nodes and", nrow(all_interactions), "edges\n")
force_gc()

# ── Build igraph object ───────────────────────────────────────────────────────
network <- safe_execute(
  graph_from_data_frame(
    d = all_interactions %>% dplyr::select(tf, target, mor),
    directed = TRUE,
    vertices = NULL
  ),
  "Failed to create network graph")

if (is.null(network)) stop("Network creation failed. Unable to continue.")
force_gc()

V(network)$type   <- ifelse(V(network)$name %in% genes_of_interest, "primary", "secondary")
E(network)$type   <- ifelse(E(network)$mor > 0, "activation", "inhibition")
E(network)$weight <- abs(E(network)$mor)
cat("Created network with", vcount(network), "nodes and", ecount(network), "edges.\n")

# ── Simplify for visualisation ────────────────────────────────────────────────
simplify_for_vis <- function(g, max_nodes = 40) {
  if (vcount(g) <= max_nodes) return(g)
  cat("Simplifying network for visualization (limiting to", max_nodes, "nodes)...\n")
  primary_nodes   <- V(g)[V(g)$type == "primary"]$name
  remaining_slots <- max_nodes - length(primary_nodes)
  if (remaining_slots <= 0)
    return(induced_subgraph(g, V(g)[V(g)$name %in% primary_nodes]))
  secondary_nodes <- V(g)[V(g)$type == "secondary"]
  if (length(secondary_nodes) == 0) return(g)
  secondary_data  <- data.frame(name   = secondary_nodes$name,
                                degree = degree(g, v = secondary_nodes))
  secondary_data  <- secondary_data[order(secondary_data$degree, decreasing = TRUE), ]
  top_secondaries <- head(secondary_data$name, remaining_slots)
  return(induced_subgraph(g, V(g)[V(g)$name %in% c(primary_nodes, top_secondaries)]))
}

# ── plot_network ──────────────────────────────────────────────────────────
plot_network <- function(g, layout_type, title) {
  g <- simplify_for_vis(g)
  cat("Calculating", layout_type, "layout...\n")
  
  tryCatch({
    # ── Layout ──────────────────────────────────────────────────────────────
    layout_coords <- if      (layout_type == "fr")       layout_with_fr(g)
    else if (layout_type == "kk")       layout_with_kk(g)
    else if (layout_type == "circle")   layout_in_circle(g)
    else if (layout_type == "graphopt") layout_with_graphopt(g)
    else                               layout_with_fr(g)
    
    # ── Per-node display attributes ──────────────────────────────────────────
    disp        <- lapply(V(g)$name, node_display)
    node_shapes <- sapply(disp, `[[`, "shape")
    node_fills  <- sapply(disp, `[[`, "fill")
    # Base label WITHOUT asterisk — asterisk drawn separately in red
    node_labels_base <- V(g)$name
    # Whether each node gets an asterisk overlay
    node_vitd   <- V(g)$name %in% vitd_direct
    
    NODE_SIZE <- 8   # uniform size for ALL nodes (ggplot size units)
    
    vertices_df <- data.frame(
      name        = V(g)$name,
      type        = V(g)$type,
      x           = layout_coords[, 1],
      y           = layout_coords[, 2],
      node_shape  = node_shapes,
      node_fill   = node_fills,
      node_label  = node_labels_base,
      node_vitd   = node_vitd,
      stringsAsFactors = FALSE
    )
    
    # ── Edges ────────────────────────────────────────────────────────────────
    el <- as_edgelist(g)
    edges_df <- data.frame(
      from       = el[, 1],
      to         = el[, 2],
      etype      = E(g)$type,
      edge_color = ifelse(E(g)$type == "activation", COL_ACT, COL_INH),
      stringsAsFactors = FALSE
    )
    
    coords_df <- vertices_df %>% dplyr::select(name, x, y)
    edges_with_coords <- edges_df %>%
      left_join(coords_df, by = c("from" = "name")) %>%
      dplyr::rename(x_from = x, y_from = y) %>%
      left_join(coords_df, by = c("to" = "name")) %>%
      dplyr::rename(x_to = x, y_to = y)
    
    # ── Circle-layout normalisation & radial label placement ─────────────────
    if (layout_type == "circle") {
      rng <- range(c(vertices_df$x, vertices_df$y))
      vertices_df$x <- (vertices_df$x - rng[1]) / diff(rng)
      vertices_df$y <- (vertices_df$y - rng[1]) / diff(rng)
      
      coords_df <- vertices_df %>% dplyr::select(name, x, y)
      edges_with_coords <- edges_df %>%
        left_join(coords_df, by = c("from" = "name")) %>%
        dplyr::rename(x_from = x, y_from = y) %>%
        left_join(coords_df, by = c("to" = "name")) %>%
        dplyr::rename(x_to = x, y_to = y)
      
      # Shorten each edge so arrow tip stops BEFORE the destination node centre.
      # NODE_SIZE ggplot units ≈ radius in data units needs estimating.
      # For a [0,1] normalised circle layout, node radius ≈ 0.028 data units
      # (empirically: ~5% of the 0.5-radius circle).
      node_r <- 0.038   # data-unit offset calibrated for NODE_SIZE=8
      dx <- edges_with_coords$x_to - edges_with_coords$x_from
      dy <- edges_with_coords$y_to - edges_with_coords$y_from
      dist_e <- sqrt(dx^2 + dy^2)
      # Pull the tip back by node_r along the edge direction
      edges_with_coords$x_to_adj <- edges_with_coords$x_to   - node_r * dx / dist_e
      edges_with_coords$y_to_adj <- edges_with_coords$y_to   - node_r * dy / dist_e
      # Also pull the tail forward slightly so it starts outside the source node
      edges_with_coords$x_from_adj <- edges_with_coords$x_from + node_r * dx / dist_e
      edges_with_coords$y_from_adj <- edges_with_coords$y_from + node_r * dy / dist_e
      
      cx <- mean(vertices_df$x); cy <- mean(vertices_df$y)
      vertices_df$angle   <- atan2(vertices_df$y - cy, vertices_df$x - cx)
      push <- 0.06
      vertices_df$label_x <- vertices_df$x + push * cos(vertices_df$angle)
      vertices_df$label_y <- vertices_df$y + push * sin(vertices_df$angle)
      vertices_df$hjust   <- ifelse(cos(vertices_df$angle) < 0, 1, 0)
      vertices_df$vjust   <- ifelse(sin(vertices_df$angle) >  0.7, 0,
                                    ifelse(sin(vertices_df$angle) < -0.7, 1, 0.5))
    }
    
    cat("Creating", layout_type, "network visualization...\n")
    
    # Layer order (bottom to top):
    #   1. Edges  — arrowheads must be visible, so edges drawn first
    #   2. Nodes  — icons cover the edge endpoints, revealing arrowheads between nodes
    #   3. Red *  — drawn on top of node icon, offset to top-right corner
    #   4. Labels — gene names outermost
    
    # ── 1. Edges ──────────────────────────────────────────────────────────────
    if (layout_type == "circle") {
      p <- ggplot() +
        geom_segment(
          data      = edges_with_coords,
          aes(x = x_from_adj, y = y_from_adj,
              xend = x_to_adj, yend = y_to_adj,
              colour = edge_color),
          arrow     = arrow(length = unit(6, "pt"), type = "closed"),
          linewidth = 0.5,
          alpha     = 0.80,
          lineend   = "butt"
        )
    } else {
      p <- ggplot() +
        geom_segment(
          data      = edges_with_coords,
          aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
              colour = edge_color),
          arrow     = arrow(length = unit(6, "pt"), type = "closed"),
          linewidth = 0.5,
          alpha     = 0.80,
          lineend   = "butt"
        )
    }
    p <- p + scale_colour_identity()
    
    # ── 2. Node icons (uniform size for every node) ───────────────────────────
    for (pch_val in unique(vertices_df$node_shape)) {
      sub <- vertices_df[vertices_df$node_shape == pch_val, ]
      p <- p + geom_point(
        data   = sub,
        aes(x = x, y = y, fill = node_fill),
        shape  = pch_val,
        size   = NODE_SIZE,
        colour = "black",
        stroke = 0.8
      )
    }
    p <- p + scale_fill_identity()
    
    # ── 3. Red asterisk ON TOP of node, offset to top-right corner ───────────
    # ast_offset pushes the * outside the filled area of the node so it is
    # always visible against any fill colour.
    vitd_nodes <- vertices_df[vertices_df$node_vitd, ]
    if (nrow(vitd_nodes) > 0) {
      ast_off <- 0.008   # small offset: * sits at icon edge, not far away
      p <- p + geom_text(
        data     = vitd_nodes,
        aes(x = x + ast_off, y = y + ast_off),
        label    = "*",
        size     = 7,
        color    = "red",
        fontface = "bold",
        vjust    = 0,
        hjust    = 0
      )
    }
    
    # ── 4. Gene name labels ───────────────────────────────────────────────────
    if (layout_type == "circle") {
      p <- p + geom_text(
        data     = vertices_df,
        aes(x = label_x, y = label_y, label = node_label,
            hjust = hjust, vjust = vjust),
        size     = 3.5,
        color    = "black",
        fontface = "bold.italic"
      )
    } else {
      p <- p + geom_text(
        data     = vertices_df,
        aes(x = x, y = y, label = node_label),
        size     = 3.5,
        vjust    = -1.2,
        color    = "black",
        fontface = "bold.italic"
      )
    }
    
    p <- p +
      labs(title = title) +
      theme_minimal() +
      theme(
        panel.grid      = element_blank(),
        axis.text       = element_blank(),
        axis.title      = element_blank(),
        axis.ticks      = element_blank(),
        legend.position = "none",
        plot.title      = element_blank()
      )
    
    if (layout_type == "circle") {
      p <- p + coord_fixed(ratio = 1, clip = "off") +
        theme(plot.margin = margin(2, 2, 2, 20))
    }
    
    force_gc()
    return(p)
    
  }, error = function(e) {
    cat("Failed to create", layout_type, "network plot:", e$message, "\n")
    return(NULL)
  })
}

# ── Manual legend builder ────────────────────────────────────────────────────
# The main plot uses scale_colour_identity()/scale_fill_identity() so ggplot
# produces no guide-box.  Build a standalone legend plot instead.
build_legend <- function() {
  items <- list(
    list(y=8, type="line", col=COL_ACT,     label="Activation edge"),
    list(y=7, type="line", col=COL_INH,     label="Inhibition (repression) edge"),
    list(y=6, type="pt", shape=24L, fill=COL_INPUT,  label="Input gene  (7 starting genes, blue triangle/circle)"),
    list(y=5, type="pt", shape=21L, fill=COL_INTERM, label="Intermediate gene  (yellow)"),
    list(y=4, type="pt", shape=21L, fill=COL_OUTPUT, label="Output gene  (pink circle)"),
    list(y=3, type="pt", shape=24L, fill=COL_INTERM, label="TF = triangle shape"),
    list(y=2, type="pt", shape=21L, fill=COL_INTERM, label="Non-TF target = circle shape"),
    list(y=1, type="pt", shape=24L, fill=COL_INPUT,  label="* = Direct Vitamin D target")
  )
  p_leg <- ggplot() + theme_void() +
    theme(plot.margin = margin(12, 12, 12, 12))
  for (it in items) {
    if (it$type == "line") {
      p_leg <- p_leg +
        annotate("segment", x=0.5, xend=0.85, y=it$y, yend=it$y,
                 colour=it$col, linewidth=2) +
        annotate("text", x=0.95, y=it$y, label=it$label,
                 hjust=0, size=3.5, fontface="bold")
    } else {
      p_leg <- p_leg +
        annotate("point", x=0.65, y=it$y,
                 shape=it$shape, fill=it$fill, colour="black",
                 size=4, stroke=0.6) +
        annotate("text", x=0.95, y=it$y, label=it$label,
                 hjust=0, size=3.5, fontface="bold")
    }
  }
  p_leg + xlim(0.4, 6.0) + ylim(0.5, 8.5)
}

# ── Save helper ───────────────────────────────────────────────────────────────
save_plots <- function(plot_obj, basename, width = 8, height = 6) {
  if (is.null(plot_obj)) {
    cat("Cannot save plot:", basename, "(plot object is NULL)\n"); return()
  }
  png_file <- file.path("png_output", paste0(basename, ".png"))
  safe_execute(
    ggsave(png_file, plot_obj, width = width, height = height, dpi = 100),
    paste("Failed to save PNG:", png_file))
  cat("Saved PNG:", png_file, "\n")
  
  svg_file <- file.path("svg_output", paste0(basename, ".svg"))
  safe_execute(
    ggsave(svg_file, plot_obj, width = width, height = height, device = svglite),
    paste("Failed to save SVG:", svg_file))
  cat("Saved SVG:", svg_file, "\n")
  
  # Build and save the manual legend
  legend_obj <- safe_execute(build_legend(),
                             paste("Failed to build legend for", basename))
  if (!is.null(legend_obj)) {
    safe_execute(
      ggsave(file.path("png_output", paste0(basename, "_legend.png")),
             legend_obj, width = 6.5, height = 4, dpi = 100),
      "Failed to save legend PNG")
    cat("Saved legend PNG:", file.path("png_output", paste0(basename, "_legend.png")), "\n")
    safe_execute(
      ggsave(file.path("svg_output", paste0(basename, "_legend.svg")),
             legend_obj, width = 6.5, height = 4, device = svglite),
      "Failed to save legend SVG")
    cat("Saved legend SVG:", file.path("svg_output", paste0(basename, "_legend.svg")), "\n")
  }
  force_gc()
}

# ── Produce all layouts ───────────────────────────────────────────────────────
cat("\nCreating network visualizations (this may take some time)...\n")

cat("\nCreating Fruchterman-Reingold layout...\n")
p1 <- safe_execute(plot_network(network, "fr", "Force-directed Layout (Fruchterman-Reingold)"),
                   "Failed to create force-directed layout")
if (!is.null(p1)) save_plots(p1, "network_fr_layout")
rm(p1); force_gc()

cat("\nCreating Kamada-Kawai layout...\n")
p2 <- safe_execute(plot_network(network, "kk", "Kamada-Kawai Layout"),
                   "Failed to create Kamada-Kawai layout")
if (!is.null(p2)) save_plots(p2, "network_kk_layout")
rm(p2); force_gc()

cat("\nCreating Circular layout...\n")
p3 <- safe_execute(plot_network(network, "circle", "Circular Layout"),
                   "Failed to create circular layout")
if (!is.null(p3)) save_plots(p3, "network_circle_layout", width = 9, height = 9)
rm(p3); force_gc()

cat("\nCreating Graph Optimization layout...\n")
p4 <- safe_execute(plot_network(network, "graphopt", "Graph Optimization Layout"),
                   "Failed to create graph optimization layout")
if (!is.null(p4)) save_plots(p4, "network_graphopt_layout")
rm(p4); force_gc()

# ── Network analysis ──────────────────────────────────────────────────────────
cat("\nCalculating network centrality measures...\n")
centrality_measures <- safe_execute({
  data.frame(
    node        = V(network)$name,
    type        = V(network)$type,
    degree      = degree(network),
    in_degree   = degree(network, mode = "in"),
    out_degree  = degree(network, mode = "out"),
    betweenness = betweenness(network),
    closeness   = closeness(network),
    eigenvector = eigen_centrality(network)$vector
  )
}, "Failed to calculate centrality measures")

if (!is.null(centrality_measures)) {
  centrality_measures <- centrality_measures[order(centrality_measures$degree,
                                                   decreasing = TRUE), ]
  cat("\n\nTop 10 nodes by degree centrality:\n")
  print(head(centrality_measures, 10))
  write.csv(centrality_measures, "node_centrality_measures.csv", row.names = FALSE)
  cat("Saved: node_centrality_measures.csv\n")
}
force_gc()

# ── Community detection ───────────────────────────────────────────────────────
cat("\nPerforming community detection...\n")
communities <- safe_execute(
  cluster_louvain(as.undirected(network)),
  "Failed to perform community detection")
if (!is.null(communities)) {
  V(network)$community <- membership(communities)
  cat("Identified", length(unique(membership(communities))), "communities\n")
}
force_gc()

# ── Adjacency matrix & heatmap (salmon inhibition) ────────────────────────────
cat("\nCreating adjacency matrix for direct interactions...\n")
direct_network <- safe_execute(
  graph_from_data_frame(
    d = direct_interactions %>%
      dplyr::filter(tf %in% genes_of_interest & target %in% genes_of_interest) %>%
      dplyr::select(tf, target, mor),
    directed = TRUE,
    vertices = genes_of_interest
  ),
  "Failed to create direct interactions network")

if (!is.null(direct_network)) {
  adj_matrix <- safe_execute(
    as_adjacency_matrix(direct_network, attr = "mor", sparse = FALSE),
    "Failed to create adjacency matrix")
  
  if (!is.null(adj_matrix)) {
    cat("Creating regulatory heatmap...\n")
    pheatmap_args <- list(
      mat             = adj_matrix,
      main            = "Regulatory Interactions Between Target Genes",
      # Salmon for inhibition (-1), white for 0, green for activation (+1)
      color           = colorRampPalette(c("#E8837A", "white", "#4DAF4A"))(100),
      breaks          = seq(-1, 1, length.out = 101),
      cluster_rows    = TRUE,
      cluster_cols    = TRUE,
      display_numbers = TRUE,
      number_color    = "black",
      fontsize_number = 12,
      fontsize_row    = 14,
      fontsize_col    = 14,
      fontsize        = 14,
      angle_col       = 45
    )
    draw_pheatmap_bolditalic <- function(args) {
      old_par <- par(font = 4, col = "black")
      on.exit(par(old_par))
      do.call(pheatmap, args)
    }
    png_heatmap <- file.path("png_output", "regulatory_heatmap.png")
    safe_execute({
      png(png_heatmap, width = 800, height = 600)
      draw_pheatmap_bolditalic(pheatmap_args)
      dev.off()
      cat("Saved:", png_heatmap, "\n")
    }, "Failed to create PNG heatmap")
    
    svg_heatmap <- file.path("svg_output", "regulatory_heatmap.svg")
    safe_execute({
      svglite(svg_heatmap, width = 8, height = 6)
      draw_pheatmap_bolditalic(pheatmap_args)
      dev.off()
      cat("Saved:", svg_heatmap, "\n")
    }, "Failed to create SVG heatmap")
  }
}
force_gc()

# ── Network statistics ────────────────────────────────────────────────────────
cat("\n\nNetwork Statistics:\n")
cat("Number of nodes:", vcount(network), "\n")
cat("Number of edges:", ecount(network), "\n")
cat("Network density:", safe_execute(edge_density(network), "Failed"), "\n")
cat("Network diameter:", safe_execute(diameter(network), "Failed"), "\n")
cat("Average path length:", safe_execute(mean_distance(network), "Failed"), "\n")
if (!is.null(communities))
  cat("Number of communities:", length(unique(V(network)$community)), "\n\n")

# ── Export ────────────────────────────────────────────────────────────────────
cat("\nExporting network data...\n")
safe_execute({
  write.csv(as_data_frame(network), "network_edges.csv", row.names = FALSE)
  cat("Saved: network_edges.csv\n")
}, "Failed to write network edges CSV")
force_gc()

cat("\nCreating direct interactions subnetwork...\n")
direct_subnet <- safe_execute(
  induced_subgraph(network, which(V(network)$name %in% genes_of_interest)),
  "Failed to create direct interactions subnetwork")
if (!is.null(direct_subnet)) {
  p5 <- safe_execute(
    plot_network(direct_subnet, "circle", "Direct Interactions Between Genes of Interest"),
    "Failed to create direct interactions network plot")
  if (!is.null(p5)) save_plots(p5, "direct_interactions_network")
  rm(p5)
}
force_gc()

cat("\nCreating ego networks for each gene...\n")
for (gene in genes_of_interest) {
  if (gene %in% V(network)$name) {
    cat("Processing ego network for", gene, "...\n")
    ego_net <- safe_execute(
      make_ego_graph(network, order = 1,
                     nodes = which(V(network)$name == gene), mode = "all")[[1]],
      paste("Failed to create ego network for", gene))
    if (!is.null(ego_net)) {
      p_ego <- safe_execute(
        plot_network(ego_net, "fr", paste0("Ego Network for ", gene)),
        paste("Failed to plot ego network for", gene))
      if (!is.null(p_ego)) { save_plots(p_ego, paste0("ego_network_", gene)); rm(p_ego) }
    }
    force_gc()
  }
}

cat("\nCreating TF-target interaction table...\n")
tf_target_table <- safe_execute({
  all_interactions %>%
    filter(tf %in% genes_of_interest | target %in% genes_of_interest) %>%
    mutate(interaction_type = ifelse(mor > 0, "Activation", "Inhibition"),
           evidence = paste0(confidence, " (", mor, ")")) %>%
    dplyr::select(TF = tf, Target = target, Interaction = interaction_type, Evidence = evidence) %>%
    arrange(TF, Target)
}, "Failed to create TF-target interaction table")
if (!is.null(tf_target_table)) {
  safe_execute({
    write.csv(tf_target_table, "tf_target_interactions.csv", row.names = FALSE)
    cat("Saved: tf_target_interactions.csv\n")
  }, "Failed to write TF-target interactions CSV")
}



### HEATMAP ##
library(igraph)
library(dplyr)
library(svglite)

# ── Gene list (all 40 genes) ──────────────────────────────────────────────────
genes_of_interest <- c(
  "ETS2", "FOS", "HIF1A", "ATF2", "CREB1", "EGR1", "ELK1", "ESR1", "ETS1",
  "FLI1", "FOXM1", "JUN", "MITF", "MYC", "NFKB1", "SP1", "STAT1", "STAT3",
  "TP53", "TWIST1", "ANGPT2", "CSF2", "FOSL1", "IL12B", "MMP3", "MMP9",
  "TERT", "BCL2L1", "CCL2", "CCND1", "CDKN1A", "FGF2", "IL1B", "NOS2",
  "NR3C1", "PLA2G4A", "PLAUR", "NR4A2", "KLF10", "JUNB"
)

# ── Load Dorothea & build adjacency matrix ────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("dorothea",   quietly = TRUE)) BiocManager::install("dorothea")
library(dorothea)
data(dorothea_hs, package = "dorothea")

regulons <- dorothea_hs %>% dplyr::filter(confidence %in% c("A", "B"))

subnetwork_interactions <- regulons %>%
  dplyr::filter(tf %in% genes_of_interest & target %in% genes_of_interest) %>%
  dplyr::select(tf, target, mor)

simplified_interactions <- subnetwork_interactions %>%
  dplyr::mutate(simplified_mor = ifelse(mor > 0, 1, -1)) %>%
  dplyr::select(tf, target, simplified_mor)

n <- length(genes_of_interest)
adjacency_matrix <- matrix(0, nrow = n, ncol = n,
                           dimnames = list(genes_of_interest, genes_of_interest))
for (i in seq_len(nrow(simplified_interactions))) {
  r <- simplified_interactions[i, ]
  adjacency_matrix[r$tf, r$target] <- r$simplified_mor
}

# ── Sort genes so that the white-space block is contiguous at bottom-right ────
# Strategy:
#   1. Genes with ANY interaction (row or column nonzero) come first,
#      sorted by out-degree DESC then in-degree DESC → the active block
#   2. Genes with NO interactions at all come last → the white block
out_deg   <- rowSums(adjacency_matrix != 0)
in_deg    <- colSums(adjacency_matrix != 0)
total_deg <- out_deg + in_deg

has_interaction <- total_deg > 0
active_genes  <- genes_of_interest[has_interaction]
silent_genes  <- genes_of_interest[!has_interaction]

# Sort active block: primary = out_deg DESC, secondary = in_deg DESC
active_order <- order(out_deg[active_genes], in_deg[active_genes], decreasing = TRUE)
active_sorted <- active_genes[active_order]

# Silent genes keep their original order (all-white rows, doesn't matter)
genes_sorted <- c(active_sorted, silent_genes)

# Reorder adjacency matrix rows AND columns identically
adjacency_matrix_sorted <- adjacency_matrix[genes_sorted, genes_sorted]
n_sorted <- length(genes_sorted)

# ── Colours ───────────────────────────────────────────────────────────────────
colors <- c("#E8837A", "white", "#4DAF4A")   # salmon | white | green
breaks <- c(-1, -0.01, 0.01, 1)

# ── Draw function ─────────────────────────────────────────────────────────────
draw_heatmap <- function() {
  # Extra-large left & bottom margins to prevent label clipping
  # mar = c(bottom, left, top, right) in lines
  par(
    mar  = c(22, 20, 8, 4),
    bty  = "n",
    xpd  = TRUE    # allow drawing outside plot region for labels
  )
  
  image(
    x      = 1:n_sorted,
    y      = 1:n_sorted,
    z      = t(adjacency_matrix_sorted)[, n_sorted:1],
    col    = colorRampPalette(colors)(3),
    breaks = breaks,
    xlab   = "",
    ylab   = "",
    axes   = FALSE,
    main   = ""
  )
  
  title(
    main      = "Gene Regulatory Interactions (+1 activation, -1 inhibition)",
    cex.main  = 3.5,
    font.main = 2,
    col.main  = "black"
  )
  
  # X-axis (columns = targets, left to right)
  axis(
    side      = 1,
    at        = 1:n_sorted,
    labels    = genes_sorted,
    las       = 2,
    cex.axis  = 3.8,
    font.axis = 4,
    col.axis  = "black",
    col       = "black",
    lwd       = 5,
    lwd.ticks = 5,
    mgp       = c(3, 2.5, 0)   # push labels further from axis to avoid clipping
  )
  
  # Y-axis (rows = TFs, top to bottom — reversed because image() plots y bottom-up)
  axis(
    side      = 2,
    at        = 1:n_sorted,
    labels    = rev(genes_sorted),
    las       = 2,
    cex.axis  = 3.8,
    font.axis = 4,
    col.axis  = "black",
    col       = "black",
    lwd       = 5,
    lwd.ticks = 5,
    mgp       = c(3, 2.5, 0)
  )
}

# ── Save outputs ──────────────────────────────────────────────────────────────
# Increase canvas size and resolution so labels have breathing room
png("subnetwork_heatmap.png", width = 10000, height = 9500, res = 150)
draw_heatmap()
dev.off()
print("Saved subnetwork_heatmap.png")

svglite("subnetwork_heatmap.svg", width = 56, height = 52)
draw_heatmap()
dev.off()
print("Saved subnetwork_heatmap.svg")

# Diagnostic: show the final gene order and degree counts
cat("\nActive genes (sorted, with interactions):\n")
cat(paste(active_sorted, collapse = " > "), "\n")
cat("\nSilent genes (no interactions, white block):\n")
cat(paste(silent_genes, collapse = ", "), "\n")
cat("\nDegree summary:\n")
print(data.frame(
  gene      = genes_sorted,
  out_deg   = out_deg[genes_sorted],
  in_deg    = in_deg[genes_sorted],
  total_deg = total_deg[genes_sorted]
))

print("Done!")






###KERNEL REDUCTION FINAL TAKE
# Kernel Reduction for Gene Regulatory Network with Dual Effects and Vitamin D Input
# Modified to aggregate effects and handle three scenarios

# Read the adjacency matrix with Vitamin D row
data <- read.csv("adj_matrix_with_VitD.csv", sep=";", check.names=FALSE)
rownames(data) <- data$Genes
data$Genes <- NULL

# Extract the Vitamin D row and store it separately
vitD_row <- data["VitD", ]
data <- data[rownames(data) != "VitD", ]

# Convert to matrix
adj_matrix <- as.matrix(data)
vitD_vector <- as.numeric(as.matrix(vitD_row))

# Store effects as lists to handle dual effects (both positive and negative)
# Each cell will contain a list of all effects between that pair
effects_list <- list()
n_genes <- nrow(adj_matrix)
genes <- rownames(adj_matrix)

# Initialize effects list
for(i in 1:n_genes) {
  effects_list[[genes[i]]] <- list()
  for(j in 1:n_genes) {
    effects_list[[genes[i]]][[genes[j]]] <- c()
    if(adj_matrix[i, j] != 0) {
      effects_list[[genes[i]]][[genes[j]]] <- c(adj_matrix[i, j])
    }
  }
}

# Store VitD effects as lists too
vitD_effects <- list()
for(i in 1:n_genes) {
  vitD_effects[[genes[i]]] <- c()
  if(vitD_vector[i] != 0) {
    vitD_effects[[genes[i]]] <- c(vitD_vector[i])
  }
}

# Function to aggregate effects (sum positive and negative separately)
aggregate_effects <- function(effects) {
  if(length(effects) == 0) return(c())
  
  pos_sum <- sum(effects[effects > 0])
  neg_sum <- sum(effects[effects < 0])
  
  result <- c()
  if(pos_sum != 0) result <- c(result, pos_sum)
  if(neg_sum != 0) result <- c(result, neg_sum)
  
  return(result)
}

# Function to calculate in-degree and out-degree
calculate_degrees <- function(effects_list, vitD_effects, active_genes) {
  in_degree <- rep(0, length(active_genes))
  out_degree <- rep(0, length(active_genes))
  names(in_degree) <- active_genes
  names(out_degree) <- active_genes
  
  for(i in 1:length(active_genes)) {
    gene_i <- active_genes[i]
    
    # Count outgoing edges
    for(j in 1:length(active_genes)) {
      gene_j <- active_genes[j]
      if(length(effects_list[[gene_i]][[gene_j]]) > 0) {
        out_degree[i] <- out_degree[i] + 1
      }
    }
    
    # Count incoming edges from other genes
    for(j in 1:length(active_genes)) {
      gene_j <- active_genes[j]
      if(length(effects_list[[gene_j]][[gene_i]]) > 0) {
        in_degree[i] <- in_degree[i] + 1
      }
    }
    
    # Add VitD input if exists
    if(length(vitD_effects[[gene_i]]) > 0) {
      in_degree[i] <- in_degree[i] + 1
    }
  }
  
  return(list(in_degree = in_degree, out_degree = out_degree))
}

# Function to perform kernel reduction
kernel_reduction <- function(effects_list, vitD_effects, genes, 
                             force_remove = NULL, force_keep = NULL) {
  active_genes <- genes
  removed <- c()
  iteration <- 1
  
  # Initial VitD targets
  initial_vitD_targets <- names(vitD_effects)[sapply(vitD_effects, function(x) length(x) > 0)]
  cat("Initial VitD targets:", paste(initial_vitD_targets, collapse=", "), "\n\n")
  
  while(TRUE) {
    # Calculate degrees
    degrees <- calculate_degrees(effects_list, vitD_effects, active_genes)
    in_degree <- degrees$in_degree
    out_degree <- degrees$out_degree
    
    # Find genes with in-degree or out-degree <= 2
    candidates <- c()
    for(gene in active_genes) {
      # Skip force_keep genes
      if(!is.null(force_keep) && gene %in% force_keep) next
      
      if(in_degree[gene] <= 2 || out_degree[gene] <= 2) {
        candidates <- c(candidates, gene)
      }
    }
    
    # Add forced removal genes if they're still active
    if(!is.null(force_remove)) {
      for(gene in force_remove) {
        if(gene %in% active_genes && !(gene %in% candidates)) {
          candidates <- c(candidates, gene)
        }
      }
    }
    
    if(length(candidates) == 0) {
      break
    }
    
    cat(paste("Iteration", iteration, ": Removing", length(candidates), "genes\n"))
    for(gene in candidates) {
      cat(paste("  Removing gene:", gene,
                " (in-degree:", in_degree[gene],
                ", out-degree:", out_degree[gene], ")\n"))
    }
    
    # Process each candidate gene
    for(gene_to_remove in candidates) {
      # Find incoming and outgoing connections
      incoming_genes <- c()
      outgoing_genes <- c()
      
      for(other_gene in active_genes) {
        if(other_gene == gene_to_remove) next
        
        if(length(effects_list[[other_gene]][[gene_to_remove]]) > 0) {
          incoming_genes <- c(incoming_genes, other_gene)
        }
        if(length(effects_list[[gene_to_remove]][[other_gene]]) > 0) {
          outgoing_genes <- c(outgoing_genes, other_gene)
        }
      }
      
      # Pass through connections from incoming to outgoing genes
      for(in_gene in incoming_genes) {
        for(out_gene in outgoing_genes) {
          if(in_gene != out_gene) {
            # Calculate all combinations of pass-through effects
            in_effects <- effects_list[[in_gene]][[gene_to_remove]]
            out_effects <- effects_list[[gene_to_remove]][[out_gene]]
            
            for(in_eff in in_effects) {
              for(out_eff in out_effects) {
                pass_through <- in_eff * out_eff
                # Append the pass-through effect to the list
                effects_list[[in_gene]][[out_gene]] <- c(effects_list[[in_gene]][[out_gene]], pass_through)
              }
            }
          }
        }
      }
      
      # Pass through VitD influence if this gene is a VitD target
      vitD_in_effects <- vitD_effects[[gene_to_remove]]
      if(length(vitD_in_effects) > 0) {
        cat(paste("    Passing VitD effect through", gene_to_remove,
                  "with VitD effects:", paste(vitD_in_effects, collapse=", "), "\n"))
        
        for(out_gene in outgoing_genes) {
          out_effects <- effects_list[[gene_to_remove]][[out_gene]]
          
          for(vitD_eff in vitD_in_effects) {
            for(out_eff in out_effects) {
              pass_through <- vitD_eff * out_eff
              # Append the pass-through VitD effect to the list
              vitD_effects[[out_gene]] <- c(vitD_effects[[out_gene]], pass_through)
              cat(paste("      VitD ->", out_gene, ": added effect", pass_through, "\n"))
            }
          }
        }
      }
      
      # Remove this gene from effects_list
      for(other_gene in active_genes) {
        effects_list[[other_gene]][[gene_to_remove]] <- c()
        effects_list[[gene_to_remove]][[other_gene]] <- c()
      }
      
      removed <- c(removed, gene_to_remove)
    }
    
    # Update active genes
    active_genes <- setdiff(active_genes, candidates)
    
    # Print current VitD targets
    current_vitD_targets <- names(vitD_effects)[sapply(vitD_effects[active_genes], function(x) length(x) > 0)]
    cat(paste("VitD targets after iteration", iteration, ":", paste(current_vitD_targets, collapse=", "), "\n\n"))
    
    iteration <- iteration + 1
  }
  
  return(list(
    effects_list = effects_list,
    vitD_effects = vitD_effects,
    removed_genes = removed,
    remaining_genes = active_genes
  ))
}

# Function to create output files
create_output_files <- function(result, prefix) {
  remaining_genes <- result$remaining_genes
  n_remaining <- length(remaining_genes)
  
  # Initialize matrices
  combined_matrix <- matrix(0, nrow=n_remaining, ncol=n_remaining)
  positive_matrix <- matrix(0, nrow=n_remaining, ncol=n_remaining)
  negative_matrix <- matrix(0, nrow=n_remaining, ncol=n_remaining)
  interactions_matrix <- matrix("0", nrow=n_remaining, ncol=n_remaining)
  
  rownames(combined_matrix) <- remaining_genes
  colnames(combined_matrix) <- remaining_genes
  rownames(positive_matrix) <- remaining_genes
  colnames(positive_matrix) <- remaining_genes
  rownames(negative_matrix) <- remaining_genes
  colnames(negative_matrix) <- remaining_genes
  rownames(interactions_matrix) <- remaining_genes
  colnames(interactions_matrix) <- remaining_genes
  
  # Fill matrices with aggregated effects
  for(i in 1:n_remaining) {
    for(j in 1:n_remaining) {
      gene_i <- remaining_genes[i]
      gene_j <- remaining_genes[j]
      
      effects <- result$effects_list[[gene_i]][[gene_j]]
      
      if(length(effects) > 0) {
        # Aggregate effects
        agg_effects <- aggregate_effects(effects)
        
        # Store aggregated effects as comma-separated
        interactions_matrix[i, j] <- paste(agg_effects, collapse=", ")
        
        # For combined matrix, sum all effects
        combined_matrix[i, j] <- sum(agg_effects)
        
        # For positive/negative matrices
        pos_effects <- agg_effects[agg_effects > 0]
        neg_effects <- agg_effects[agg_effects < 0]
        
        if(length(pos_effects) > 0) {
          positive_matrix[i, j] <- sum(pos_effects)
        }
        if(length(neg_effects) > 0) {
          negative_matrix[i, j] <- sum(neg_effects)
        }
      }
    }
  }
  
  # Create VitD rows with aggregated effects
  vitD_row_combined <- rep(0, n_remaining)
  vitD_row_positive <- rep(0, n_remaining)
  vitD_row_negative <- rep(0, n_remaining)
  vitD_row_interactions <- rep("0", n_remaining)
  
  for(i in 1:n_remaining) {
    gene <- remaining_genes[i]
    effects <- result$vitD_effects[[gene]]
    
    if(length(effects) > 0) {
      # Aggregate effects
      agg_effects <- aggregate_effects(effects)
      
      vitD_row_interactions[i] <- paste(agg_effects, collapse=", ")
      vitD_row_combined[i] <- sum(agg_effects)
      
      pos_effects <- agg_effects[agg_effects > 0]
      neg_effects <- agg_effects[agg_effects < 0]
      
      if(length(pos_effects) > 0) {
        vitD_row_positive[i] <- sum(pos_effects)
      }
      if(length(neg_effects) > 0) {
        vitD_row_negative[i] <- sum(neg_effects)
      }
    }
  }
  
  # Append VitD row to matrices
  combined_with_vitD <- rbind(combined_matrix, "VitD" = vitD_row_combined)
  positive_with_vitD <- rbind(positive_matrix, "VitD" = vitD_row_positive)
  negative_with_vitD <- rbind(negative_matrix, "VitD" = vitD_row_negative)
  interactions_with_vitD <- rbind(interactions_matrix, "VitD" = vitD_row_interactions)
  
  # Write CSV files
  combined_df <- data.frame(Genes = c(remaining_genes, "VitD"), combined_with_vitD, check.names=FALSE)
  write.csv(combined_df, paste0(prefix, "_combined_matrix_with_vitD.csv"), row.names=FALSE)
  
  positive_df <- data.frame(Genes = c(remaining_genes, "VitD"), positive_with_vitD, check.names=FALSE)
  write.csv(positive_df, paste0(prefix, "_positive_matrix_with_vitD.csv"), row.names=FALSE)
  
  negative_df <- data.frame(Genes = c(remaining_genes, "VitD"), negative_with_vitD, check.names=FALSE)
  write.csv(negative_df, paste0(prefix, "_negative_matrix_with_vitD.csv"), row.names=FALSE)
  
  interactions_df <- data.frame(Genes = c(remaining_genes, "VitD"), interactions_with_vitD, check.names=FALSE)
  write.csv(interactions_df, paste0(prefix, "_interactions_matrix_with_vitD.csv"), row.names=FALSE)
  
  # Generate equations
  generate_equations(interactions_with_vitD, prefix)
  
  cat(paste0("\n", prefix, " files saved successfully!\n\n"))
}

# Function to create Hill-type equations
generate_equations <- function(interactions_with_vitD, prefix) {
  # Transpose matrix
  transposed_matrix <- t(interactions_with_vitD)
  write.csv(transposed_matrix, paste0(prefix, "_transposed_matrix_with_vitD.csv"))
  
  # Save as Excel file
  library(openxlsx)
  write.xlsx(transposed_matrix, paste0(prefix, "_transposed_matrix_with_vitD.xlsx"))
  
  # Create equations
  equations <- list()
  for (i in 1:nrow(transposed_matrix)) {
    gene_name <- rownames(transposed_matrix)[i]
    regulators_row <- transposed_matrix[i, ]
    equations[[gene_name]] <- create_gene_equation(gene_name, regulators_row)
  }
  
  # Save equations
  writeLines(unlist(equations), paste0(prefix, "_hill_equations_with_VitD.txt"))
  
  cat(paste0("\nGenerated Hill-type Equations for ", prefix, ":\n"))
  for (gene in names(equations)) {
    cat(equations[[gene]], "\n")
  }
}

# Function to create Hill-type equations
create_gene_equation <- function(gene_name, regulators_row) {
  activator_terms <- c()
  repressor_terms <- c()
  
  for (regulator_name in names(regulators_row)) {
    effects <- as.character(regulators_row[[regulator_name]])
    if (is.na(effects) || effects == "0") next
    
    # Split comma-separated effects
    effects_list <- unlist(strsplit(effects, ","))
    
    for (effect in effects_list) {
      effect <- trimws(effect)
      if (effect == "") next
      
      effect_num <- suppressWarnings(as.numeric(effect))
      if (is.na(effect_num) || effect_num == 0) next
      
      if (effect_num > 0) {
        activator_terms <- c(activator_terms,
                             paste0("γ_", regulator_name, "_", gene_name,
                                    " * ", regulator_name, "^n / (K^n + ", regulator_name, "^n)"))
      } else if (effect_num < 0) {
        repressor_terms <- c(repressor_terms,
                             paste0("γ_", regulator_name, "_", gene_name,
                                    " / (1 + (", regulator_name, " / K)^n)"))
      }
    }
  }
  
  # Combine activator terms
  activator_sum <- if (length(activator_terms) > 0) {
    paste(c(paste0("α_", gene_name), activator_terms), collapse = " + ")
  } else {
    paste0("α_", gene_name)
  }
  
  # Combine repressor terms
  repressor_product <- if (length(repressor_terms) > 0) {
    paste0("(", paste(repressor_terms, collapse = " * "), ")")
  } else {
    "1"
  }
  
  # Construct final equation
  equation <- paste0("d", gene_name, "/dt = (", activator_sum, ") * ", repressor_product,
                     " - δ_", gene_name, " * ", gene_name)
  
  return(equation)
}

# ============================================================================
# PART 1: Standard kernel reduction
# ============================================================================
cat("\n========================================\n")
cat("PART 1: STANDARD KERNEL REDUCTION\n")
cat("========================================\n\n")

result1 <- kernel_reduction(effects_list, vitD_effects, genes)
cat("\nRemoved genes:", paste(result1$removed_genes, collapse=", "), "\n")
cat("Remaining genes:", paste(result1$remaining_genes, collapse=", "), "\n")
create_output_files(result1, "part1_standard")

# ============================================================================
# PART 2: Forced removal of ESR1
# ============================================================================
cat("\n========================================\n")
cat("PART 2: FORCED REMOVAL OF ESR1\n")
cat("========================================\n\n")

# Re-initialize effects lists for part 2
effects_list2 <- list()
for(i in 1:n_genes) {
  effects_list2[[genes[i]]] <- list()
  for(j in 1:n_genes) {
    effects_list2[[genes[i]]][[genes[j]]] <- c()
    if(adj_matrix[i, j] != 0) {
      effects_list2[[genes[i]]][[genes[j]]] <- c(adj_matrix[i, j])
    }
  }
}

vitD_effects2 <- list()
for(i in 1:n_genes) {
  vitD_effects2[[genes[i]]] <- c()
  if(vitD_vector[i] != 0) {
    vitD_effects2[[genes[i]]] <- c(vitD_vector[i])
  }
}

result2 <- kernel_reduction(effects_list2, vitD_effects2, genes, force_remove = "ESR1")
cat("\nRemoved genes:", paste(result2$removed_genes, collapse=", "), "\n")
cat("Remaining genes:", paste(result2$remaining_genes, collapse=", "), "\n")
create_output_files(result2, "part2_remove_ESR1")







### NETWORK FIGURES
# Required packages
required_packages <- c(
  "igraph", "ggplot2", "dplyr", "readxl",
  "gridExtra", "cowplot", "ggrepel", "tidyr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

if (!requireNamespace("dorothea", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("dorothea")
}
library(dorothea)

################################################ 
# Load TF list
################################################ 
tf_list_data <- read_excel("TF_list.xlsx")
tf_list <- unique(tf_list_data$hgnc_symbol)

################################################ 
# Define direct Vitamin D targets
################################################ 
vitd_targets <- c(
  "ETS2", "FOS", "HIF1A", "ETS1", "SP1", "CDKN1A", "NR3C1", "PLAUR", 
  "JUNB", "CREB1", "ELK1", "FLI1", "MYC", "STAT1", "STAT3", "BCL2L1")

################################################ 
# Gene sets
################################################ 
starting_genes <- c("JUNB","ETS2","FOS","NR4A2","KLF10","CSRNP1","HIF1A")

intermediate_genes <- c(
  "ETS2","FOS","HIF1A","ATF2","CREB1","EGR1","ELK1","ESR1","ETS1","FLI1",
  "FOXM1","JUN","MITF","MYC","NFKB1","SP1","STAT1","STAT3","TP53","TWIST1",
  "ANGPT2","CSF2","FOSL1","IL12B","MMP3","MMP9","TERT","BCL2L1","CCL2","CCND1",
  "CDKN1A","FGF2","IL1B","NOS2","NR3C1","PLA2G4A","PLAUR","NR4A2","KLF10","CSRNP1","JUNB")

stage3_genes <- c(
  "ETS2","FOS","HIF1A","ETS1","SP1","CDKN1A","NR3C1","PLAUR","JUNB","CREB1",
  "ELK1","FLI1","MYC","STAT1","STAT3","BCL2L1","EGR1","ESR1","MITF","CSF2",
  "FOSL1","IL12B","MMP3","MMP9","TERT","CCND1","IL1B","NOS2","ATF2","FOXM1",
  "JUN","NFKB1","TP53","TWIST1","ANGPT2","CCL2","FGF2","PLA2G4A","NR4A2")

final_genes <- c(
  "ETS2","FOS","HIF1A","EGR1","ESR1","ETS1","JUN","MYC","NFKB1",
  "STAT1","STAT3","TP53")

custom_genes_stage6 <- c(
  "ETS2","FOS","HIF1A","EGR1","ETS1","JUN",
  "MYC","NFKB1","STAT1","STAT3","TP53")

################################################ 
# Load adjacency matrices
################################################ 
stage2_adj <- read_excel("processed_adj_matrix.xlsx")
stage4_adj <- read_excel("part1_standard_transposed_matrix_with_vitD.xlsx")
stage5_adj <- read_excel("part2_remove_ESR1_transposed_matrix_with_vitD.xlsx")

################################################ 
# Helper functions
################################################ 
count_edges <- function(value_str) {
  if (is.na(value_str) || value_str == "" || value_str == "0") return(0)
  values <- strsplit(as.character(value_str), ",")[[1]]
  values <- trimws(values)
  values <- as.numeric(values[values != "0" & values != ""])
  values <- values[!is.na(values)]
  return(length(values))
}

extract_edges <- function(value_str) {
  if (is.na(value_str) || value_str == "" || value_str == "0") return(NULL)
  values <- strsplit(as.character(value_str), ",")[[1]]
  values <- trimws(values)
  values <- as.numeric(values)
  values <- values[!is.na(values) & values != 0]
  return(values)
}

calculate_degrees <- function(adj_matrix, gene_col = "Genes", transposed = FALSE) {
  genes <- adj_matrix[[gene_col]]
  mat_values <- adj_matrix[, -1]
  degrees <- data.frame(gene = genes, in_degree = 0, out_degree = 0)
  for (i in 1:length(genes)) {
    out_deg <- 0
    for (j in 1:ncol(mat_values)) out_deg <- out_deg + count_edges(mat_values[i, j])
    in_deg <- 0
    for (k in 1:nrow(mat_values)) in_deg <- in_deg + count_edges(mat_values[k, i])
    if (transposed) {
      degrees$in_degree[i] <- out_deg
      degrees$out_degree[i] <- in_deg
    } else {
      degrees$in_degree[i] <- in_deg
      degrees$out_degree[i] <- out_deg
    }
  }
  return(degrees)
}

assign_layers_from_degrees <- function(gene_list, degrees_df, stage_name = "") {
  if (stage_name == "Stage 1") {
    return(list(input = gene_list, intermediate = character(0), output = character(0)))
  }
  deg_subset <- degrees_df[degrees_df$gene %in% gene_list, ]
  input_nodes <- deg_subset$gene[deg_subset$in_degree == 0]
  force_inputs <- c("KLF10", "CSRNP1")
  force_inputs <- force_inputs[force_inputs %in% gene_list]
  input_nodes <- unique(c(input_nodes, force_inputs))
  output_nodes <- deg_subset$gene[deg_subset$out_degree == 0]
  intermediate_nodes <- setdiff(gene_list, c(input_nodes, output_nodes))
  list(input = input_nodes, intermediate = intermediate_nodes, output = output_nodes)
}

extract_edges_from_adj <- function(adj_matrix, gene_col = "Genes", transposed = FALSE) {
  genes <- adj_matrix[[gene_col]]
  mat_values <- adj_matrix[, -1]
  edges_list <- list()
  edge_idx <- 1
  for (i in 1:length(genes)) {
    from_gene <- genes[i]
    for (j in 1:ncol(mat_values)) {
      to_gene <- names(mat_values)[j]
      edge_values <- extract_edges(mat_values[i, j])
      if (!is.null(edge_values) && length(edge_values) > 0) {
        for (val in edge_values) {
          if (transposed) {
            edges_list[[edge_idx]] <- data.frame(from = to_gene, to = from_gene, mor = val)
          } else {
            edges_list[[edge_idx]] <- data.frame(from = from_gene, to = to_gene, mor = val)
          }
          edge_idx <- edge_idx + 1
        }
      }
    }
  }
  if (length(edges_list) > 0) return(bind_rows(edges_list))
  else return(data.frame(from = character(), to = character(), mor = numeric()))
}

################################################ 
# CREATE SEPARATE LEGEND FIGURE - COMPACT & PROFESSIONAL
################################################ 
create_legend_figure <- function() {
  
  # Uniform font size for all text - small and professional
  uniform_text_size <- 5.5
  
  # Build a clean, compact legend
  p_legend <- ggplot() +
    
    # --- Node type legend (shapes) ---
    # Blue circle = Input
    geom_point(aes(x = 1, y = 4), shape = 21, fill = "#3498DB", color = "black", size = 8, stroke = 1.2) +
    annotate("text", x = 1.5, y = 4, label = "Input Gene", hjust = 0, size = uniform_text_size, color = "black") +
    
    # Yellow circle = Intermediate
    geom_point(aes(x = 1, y = 3.3), shape = 21, fill = "#F1C40F", color = "black", size = 8, stroke = 1.2) +
    annotate("text", x = 1.5, y = 3.3, label = "Intermediate Gene", hjust = 0, size = uniform_text_size, color = "black") +
    
    # Pink circle = Output
    geom_point(aes(x = 1, y = 2.6), shape = 21, fill = "#FF69B4", color = "black", size = 8, stroke = 1.2) +
    annotate("text", x = 1.5, y = 2.6, label = "Output Gene", hjust = 0, size = uniform_text_size, color = "black") +
    
    # Blue triangle = TF
    geom_point(aes(x = 1, y = 1.9), shape = 24, fill = "#3498DB", color = "black", size = 8, stroke = 1.2) +
    annotate("text", x = 1.5, y = 1.9, label = "Transcription Factor (TF) – triangle shape", hjust = 0, size = uniform_text_size, color = "black") +
    
    # Red asterisk = VitD target
    annotate("text", x = 1, y = 1.2, label = "*", size = 12, fontface = "bold", color = "red", vjust = 0.5) +
    annotate("text", x = 1.5, y = 1.2, label = "Direct Vitamin D Target", hjust = 0, size = uniform_text_size, color = "black") +
    
    # --- Edge type legend (arrows) ---
    # Activation
    annotate("segment",
             x = 0.95, xend = 1.35, y = 0.4, yend = 0.4,
             color = "#27AE60", linewidth = 1.5,
             arrow = arrow(length = unit(5, "mm"), type = "closed")) +
    annotate("text", x = 1.5, y = 0.4, label = "Activation Edge", hjust = 0, size = uniform_text_size, color = "black") +
    
    # Repression
    annotate("segment",
             x = 0.95, xend = 1.35, y = -0.3, yend = -0.3,
             color = "#E74C3C", linewidth = 1.5,
             arrow = arrow(length = unit(5, "mm"), type = "closed")) +
    annotate("text", x = 1.5, y = -0.3, label = "Repression Edge", hjust = 0, size = uniform_text_size, color = "black") +
    
    xlim(0.8, 5.5) +
    ylim(-0.7, 4.5) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = "gray60", linewidth = 1),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  return(p_legend)
}

################################################ 
# Main graph function - INCREASED FONTS & SPACING
################################################ 
create_layered_graph <- function(gene_list, stage_name, title_line1, title_line2,
                                 degrees_df, edges_df, curved = FALSE) {
  
  interactions <- edges_df %>%
    filter(from %in% gene_list & to %in% gene_list)
  
  if (nrow(interactions) == 0) {
    cat(paste("Warning: No interactions found for", stage_name, "\n"))
    return(NULL)
  }
  
  layers <- assign_layers_from_degrees(gene_list, degrees_df, stage_name)
  
  cat(paste("\n", stage_name, "- Input:", length(layers$input),
            "| Intermediate:", length(layers$intermediate),
            "| Output:", length(layers$output), "\n"))
  
  nodes <- data.frame(name = gene_list)
  nodes$layer <- ifelse(nodes$name %in% layers$input, "Input",
                        ifelse(nodes$name %in% layers$output, "Output", "Intermediate"))
  nodes$is_tf <- nodes$name %in% tf_list
  nodes$is_vitd_target <- nodes$name %in% vitd_targets
  nodes$original_input <- nodes$name %in% starting_genes
  
  nodes$node_color <- ifelse(nodes$original_input, "#3498DB",
                             ifelse(nodes$layer == "Input", "#3498DB",
                                    ifelse(nodes$layer == "Output", "#FF69B4", "#F1C40F")))
  nodes$label <- nodes$name
  
  n_tfs <- sum(nodes$is_tf)
  n_vitd <- sum(nodes$is_vitd_target)
  
  layer_y <- c("Input" = 100, "Intermediate" = 50, "Output" = 0)
  nodes$y <- layer_y[nodes$layer]
  
  N <- length(gene_list)
  # Increased horizontal spacing significantly for larger labels
  h_spacing <- max(10, 160 / sqrt(N))
  
  for (L in unique(nodes$layer)) {
    subset_names <- nodes$name[nodes$layer == L]
    if (length(subset_names) > 0) {
      nodes$x[nodes$layer == L] <- seq(0, by = h_spacing, length.out = length(subset_names))
    }
  }
  
  # Increased padding significantly for much larger labels
  x_min <- min(nodes$x)
  x_max <- max(nodes$x)
  x_pad <- h_spacing * 2.5   # increased from 2.0
  x_limits <- c(x_min - x_pad, x_max + x_pad)
  
  edges <- interactions %>%
    left_join(nodes, by = c("from" = "name")) %>%
    rename(x_from = x, y_from = y, layer_from = layer, is_tf_from = is_tf) %>%
    left_join(nodes, by = c("to" = "name")) %>%
    rename(x_to = x, y_to = y, layer_to = layer, is_tf_to = is_tf)
  
  edges$edge_color <- ifelse(edges$mor > 0, "#27AE60", "#E74C3C")
  
  edges <- edges %>%
    group_by(from, to) %>%
    mutate(edge_id = row_number(), n_edges = n()) %>%
    ungroup()
  
  has_all_layers <- length(layers$input) > 0 &&
    length(layers$intermediate) > 0 &&
    length(layers$output) > 0
  
  p <- ggplot() +
    # Expanded coord limits significantly to prevent clipping with much larger labels
    coord_cartesian(xlim = x_limits, ylim = c(-25, 125), clip = "off") +
    
    # Separator dashed lines
    {
      if (has_all_layers) {
        list(
          geom_hline(yintercept = 75, linetype = "dashed", color = "gray40", linewidth = 1.2),
          geom_hline(yintercept = 25, linetype = "dashed", color = "gray40", linewidth = 1.2)
        )
      } else list()
    } +
    
    # Edges
    {
      if (curved) {
        list(
          geom_curve(
            data = edges,
            aes(x = x_from, y = y_from, xend = x_to, yend = y_to, color = edge_color),
            curvature = 0.3,
            arrow = arrow(length = unit(5, "mm"), type = "closed"),
            linewidth = 1.0, alpha = 0.4
          )
        )
      } else {
        list(
          geom_segment(
            data = edges %>% mutate(
              x_from_offset = x_from + ifelse(n_edges > 1, (edge_id - 1.5) * 0.5, 0),
              x_to_offset   = x_to   + ifelse(n_edges > 1, (edge_id - 1.5) * 0.5, 0)
            ),
            aes(x = x_from_offset, y = y_from, xend = x_to_offset, yend = y_to, color = edge_color),
            arrow = arrow(length = unit(5, "mm"), type = "closed"),
            linewidth = 1.0, alpha = 0.4
          )
        )
      }
    } +
    
    scale_color_identity() +
    scale_shape_manual(values = c("TF" = 24, "Target" = 21)) +
    scale_fill_identity() +
    
    # Nodes - MUCH LARGER
    geom_point(
      data = nodes,
      aes(x = x, y = y, fill = node_color, shape = ifelse(is_tf, "TF", "Target")),
      size = 9, color = "black", stroke = 2.0
    ) +
    
    # Asterisk for VitD targets - MUCH LARGER
    geom_text(
      data = nodes %>% filter(is_vitd_target),
      aes(x = x, y = y, label = "*"),
      size = 20, fontface = "bold", color = "red", vjust = 0.5
    ) +
    
    # Gene labels - MUCH LARGER font for combined figure visibility
    geom_text_repel(
      data = nodes,
      aes(x = x, y = y, label = label),
      size = 9,                   # increased from 7
      max.overlaps = Inf,
      fontface = "bold",
      color = "black",
      box.padding   = 2.0,        # increased from 1.5
      point.padding = 5.0,        # increased from 4.0
      min.segment.length = 0,
      segment.color = "gray50",
      segment.size  = 0.6
    ) +
    
    # Title - MUCH LARGER fonts for combined figure
    ggtitle(
      label    = title_line1,
      subtitle = paste0(title_line2,
                        " | TFs: ", n_tfs,
                        " | VitD targets: ", n_vitd)
    ) +
    
    theme_minimal(base_size = 18) +
    theme(
      panel.grid   = element_blank(),
      axis.text    = element_blank(),
      axis.title   = element_blank(),
      axis.ticks   = element_blank(),
      plot.title   = element_text(size = 32, face = "bold", color = "black",
                                  hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 28, face = "bold", color = "black",
                                   hjust = 0.5, margin = margin(b = 15)),
      legend.position = "none",
      plot.margin = margin(20, 20, 20, 20)
    )
  
  return(p)
}

################################################ 
# Output folders
################################################ 
dir.create("layered_networks_vitd", showWarnings = FALSE)
dir.create("layered_networks_vitd/straight", showWarnings = FALSE)
dir.create("layered_networks_vitd/curved", showWarnings = FALSE)
dir.create("layered_networks_vitd/combined", showWarnings = FALSE)

################################################ 
# Calculate degrees
################################################ 
cat("\n===== Calculating degrees from adjacency matrices =====\n")

degrees_stage2 <- calculate_degrees(stage2_adj, transposed = FALSE)
edges_stage2   <- extract_edges_from_adj(stage2_adj, transposed = FALSE)

degrees_stage3 <- degrees_stage2; edges_stage3 <- edges_stage2
degrees_stage4 <- degrees_stage2; edges_stage4 <- edges_stage2

degrees_stage5 <- calculate_degrees(stage4_adj, transposed = TRUE)
edges_stage5   <- extract_edges_from_adj(stage4_adj, transposed = TRUE)

degrees_stage6 <- calculate_degrees(stage5_adj, transposed = TRUE)
edges_stage6   <- extract_edges_from_adj(stage5_adj, transposed = TRUE)

degrees_stage1 <- degrees_stage2; edges_stage1 <- edges_stage2

################################################ 
# Stage definitions
################################################ 
stages <- list(
  list(genes = starting_genes,   name = "Stage 1",
       title_line1 = "Stage 1 - 7 Starting Genes",
       title_line2 = "Nodes: 7 | Edges: XXX | Triangles = TFs, Circles = Targets",
       degrees = degrees_stage1, edges = edges_stage1),
  
  list(genes = intermediate_genes, name = "Stage 2",
       title_line1 = "Stage 2 - 41 Genes",
       title_line2 = "Nodes: 41 | Edges: XXX | Triangles = TFs, Circles = Targets",
       degrees = degrees_stage2, edges = edges_stage2),
  
  list(genes = stage3_genes, name = "Stage 3",
       title_line1 = "Stage 3 - 39 Genes",
       title_line2 = "Nodes: 39 | Edges: XXX | Triangles = TFs, Circles = Targets",
       degrees = degrees_stage3, edges = edges_stage3),
  
  list(genes = stage3_genes[stage3_genes %in% tf_list], name = "Stage 4",
       title_line1 = "Stage 4 - TFs Only",
       title_line2 = "Nodes: XXX | Edges: XXX | Triangles = TFs, Circles = Targets",
       degrees = degrees_stage4, edges = edges_stage4),
  
  list(genes = final_genes, name = "Stage 5",
       title_line1 = "Stage 5 - Final Genes",
       title_line2 = "Nodes: 13 | Edges: XXX | Triangles = TFs, Circles = Targets",
       degrees = degrees_stage5, edges = edges_stage5),
  
  list(genes = custom_genes_stage6, name = "Stage 6",
       title_line1 = "Stage 5 - 11 Genes",
       title_line2 = "Nodes: 12 | Edges: XXX | Triangles = TFs, Circles = Targets",
       degrees = degrees_stage6, edges = edges_stage6))

all_straight <- list()
all_curved   <- list()

################################################ 
# Generate plots - INCREASED dimensions
################################################ 
for (i in seq_along(stages)) {
  st <- stages[[i]]
  cat(paste("\n===== Processing", st$name, "=====\n"))
  
  edges_for_stage <- st$edges %>% filter(from %in% st$genes & to %in% st$genes)
  n_edges <- nrow(edges_for_stage)
  n_nodes <- length(st$genes)
  
  title_line2_updated <- gsub("Nodes: XXX", paste0("Nodes: ", n_nodes), st$title_line2)
  title_line2_updated <- gsub("Edges: XXX",  paste0("Edges: ",  n_edges), title_line2_updated)
  
  p1 <- create_layered_graph(st$genes, st$name, st$title_line1, title_line2_updated,
                             st$degrees, st$edges, curved = FALSE)
  p2 <- create_layered_graph(st$genes, st$name, st$title_line1, title_line2_updated,
                             st$degrees, st$edges, curved = TRUE)
  
  # Increased plot dimensions significantly to accommodate much larger text
  plot_width  <- max(18, min(36, length(st$genes) / 1.5))
  plot_height <- 14
  
  if (!is.null(p1)) {
    ggsave(paste0("layered_networks_vitd/straight/", st$name, "_straight.png"),
           p1, width = plot_width, height = plot_height, dpi = 300, bg = "white")
    ggsave(paste0("layered_networks_vitd/straight/", st$name, "_straight.svg"),
           p1, width = plot_width, height = plot_height, bg = "white")
    all_straight[[i]] <- p1
  }
  
  if (!is.null(p2)) {
    ggsave(paste0("layered_networks_vitd/curved/", st$name, "_curved.png"),
           p2, width = plot_width, height = plot_height, dpi = 300, bg = "white")
    ggsave(paste0("layered_networks_vitd/curved/", st$name, "_curved.svg"),
           p2, width = plot_width, height = plot_height, bg = "white")
    all_curved[[i]] <- p2
  }
}

################################################ 
# SAVE SEPARATE LEGEND FIGURE
################################################ 
cat("\n===== Creating Separate Legend Figure =====\n")

p_legend <- create_legend_figure()

ggsave("layered_networks_vitd/combined/legend.png",
       p_legend, width = 6, height = 4, dpi = 300, bg = "white")
ggsave("layered_networks_vitd/combined/legend.svg",
       p_legend, width = 6, height = 4, bg = "white")

cat("Legend saved: layered_networks_vitd/combined/legend.png\n")

################################################ 
# COMBINED PLOTS: 2x2 top + 1 centered bottom
################################################ 
cat("\n===== Creating Combined Figures =====\n")

make_combined <- function(plot_list_4, plot_last,
                          labels_4 = c("A","B","C","D"), label_last = "E") {
  top_grid <- plot_grid(
    plotlist = plot_list_4, ncol = 2,
    labels = labels_4, label_size = 40
  )
  empty_plot  <- ggplot() + theme_void()
  bottom_grid <- plot_grid(
    empty_plot, plot_last, empty_plot,
    ncol = 3, rel_widths = c(0.5, 1, 0.5),
    labels = c("", label_last, ""), label_size = 40
  )
  plot_grid(top_grid, bottom_grid, ncol = 1, rel_heights = c(2, 1))
}

combined_s <- make_combined(
  list(all_straight[[1]], all_straight[[2]], all_straight[[3]], all_straight[[4]]),
  all_straight[[6]]
)
combined_c <- make_combined(
  list(all_curved[[1]], all_curved[[2]], all_curved[[3]], all_curved[[4]]),
  all_curved[[6]]
)

# Increased combined figure dimensions for publication quality
ggsave("layered_networks_vitd/combined/all_stages_straight.png",
       combined_s, width = 36, height = 32, dpi = 300, bg = "white")
ggsave("layered_networks_vitd/combined/all_stages_straight.svg",
       combined_s, width = 36, height = 32, bg = "white")

ggsave("layered_networks_vitd/combined/all_stages_curved.png",
       combined_c, width = 36, height = 32, dpi = 300, bg = "white")
ggsave("layered_networks_vitd/combined/all_stages_curved.svg",
       combined_c, width = 36, height = 32, bg = "white")



