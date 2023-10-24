#!/usr/bin/env Rscript

# Load required libraries
library("gplots")       # Provides plotting utilities
library("here")         # Helps with file paths
library("RColorBrewer") # Color palettes
library("flipr")        # A flexible permutation framework for making inferences 
library("reshape2")     # Reshaping data
library("ggplot2")      # Primary plotting library
library("colorspace")   # Color utilities
library("pals")         # More color palettes
library("cowplot")      # ggplot themes and utilities
library("scales")       # Scale functions for graphics
library("sommer")       # Structural multivariate-univariate linear mixed model solver
library("dplyr")        # Data manipulation
library("corrplot")     # Visualize a correlation matrix
library("data.table")   # Enhanced data frame
library("glue")         # Interpreted string literals



# Set options for the R session
options(warn = -1)
gc()

# Clear the R environment
rm(list = ls(all.names = TRUE))

# Define the current working directory
cwd <- here() 

# Check for input file in command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args) == 1) {
  args[2] <- "out.txt"
}

# Set the encoding to UTF-8
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

# Read data
df <- fread(paste0("averages_and_counts_Mode", args[3], ".csv"))


# Rescale columns
rescale_column <- function(df, column_name, new_name, to_range) {
  local_min <- min(df[[column_name]])
  local_max <- max(df[[column_name]])
  rescaled_data <- scales::rescale(df[[column_name]], to = to_range, from = c(local_min, local_max))
  df[[new_name]] <- rescaled_data
  return(df)
}




df <- df %>% 
  rescale_column("Average_Cosine_Similarity", "Scaled_Avg_Cosine_Similarity", c(-1, 1)) %>%
  rescale_column("Average_Sequence_Similarity", "Scaled_Avg_Sequence_Similarity", c(-1, 1)) %>%
  rescale_column("Average_Frobenius_Similarity", "Scaled_Avg_frobenius_similarity", c(0, 1)) %>%
  rescale_column("Average_Weighted_Sum", "Scaled_weighted_sum", c(0, 1))


# Set up common ggplot theme for the plots
custom_theme <- function() {
  theme_bw() +
  theme(
    panel.grid.major.x = element_line(color = "gray70", linetype = "dashed", linewidth = 0.5),
    panel.grid.minor.x = element_line(color = "#d3d3d3"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 2),
    axis.ticks = element_line(linewidth = 1.5), 
    axis.text = element_text(size = 21),
    axis.title.x = element_text(size = 21, margin = margin(t = 12)),
    axis.title.y = element_text(size = 21, margin = margin(r = 12)),
    plot.title = element_text(size = 25)
  )
}

# Function to create a ggplot
create_plot <- function(df, fill_var, title) {
  
  # Define the base plot
  p <- ggplot(df, aes(x=Tile_Center_Index, y=Tile_Length)) + 
    coord_fixed() + 
    geom_tile(aes(fill = .data[[fill_var]]), width = 1.00, height = 1.00)

  # Conditionally set the color scale
  if (fill_var == "Scaled_Avg_frobenius_similarity") {
    p <- p + scale_fill_gradientn(
      colors = brewer.spectral(500), # Directly use the palette without reversing
      oob = scales::oob_squish
    )
  } else {
    p <- p + scale_fill_gradientn(
      colors = rev(brewer.spectral(500)), # Reverse the palette
      oob = scales::oob_squish
    )
  }
  
  # Add the remaining plot components
  p <- p + 
    scale_y_reverse(breaks = scales::pretty_breaks(n = 6), sec.axis = sec_axis(~., breaks = NULL)) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 20), sec.axis = sec_axis(~., breaks = NULL)) + 
    labs(title = title, x = "Tile Center Index", y = "Tile Length") +
    guides(
      fill = guide_colorbar(
        barwidth = 2.0,
        barheight = 35,
        label = TRUE,
        label.position = "right",
        label.theme = element_text(size = 20), 
        ticks = TRUE,
        ticks.linewidth = 2.0,
        ticks.colour = "black",
        title = NULL,
        frame.colour = "black",
        frame.linewidth = 1.3
      )
    ) +
    custom_theme()
  
  return(p)
}



# Define plot dimensions
width <- 30
height <- 12 

# Create the plots
all_avgs_and_counts <- create_plot(df, "Scaled_Avg_Cosine_Similarity", paste0("Average Cosine Similarity for ", toupper(args[1]), " ", args[2], " - Top 6 ANM Modes"))
seqSim_and_counts <- create_plot(df, "Scaled_Avg_Sequence_Similarity", paste0("Average Sequence Similarity for ", toupper(args[1]), " ", args[2]))
FrobSim_and_counts <- create_plot(df, "Scaled_Avg_frobenius_similarity", paste0("Average Frobenius Similarity for ", toupper(args[1]), " ", args[2], " - Top 6 ANM Modes"))

# Save the plots
ggsave(paste0("AvgCounts_Mode", args[3], "Cos-SeqSim-Frob.pdf"), all_avgs_and_counts, width = width, height = height)
ggsave(paste0("AvgCounts_Mode", args[3], "SeqSim.pdf"), seqSim_and_counts, width = width, height = height)
ggsave(paste0("AvgCounts_Mode", args[3], "FrobSim.pdf"), FrobSim_and_counts, width = width, height = height)
