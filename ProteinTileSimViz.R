#!/usr/bin/env Rscript

# Library imports
libraries <- c("gplots", "here", "RColorBrewer", "flipr", "reshape2", 
               "ggplot2", "colorspace", "pals", "cowplot", "scales", 
               "sommer", "dplyr", "corrplot", "data.table", "ggExtra")
invisible(lapply(libraries, library, character.only = TRUE))

# Set options for the R session and clear the environment
options(warn = -1)
gc()
rm(list = ls(all.names = TRUE))

# Define the current working directory and check for input file in command line arguments
cwd <- here() 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args) == 1) {
  args[2] <- "out.txt"
}

# Set the encoding to UTF-8 and read in the data
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")
df <- fread(paste0("cosine_similarities_Mode", args[3], "_updated_ranked_LocalAln.csv"))

# Define global components that are reused
breaks <- seq(min(df$tile_length) - 0.01, max(df$tile_length), length.out = 5)
colors <- brewer.spectral(length(breaks) - 1)
df$colors <- cut(df$tile_length, breaks = breaks, labels = colors)
df$colors <- as.factor(df$colors)
width <- 18  # Width in inches
height <- 12  # Height in inches

# Default ggplot theme
custom_theme <- theme_bw() +
  theme(panel.grid.major.x = element_line(color = "gray70", linetype = "dashed", linewidth = 0.5),
        panel.grid.minor.x = element_line(color = "#d3d3d3"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 2),
        axis.ticks = element_line(linewidth = 1.5), 
        axis.text = element_text(size = 21),
        axis.title.x = element_text(size = 21, margin = margin(t = 12)),
        plot.margin = margin(1, 1, 2, 1, "cm"),
        plot.title = element_text(size = 25))




# Modularized plotting functions

plotSequenceVsTileLength <- function(df, args, breaks, colors, width, height, custom_theme) {
  filename <- paste0("seqSim_against_TileLen_Mode", args[3], ".pdf")
  pdf(filename, width = width, height = height)
  
  plot <- ggplot(df, aes(x=tile_length, y=sequence_similarity)) + 
    stat_density_2d(geom = "polygon", contour = TRUE,
                    aes(fill = after_stat(level)),
                    bins = 7) +
    scale_fill_distiller(palette = "Spectral", direction = -1) + 
    scale_y_reverse(expand = c(0, 0), 
                    breaks = scales::pretty_breaks(n = 10), 
                    sec.axis = sec_axis(~., breaks = NULL)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = scales::pretty_breaks(n = 20), 
                       sec.axis = sec_axis(~., breaks = NULL)) + 
    labs(title = paste0("Sequence Similarity vs. Tile Length: 2FDN Against ", toupper(args[2]), " - Top 6 ANM Modes"),
         x = "Tile Length",
         y = "Sequence Similarity")  +
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
          frame.linewidth = 1.3,
          )
      ) +
    custom_theme
  
  print(plot)
  dev.off()
}

plotBarSequenceVsTileLength <- function(df, args, breaks, colors, width, height, custom_theme) {
  filename <- paste0("seqSim_against_TileLen_BarPlt_Mode", args[3], ".pdf")
  pdf(filename, width = width, height = height)
  
  plot <- ggplot(df, aes(x=tile_length, y=sequence_similarity)) +
    geom_bar(stat='identity', fill='blue') +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 10), 
      sec.axis = sec_axis(~., breaks = NULL)) + 
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 20), 
      sec.axis = sec_axis(~., breaks = NULL)) + 
    labs(title = paste0("Sequence Similarity vs. Tile Length: 2FDN Against ", toupper(args[2]), " - Top 6 ANM Modes"),
         x = "Tile Length",
         y = "Sequence Similarity")  +
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
          frame.linewidth = 1.3,
          )
      ) +
    custom_theme
  
  print(plot)
  dev.off()
}

plotSequenceVsCenterIndex <- function(df, args, breaks, colors, width, height, custom_theme) {
  filename <- paste0("seqSim_against_CentIdx1_Tst_Mode", args[3], ".tiff")
  tiff(filename, width = width, height = height, units = "in", res = 300)
  
  plot <- ggplot(df, aes(x=center_index_1, y=sequence_similarity)) +
    geom_point(aes(fill = colors), shape = 21, colour = "white", stroke = 0.5, alpha = 6/10, size = 3.00) +
    scale_fill_manual(values = colors, labels = c("6-18", "18-30", "30-42", "42-54")) +
    scale_y_reverse(expand = c(0, 0), 
                    breaks = scales::pretty_breaks(n = 10), 
                    sec.axis = sec_axis(~., breaks = NULL)) + 
    scale_x_continuous(expand = c(0, 0), 
                       breaks = scales::pretty_breaks(n = 20), 
                       sec.axis = sec_axis(~., breaks = NULL)) + 
    labs(title = paste0("Sequence Similarity vs. Tile Length: 2FDN Against ", toupper(args[2]), " - Top 6 ANM Modes"),
         x = "Tile Center Index (from 5T5I)",
         y = "Sequence Similarity",
         fill = "Tile Lengths")  +
    guides(
        fill = guide_legend(
          title = "Tile Lengths",
          label.position = "right",
          label.theme = element_text(size = 15), 
          title.theme = element_text(size = 18),
          keywidth = 1.0,
          keyheight = 1.0,
          default.unit = "cm",
          override.aes = list(size = 8)
          )
      ) +
    custom_theme
  
  print(plot)
  dev.off()
}


# Invoking the functions
plotSequenceVsTileLength(df, args, breaks, colors, width, height, custom_theme)
plotBarSequenceVsTileLength(df, args, breaks, colors, width, height, custom_theme)
plotSequenceVsCenterIndex(df, args, breaks, colors, width, height, custom_theme)