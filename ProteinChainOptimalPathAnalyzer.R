#!/usr/bin/env Rscript

# Libraries
library(ggplot2)
library(dplyr)
library(here)
library(RColorBrewer)
library(pals)
library(glue)

# --------------------------
# Utility Functions
# --------------------------

# Extract value by flag from command line arguments
get_value_by_flag <- function(flag) {
  all_args <- commandArgs()
  arg_with_equal <- all_args[grep(paste0(flag, "="), all_args)]
  if (length(arg_with_equal) > 0) {
    return(sub(paste0(flag, "="), "", arg_with_equal))
  }
  short_flag_index <- grep(paste0("-", substr(flag, 3, 3)), all_args)
  if (length(short_flag_index) > 0) {
    return(all_args[short_flag_index + 1])
  }
  return(NULL)
}

# Read and preprocess data from files
read_and_preprocess_data <- function() {
  path_sequences <- read.csv(paste0(here(), "/path_sequences.csv"))
  path_sequences$Label <- as.factor(path_sequences$Label)
  path_sequences$Best_Path <- as.logical(path_sequences$Best_Path)

  all_labeled_points <- read.csv(paste0(here(), "/all_labeled_points.csv"))
  all_labeled_points$Label <- as.factor(all_labeled_points$Label)

  # Create a dataframe with starting and ending coordinates for each segment
  segments_df <- path_sequences %>%
    group_by(Path_ID) %>%
    mutate(
      xend = lead(Tile_Center_Index),
      yend = lead(Tile_Length)
    ) %>%
    ungroup() %>%
    filter(!is.na(xend))

  # Automatically create a color palette
  labels <- unique(path_sequences$Label)
  colors <- pals::turbo(length(labels)) # Generate a color palette using pals palettes

  # Map the colors to the labels A to the last letter
  label_colors <- setNames(colors, unique(path_sequences$Label))

  list(path_sequences = path_sequences, all_labeled_points = all_labeled_points,
       segments_df = segments_df, label_colors = label_colors)
}

# Identify and print boundaries based on differences in Tile_Center_Index
identify_and_plot_indices <- function(df) {
  boundaries <- c(min(df$Tile_Center_Index) - 1) # Start with min value

  for(tile_length in unique(df$Tile_Length)) {
    group <- df[df$Tile_Length == tile_length,]
    tile_center_indices <- group$Tile_Center_Index

    for(i in 2:(length(tile_center_indices) - 1)) { # Start from index 2 to check previous difference
      current_difference <- tile_center_indices[i + 1] - tile_center_indices[i]
      previous_difference <- tile_center_indices[i] - tile_center_indices[i - 1]

      # Check if the current difference is significantly greater than the previous difference
      if (current_difference >= 2 * previous_difference && tile_length != 6) {
        midpoint <- (tile_center_indices[i] + tile_center_indices[i + 1]) / 2

        # Check if the midpoint is not already in the boundaries
        if (!(midpoint %in% boundaries)) {
          boundaries <- c(boundaries, midpoint) # Add the midpoint to boundaries
        }
      }
    }
  }

  # Add the maximum Tile_Center_Index value as the last boundary
  boundaries <- c(boundaries, max(df$Tile_Center_Index) + 1)

  # Print the accumulated boundaries
  for (i in 1:(length(boundaries) - 1)) {
    cat(paste("Boundary between region", i, "and region", i + 1, "is between Tile Center Indices",
              boundaries[i], "and", boundaries[i + 1], "\n"))
  }

  return(boundaries)
}


report_tile_indices <- function(file_path, sequence_file_path) {
  
  # Read the CSV file
  path_sequences <- read.csv(file_path)

  # Convert Best_Path to logical
  path_sequences$Best_Path <- as.logical(path_sequences$Best_Path)

  # Extract the maximum tile length
  max_tile_length <- max(path_sequences$Tile_Length)

  # Define points of interest
  half_max_tile_length <- 30
  
  # Tolerance for floating-point comparison
  tolerance <- 1e-5

  # Retrieve sequence based on tile length and center index
  get_sequence <- function(tile_length, tile_center_index) {
    tile_center_index_formatted <- formatC(tile_center_index, format = "f", digits = 1)
    awk_command <- sprintf("gawk -F ',' '{ if ($5 == \"%s\" && $9 == \"%s\") print $1}' %s", 
                           tile_center_index_formatted, tile_length, sequence_file_path)
    system(awk_command, intern = TRUE)
  }

  # Identify overlaps between sequences
  find_overlap <- function(seq1, seq2) {
    for (i in 1:nchar(seq1)) {
      if (startsWith(seq2, substring(seq1, i))) {
        return(substring(seq1, i))
      }
    }
    return(NULL)
  }

  # Mix two colors together
  mix_colors <- function(color1, color2) {
    mixed_rgb <- (col2rgb(color1) + col2rgb(color2)) / 2
    return(rgb(mixed_rgb[1], mixed_rgb[2], mixed_rgb[3], maxColorValue = 255))
  }

  # Print sequences and handle overlaps
  print_sequences <- function(df, tile_length, colors) {
    sequences <- sapply(df$Tile_Center_Index, function(index) unique(get_sequence(tile_length, index)))
    seen_overlaps <- character(0)
    overlaps <- character(0)
    overlap_colors <- character(0)
    overlap_threshold <- 2  # Minimum length for overlaps
    
    for (i in seq_len(nrow(df))) {
      for (j in seq_len(nrow(df))) {
        if (i != j) {
          overlap <- find_overlap(sequences[i], sequences[j])
          if (!is.null(overlap) && !(overlap %in% seen_overlaps) && nchar(overlap) >= overlap_threshold) {
            seen_overlaps <- c(seen_overlaps, overlap)
            overlaps <- c(overlaps, overlap)
            overlap_colors <- c(overlap_colors, mix_colors(colors[i], colors[j]))
          }
        }
      }
    }
    
    # Print sequences and overlaps
    for (i in order(nchar(sequences), decreasing = TRUE)) {
      cat("color", gsub("#", "0x", colors[i]), ", ( ps. ", sequences[i], ")\n")
    }
    for (i in seq_along(overlaps)) {
      cat("color", gsub("#", "0x", overlap_colors[i]), ", ( ps. ", overlaps[i], ")\n")
    }
  }

  # Retrieve color palette
  get_colors <- function(n) {
    # Sample color palette, can be replaced or extended
    colors <- c("#E88180", "#D9C47C", "#BE6EDB", "#58DB80", "#6571E0", 
                "#F4A460", "#8A2BE2", "#20B2AA", "#FF4500", "#FFD700")
    
    if (length(colors) < n) {
      colors <- colorRampPalette(rainbow(n))(n)
    }
    return(colors)
  }

  # Filter and print sequences based on Tile_Length and Best_Path
  half_max_tile_info <- path_sequences[abs(path_sequences$Tile_Length - half_max_tile_length) < tolerance & 
                                       path_sequences$Best_Path,]
  max_tile_info <- path_sequences[abs(path_sequences$Tile_Length - max_tile_length) < tolerance & 
                                  path_sequences$Best_Path,]

  cat("Half max tile info:\n")
  print(half_max_tile_info)
  print_sequences(half_max_tile_info, half_max_tile_length, get_colors(nrow(half_max_tile_info) + 1))

  cat("\nMax tile info:\n")
  print(max_tile_info)
  print_sequences(max_tile_info, max_tile_length, get_colors(nrow(max_tile_info)))
}


for_nine_islands <- c("#E88180", "#D9C47C", "#BE6EDB", "#58DB80", "#6571E0",
                      "#F4A460", "#8A2BE2", "#20B2AA", "#FF4500", "#FFD700")


# Plot the data with ggplot2
plot_data <- function(path_sequences, all_labeled_points, segments_df, label_colors, boundaries) {
  plot <- ggplot(all_labeled_points, aes(x = Tile_Center_Index, y = Tile_Length, color = Label)) +
    geom_point(size = 3.5, shape = 15, stroke = 1.0) +
    geom_segment(data = segments_df, aes(xend = xend, yend = yend),
                 linewidth = 0.95, alpha = 0.25, color = "red") +
    geom_point(data = path_sequences, aes(x = Tile_Center_Index, y = Tile_Length),
               size = 2.75, shape = 21, alpha = 0.25, fill = "red", color = "white") +
    geom_vline(xintercept = boundaries, color = 'black') +
    scale_color_manual(values = for_nine_islands, name = "Tile Island IDs") +
    labs(title = glue("Optimal Path Selection for {toupper(pdb_id)} Chain {chain_id} Using Dynamic Programming"),
         x = 'Tile Center Index', y = 'Tile Length') +
    guides(color = guide_legend(override.aes = list(size = 14))) + 
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 2),
          axis.ticks = element_line(linewidth = 1.5),
          axis.text = element_text(size = 25),
          axis.title.x = element_text(size = 27, margin = margin(t = 12)),
          axis.title.y = element_text(size = 27, margin = margin(r = 12)),
          plot.title = element_text(size = 28),
          legend.title = element_text(size = rel(2.0)),
          legend.text = element_text(size = 24),
          legend.key.size = unit(1.5, "lines")) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_y_reverse()

  return(plot)
}




# --------------------------
# Main Execution
# --------------------------

# Get command line arguments
pdb_id <- get_value_by_flag("--pdb-id")
chain_id <- get_value_by_flag("--chain-id")

# Read and preprocess data
data <- read_and_preprocess_data()

# Extract relevant data from the list
path_sequences <- data$path_sequences
all_labeled_points <- data$all_labeled_points
segments_df <- data$segments_df
label_colors <- data$label_colors

# Identify boundaries
boundaries <- identify_and_plot_indices(data$all_labeled_points)

# Plot data and save the plot
plot <- plot_data(path_sequences, all_labeled_points, segments_df, label_colors, boundaries)
ggsave(glue("dynamicProgrammingResults_{pdb_id}_chain{chain_id}.pdf"), plot, width = 20, height = 10, dpi = 320, units = "in")

# Report tile indices
report_tile_indices("path_sequences.csv", "final_scan_ModeAll_withFalse.csv")
