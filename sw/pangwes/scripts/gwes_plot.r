#!/usr/bin/env Rscript

# Help output
print_usage <- function() {
  cat("Usage: ./gwes_plot.r -i input_links_file -n number_of_assemblies [other options]\n")
  cat("Options:\n")
  cat("  -i, --input            Input links file (required)\n")
  cat("  -n  --num-assemblies   Number of assemblies in the analysis (required)\n")
  cat("  -o, --output           Output plot path (default: ./pangwes_plot.png)\n")
  cat("  -l  --ld-dist          LD distance cutoff (default: 50000)\n")
  cat("  --no-deps              Generate a basic plot without any R package dependencies (Much slower!)\n")
  cat("  -h  --help             Display this help message\n")
}
# Function to parse command-line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Initialize a list to store options
  options <- list(input = NULL, output = "pangwes_plot.png", number = NULL, basic_plot_only = FALSE, ld_dist = 50000)
  
  # Loop through the arguments and assign values to options
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "-i" || args[i] == "--input") {
      options$input <- args[i + 1]
      i <- i + 1
    } else if (args[i] == "-o" || args[i] == "--output") {
      options$output <- args[i + 1]
      i <- i + 1
    } else if (args[i] == "-n" || args[i] == "--number") {
      options$number <- as.integer(args[i + 1])
      if(is.na(options$number)) {
        stop(paste("Number of assemblies must be an integer"), call. = FALSE)
      }
      i <- i + 1
    } else if (args[i] == "--no-deps") {
      options$basic_plot_only <- TRUE
    } else if (args[i] == "-l" | args[i] == "--ld-dist") {
      options$ld_dist <- as.integer(args[i + 1])
    } else if (args[i] == "-h" | args[i] == "--help") {
      print_usage()
      quit(status = 0)
    } else {
      stop(paste("Unknown input option:", args[i],". Use --help for options."), call. = FALSE)
    }
    i <- i + 1
  }
  
  # Check if input and output are provided
  if (is.null(options$input) || is.null(options$number)) {
    # print_usage()
    stop("Input links file and the number of assemblies must be provided, use --help for options.", call. = FALSE)
  }
  
  return(options)
}
# parse the 9 column input (hopefully there are not other versions)
parse_input <- function(input) {
  cat("Parsing input...\n")
  # The input must have these 9 columns
  # # Column 1:     Unitig 1 (#1)
  # # Column 2:     Unitig 2 (#2)
  # # Column 3:     Average graph distance (#3)
  # # Column 4:     SpydrPick ARACNE flag (#4)
  # # Column 5:     SpydrPick Mutual Information score (#5)
  # # Column 6:     Graph distance count (#6)
  # # Column 7:     Graph distances' M2 value (sum of squares of differences from the current mean) (#19)
  # # Column 8:     Graph distances' minimum value (#20)
  # # Column 9:     Graph distances' maximum value (#21)
  if(ncol(input) != 9) stop("Input links file must have 9 columns", call. = FALSE)
  if(nrow(input) < 100) stop("Input has less than 100 rows, cannot generate plot!", call. = FALSE)
  colnames(input) = c("U1", "U2", "sep", "ARACNE", "MI", "count", "d_M2", "d_min", "d_max")
  
  # Check for connected links
  n_K_disconnected <- length(which(input$sep == -1))
  if(length(n_K_disconnected) == 0) n_K_disconnected = 0
  
  cat("Filtering connected links...")
  connected <- input[input$sep != -1, ]
  n_K_connected <- 0
  if (nrow(connected) > 0) { n_K_connected <- floor(nrow(connected) / 1000) } else stop("Cannot generate plot, no connected links found in the input links file!")
  connected$rsd = (sqrt(connected$d_M2 / (num_assemblies - 1)))/connected$sep
  
  # Filtering
  idx_to_keep = which( (connected$count >= ceiling(0.05 * num_assemblies) & connected$rsd <= 1) )
  if(length(idx_to_keep) < 100) stop("Cannot generate plot, less than 100 links remain after filtering!")
  
  filtered_data <- connected[ idx_to_keep , ]
  # Compute outlier thresholds
  q <- quantile(filtered_data$MI, c(0.25, 0.75))
  
  cat("Computing outlier thresholds...")
  outlier_thresh <- q[[2]] + 1.5 * (q[[2]]-q[[1]])
  extreme_outlier_thresh  <- q[[2]] + 3 * (q[[2]]-q[[1]])
  
  return(list(filtered_data = filtered_data, n_K_disconnected = n_K_disconnected, n_K_connected = n_K_connected,
              outlier_thresh = outlier_thresh, extreme_outlier_thresh = extreme_outlier_thresh))
}
# Check and install packages
check_and_install_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      cat(paste("Package", pkg, "is not installed.\n"))
      answer <- readline(prompt = "Do you want to install it? (yes/no): ")
      if (tolower(answer) %in% c("yes", "y")) {
        install.packages(pkg, repos = "http://cran.us.r-project.org")
        library(pkg, character.only = TRUE)
      } else {
        stop(paste("Package", pkg, "is required. Exiting script."))
      }
    }
  }
}
# Parse the command-line options
opts <- parse_args()

# Script 
# Obvious sanity checks
input_path <- opts$input
if(!file.exists(input_path)) {
  stop("Input file does not exist", call. = FALSE)
}

output_path <- opts$output
if(file.exists(output_path)) {
  stop("Plot already exists, please delete before running again", call. = FALSE)
}

# Other inputs
num_assemblies <- opts$number
basic_plot_only <- opts$basic_plot_only

# Echo primary parameters
cat("Input file:", input_path, "\n")
cat("Output file:", output_path, "\n")
cat("Number of assemblies:", num_assemblies, "\n")


## Let's install the packages 
if(!basic_plot_only){
  # List of required packages
  required_packages <- c("data.table", "ggplot2", "ggrastr", "ggthemes")
  # Check and install required packages
  check_and_install_packages(required_packages)
}

## Read the input
cat("Reading input file...")
time_reading_start <- proc.time()
if(basic_plot_only){
  input <- read.csv(input_path, header = FALSE, sep = " ") # May take a few minutes.
} else {
  library(data.table)
  library(ggplot2)
  library(ggrastr)
  library(ggthemes)
  input <- data.table::fread(input_path, header = FALSE, sep = " ") # May take a few minutes.
}
time_reading_end <- proc.time()
time_reading <- (time_reading_end - time_reading_start)[[3]]
cat(paste0("Done in ", time_reading, "s\n"))

input_list <- parse_input(input) # perform some sanity checks

connected = input_list$filtered_data # This is now the filtered data ready for plotting
outlier_threshold = input_list$outlier_thresh
extreme_outlier_threshold = input_list$extreme_outlier_thresh

# Fixed plot parameters
plot_width <- 2400
plot_height <- 1200
color_ld <- "black"
color_outlier <- rgb(165, 0, 38, maxColorValue = 255)
color_extreme_outlier <- rgb(215, 48, 39, maxColorValue = 255)
color_direct <- rgb(0, 115, 190, maxColorValue = 255)
color_indirect <- rgb(192, 192, 192, maxColorValue = 255)

## Generate plot
cat("Generating plot... \n")
if(basic_plot_only){
  # Fixed parameters
  plot_width <- 1920
  plot_height <- 1080
  plot_pointsize <- 16
  
 
  plot_symbol <- 19 # 19 - solid circle.
  cex_direct <- 0.2
  cex_indirect <- 0.1
  cex_legend <- 1.2
  
  # Estimated parameters
  min_mi <- min(connected[, 5], na.rm = TRUE)
  max_mi <- max(connected[, 5], na.rm = TRUE)
  max_distance <- max(connected[, 3], na.rm = TRUE)
  exponent <- round(log10(max_distance)) - 1
  
  time_plotting_start <- proc.time()
  png(output_path, width = plot_width, height = plot_height, pointsize = plot_pointsize)
  plot(connected[!connected[, 4], 3], connected[!connected[, 4], 5], col = color_indirect, type = "p", pch = plot_symbol, cex = cex_indirect,
       xlim = c(0, max_distance), ylim = c(min_mi, max_mi), xaxs = "i", yaxs = "i",
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  lines(connected[as.logical(connected[, 4]), 3], connected[as.logical(connected[, 4]), 5], col = color_direct, type = "p", pch = plot_symbol, cex = cex_direct)
  axis(1, at = seq(0, max_distance, 10^exponent), tick = FALSE, labels = seq(0, max_distance / 10^exponent), line = -0.8)
  axis(2, at = seq(0.05, 1, 0.05), labels = FALSE, tcl = -0.5)
  axis(2, at = seq(0.1, 1, 0.1), labels = seq(0.1, 1, 0.1), las = 1, tcl = -0.5)
  title(xlab = "Distance between unitigs (bp)", line = 1.2)
  title(xlab = substitute(x10^exp, list(exp = exponent)), line = 1.4, adj = 1)
  title(ylab = "Mutual Information", line = 2.5)
  
  
  # Outlier thresholds.
  if (outlier_threshold > 0) {
    segments(0, outlier_threshold, max_distance, outlier_threshold, col = color_outlier, lty = 2, lwd = 2)
    text(0, outlier_threshold, "*", col = color_outlier, pos = 2, offset = 0.2, cex = 1, xpd = NA)
  }
  if (extreme_outlier_threshold > 0) {
    segments(0, extreme_outlier_threshold, max_distance, extreme_outlier_threshold, col = color_extreme_outlier, lty = 2, lwd = 2)
    text(0, extreme_outlier_threshold, "**", col = color_extreme_outlier, pos = 2, offset = 0.2, cex = 1, xpd = NA)
  }
  
  # LD Distance
  if (opts$ld_dist > 0) { segments(opts$ld_dist, min_mi, opts$ld_dist, 1, col = color_ld, lty = 2, lwd = 2) }
  
  legend(x = max_distance, y = max_mi, cex = cex_legend, pch = plot_symbol, bty = "n", xjust = 1.2, yjust = 1.2,
         c("Indirect", "Direct"), col = c(color_indirect, color_direct))
  
  text(0.90 * max_distance, min_mi + 0.8 * (max_mi - min_mi), paste0("#Connected: ", input_list$n_K_connected, "K"), cex = 0.8, adj = 0)
  text(0.90 * max_distance, min_mi + 0.78 * (max_mi - min_mi), paste0("#Disconnected: ", input_list$n_K_disconnected, "K"), cex = 0.8, adj = 0)
  
  dev.off()
  
} else {
  
  time_plotting_start <- proc.time()
  p1 = ggplot() + #(filtered_data, aes(x=V3, y=V5)) +
    ggrastr::rasterise( geom_hex(data = connected[which(connected$ARACNE != 1), ], aes(x = sep, y = MI), bins=500, fill="#c0c0c0"), dpi=1000) +
    ggrastr::rasterise( geom_hex(data = connected[which(connected$ARACNE == 1), ], aes(x = sep, y = MI), bins=500, fill="#3182bd"), dpi=1000) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 10)) +
    geom_hline(yintercept = outlier_threshold, col='#a50026', linetype='dashed', linewidth = 0.3) +
    geom_hline(yintercept = extreme_outlier_threshold, col='#d73027', linetype='dashed', linewidth = 0.3) +
    geom_vline(xintercept = c(opts$ld_dist), col='black', linetype='dashed', linewidth = 0.3) +
    theme_clean(base_size = 8) +
    theme(legend.position = "None") +
    xlab("Graph distance") +
    ylab("Mutual information")
  
  ggsave(filename = output_path, plot = p1, device = "png", width = plot_width, height = plot_height, units = "px")
  
}

time_plotting_end <- proc.time()
time_plotting <- (time_plotting_end - time_plotting_start)[[3]]
cat(paste0("All done! Time plotting: ", time_plotting, "s\n"))