## Args
args <- commandArgs(trailingOnly = TRUE)
tree_in <- args[1]
figfile_out <- args[2]
if (length(args >= 3)) outgroup <- args[3] else outgroup <- NA

# tree_in <- "results/iqtree/COI.treefile"
# outgroup <- "KM638737.1"

## Report
cat("Nr args:        ", length(args), "\n")
cat("Input tree:     ", tree_in, "\n")
cat("Output figure:  ", figfile_out, "\n")
cat("Outgroup:       ", outgroup, "\n")

## Packages
packages <- c("tidyverse", "here", "ape", "ggtree", "argparse")
pacman::p_load(char = packages, install = TRUE, repos = "http://cran.us.r-project.org")

## Read the tree file
tree <- read.tree(tree_in)
if (!is.na(outgroup)) tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)

## Annot
focal_IDs <- paste0("Wooster", 1:6, "_Ohio_USA_2021")
annot <- data.frame(tiplab = tree$tip.label) %>%
    mutate(
        plotlab = gsub("^_R_", "", tiplab),
        is_query = ifelse(plotlab %in% focal_IDs, TRUE, FALSE)
    )

## Make the tree
p <- ggtree(tree) %<+% annot +
    geom_tiplab(aes(label = plotlab, color = is_query)) +
    scale_color_manual(values = c("black", "red"), guide = "none") +
    guides(size = "none", color = "none") +
    theme(plot.margin = margin(0.2, 7, 0.2, 0.2, "cm")) +
    coord_cartesian(clip = "off")

ggsave(figfile_out, p, width = 14, height = 7)