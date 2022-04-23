## Args
# args <- commandArgs(trailingOnly = TRUE)
# tree_in <- args[1]
# figfile_out <- args[2]
# if (length(args >= 3)) outgroup <- args[3] else outgroup <- NA

iqtree_dir <- "results/iqtree/woutgroups"
tree_in <- file.path(iqtree_dir, "COI.treefile")
figfile_out_png <- file.path(iqtree_dir, "COI.png")
figfile_out_svg <- file.path(iqtree_dir, "COI.svg")
outgroup <- "KX421869.1"

## Report
#cat("Nr args:        ", length(args), "\n")
cat("Input tree:     ", tree_in, "\n")
cat("Output figure:  ", figfile_out, "\n")
cat("Outgroup:       ", outgroup, "\n")

## Packages
if (!require(pacman)) install.packages("pacman")
if (!require(BiocManager)) install.packages("BiocManager")
packages <- c("tidyverse", "here", "ape", "ggtree", "argparse")
pacman::p_load(char = packages, install = TRUE, repos = "http://cran.us.r-project.org")

## Read the tree file
tree <- read.tree(tree_in)
if (!is.na(outgroup)) tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)

## Annot
focal_IDs <- paste0("Wooster", 1:6, "_Ohio_USA_2021")
annot <- data.frame(tiplab = tree$tip.label) %>%
  mutate(plotlab = sub("^_R_", "", tiplab),
         plotlab = sub("HQ944971.1", "HQ944971.1_Ochlerotatus canadensis", plotlab),
         plotlab = sub("KX421869.1", "KX421869.1_Lebertia maderigena", plotlab),
         plotlab = sub("JX629050.1", "JX629050.1_Torrenticola lundbladi", plotlab),
         plotlab = sub("GQ254798.1", "GQ254798.1_Iowa_USA_2008", plotlab),
         plotlab = sub("GQ254793.1_Iowa_USA_2007", "GQ254793.1_Illinois_USA_2007", plotlab),
         is_query = ifelse(plotlab %in% focal_IDs, TRUE, FALSE),
         plotlab = gsub("_", " - ", plotlab))

## Remove the distant outgroups
to_drop <- c("KX421869.1", "JX629050.1")
tree <- drop.tip(tree, to_drop)

## Make the tree
p <- ggtree(tree, size = 1) %<+% annot +
  geom_tiplab(aes(label = plotlab, color = is_query), size = 4) +
  geom_text2(aes(label = label,
                 subset = !is.na(as.numeric(label)) & as.numeric(label) > 95),
             nudge_x = -0.0027, nudge_y = 0.17,
             size = 4, fontface = "italic", color = "grey50") +
  geom_treescale(x = 0.1, y = 0, color = "grey50") +
  geom_rootedge(rootedge = 0.005) +
  scale_color_manual(values = c("black", "red"), guide = "none") +
  guides(size = "none", color = "none") +
  theme(plot.margin = margin(0.2, 7, 0.2, 0.2, "cm")) +
  coord_cartesian(clip = "off")
p
ggsave(figfile_out_png, p, width = 14, height = 7)
ggsave(figfile_out_svg, p, width = 14, height = 7)
