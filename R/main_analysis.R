# ================================
## BINF*6210 – Assignment 4
## Theme 3: Geography & Evolutionary Diversification
## Student: Indhu Ayyappan
## Project: Geophylogeny of genus Danio (COI + GBIF)
##
## Main question:
## Do closely related Danio species occur in the same
## geographic regions or different regions?
##
## Optional sub-question:
## Is geographic distance between species ranges
## correlated with phylogenetic distance?
# ================================

#### 1. SETUP ----

# Load required packages
#uncomment packages that needs to be installed
#BiocManager::install("ggtree")
library(tidyverse)
library(ape)
library(phangorn)
library(DECIPHER)
library(rgbif)
library(rentrez)
library(stringr)
library(ggtree)
library(ggplot2)

# Folder paths (relative)
dir_raw   <- "data_raw"
dir_clean <- "data_clean"
dir_figs  <- "figs"

#### 2. DATA ACQUISITION – SEQUENCES (NCBI) ----
query <- "Danio[Organism] AND (COI OR CO1)"

seq_search <- entrez_search(
  db = "nucleotide",
  term = query,
  retmax = 500
)
#### 2B. FETCH FASTA SEQUENCES ----

# fetch all sequences in FASTA format
seq_fasta <- entrez_fetch(
  db = "nucleotide",
  id = seq_search$ids,
  rettype = "fasta",
  retmode = "text"
)
# save raw fasta file into data_raw
fasta_file <- file.path(dir_raw, "danio_coi_raw.fasta")
write(seq_fasta, file = fasta_file)

#### 3. SEQUENCE QC & CLEANING ----
# In this section, I perform basic quality control:
# - extract species names
# - filter by expected COI length
# - retain the longest sequence per species
# - prepare a clean FASTA file for alignment

# read fasta into R for QC
danio_raw <- ape::read.FASTA(fasta_file)
length(danio_raw)          
head(names(danio_raw), 10)
#### 3A. BASIC QC: SEQUENCE LENGTHS ----
seq_lengths <- sapply(danio_raw, length)
summary(seq_lengths)
#### 3B. EXTRACT SPECIES NAMES ----
headers <- names(danio_raw)
species_tbl <- tibble(
  header  = headers,
  genus   = word(headers, 2),
  species = word(headers, 3)
) %>%
  mutate(
    species_name = paste(genus, species)  
  )

head(species_tbl, 10)
#### 3C. ADD SEQUENCE LENGTHS ----
length_tbl <- tibble(
  header = names(danio_raw),
  length_bp = seq_lengths
)
# join with species table
seq_info <- left_join(species_tbl, length_tbl, by = "header")

head(seq_info, 10)
#### 3D. FILTER BY LENGTH (KEEP COI-SIZE) ----
seq_filt <- seq_info %>%
  filter(
    length_bp >= 500,
    length_bp <= 800
  )

nrow(seq_filt)
length(unique(seq_filt$species_name))
head(seq_filt, 10)
#### 3E. KEEP ONE SEQUENCE PER SPECIES (LONGEST) ----
seq_best <- seq_filt %>%
  arrange(species_name, desc(length_bp)) %>%
  distinct(species_name, .keep_all = TRUE)

nrow(seq_best)          
head(seq_best, 10)
#### 3F. SUBSET FASTA TO CLEANED HEADERS ----

clean_headers <- seq_best$header

danio_clean <- danio_raw[clean_headers]

length(danio_clean)
names(danio_clean)[1:5]
#### 3G. SAVE CLEAN FASTA ----

clean_fasta_file <- file.path(dir_clean, "danio_coi_clean.fasta")
write.FASTA(danio_clean, clean_fasta_file)
#### 4. ALIGNMENT (DECIPHER) ----

clean_fasta <- file.path(dir_clean, "danio_coi_clean.fasta")
danio_clean <- readDNAStringSet(clean_fasta)
danio_clean

#### 4B. ALIGN WITH DECIPHER ----

alignment <- AlignSeqs(danio_clean)
alignment

# use genus + species as labels, and abbreviate Danio -> D.
short_labels <- str_extract(
  names(alignment),
  "Danio [A-Za-z'.]+|Microrasbora [A-Za-z'.]+"
)

# if the pattern fails for any record, fall back to the original name
short_labels <- ifelse(is.na(short_labels), names(alignment), short_labels)

#Replace full genus names ("Danio species") with abbreviated labels ("D. species") for cleaner plotting
short_labels <- sub("^Danio ", "D. ", short_labels)

names(alignment) <- short_labels

#### 4C. CONVERT ALIGNMENT TO phyDat (phangorn) ----

alignment_phydat <- phyDat(as.matrix(alignment), type = "DNA")
alignment_phydat

#### 4D. DISTANCE MATRIX (JC69) ----

dist_coi <- dist.ml(alignment_phydat, model = "JC69")
dist_coi

#### 5. PHYLOGENETIC TREE – NEIGHBOUR JOINING ----

#5A. Build an initial NJ tree from the JC69 distance matrix.
#This gives us an unrooted tree based only on COI pairwise distances.
tree_nj <- nj(dist_coi)

#5B. Root the tree using Microrasbora as the outgroup.
#Biologically, Microrasbora is outside the Danio clade, so this gives a direction to the tree (root at Microrasbora branch).
out_lbl <- grep("Microrasbora", tree_nj$tip.label, value = TRUE)

tree_rooted <- root(tree_nj,
                    outgroup = out_lbl,
                    resolve.root = TRUE)

#5C. Ladderize to "tidy" the topology.
#Ladderizing just reorders the nodes so the tree plots in a cleaner, more readable way (no effect on branch lengths or relationships).
tree_rooted <- ladderize(tree_rooted)

#5D. Re-compute branch lengths for plotting.
#NJ + rooting can sometimes introduce tiny negative branches, which ggtree complains about. Here I apply a Grafen transform to get a clean, strictly non-negative set of edge lengths for plotting.
tree_plot <- compute.brlen(tree_rooted, method = "Grafen")

# 5E. Calculate overall tree height.
#node.depth.edgelength() returns the cumulative root-to-tip distance. After Grafen transformation this value = 1, but we keep it dynamic so the code remains generalizable.
max_d <- max(node.depth.edgelength(tree_plot))

#### 5F. BASE R QUICK CHECK ----
#Quick sanity check of the rooted tree before making the nicer plot.
par(mar = c(1, 1, 1, 6))
plot(tree_plot, cex = 0.9, no.margin = TRUE)
add.scale.bar()

#### 5G. PUBLICATION-STYLE TREE WITH GGTREE ----

#Now I use ggtree so that the labels are aligned and the x-axis is scaled in substitutions per site. The xlim() gives extra room for the dotted label guides.
p_tree <- ggtree(tree_plot) +
  xlim(0, max_d * 1.2) +                    
  geom_tiplab(align = TRUE,
              linetype = "dotted",
              size = 3,
              hjust = 0) +
  theme_tree2() +                           
  xlab("Relative divergence (scaled branch lengths)")

p_tree
#### 6. MATCH PHYLOGENY + GEOGRAPHY ----
# (placeholder)

#### 7. VISUALIZATIONS (GEOPHYLOGENY, MAPS, PLOTS) ----
# (placeholder)

#### 8. SAVE CLEAN DATA & FIGURES ----
# (placeholder)