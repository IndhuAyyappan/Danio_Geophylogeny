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
#install.packages(c("maps", "patchwork"))
# install.packages("geosphere")
library(tidyverse)
library(ape)
library(phangorn)
library(DECIPHER)
library(rgbif)
library(rentrez)
library(stringr)
library(ggtree)
library(ggplot2)
library(maps)     
library(patchwork)
library(viridis)
library(dplyr)
library(geosphere)

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
  geom_tiplab(
    size  = 2.4,
    hjust = 0
  ) +
  xlim(0, max_d * 1.35) +   # a bit of extra room for labels
  theme_tree2() +
  xlab("Relative divergence (scaled branch lengths)")

p_tree
#### 6. GBIF OCCURRENCE DATA ----

# In this section I pull occurrence records from GBIF for the same Danio
# species that I used in the COI phylogeny. The goal is to later link
# range locations to the tips of the tree.

#### 6A. BUILD SPECIES LIST FROM PHYLOGENY ----

# I use the cleaned QC table (seq_best) to get the full species names
# that actually made it into the phylogeny (1 sequence per species).
danio_species <- seq_best %>%
  dplyr::pull(species_name) %>%
  unique() %>%
  sort()

danio_species
length(danio_species)   # number of Danio species represented

#### 6B. DOWNLOAD GBIF RECORDS PER SPECIES ----

# For each Danio species, I call rgbif::occ_data() and request up to
# 500 records with coordinates. This is usually plenty for an assignment-scale
# project and keeps the download manageable.
gbif_list <- purrr::map(
  danio_species,
  ~{
    message("Downloading GBIF records for: ", .x)
    rgbif::occ_data(
      scientificName = .x,
      hasCoordinate  = TRUE,
      limit          = 500
    )$data
  }
)

# Combine all species into one data frame and keep track of the name.
gbif_raw <- dplyr::bind_rows(gbif_list, .id = "species_index") %>%
  dplyr::mutate(
    species_name = scientificName   # standardise column name
  )

dim(gbif_raw)
dplyr::count(gbif_raw, species_name, name = "n_records_raw")

#### 6C. BASIC GBIF CLEANING ----

# Here I do a light but sensible cleaning step so the downstream maps
# aren’t dominated by obviously bad points.
gbif_clean <- gbif_raw %>%
  # keep only records with valid lat/lon
  dplyr::filter(!is.na(decimalLongitude),
                !is.na(decimalLatitude)) %>%
  # drop the classic (0, 0) "in the ocean" errors
  dplyr::filter(!(decimalLongitude == 0 & decimalLatitude == 0)) %>%
  # remove records with very high coordinate uncertainty (> 50 km)
  dplyr::filter(
    is.na(coordinateUncertaintyInMeters) |
      coordinateUncertaintyInMeters <= 50000
  ) %>%
  # keep a few sensible basisOfRecord types
  dplyr::filter(
    basisOfRecord %in% c(
      "HUMAN_OBSERVATION",
      "OBSERVATION",
      "MATERIAL_SAMPLE",
      "PRESERVED_SPECIMEN"
    )
  ) %>%
  # keep the columns that will be useful for mapping and summaries
  dplyr::transmute(
    species_name,
    lon  = decimalLongitude,
    lat  = decimalLatitude,
    countryCode,
    year,
    basisOfRecord,
    coord_uncert_m = coordinateUncertaintyInMeters
  )

dim(gbif_clean)

# Quick summary: how many cleaned records per species?
gbif_summary <- gbif_clean %>%
  dplyr::count(species_name, name = "n_records_clean") %>%
  dplyr::arrange(dplyr::desc(n_records_clean))

gbif_summary

#### 6D. SAVE CLEANED GBIF DATA ----

# I save both the cleaned table (all points) and the per-species summary
# so that I can reuse them later for mapping and for the write-up.
gbif_clean_file   <- file.path(dir_clean, "danio_gbif_clean.csv")
gbif_summary_file <- file.path(dir_clean, "danio_gbif_summary.csv")

readr::write_csv(gbif_clean,   gbif_clean_file)
readr::write_csv(gbif_summary, gbif_summary_file)

#### 7. VISUALIZATIONS (GEOPHYLOGENY, MAPS, SAMPLING) ----

### 7A. Make GBIF names line up with tree tip labels ----
# Problem we saw:
#   - Tree tips look like: "D. aesculapii", "Microrasbora erythromicron"
#   - GBIF species_name look like: "Danio aesculapii Kullander & Fang, 2009"
#   - So our earlier tip_label did NOT match tree_tips and everything was filtered out.
#
# Here I:
#   1. Strip author + year off the GBIF names (keep only genus + species)
#   2. Abbreviate "Danio <species>" -> "D. <species>" (to match the phylogeny)
#   3. Keep Microrasbora as full genus (matches tree labels)
#   4. Handle the Brachydanio synonym manually.
#   5. Drop BOLD: pseudo-taxa.

gbif_named <- gbif_clean %>%
  # remove weird BOLD “species” rows
  dplyr::filter(!stringr::str_starts(species_name, "BOLD:")) %>%
  # pull out just genus + species (first two words)
  dplyr::mutate(
    genus   = stringr::word(species_name, 1),
    sp      = stringr::word(species_name, 2),
    short_name = paste(genus, sp)
  ) %>%
  # create a tip_label that is formatted like the tree labels
  dplyr::mutate(
    tip_label = dplyr::case_when(
      genus == "Danio"        ~ paste("D.", sp),
      genus == "Microrasbora" ~ paste(genus, sp),
      genus == "Brachydanio"  ~ "D. albolineatus",  # known synonym
      TRUE                    ~ short_name          # fallback
    )
  )

# Tip labels from the phylogeny
tree_tips <- tree_plot$tip.label

# Quick sanity check: which tip labels are shared?
intersect(tree_tips, unique(gbif_named$tip_label))

# Keep only GBIF records for species that actually appear in the tree
gbif_phylo <- gbif_named %>%
  dplyr::filter(tip_label %in% tree_tips)

# Sanity: how many records per tip?
gbif_phylo %>%
  dplyr::count(tip_label, name = "n_records") %>%
  dplyr::arrange(dplyr::desc(n_records))
# One mean coordinate per species (for labelled centroids)
gbif_centroids <- gbif_phylo %>%
  dplyr::group_by(tip_label) %>%
  dplyr::summarise(
    lon = mean(lon, na.rm = TRUE),
    lat = mean(lat, na.rm = TRUE),
    .groups = "drop"
  )

### 7B. Base map for South / Southeast Asia ----
# I draw a simple world map and crop it to the region where Danio occurs.

world_df <- map_data("world")

asia_bb <- world_df %>%
  dplyr::filter(long > 60, long < 130,
                lat  > -10, lat < 40)


### 7C. Occurrence map for Danio species used in the phylogeny ----
# Each point = a cleaned GBIF record, coloured by species (tree tip label).
# NOTE: use lon/lat columns here (decimalLongitude/Latitude no longer exist).

p_map <- ggplot() +
  geom_polygon(
    data = asia_bb,
    aes(x = long, y = lat, group = group),
    fill   = "grey96",
    colour = "grey80",
    linewidth = 0.2
  ) +
  geom_point(
    data = gbif_phylo,
    aes(x = lon, y = lat),
    colour = "grey75",
    alpha  = 0.4,
    size   = 1
  ) +
  geom_point(
    data = gbif_centroids,
    aes(x = lon, y = lat, colour = tip_label),
    size  = 3,
    alpha = 0.9
  ) +
  coord_quickmap(xlim = c(70, 125),
                 ylim = c(-5, 35)) +
  scale_colour_viridis_d(
    option = "D",
    name   = "Species"
  ) +
  guides(
    colour = guide_legend(
      ncol = 2,
      override.aes = list(size = 4, alpha = 1)
    )
  ) +
  labs(
    title = "GBIF occurrences for Danio species used in the COI phylogeny",
    x     = "Longitude",
    y     = "Latitude"
  ) +
  theme_bw() +
  theme(
    plot.title        = element_text(face = "bold", hjust = 0),
    legend.position   = "right",
    legend.title      = element_text(face = "bold"),
    legend.key.height = unit(0.45, "lines")
  )

### 7D. Geophylogeny: tree + map in one figure ----
# I put the Grafen-scaled COI tree on top and the occurrence map below.
# This keeps both readable and still shows “phylogeny + geography”.

p_tree_clean <- p_tree +
  labs(title = "COI phylogeny for Danio (Grafen-scaled)") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0)
  )

fig_geophylo <- p_tree_clean / p_map +
  plot_layout(heights = c(1, 1.2))


### 7E. Sampling intensity barplot ----
# Now I summarise how many cleaned GBIF records each phylogeny species has.
# I summarise directly from gbif_phylo so that:
#   - names already match the tree
#   - we use the correct count column name.

gbif_summary_phylo <- gbif_phylo %>%
  dplyr::count(tip_label, name = "n_records_clean") %>%
  dplyr::arrange(dplyr::desc(n_records_clean))

fig_sampling <- gbif_summary_phylo %>%
  mutate(tip_label = forcats::fct_reorder(tip_label, n_records_clean)) %>%
  ggplot(aes(x = n_records_clean, y = tip_label)) +
  geom_col(width = 0.6, fill = "grey30") +
  labs(
    title = "Sampling intensity for Danio species in COI phylogeny",
    x     = "Number of cleaned GBIF records",
    y     = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0)
  )

p_tree_clean <- p_tree +
  labs(title = "COI phylogeny for Danio (Grafen-scaled)") +
  theme(plot.title = element_text(face = "bold", hjust = 0))

fig_geophylo <- p_tree_clean / p_map + plot_layout(heights = c(1, 1.2))

fig_geophylo
fig_sampling