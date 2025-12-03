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
#install.packages("paletteer")
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
p_tree_standalone <- ggtree(tree_plot) +
  geom_tiplab(
    size  = 2.8,
    hjust = 0
  ) +
  xlim(0, max_d * 1.30) +     # room for labels
  theme_tree() +            
  labs(
    title = "COI phylogeny for Danio "
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
  )

p_tree_standalone
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

gbif_named <- gbif_clean %>%
  # remove BOLD “species” rows
  dplyr::filter(!stringr::str_starts(species_name, "BOLD:")) %>%
  # pull out just genus + species (first two words)
  dplyr::mutate(
    genus      = stringr::word(species_name, 1),
    sp         = stringr::word(species_name, 2),
    short_name = paste(genus, sp)
  ) %>%
  # create a tip_label that is formatted like the tree labels
  dplyr::mutate(
    tip_label = dplyr::case_when(
      # special case: Danio erythromicron actually belongs in Microrasbora
      genus == "Danio"       & sp == "erythromicron" ~ "Microrasbora erythromicron",
      genus == "Danio"                               ~ paste("D.", sp),
      genus == "Microrasbora"                        ~ paste(genus, sp),
      genus == "Brachydanio"                         ~ "D. albolineatus",
      TRUE                                           ~ short_name
    )
  )

# Tip labels from the phylogeny (in tree order)
tree_tips <- tree_plot$tip.label

# Keep only GBIF records for species that actually appear in the tree
# and lock the factor levels to the tree tip order
gbif_phylo <- gbif_named %>%
  dplyr::filter(tip_label %in% tree_tips) %>%
  dplyr::mutate(
    tip_label = factor(tip_label, levels = tree_tips)
  )

## Create clade assignments (broad COI clades from the tree)
clade_df <- tibble(
  tip_label = tree_tips
) %>%
  dplyr::mutate(
    clade = dplyr::case_when(
      tip_label %in% c("D. roseus", "D. albolineatus", "D. tweediei") ~ "Roseus clade",
      tip_label %in% c("D. kerri", "D. sp.", "D. kyathit", "D. aff.",
                       "D. aesculapii", "D. tinwini", "D. nigrofasciatus") ~ "Kerri–kyathit clade",
      tip_label %in% c("D. catenatus", "D. annulosus", "D. sysphigmatus",
                       "D. dangila", "D. meghalayensis", "D. assamila",
                       "D. rerio", "D. cf.") ~ "Annulatus clade",
      tip_label %in% c("D. choprae", "D. choprai", "D. feegrarei",
                       "D. htammanthin us", "D. flagrans") ~ "Choprae clade",
      tip_label %in% c("D. margaritatus", "D. erythromicron") ~ "Margaritatus clade",
      tip_label == "Microrasbora erythromicron" ~ "Outgroup",
      TRUE ~ "Other"
    )
  )

## Join clade labels to GBIF data
gbif_phylo <- gbif_phylo %>%
  dplyr::left_join(clade_df, by = "tip_label")

# One mean coordinate per species (for labelled centroids)
gbif_centroids <- gbif_phylo %>%
  dplyr::group_by(tip_label, clade) %>%
  dplyr::summarise(
    lon = mean(lon, na.rm = TRUE),
    lat = mean(lat, na.rm = TRUE),
    .groups = "drop"
  )

### 7B. Base map for South / Southeast Asia ----

world_df <- map_data("world")

asia_bb <- world_df %>%
  dplyr::filter(long > 60, long < 130,
                lat  > -10, lat < 40)

### 7C. Occurrence map for Danio species used in the phylogeny ----

p_map <- ggplot() +
  geom_polygon(
    data = asia_bb,
    aes(x = long, y = lat, group = group),
    fill   = "grey96",
    colour = "grey80",
    linewidth = 0.2
  ) +
  # all cleaned GBIF points (faint background in grey)
  geom_point(
    data = gbif_phylo,
    aes(x = lon, y = lat),
    colour = "grey20",
    alpha  = 0.25,
    size   = 0.5
  ) +
  # species centroids (coloured by COI clade)
  geom_point(
    data   = gbif_centroids %>% left_join(clade_df),
    aes(x = lon, y = lat, fill = clade),
    shape  = 21,        # filled circle with border
    colour = "black",   # outline
    size   = 3,
    alpha  = 0.6
  ) +
  scale_fill_manual(
    name   = "Major Danio clades",
    values = c(
      "Annulatus clade"      = "#1b9e77",
      "Choprae clade"        = "#7570b3",
      "Kerri–kyathit clade"  = "#d95f02",
      "Margaritatus clade"   = "#e7298a",
      "Roseus clade"         = "#1f78b4",
      "Outgroup"             = "#000000",
      "Other"                = "#666666"
    )
  ) +
  coord_quickmap(xlim = c(70, 125),
                 ylim = c(-5, 35)) +
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

p_tree_clean <- ggtree(tree_plot) +
  geom_tiplab(
    size  = 2.4,
    hjust = 0
  ) +
  xlim(0, max_d * 2.2) +   # extra room for labels
  labs(
    title = "COI phylogeny for Danio"
  ) +
  theme_tree() +           # removes x-axis line and ticks
  theme(
    plot.title = element_text(face = "bold", hjust = 0)
  )

fig_geophylo <- p_tree_clean / p_map +
  plot_layout(heights = c(1, 1.2))

### 7E. Sampling intensity barplot ----

gbif_summary_phylo <- gbif_phylo %>%
  dplyr::count(tip_label, name = "n_records_clean") %>%
  dplyr::arrange(dplyr::desc(n_records_clean))

fig_sampling <- gbif_summary_phylo %>%
  dplyr::mutate(
    tip_label = forcats::fct_reorder(tip_label, n_records_clean)
  ) %>%
  ggplot(aes(x = n_records_clean, y = tip_label)) +
  geom_col(width = 0.6, fill = "#1f78b4") +
  labs(
    title = "Sampling intensity for Danio species in COI phylogeny",
    x     = "Number of cleaned GBIF records",
    y     = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0)
  )

fig_geophylo
fig_sampling
### 8. Sister-species geographic separation ----
## 8A. Use existing tree for sister-pair search ----
tree <- tree_plot

# Helper: extract sister pairs (both children are tips)
get_sister_pairs <- function(tree) {
  pairs <- list()
  count <- 1
  
  for (node in (Ntip(tree) + 1):(tree$Nnode + Ntip(tree))) {
    children <- tree$edge[tree$edge[, 1] == node, 2]
    if (all(children <= Ntip(tree))) {
      pairs[[count]] <- tree$tip.label[children]
      count <- count + 1
    }
  }
  pairs
}

sister_pairs <- get_sister_pairs(tree)
sister_pairs  

## 8B. Centroids for each species used in the phylogeny ----
centroids <- gbif_phylo %>%
  dplyr::group_by(tip_label) %>%
  dplyr::summarise(
    lon = mean(lon, na.rm = TRUE),
    lat = mean(lat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(centroids) <- centroids$tip_label

## 8C. Keep only sister pairs where BOTH species have centroids ----
usable_sisters <- purrr::keep(
  sister_pairs,
  ~ all(.x %in% rownames(centroids))
)

usable_sisters

## 8D. Geographic distance (km) between sister centroids ----
sister_geo_dist <- purrr::map_dfr(
  seq_along(usable_sisters),
  function(i) {
    sp <- usable_sisters[[i]]
    sp1 <- sp[1]
    sp2 <- sp[2]
    
    d_km <- geosphere::distHaversine(
      centroids[sp1, c("lon", "lat")],
      centroids[sp2, c("lon", "lat")]
    ) / 1000  # metres → km
    
    data.frame(
      pair_id = paste0("Pair ", i, ": ", sp1, " – ", sp2),
      sp1     = sp1,
      sp2     = sp2,
      dist_km = d_km,
      stringsAsFactors = FALSE
    )
  }
)

sister_geo_dist <- sister_geo_dist %>%
  mutate(
    pair_label = paste0(
      "Pair ", row_number(), ": ",
      sp1, " – ", sp2,
      " (", round(dist_km), " km)"
    )
  )
sister_geo_dist
## 8E. Build line segments between sister centroids ----
sister_lines <- purrr::map_dfr(
  seq_along(usable_sisters),
  function(i) {
    sp <- usable_sisters[[i]]
    
    this_label <- sister_geo_dist$pair_label[i]
    
    data.frame(
      pair_label = this_label,
      sp1        = sp[1],
      sp2        = sp[2],
      lon1       = centroids[sp[1], "lon"],
      lat1       = centroids[sp[1], "lat"],
      lon2       = centroids[sp[2], "lon"],
      lat2       = centroids[sp[2], "lat"],
      stringsAsFactors = FALSE
    )
  }
)
sister_lines

## 8F. Map: sister-pair connections on top of all Danio points ----
fig_sister_map <- ggplot() +
  geom_polygon(data = asia_bb,
               aes(x = long, y = lat, group = group),
               fill   = "grey95",
               colour = "grey80",
               linewidth = 0.2) +
  geom_point(data = gbif_phylo,
             aes(x = lon, y = lat),
             colour = "grey80",
             alpha  = 0.4,
             size   = 0.7) +
  geom_segment(data = sister_lines,
               aes(x = lon1, y = lat1,
                   xend = lon2, yend = lat2,
                   colour = pair_label),
               linewidth = 1.1) +
  geom_point(data = sister_lines,
             aes(x = lon1, y = lat1, colour = pair_label),
             size = 2.3) +
  geom_point(data = sister_lines,
             aes(x = lon2, y = lat2, colour = pair_label),
             size = 2.3) +
  coord_quickmap(xlim = c(70, 125),
                 ylim = c(-5, 35)) +
  scale_colour_brewer(palette = "Dark2", name = "Sister pairs") +
  labs(
    title = "Geographic separation between sister Danio species",
    x     = "Longitude",
    y     = "Latitude"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0),
    legend.position = "right"
  )

fig_sister_map
