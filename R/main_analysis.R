##***************************
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
##***************************

#### 1. SETUP ----

# Load required packages
library(tidyverse)
library(ape)
library(phangorn)
library(DECIPHER)
library(rgbif)
library(rentrez)

# Folder paths (relative)
dir_raw   <- "data_raw"
dir_clean <- "data_clean"
dir_figs  <- "figs"

#### 2. DATA ACQUISITION – SEQUENCES (NCBI) ----
# (placeholder – we'll fill this section together)
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

#### 3. SEQUENCE QC & ALIGNMENT ----
# read fasta into R for QC
danio_raw <- read.FASTA(fasta_file)

length(danio_raw)          
head(names(danio_raw), 10)

#### 4. PHYLOGENY (MODEL TEST + ML TREE) ----
# (placeholder)

#### 5. DATA ACQUISITION – GEOGRAPHY (GBIF) ----
# (placeholder)

#### 6. MATCH PHYLOGENY + GEOGRAPHY ----
# (placeholder)

#### 7. VISUALIZATIONS (GEOPHYLOGENY, MAPS, PLOTS) ----
# (placeholder)

#### 8. SAVE CLEAN DATA & FIGURES ----
# (placeholder)