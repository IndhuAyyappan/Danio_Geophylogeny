# Assignment 4 – Geography & Evolutionary Diversification in *Danio*

**Student:** Indhu Ayyappan\
**Course:** BINF\*6210 – Software Tools\
**Project Theme:** Theme 3 – Geography & Evolutionary Diversification

**Objective:**\
Investigate whether closely related species within the genus *Danio* occupy the same geographic regions or different regions, using COI sequences (NCBI) and species occurrences (GBIF).

------------------------------------------------------------------------

## Current Status

### Completed so far

-   Project structure initialized\
-   COI sequence dataset downloaded from NCBI (via `rentrez`)\
-   QC performed (sequence lengths, species name parsing, duplicates removed)\
-   One representative sequence per species retained\
-   DECIPHER multiple sequence alignment completed\
-   Distance matrix calculated using JC69 substitution model\
-   Neighbour-Joining tree constructed and rooted using *Microrasbora* outgroup\
-   Clean, publication-ready phylogeny generated (using `ggtree`)

### Next Steps

1.  Download species occurrence records from GBIF\
2.  Clean GBIF data (coordinates, species matching, outlier removal)\
3.  Match GBIF species names to phylogeny tip labels\
4.  Generate geophylogeny visualizations (tree + mapped range points)\
5.  Examine concordance between phylogenetic distances and spatial distances\
6.  Write interpretations (sympatric vs allopatric patterns)

------------------------------------------------------------------------
