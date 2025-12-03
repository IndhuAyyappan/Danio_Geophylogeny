Assignment 4 – Geography & Evolutionary Diversification in Danio

Student: Indhu Ayyappan Course: BINF\*6210 – Software Tools Project Theme: Theme 3 – Geography & Evolutionary Diversification

Objective: Investigate whether closely related species within the genus Danio occupy the same geographic regions or different regions, using COI sequences (NCBI) and species occurrences (GBIF).

⸻

Current Status

Completed so far • Project structure initialized\
• COI sequence dataset downloaded from NCBI (rentrez)\
• QC performed (sequence lengths, species parsing, duplicates removed)\
• One representative sequence per species retained\
• Multiple sequence alignment completed using DECIPHER\
• Distance matrix calculated using JC69 substitution model\
• Neighbour–Joining tree constructed\
• Tree successfully rooted using Microrasbora erythromicron as outgroup\
• Tree ladderized and Grafen-scaled for clearer visualization\
• Final publication-ready phylogeny generated using ggtree\
• GBIF occurrence data downloaded and cleaned (coordinate QC, species name filtering, deduplication)\
• GBIF species names reconciled to phylogeny tip labels using a custom mapping (Danio → D. abbreviations, Microrasbora handled separately)\
• Clades assigned from the COI phylogeny (Annulatus, Kerri–kyathit, Choprae, Margaritatus, Roseus, Outgroup, Other)\
• Species-centroid coordinates calculated for map visualization\
• Asia basemap generated using map_data()\
• Final geophylogeny figure produced (COI phylogeny + geographic occurrences)\
• Final sampling-intensity barplot produced (clean GBIF counts per phylogeny species)

⸻

Visualizations Produced

1.  COI Phylogeny for Danio • Neighbour–Joining tree • Rooted using Microrasbora erythromicron • Ladderized and Grafen-scaled • Clean tip labels • High-resolution layout using ggtree

2.  Geographic Occurrence Map • Cleaned GBIF points shown in faint grey (background sampling density) • One centroid point per species (coloured by major clade) • Custom, colour-blind-safe clade palette • Legend ordered + formatted cleanly • High-resolution basemap (South + Southeast Asia)

3.  Geophylogeny (Tree + Map) • Combined using patchwork • Phylogeny on top, map below • Consistent styling and font size • Publication-quality layout

4.  Sampling Intensity Barplot • Cleaned GBIF records per phylogeny species • Species ordered by count • Useful for comparing data availability and spatial sampling bias

⸻

Next Steps 1. Examine spatial patterns for evidence of: • sympatry (co-occurring closely related species) • allopatry (geographic separation along phylogeny) 2. Quantify phylogenetic vs. geographic distances (e.g., Mantel test)\
3. Interpret biogeographical patterns at the clade level\
4. Write final report section discussing diversification and spatial structuring

⸻
