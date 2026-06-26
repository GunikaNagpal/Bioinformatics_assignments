# ============================================================
# Phylogenetic Analysis of BRCA1 Gene — R Script
# Assignment | Bioinformatics | Skill-Up Internship, BioTecnika
# Species: Human, Mouse, Rat, Dog (4 species)
# Distance matrix and Newick trees from Python (Biopython)
# ============================================================

# ---- 1. Install and load packages ----
if (!requireNamespace("ape",      quietly = TRUE)) install.packages("ape")
if (!requireNamespace("phangorn", quietly = TRUE)) install.packages("phangorn")

library(ape)
library(phangorn)
cat("Packages loaded.\n")

# ---- 2. Load Newick trees from Python output ----
upgma_newick <- "(Dog:0.41807,((Rat:0.09365,Mouse:0.09365)Inner1:0.12003,Human:0.21368)Inner2:0.20439)Inner3:0.00000;"
nj_newick    <- "(Dog:0.62247,(Rat:0.09148,Mouse:0.09582)Inner1:0.10827,Human:0.22544)Inner2:0.00000;"

upgma_tree <- read.tree(text = upgma_newick)
nj_tree    <- read.tree(text = nj_newick)
cat("Trees loaded.\n")

# ---- 3. Distance matrix from Python output ----
species <- c("Human", "Mouse", "Rat", "Dog")

dist_values <- matrix(c(
  0.0000, 0.4320, 0.4228, 0.8479,
  0.4320, 0.0000, 0.1873, 0.8241,
  0.4228, 0.1873, 0.0000, 0.8247,
  0.8479, 0.8241, 0.8247, 0.0000
), nrow = 4, byrow = TRUE,
dimnames = list(species, species))

cat("\nPairwise Distance Matrix:\n")
print(round(dist_values, 4))

# ---- 4. Plot UPGMA Tree ----
png("brca1_upgma_R.png", width = 900, height = 600, res = 120)
plot(upgma_tree,
     main       = "UPGMA Phylogenetic Tree — BRCA1",
     cex        = 1.4,
     tip.color  = c("red", "blue", "blue", "darkgreen"),
     edge.width = 2.5,
     edge.color = "gray30",
     font       = 2)
add.scale.bar(cex = 0.9, col = "black", lwd = 2)
axisPhylo()
legend("topleft",
       legend = c("Primates (Human)",
                  "Rodents (Mouse, Rat)",
                  "Carnivores (Dog)"),
       col    = c("red", "blue", "darkgreen"),
       lwd    = 3, bty = "n", cex = 0.9)
dev.off()
cat("UPGMA tree saved.\n")

# ---- 5. Plot NJ Tree ----
png("brca1_nj_R.png", width = 900, height = 600, res = 120)
plot(nj_tree,
     main       = "Neighbor-Joining Tree — BRCA1",
     cex        = 1.4,
     tip.color  = c("darkgreen", "blue", "blue", "red"),
     edge.width = 2.5,
     edge.color = "gray30",
     font       = 2)
add.scale.bar(cex = 0.9, col = "black", lwd = 2)
axisPhylo()
legend("topleft",
       legend = c("Primates (Human)",
                  "Rodents (Mouse, Rat)",
                  "Carnivores (Dog)"),
       col    = c("red", "blue", "darkgreen"),
       lwd    = 3, bty = "n", cex = 0.9)
dev.off()
cat("NJ tree saved.\n")

# ---- 6. Fan/Circular tree ----
png("brca1_fan_R.png", width = 800, height = 800, res = 120)
plot(nj_tree,
     type       = "fan",
     main       = "BRCA1 Phylogenetic Tree — Fan Layout (NJ)",
     cex        = 1.5,
     tip.color  = c("darkgreen", "blue", "blue", "red"),
     edge.width = 3,
     edge.color = "steelblue",
     font       = 2)
dev.off()
cat("Fan tree saved.\n")

# ---- 7. Save Newick files ----
write.tree(upgma_tree, file = "brca1_upgma_from_R.nwk")
write.tree(nj_tree,    file = "brca1_nj_from_R.nwk")
cat("Newick files written.\n")

# ---- 8. Interpretation ----
cat("\n=== INTERPRETATION ===\n")
cat("Distance matrix observations:\n")
cat("  Human-Mouse distance: 0.4320\n")
cat("  Human-Rat distance:   0.4228\n")
cat("  Mouse-Rat distance:   0.1873 (closest pair — rodent monophyly confirmed)\n")
cat("  Human-Dog distance:   0.8479 (most distant)\n")
cat("  Mouse-Dog distance:   0.8241\n")
cat("  Rat-Dog distance:     0.8247\n")
cat("\nTree topology (both UPGMA and NJ agree):\n")
cat("  - Dog is correctly the outgroup (longest branch, most distant)\n")
cat("  - Mouse and Rat correctly cluster as sister taxa\n")
cat("  - Human is sister to the Mouse-Rat clade\n")
cat("  - Topology is biologically consistent with known mammalian phylogeny\n")
cat("\nConclusion:\n")
cat("  BRCA1 sequence-based phylogeny correctly recovers mammalian relationships.\n")
cat("  Rodents (Mouse, Rat) diverged ~87 Mya from primates.\n")
cat("  Dog (Carnivora) is the most distantly related species in this analysis.\n")
