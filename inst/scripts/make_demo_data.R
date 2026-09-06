# Generates the package's demo datasets (demo_z, demo_pathways).
# Fully synthetic and seeded; rerun to reproduce data/*.rda exactly.
set.seed(20260905)

genes <- sprintf("g%03d", 1:400)
demo_pathways <- list(
    PATHWAY_A = genes[1:30],    PATHWAY_B = genes[51:80],
    PATHWAY_C = genes[101:130], PATHWAY_D = genes[151:180],
    PATHWAY_E = genes[201:230], PATHWAY_F = genes[251:280])

# perturb every member of A, C, E plus 30 pathway-free genes
perts <- c(demo_pathways$PATHWAY_A, demo_pathways$PATHWAY_C,
    demo_pathways$PATHWAY_E, genes[301:330])
demo_z <- matrix(stats::rnorm(length(perts) * length(genes)),
    length(perts), length(genes),
    dimnames = list(perts, genes))

# planted regulation (CRISPRi orientation: knocking down an activator
# LOWERS the target pathway's genes, z < 0):
#   A activates B  -> perturbing A members: z(B genes) ~ -1.5
#   C suppresses D -> perturbing C members: z(D genes) ~ +1.5
demo_z[demo_pathways$PATHWAY_A, demo_pathways$PATHWAY_B] <-
    demo_z[demo_pathways$PATHWAY_A, demo_pathways$PATHWAY_B] - 1.5
demo_z[demo_pathways$PATHWAY_C, demo_pathways$PATHWAY_D] <-
    demo_z[demo_pathways$PATHWAY_C, demo_pathways$PATHWAY_D] + 1.5

# strong on-target knockdown (exercises the tautology exclusion)
for (p in perts) demo_z[p, p] <- -8

save(demo_pathways, file = "data/demo_pathways.rda", compress = "xz")
save(demo_z, file = "data/demo_z.rda", compress = "xz")
