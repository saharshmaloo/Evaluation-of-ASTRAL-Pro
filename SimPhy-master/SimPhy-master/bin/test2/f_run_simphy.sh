#!/bin/bash

# exit 0

# Simpy Parameters:
# -----------------
# -rs : Number of species tree replicates
# -rl : Number of locus trees per species tree
# -rg : Number of gene trees per locus tree
# -gb : Genome-wide (i.e., sampled for each species tree) duplication rate
#       parameter (to use with LB)
# -gd : Genome-wide (i.e., sampled for each species tree) loss parameter
#       (to use with LD)
# -gg : Genome-wide (i.e., sampled for each species tree) gene conversion
#       parameter (to use with LG)
# -gp : Genome-wide (i.e., sampled for each species tree) gene-by-lineage-
#       specific parameter (to use with HG)
#  -s : Species tree (branch lengths in generations)
# -si : Number of individuals per species
# -sp : Tree-wide effective population size
# -su : Tree-wide substitution rate
# -sg : Tree-wide generation time
# -lb : Duplication rate (events/generation)
# -ld : Loss rate (events/generation)
# -lg : Gene conversion rate (events/generation)
# -hs : Species-specific branch rate heterogeneity modifiers
# -hl : Gene-family-specific branch rate heterogeneity modifiers
# -hh : Gene-by-lineage-specific locus tree parameter
# -hg : Gene-by-lineage-specific rate heterogeneity modifiers
#  -o : Common output prefix­name (for folders and files).
# -ot : Determines whether the species and locus tree branches are written in
#       number of generations (0) or time units (1)
# -om : Activates the tree mapping output
# -od : Activates the SQLite database output
# -op : Activates logging of sampled options
# -oc : Activates loggin of the original command line parameters, etc.
# -ol : Activates the output of trees with internal nodes labeled by its
#       post-order id starting from 0
#  -v : Verbosity
# -cs : Random number generator seed


# Notes:
# ------
# 1. Most parameters (-rs, -s, -sp, -su, -lb, -ld) are taken from Du, Hahn, and
#    Nakhleh (biorxiv 2019). Note that because the genome-wide parameters (-gb
#    and -gd) are not set, the duplication and loss rates are the same for each
#    replicate species tree.
# 2. Unlike Du, Hahn, & Nakhleh (biorxiv 2019), we did NOT enable gene
#    conversion; they set gene conversion (-lg) to 50X the mutation rate.
# 3. We did not enable species-specific (-hs) branch rate heterogeneity
#    modifiers. This option sets the alpha parameter for a Gamma distribution
#    with mean 1. Multipliers for each branch in the species tree are sampled
#    from this distribution. Note that if this option is a fixed number, then
#    alpha is the same across all replicate species trees, and otherwise,
#    alpha is sampled for each species tree. 
# 4. We did not enable gene-family-specific (-hl) branch rate heterogeneity
#    modifiers. This option sets the alpha parameter for a Gamma distribution
#    with mean 1. Multipliers for each locus tree are sampled from this
#    distribution. Note that if this option is a fixed number, then alpha is
#    the same across all replicate species trees, and otherwise, alpha is
#    sampled for each species tree.
# 5. We wanted gene trees to deviate from a molecular clock, so we did enable
#    we did enable gene-by-lineage-specific rate heterogeneity modifiers. This
#    requires considering three parameters (-gh, -hg, and -hh). The genome-wide
#    parameter (-gh) defines a distribution that is then sampled for each
#    replicate species tree. The locus tree parameter (-hh) defines a
#    distribution that is then sampled for each locus tree, and the branch rate
#    heterogeneity parameter (-hg) defines alpha for a Gamma distribution that
#    is then sampled to get a multipler for each branch in the gene tree. Here
#    we only set the -hg parameter to LN:1.5,1, which is the same value used by
#    Zhang et al., (BMC Bioinformatics 2018).
# 6. We used the random seed (-cs) from Mirarab et al., (Bioinformatics 2015).
# 7. For the simulation with 25, 50, and 100 taxa, we used the same species tree
#    height as the fungi species tree (1.8 x 10^9 generations). We found
#    (through visualization of 16-taxon species trees) that speciation rates of
#    1.8 x 10^-10 and 1.8 x 10^-8 corresponded to deep and recent speciation,
#    respectively. We decided to use the intermediate value of 1.8 x 10^-9.

base="."
data="$base/data"
sfwr="$base/software"

ntaxs=( 100 )
rates=( 0.0 0.0000000001 0.0000000002 0.0000000005 )  # Events per generation
sizes=( 10000000 50000000 )

for ntax in ${ntaxs[@]}; do
    for rate in ${rates[@]}; do
        for size in ${sizes[@]}; do
            ../simphy \
                -rs 10 \
                -rl F:1000 \
                -rg 1 \
                -sb F:0.0000000018 \
                -st F:1800000337.5 \
                -sl F:$ntax \
                -si F:1 \
                -sp F:$size \
                -su F:0.0000000004 \
                -lb F:$rate \
                -ld F:lb \
                -hg LN:1.5,1 \
                -o $data/ntaxa-$ntax.dlrate-$rate.psize-$size \
                -ot 0 \
                -om 1 \
                -od 1 \
                -op 1 \
                -oc 1 \
                -ol 1 \
                -v 3 \
                -cs 293745 &> $data/ntaxa-$ntax.dlrate-$rate.psize-$size.txt
        done
    done
done
