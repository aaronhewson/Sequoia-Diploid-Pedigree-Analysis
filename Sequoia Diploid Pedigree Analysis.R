# Pedigree analysis using Sequoia for Jim Dunckley Heritage Orchard and Plant & Food Research samples SNP genotyped in 2025 for Aaron Hewson's MSc project.
# The majority of samples were identified as known cultivars by duplicate matching, this is an attempt to identify parents of remaining samples.

# Prior to running this script, all samples and reference accessions were combined into one PLINK .ped/.map file pair, with ACTG allele codes, 
# using a list from Nick Howard's to identify the high-quality SNPs common across the apple 20K, 50K, and 480K SNP array platforms.


# Load Packages -----------------------------------------------------------

#install.packages("sequoia")

library(sequoia)


# Set Working Directory ---------------------------------------------------

setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/Sequoia Diploid Pedigree Analysis")


# Load Input Data for Sequoia ---------------------------------------------
# All diacritical marks and apostrophes were removed from sample names before loading!

#Load genotypes (PLINK .ped file)
genoIn <- as.matrix(read.table("Inputs/All_Samples.ped", header = FALSE))
geno <- GenoConvert(InData = genoIn, InFormat = "ped", Missing = "0", quiet = FALSE)

#Load life history (Sex and BirthYear given as "NA" = unknown)

#Load pedigree (Mum and Dad given as "NA" = unknown)
