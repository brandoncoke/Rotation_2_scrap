if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
if("openPrimeR" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("openPrimeR")}
library(openPrimeR)
hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
file="J:\\gene.FASTA"
seqdf= read_templates(file, hdr.structure, delim = "|", id.column = "GROUP")
# Show the loaded metadata for the first template
c(seqdf$Accession[1], seqdf$Group[1], seqdf$Species[1], seqdf$Function[1])
template.df.uni <- assign_binding_regions(seqdf, fw = c(1,50), rev = c(1,40))
list.files(system.file("extdata", "settings", package = "openPrimeR"), pattern = "*\\.xml")
settings.xml <- system.file("extdata", "settings", 
                            "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
settings <- read_settings(settings.xml)
cvg_constraints(settings) <- list()
conOptions(settings)$allowed_mismatches <- 3
design.settings <- settings
constraints(design.settings) <- constraints(design.settings)[!grepl(
  "gc_clamp", names(constraints(design.settings)))]
constraints(design.settings)[["primer_length"]] <- c("min" = 25, "max" = 25)
out.file <- tempfile("settings", fileext = ".xml")
write_settings(design.settings, out.file)#
template.df <- assign_binding_regions(seqdf[1,], fw = file, rev = NULL)
optimal.primers <- design_primers(template.df[1:2,], mode.directionality = "fw",
                                  settings = design.settings)
out.file <- tempfile("my_primers", fileext = ".fasta")
write_primers(optimal.primers$opti, out.file)
# Define the FASTA primer file to load
primer.location <- system.file("extdata", "IMGT_data", "primers", "IGHV", 
                               "Ippolito2012.fasta", package = "openPrimeR")
# Load the primers
primer.df <- read_primers(primer.location, fw.id = "_fw")


