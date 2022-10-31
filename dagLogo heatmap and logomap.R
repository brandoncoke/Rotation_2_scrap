################################################################################
#Package install
################################################################################
if("BiocManager" %in% rownames(installed.packages()) == FALSE){
  install.packages("BiocManager")}
library(Biobase)
if("dagLogo" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("dagLogo")}
library(dagLogo)
if("biomaRt" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("biomaRt")}
library(biomaRt)
if("UniProt.ws" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("UniProt.ws")}
library(UniProt.ws)

################################################################################
#Package install
################################################################################
Humangenes = useDataset(mart = useMart("ensembl"),
                        dataset = "hsapiens_gene_ensembl")
sequences_and_genes <- read.csv("~\\Rotation 2\\Raw data\\data for dagLogo.csv")
anchorPos=as.character(sequences_and_genes$Sequences)
IDs = as.character(sequences_and_genes$ENTREZID)
mart <- useMart("ensembl")
range=20
humanmart <-
  useDataset(mart = mart, dataset = "hsapiens_gene_ensembl")
seq <- fetchSequence(IDs = IDs,
                     anchorAA = "*",
                     anchorPos = anchorPos,
                     mart = humanmart,
                     upstreamOffset = range,
                     downstreamOffset = range)
head(seq@peptides)

proteome <- dagLogo::prepareProteome("UniProt", species = "Homo sapiens")

bg_fisher <- buildBackgroundModel(seq, background = "wholeProteome", 
                                  proteome = proteome, testType = "fisher")
bg_ztest <- buildBackgroundModel(seq, background = "wholeProteome", 
                                 proteome = proteome, testType = "ztest")

t0 <- testDAU(seq, dagBackground = bg_ztest)

dagHeatmap(t0,type="statistics") 
dagLogo(testDAUresults= t0,
        #,title="Addjacent residues to phosphorylation sites",type="diff",
        fontface=1,fontsize = 10) 

#Plk1 consensus L-hy-D/E/N(Q)-X-S/T*-L-hy hy=hydrophobic minimum (D/E)X(S/T)
#Or (E/D)-X-(S/T)*-(I/L/V/M)-X-(E)
#CK1 (S/T)*-X-X-S/T
#Cdk2 
#Auroa B consensus motif (K/R)X(S/T)(I/L/V)
#Auroa A R-R-X(S/T)*.
#CK1-alpha D/E-X-X-S/T
#Cdk1 optimal S/T-P-X-R/K or S/T-P as minimum
#Cdk1 non canonical- S-A/G-X-R/K, S-X-R/K and P-X-S-X-(R/K)-R/K-(R)
#Core Nek hy-x-x-S/T-[non-P], acidic. Generally S over T exception Nek10
#Nek10 goes for S/T-hy very high preference
#Nek1 & 4 K/R 2+ & R -1
#Nek 5 and 8 K/R 2+ less about R -1
#Nek 6,7,9 acid -5 and -4 partially -2. Y -1 also improtant!  similar to Plk1
################################################################################
#Startifying based on localisatiion mitochondrial first
################################################################################
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
go=c("GO:0005739") #Mitochondrial proteins
mitochondrial_genes= getBM(attributes= "hgnc_symbol",filters=c("go"),values=list(go), mart=ensembl)
mitochondrial_genes=mitochondrial_genes$hgnc_symbol
sequences_and_genes$Is_mitochondrial= FALSE
sequences_and_genes$Is_mitochondrial[sequences_and_genes$Gene_name %in%
                                       mitochondrial_genes] = TRUE
sequences_and_genes_mitochondrial_only= 
  sequences_and_genes[sequences_and_genes$Is_mitochondrial == TRUE,]

anchorPos=as.character(sequences_and_genes_mitochondrial_only$Sequences)
IDs = as.character(sequences_and_genes_mitochondrial_only$ENTREZID[samples])
mart <- useMart("ensembl")
humanmart <-
  useDataset(mart = mart, dataset = "hsapiens_gene_ensembl")
seq <- fetchSequence(IDs = IDs,
                     anchorAA = "*",
                     anchorPos = anchorPos,
                     mart = humanmart,
                     upstreamOffset = range,
                     downstreamOffset = range)
head(seq@peptides)

proteome <- dagLogo::prepareProteome("UniProt", species = "Homo sapiens")

bg_fisher <- buildBackgroundModel(seq, background = "wholeProteome", 
                                  proteome = proteome, testType = "fisher")
bg_ztest <- buildBackgroundModel(seq, background = "wholeProteome", 
                                 proteome = proteome, testType = "ztest")

t0 <- testDAU(seq, dagBackground = bg_ztest)

dagHeatmap(t0) 
dagLogo(t0) 

################################################################################
#Proteins associated with ribosome and assessing the phorylation sites
################################################################################
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
go="GO:0005840" #ribosomal proteins
ribosomal_genes= getBM(attributes= "hgnc_symbol",filters=c("go"),values=list(go), mart=ensembl)
ribosomal_genes=ribosomal_genes$hgnc_symbol
sequences_and_genes$Is_ribosomal= FALSE
sequences_and_genes$Is_ribosomal[sequences_and_genes$Gene_name %in%
                                       ribosomal_genes] = TRUE
sequences_and_genes_ribosomal_only= 
  sequences_and_genes[sequences_and_genes$Is_ribosomal == TRUE,]

anchorPos=as.character(sequences_and_genes_ribosomal_only$Sequences)
IDs = as.character(sequences_and_genes_ribosomal_only$ENTREZID)
seq <- fetchSequence(IDs = IDs,
                     anchorAA = "*",
                     anchorPos = anchorPos,
                     mart = humanmart,
                     upstreamOffset = range,
                     downstreamOffset = range)
head(seq@peptides)

bg_fisher <- buildBackgroundModel(seq, background = "wholeProteome", 
                                  proteome = proteome, testType = "fisher")
bg_ztest <- buildBackgroundModel(seq, background = "wholeProteome", 
                                 proteome = proteome, testType = "ztest")

t0 <- testDAU(seq, dagBackground = bg_ztest)

dagHeatmap(t0) 
dagLogo(t0) 

################################################################################
#Proteins associated with cytoskeleton and assessing the phorylation sites
################################################################################
go="GO:0005856" #cytoskeletal GO term
cytoskeletal_genes= getBM(attributes= "hgnc_symbol",filters=c("go"),values=list(go), mart=ensembl)
cytoskeletal_genes=cytoskeletal_genes$hgnc_symbol
sequences_and_genes$Is_cytoskeletal= FALSE
sequences_and_genes$Is_cytoskeletal[sequences_and_genes$Gene_name %in%
                                       cytoskeletal_genes] = TRUE
sequences_and_genes_cytoskeletal_only= 
  sequences_and_genes[sequences_and_genes$Is_cytoskeletal == TRUE,]
set.seed(0451)
samples= sample(nrow(sequences_and_genes),1000) # Need to do it in 1000 chunks, server too feeble

anchorPos=as.character(sequences_and_genes_cytoskeletal_only$Sequences)
IDs = as.character(sequences_and_genes_cytoskeletal_only$ENTREZID)
seq <- fetchSequence(IDs = IDs,
                     anchorAA = "*",
                     anchorPos = anchorPos,
                     mart = humanmart,
                     upstreamOffset = range,
                     downstreamOffset = range)
head(seq@peptides)

bg_fisher <- buildBackgroundModel(seq, background = "wholeProteome", 
                                  proteome = proteome, testType = "fisher")
bg_ztest <- buildBackgroundModel(seq, background = "wholeProteome", 
                                 proteome = proteome, testType = "ztest")

t0 <- testDAU(seq, dagBackground = bg_ztest)

dagHeatmap(t0) 
dagLogo(t0)

################################################################################
#Proteins associated with cell cycle and assessing the phorylation sites
################################################################################
go="GO:0007049" #Cell cycle term
cytoskeletal_genes= getBM(attributes= "hgnc_symbol",filters=c("go"),values=list(go), mart=ensembl)
cytoskeletal_genes=cytoskeletal_genes$hgnc_symbol
sequences_and_genes$Is_cytoskeletal= FALSE
sequences_and_genes$Is_cytoskeletal[sequences_and_genes$Gene_name %in%
                                      cytoskeletal_genes] = TRUE
sequences_and_genes_cytoskeletal_only= 
  sequences_and_genes[sequences_and_genes$Is_cytoskeletal == TRUE,]
set.seed(0451)
samples= sample(nrow(sequences_and_genes),1000) # Need to do it in 1000 chunks, server too feeble

anchorPos=as.character(sequences_and_genes_cytoskeletal_only$Sequences)
IDs = as.character(sequences_and_genes_cytoskeletal_only$ENTREZID)
seq <- fetchSequence(IDs = IDs,
                     anchorAA = "*",
                     anchorPos = anchorPos,
                     mart = humanmart,
                     upstreamOffset = range,
                     downstreamOffset = range)
head(seq@peptides)

bg_fisher <- buildBackgroundModel(seq, background = "wholeProteome", 
                                  proteome = proteome, testType = "fisher")
bg_ztest <- buildBackgroundModel(seq, background = "wholeProteome", 
                                 proteome = proteome, testType = "ztest")

t0 <- testDAU(seq, dagBackground = bg_ztest)

dagHeatmap(t0) 
dagLogo(t0)









