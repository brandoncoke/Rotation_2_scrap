################################################################################
#Package install
################################################################################
if("BiocManager" %in% rownames(installed.packages()) == FALSE){
  install.packages("BiocManager")}
library(Biobase)
if("biomaRt" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("biomaRt")}
library(biomaRt)
################################################################################
#Sequence search
################################################################################
sequences_and_genes <- read.csv("~\\Rotation 2\\Raw data\\data for dagLogo.csv")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
IDs = as.character(sequences_and_genes$ENTREZID)
Sequences1 = getSequence(id=as.numeric(IDs[1:500]),
                        type="entrezgene_id",
                        seqType="peptide", 
                        mart=ensembl)
Sequences2 = getSequence(id=as.numeric(IDs[501:1000]),
                         type="entrezgene_id",
                         seqType="peptide", 
                         mart=ensembl)
Sequences3 = getSequence(id=as.numeric(IDs[1001:1500]),
                         type="entrezgene_id",
                         seqType="peptide", 
                         mart=ensembl)
Sequences4 = getSequence(id=as.numeric(IDs[1501:2000]),
                         type="entrezgene_id",
                         seqType="peptide", 
                         mart=ensembl)
Sequences5 = getSequence(id=as.numeric(IDs[2001:2500]),
                         type="entrezgene_id",
                         seqType="peptide", 
                         mart=ensembl)
Sequences6 = getSequence(id=as.numeric(IDs[2501:3000]),
                         type="entrezgene_id",
                         seqType="peptide", 
                         mart=ensembl)
Sequences7 = getSequence(id=as.numeric(IDs[3001:3500]),
                         type="entrezgene_id",
                         seqType="peptide", 
                         mart=ensembl)
Sequences8 = getSequence(id=as.numeric(IDs[3501:4000]),
                         type="entrezgene_id",
                         seqType="peptide", 
                         mart=ensembl)
Sequences9 = getSequence(id=as.numeric(IDs[4001:4500]),
                         type="entrezgene_id",
                         seqType="peptide", 
                         mart=ensembl)
Sequences10 = getSequence(id=as.numeric(IDs[4501:4868]),
                         type="entrezgene_id",
                         seqType="peptide", 
                         mart=ensembl)

Sequences=rbind(Sequences1,Sequences2,Sequences3,Sequences4,
            Sequences5,Sequences6,Sequences7,Sequences8,Sequences9,Sequences10)

save.image("~/Sequences of interest.RData") #Saving progress so far into working directory.
rm(Sequences1,Sequences2,Sequences3,Sequences4,
    Sequences5,Sequences6,Sequences7,Sequences8,Sequences9,Sequences10)



Sequences$Contains_S_CDK_motif=FALSE
Sequences$Contains_S_CDK_motif[grepl("R.L",Sequences$peptide)]= #As | aint working in regular expression
  TRUE
Sequences$Contains_S_CDK_motif[grepl("K.L",Sequences$peptide)]=
  TRUE
Sequences$Contains_S_CDK_motif[grepl("R.F",Sequences$peptide)]= 
  TRUE
Sequences$Contains_S_CDK_motif[grepl("K.F",Sequences$peptide)]=
  TRUE
Sequences$Contains_G2_CDK_motif=FALSE
Sequences$Contains_G2_CDK_motif[grepl("P..P.F",Sequences$peptide)]= #As | aint working in regular expression
  TRUE
Sequences$Contains_M_CDK_motif=FALSE
Sequences$Contains_M_CDK_motif[grepl("E.L.F",Sequences$peptide)]= #As | aint working in regular expression
  TRUE
Sequences$Contains_common_CDK_motif=FALSE
Sequences$Contains_common_CDK_motif[grepl("R.L.F",Sequences$peptide)]= #As | aint working in regular expression
  TRUE
Sequences$Contains_common_CDK_motif[grepl("K.L.F",Sequences$peptide)]= #As | aint working in regular expression
  TRUE
piedata=data.frame(entrezgene_id=unique(Sequences$entrezgene_id),Contains_common_CDK_motif=F,
                   Contains_G2_CDK_motif=F,
                   Contains_S_CDK_motif=F,
                   Contains_M_CDK_motif=F)
c=1
cycles=nrow(piedata) + 1
while(c<cycles){
  tempframe= subset(Sequences,entrezgene_id ==piedata$entrezgene_id[c])
  if(any(tempframe$Contains_common_CDK_motif)){
    piedata$Contains_common_CDK_motif[c]= TRUE
  }
  if(any(tempframe$Contains_S_CDK_motif)){
    piedata$Contains_S_CDK_motif[c]= TRUE
  }
  if(any(tempframe$Contains_G2_CDK_motif)){
    piedata$Contains_G2_CDK_motif[c]= TRUE
  }
  if(any(tempframe$Contains_M_CDK_motif)){
    piedata$Contains_M_CDK_motif[c]= TRUE
  }
  c=c+1
}

#Ugly tho. Venn diagram without tidyverse.
library(VennDiagram)
categories=c("Contains common CDK motif","Contains S phase specific CDK motif",
             "Contains G2 phase specific CDK motif",
             "Contains M phase specific CDK motif")
a1=subset(piedata,Contains_common_CDK_motif == TRUE)
a2=subset(piedata,Contains_S_CDK_motif == TRUE)
a3=subset(piedata,Contains_G2_CDK_motif == TRUE)
a4=subset(piedata,Contains_M_CDK_motif == TRUE)

area1= nrow(a1)
area2= nrow(a2)
area3= nrow(a3)
area4= nrow(a4)

n12= nrow(subset(a1,Contains_S_CDK_motif == TRUE))

n13= nrow(subset(a1,Contains_G2_CDK_motif == TRUE))

n14= nrow(subset(a1,Contains_M_CDK_motif == TRUE))

n23= nrow(subset(a2,Contains_G2_CDK_motif == TRUE))

n24= nrow(subset(a2,Contains_M_CDK_motif == TRUE))

n34= nrow(subset(a3,Contains_M_CDK_motif == TRUE))


n123= nrow(subset(a1,Contains_S_CDK_motif == TRUE &
                    Contains_G2_CDK_motif == TRUE))

n124= nrow(subset(a1,Contains_S_CDK_motif == TRUE &
                    Contains_M_CDK_motif == TRUE))

n134= nrow(subset(a1,Contains_G2_CDK_motif == TRUE &
                    Contains_M_CDK_motif == TRUE))

n234= nrow(subset(a2,Contains_G2_CDK_motif == TRUE &
                    Contains_M_CDK_motif == TRUE))

n1234= nrow(subset(a1,Contains_S_CDK_motif == TRUE &
                     Contains_G2_CDK_motif == TRUE &
                     Contains_M_CDK_motif == TRUE))           

venn.plot <- draw.quad.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  area4 = area4,
  n12 = n12,
  n13 = n13,
  n14 = n14,
  n23 = n23,
  n24 = n24,
  n34 = n34,
  n123 = n123,
  n124 = n124,
  n134 = n134,
  n234 = n234,
  n1234 = n1234,
  category = categories, #Defined above. Needs to be a list.
  fill = c("red", "blue","green","yellow"), #Colours of the areas.
  cex=4, #Changes size of the venn text.
  lty = 1, #Creates dashes with 2
  cat.cex = 0.2,
  fontfamily = "sans", #Defaults to serif. Could use Times, AvantGarde.
  cat.fontfamily = "sans",
  print.mode = "raw",
  sigdigs = 3) #Defaults to "raw"


categories=c("Contains common CDK motif",
             "Contains G2 phase specific CDK motif",
             "Contains M phase specific CDK motif")

a1=subset(piedata,Contains_common_CDK_motif == TRUE)
a2=subset(piedata,Contains_G2_CDK_motif == TRUE)
a3=subset(piedata,Contains_M_CDK_motif == TRUE)

area1= nrow(a1)
area2= nrow(a2)
area3= nrow(a3)

n12= nrow(subset(a1,Contains_G2_CDK_motif == TRUE))

n13= nrow(subset(a1,Contains_M_CDK_motif == TRUE))

n23= nrow(subset(a2,Contains_M_CDK_motif == TRUE))


n123= nrow(subset(a1,Contains_G2_CDK_motif == TRUE &
                    Contains_M_CDK_motif == TRUE))


venn.plot <- draw.triple.venn(area1 = area1,
                              area2 = area2,
                              area3 = area3,
                              n12= n12,
                              n13 = n13,
                              n23= n23,
                              n123=n123,
                              category = categories, #Defined above. Needs to be a list.
                              fill = c("red", "blue","green"), #Colours of the areas.
                              cex=4, #Changes size of the venn text.
                              lty = 1, #Creates dashes with 2
                              cat.cex = 0.2,
                              fontfamily = "sans", #Defaults to serif. Could use Times, AvantGarde.
                              cat.fontfamily = "sans",
                              print.mode = "raw",
                              sigdigs = 3) #Defaults to "raw"

