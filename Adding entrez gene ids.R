################################################################################
#Data prep
################################################################################
library(readxl)

library(org.Hs.eg.db)
hs <- org.Hs.eg.db
my.symbols <- unique(Sequences$hgnc_symbol)
ENTREZIDs= select(hs,
       keys = my.symbols,
       columns = c('UNIPROT', 'SYMBOL', 'ENTREZID'),
       keytype = "SYMBOL") 
ENTREZIDs= na.omit(ENTREZIDs)
AllPhosposites_no$ENTREZID= NA
i=1
cycle=length(ENTREZIDs$SYMBOL) +1
while(i < cycle){
  AllPhosposites_no$ENTREZID[which(
    ENTREZIDs$SYMBOL[i] == AllPhosposites_no$Gene.names)] = ENTREZIDs$ENTREZID[i]
  i=i+1
  }
AllPhosposites_no= subset(AllPhosposites_no,!is.na(ENTREZID))
AllPhosposites_no= AllPhosposites_no[,c(4,5,7,49)]
astring=AllPhosposites_no$Phospho..STY..Probabilities
astring= gsub("0.9","Z",astring)
astring= gsub("0.8","Z",astring)
astring= gsub("0.7","Z",astring)
astring= gsub("0.6","Z",astring)
astring= gsub("0.5","Z",astring)
astring= gsub('[[:punct:] ]+','',astring)
astring= gsub("S1","SZ",astring)
astring= gsub("T1","TZ",astring)
astring= gsub('[[:digit:]]+', '', astring)
library(stringr)
astring_counts= str_count(astring, "Z")
astring=gsub("(Z)\\1{1,}", "\\1", astring)
#astring= gsub("Z","*",astring)
astring= gsub("Z","",astring) #COMMENT ME OUT FOR DAG LOGO PREP!!!
data_for_dagLogo=data.frame(Gene_name= AllPhosposites_no$Gene.names, 
                          #ENTREZID= AllPhosposites_no$ENTREZID,
                          Sequences= astring,
                          string_counts= astring_counts,
                          Motif_start= AllPhosposites_no$Positions.within.proteins)
#data_for_dagLogo=subset(data_for_dagLogo,string_counts==1)
write.csv(data_for_dagLogo,"K:\\data for dagLogo.csv")