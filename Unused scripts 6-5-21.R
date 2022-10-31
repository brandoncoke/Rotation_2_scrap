get_longest_sequence<- function(Gene_identifiers){
  output_sequence_frame= data.frame(peptide=NA, hgnc_symbol=NA)
  c=1
  while(c <  length(Gene_identifiers)+1){
    Sequence= getSequence(id=as.character(Gene_identifiers[c]),
                type="hgnc_symbol",
                seqType="peptide", 
                mart=mart)
    longest_sequence= Sequence$peptide[max(length(Sequence$peptide))]
    longest_sequence=longest_sequence[1]
    temp_sequence_frame=data.frame(peptide= longest_sequence, hgnc_symbol= Gene_identifiers[c])
    output_sequence_frame= rbind(output_sequence_frame,temp_sequence_frame)
    c=c+1}
  output_sequence_frame}

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
IDs= NNdataframe$Substrate
IDs = as.character(unique(IDs))
#Unused obtaining longest sequence.
get_longest_sequence= function(gene_symbol){
  Sequence_dataframe= getSequence(id=as.character(gene_symbol),
                                  type="hgnc_symbol",
                                  seqType="peptide", 
                                  mart=mart)
  Sequences= Sequence_dataframe$peptide
  Longest_sequence= Sequences[max(length(Sequences))]
  Longest_sequence[1]}
Sequences= lapply(IDs,get_longest_sequence)
Sequences= as.character(Sequences)
Sequence_dataframe= data.frame(Gene_symbol=IDs[1:10] ,Sequences)
  
  
  