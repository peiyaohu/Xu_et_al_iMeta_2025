##########################
#  Purpose: filter phyloseq obj.
#  Date: 20240703
#  By: hupeiyao

# Parameter:
# @ps: raw ps to be filterd
# @prev_N: The threshold for an OTU to occur in at least prev_N samples.
# @outdir: the directory where the `discard_log` will be written.

# process:
# 1) Prevalence: less than or equal to 20
# 2) Annotated as non-bacteria: `kingdom != "Bacteria"` or `Family == "Mitochondria"`

# output:
# 1) ps.filter: filterd phyloseq obj
# 2) discard_log: a dataframe containing reads_statistics



##########################
ps_filter<- function(ps, prev_n= 20, outdir="Table/"){
  dir.create("Table")
  # 1) prevalence
  index1<- rowSums(otu_table(ps)>0)>prev_n
  ps1<- prune_taxa(!index1, ps)
  
  # 2) kingdom != Bacteria
  index2<- data.frame(tax_table(ps))$Kingdom=="Bacteria"
  ps2<- prune_taxa(!index2, ps)
  
  # 3) Family == "Mitochondria
  idx1<- grep("Mito",data.frame(tax_table(ps))$Family)
  idx2<- rownames(data.frame(tax_table(ps))[idx1,])
  index3<- ! rownames(otu_table(ps)) %in% idx2
  ps3<- prune_taxa(!index3, ps)
  
  # DISCARD:
  ps.filter<- prune_taxa(index1, ps)
  ps.filter<- prune_taxa(data.frame(tax_table(ps.filter))$Kingdom=="Bacteria", ps.filter)
  ps.filter<- prune_taxa(!data.frame(tax_table(ps.filter))$Family %in% "Mitochondria", ps.filter)

  # process_log
  discard_log<- data.frame(
    Raw= sample_sums(ps),
    LowPrevalence=sample_sums(ps1) ,
    Kingdom=sample_sums(ps2),
    Mitochondria=sample_sums(ps3),
    Retained=sample_sums(ps.filter),
    percentage= round(sample_sums(ps.filter) / sample_sums(ps) *100,2)) 
  cat("Note: The integers in the table represent read counts. The `LowPrevalence`, `Kingdom`,and `Mitochondria` statistics indicate the number of reads in the original input that do not meet the specified conditions (these categories may overlap).")
  cat(paste("Discard_log was written to", outdir))
  write.csv(discard_log, paste0(outdir,"Table0.ps.filter_log.csv"))
  
  # output
  return(list(ps.filter = ps.filter, 
              discard_log = discard_log))
}

