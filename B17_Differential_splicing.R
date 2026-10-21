library(leafcutter)
library(doMC)
library(dplyr)
library(GenomicRanges)

num_threads = 36
registerDoMC(num_threads)

# Input is counts.gz from jc_txt files
Counts=read.table(
  "/.../", 
  header=T, 
  check.names = F
)
library(stringr)

colnames(Counts) <- 
  str_replace_all(
    str_replace_all(colnames(Counts),"WNN_L15_", ""), 
    "_JC.txt", 
    ""
  )

write.table(
  data.frame(
    sample=colnames(Counts)
  ), 
  "./Sample_Groups_CZ.txt",
    sep="\t", 
    quote=FALSE, 
  row.names = FALSE  
)

meta=read.table(
  "./Sample_Groups_CZ.txt", 
  header=T, 
  stringsAsFactors = F
)

counts = Counts 

confounders=NULL
s=Inf 
timeout = 30
init = "smart"
seed = 1042 

i = 5 

min_samples_per_group = c(5)
min_coverage = c(20)

group_names=unique(meta$group) 
meta$group=factor(meta$group, group_names)
numeric_x=as.numeric(meta$group)-1
DS_Results_Dox_NoCov = list()

for(g in min_samples_per_group){
  for(c in min_coverage){
    
    print(Sys.time())
    message(
      paste0(
        'Working on Analysis: \n"', 
        g, 
        '" minimal samples per group, ', 
        'minimal reads per sample: ', 
        c
      )
    )
    results <- differential_splicing(
      counts, 
      numeric_x, 
      confounders=confounders, 
      max_cluster_size=s, 
      min_samples_per_intron=i, 
      min_samples_per_group=g, 
      min_coverage=c, 
      timeout=timeout, 
      init=init, 
      seed=seed
    ) 
    
    cluster_table = cluster_results_table(results)
    cluster_table$cluster = add_chr(cluster_table$cluster)
    
    DS_Results_Dox_NoCov[[
      paste0(
        "Min_", 
        g, 
        "Samples", 
        "_", 
        "Min_", 
        c, 
        "_Reads"
      )
    ]] <- list(
      "Min_Samples" = g,
      "Min_Reads" = c, 
      "Results" = results, 
      "Cluster_table" = cluster_table, 
      "N_Sign" <- table(
        cluster_table$p.adjust<0.05
      )["TRUE"]
    )
  }
}

saveRDS(
  DS_Results_Dox_NoCov, 
  "./Leafcutter_Differential_Splicing_Results.rds",
  compress=FALSE
)

rm(results,cluster_table,c,g,group_names,N_Sign, numeric_x)

library(dplyr)
library(leafcutter)
library(leafviz)
library(doMC)
num_threads = 36
registerDoMC(num_threads)

Results <- DS_Results_Dox_NoCov$Min_5Samples_Min_20_Reads

Sample_data <- meta


# Visualize statistics

table(Results$Cluster_table$status)

hist(
  Results$Cluster_table$p, 
  br=100
)

hist(
  Results$Cluster_table$p.adjust, 
  br=100
)

table(
  Results$Cluster_table$p.adjust < 0.05
)

# Prepare data for Leafviz visualization 

results = Results$Results
cluster_table =Results$Cluster_table
exon_file = "./gencode.v49.exons_generated_by_leafcutter.txt" 
counts = read.table(
  "./Leafcutter/leafcutter_perind_numers.counts.gz", 
  header=T, 
  check.names = F
) 

# Add gene names to clusters if an exon file is available
if (!is.null(exon_file)) {
  if (file.exists(exon_file)) {
    tryCatch( {
      exons_table     = read.table(exon_file, header=T, stringsAsFactors = F)
      intron_meta     = get_intron_meta(rownames(counts))
      exons_table$chr = add_chr(exons_table$chr)
      intron_meta$chr = add_chr(intron_meta$chr)
      clu_gene_map    = map_clusters_to_genes(intron_meta, exons_table)
      cluster_table   = merge(cluster_table, clu_gene_map, by.x="cluster", by.y="clu", all.x=TRUE)
    }, error=function(err) warning(as.character(err)) ) 
  } else warning("File ",opt$exon_file," does not exist")
} else cat("No exon_file provided.\n") 

# Write cluster table 
write.table(
  cluster_table, 
  "./DS_cluster_significance.txt", 
  quote=F, 
  sep="\t", 
  row.names = F
)

# Write effect size table
effect_size_table = leaf_cutter_effect_sizes(results)

meta <- Sample_data %>% 
  select(
    "sample"=sample, 
    "group"=group
  ) 
counts=Counts[,meta$sample]

group_names=unique(meta$group)
colnames(effect_size_table)[3:4] = group_names
effect_size_table$intron = add_chr(effect_size_table$intron) 

write.table(
  effect_size_table, 
  "./DS_effect_sizes.txt",  
  quote=F, 
  col.names = T, 
  row.names = F, 
  sep="\t"
) 

cluster_table_sign <- cluster_table %>% filter(p.adjust<0.05) 

write.table(
  cluster_table_sign, 
  "./WICHTIG_DS_Signifificant_Differential_Splicing_Events.txt", 
  sep="\t", 
  quote = FALSE
)

# new section (10.10.26)
# preprocessed text files
setwd("...")

clusters <- read.delim("DS_cluster_significance.txt", stringsAsFactors = FALSE)
effects  <- read.delim("DS_effect_sizes.txt", stringsAsFactors = FALSE)
sig      <- read.table("WICHTIG_DS_Signifificant_Differential_Splicing_Events.txt",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)

head(clusters); head(effects); head(sig)
View(sig)


effects <- read.delim("DS_effect_sizes.txt", stringsAsFactors = FALSE)

# split intron name chr:start:end:cluster
p <- do.call(rbind, strsplit(effects$intron, ":"))
effects$chr     <- p[, 1]
effects$start   <- as.numeric(p[, 2])
effects$end     <- as.numeric(p[, 3])
effects$cluster <- paste0(p[, 1], ":", p[, 4])

get_ce <- function(cl, min_dpsi = 0.05) {
  x <- effects[effects$cluster == cl, ]
  up <- x[x$deltapsi > min_dpsi, ]
  out <- list()
  for (i in seq_len(nrow(up))) for (j in seq_len(nrow(up))) {
    A <- up[i, ]; B <- up[j, ]
    if (A$end < B$start) {    
      out[[length(out) + 1]] <- data.frame(
        cluster = cl, chr = A$chr,
        ce_start = A$end, ce_end = B$start,
        ce_length = B$start - A$end,
        dpsi_A = A$deltapsi, dpsi_B = B$deltapsi,
        dpsi_canonical = if (nrow(canon)) canon$deltapsi[1] else NA)
    }
  }
  do.call(rbind, out)
}

# STMN2
x <- effects[effects$cluster == "chr8:clu_12980_-", ]
x[order(x$start), ]
colnames(effects)

x <- effects[effects$cluster == "chr8:clu_12980_-", ]
canon_end <- max(x$end)
ce_start_candidates <- subset(x, end < canon_end & deltapsi < -0.01)
unique(ce_start_candidates$end) 

library(GenomicRanges)
ex <- read.table("exon_file_STMN2.txt", col.names = c("chr","start","end","strand","gene"))
ce <- GRanges("chr8", IRanges(79616822, 79617049)) 
ex_gr <- GRanges(ex$chr, IRanges(ex$start, ex$end))
findOverlaps(ce, ex_gr)

ex <- read.table("./gencode.v49.exons_generated_by_leafcutter.txt",
                 header = TRUE, stringsAsFactors = FALSE)
colnames(ex)

ex_stmn2 <- ex[ex$gene_name == "STMN2", ]
ex_stmn2$chr <- ifelse(grepl("^chr", ex_stmn2$chr), ex_stmn2$chr, paste0("chr", ex_stmn2$chr))

ex_gr <- GRanges(ex_stmn2$chr, IRanges(ex_stmn2$start, ex_stmn2$end))

# test the junction-derived CE start
ce_start <- GRanges("chr8", IRanges(79616822, 79616822))
findOverlaps(ce_start, ex_gr)  

#------------------
ex <- read.table("./gencode.v49.exons_generated_by_leafcutter.txt",
                 header = TRUE, stringsAsFactors = FALSE)
ex$chr <- ifelse(grepl("^chr", ex$chr), ex$chr, paste0("chr", ex$chr))

ex_ends   <- unique(paste(ex$chr, ex$end))
ex_starts <- unique(paste(ex$chr, ex$start))

effects$donor_annot <- paste(effects$chr, effects$start) %in% ex_ends
effects$acc_annot   <- paste(effects$chr, effects$end)   %in% ex_starts

effects$type <- ifelse(effects$donor_annot & effects$acc_annot, "annotated_ends",
                       ifelse(!effects$donor_annot & !effects$acc_annot, "both_novel",
                              ifelse(!effects$acc_annot, "novel_acceptor", "novel_donor")))


# KALRN
ex$chr <- ifelse(grepl("^chr", ex$chr), ex$chr, paste0("chr", ex$chr))
ex_k <- ex[ex$gene_name == "KALRN", ]

cl_k <- clusters$cluster[which(grepl("KALRN", clusters$genes) & clusters$p.adjust < 0.05)]
cl_k

res <- lapply(cl_k, function(cl) {
  x <- effects[effects$cluster == cl, ]
  x$donor_annot <- x$start %in% ex_k$end
  x$acc_annot   <- x$end   %in% ex_k$start
  x$type <- ifelse(x$donor_annot & x$acc_annot, "annotated_ends",
                   ifelse(!x$donor_annot & !x$acc_annot, "both_novel",
                          ifelse(!x$acc_annot, "novel_acceptor", "novel_donor")))
  x[order(x$start), c("cluster", "intron", "X1", "X2", "deltapsi",
                      "donor_annot", "acc_annot", "type")]
})
res <- do.call(rbind, res)
nrow(res)
subset(res, type != "annotated_ends")


# MARK3
cl <- "chr14:clu_17839_-"   
grep("clu_17839", unique(effects$cluster), value = TRUE)

ex$chr <- ifelse(grepl("^chr", ex$chr), ex$chr, paste0("chr", ex$chr))
ex_m <- ex[ex$gene_name == "MARK3", ]

x <- effects[effects$cluster == cl, ]
nrow(x)   

x$donor_annot <- x$start %in% ex_m$end
x$acc_annot   <- x$end   %in% ex_m$start
x$type <- ifelse(x$donor_annot & x$acc_annot, "annotated_ends",
                 ifelse(!x$donor_annot & !x$acc_annot, "both_novel",
                        ifelse(!x$acc_annot, "novel_acceptor", "novel_donor")))

x <- x[order(x$start), c("intron", "X1", "X2", "deltapsi", "type")]
x

subset(x, type != "annotated_ends" & abs(deltapsi) > 0.01)

# TANC1
cl <- "chr2:clu_15426_-" 
grep("clu_15426", unique(effects$cluster), value = TRUE)

ex$chr <- ifelse(grepl("^chr", ex$chr), ex$chr, paste0("chr", ex$chr))
ex_m <- ex[ex$gene_name == "TANC1", ]

x <- effects[effects$cluster == cl, ]
nrow(x) 

x$donor_annot <- x$start %in% ex_m$end
x$acc_annot   <- x$end   %in% ex_m$start
x$type <- ifelse(x$donor_annot & x$acc_annot, "annotated_ends",
                 ifelse(!x$donor_annot & !x$acc_annot, "both_novel",
                        ifelse(!x$acc_annot, "novel_acceptor", "novel_donor")))

x <- x[order(x$start), c("intron", "X1", "X2", "deltapsi", "type")]
x

subset(x, type != "annotated_ends" & abs(deltapsi) > 0.01)
