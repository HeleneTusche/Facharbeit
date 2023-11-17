#Generates read tables and heatmaps for read lists
#Berlin enterovirus
#Requires the following packages: dplyr, tibble, tidyr, purrr and stringr from tidyverse, as well as pheatmap

library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(stringr)
library(pheatmap)


the_wd <- "/fast/AG_Landthaler/emanuel/wastewater/Berlin_enterovirus/joint_mapping/"
paired_filelist <- "/fast/AG_Landthaler/emanuel/wastewater/Berlin_enterovirus/joint_mapping/RNA_filelist.txt"
single_filelist <- ""

accs_names_list <- "/fast/AG_Landthaler/emanuel/wastewater/Berlin_enterovirus/acc_taxid_name.txt"
readcountvector <- "/fast/AG_Landthaler/emanuel/wastewater/ww_RNA_all/readcountvector.txt"

#We do two filters: keep accessions that are on average higher than min_average, and/or in min_nos number of samples higher than min_sample
#For min_nos, something like 2-4 seems reasonable, and for min_sample something like 3-6 times the value of min_average
min_average <- 0.1
min_nos <- 3
min_sample <- 0.5
#After the initial filtering, we remove those accessions, that have more than overlap_threshold overlap with accessions that are higher in the globalrank
#The idea being - if the overlap with a more abundant accession number is high, then likely the reads rather belong to something closer to the more abundant accession
#Still, do not set overlap_threshold too low, in order not to loose to many accessions at this point
overlap_threshold = 50

samples_to_remove <- c("mix_1", "mix_2", "BIOM", "PCRinh", "captu", "cer_")

prefix <- "/fast/AG_Landthaler/emanuel/wastewater/Berlin_enterovirus/Berlin_enterovirus"

reads_suffix <- "_joint_sequences_reads.txt"

filtered_list_file <- "/fast/AG_Landthaler/emanuel/wastewater/Berlin_enterovirus/filtered_accs_list.txt"

liste = c(1,2,3,4,5)
for (i in liste){
  print(paste0("i ist = ", i))
  j = 2 + i
  print(paste0("j ist = ", j))
}

setwd(the_wd)
paired_filelist <- try(read.table(paired_filelist, header=FALSE)$V1, silent = TRUE)
if (inherits(paired_filelist, "try-error")) {
  paired_filelist <- c()
  print("paired end file list not found/not specified")
}
single_filelist <- try(read.table(single_filelist, header=FALSE)$V1, silent = TRUE)
if (inherits(single_filelist, "try-error")) {
  single_filelist <- c()
  print("single end file list not found/not specified")
}
acc_names <- read.table(accs_names_list, header=FALSE, sep="\t")
colnames(acc_names) <- c("acc", "taxid", "name")


#Prepare acc_names
acc_names$taxid <- NULL
acc_names$name <- paste(acc_names$acc, acc_names$name, sep=" ")

#Do the counting
#Note: this code assumes that there are no multi-mappers in the data!
all_files_list <- list()
expr <- list()
for (the_file in paired_filelist) {
  #Try reading the file. If the file is empty, and error message is shown
  reads <- try(read.table(paste0(the_file, reads_suffix), header=FALSE), silent = TRUE)
  if (inherits(reads, "try-error")) {
    reads <- data.frame(matrix(nrow = 0, ncol = 3)) 
    colnames(reads) <- c("read", "acc", "sample")
    print(paste0("file ", the_file, " not found or file is empty"))
  }
  else {
    reads$sample <- the_file
    colnames(reads) <- c("read", "acc", "sample")
  }
  
  #If a read identifiers occurs twice for the at least one accession number, then likely both read 1 and read 2 have mapped
  #Be aware that this could also be a multimapping of either read 1 or read 2, which would look the same here
  #Also, if a read identifier occurs three times, this is not taken into account
  paired_reads <- reads %>% group_by(read, acc) %>% tally() %>% filter(n==2) %>% pull(var="read") %>% unique()
  single_reads <- setdiff(unique(reads$read), paired_reads)
  paired_reads <- subset(reads, read %in% paired_reads)
  single_reads <- subset(reads, read %in% single_reads)

  #Add the paired_ranks information to the paired_reads information. This gives the rank of accession numbers within the current sample.
  #So, the accession numbers with most paired reads in this sample is on rank 1, etc.
  #For a given read, only those accs participate in the ranking that have the read identifier twice, i.e. for which both reads mapped
  paired_ranks <- paired_reads %>% group_by(read, acc) %>% tally() %>% filter(n==2) %>% ungroup() %>%
    group_by(acc) %>% tally() %>% dplyr::arrange(desc(n)) %>% rownames_to_column(var="order") %>% select(-n)
  paired_ranks$order <- as.numeric(paired_ranks$order)
  paired_reads <- merge(paired_reads, paired_ranks, by="acc")
  #Every reads is counted only once, for the accession number with the best rank within the sample
  paired_counts <- paired_reads %>% group_by(read) %>% slice_min(order_by = order, n=1, with_ties = FALSE) %>% group_by(acc) %>% tally(name = "count") %>% dplyr::arrange(desc(count))
  #Because it's paired-end and both reads mapped, we count them twice
  paired_counts$count = paired_counts$count*2
  
  #Now the same for single reads
  single_ranks <- single_reads %>% group_by(acc) %>% tally() %>% dplyr::arrange(desc(n)) %>% rownames_to_column(var="order") %>% select(-n)
  single_ranks$order <- as.numeric(single_ranks$order)
  single_reads <- merge(single_reads, single_ranks, by="acc")
  single_counts <- single_reads %>% group_by(read) %>% slice_min(order_by = order, n=1, with_ties = FALSE) %>% group_by(acc) %>% tally(name = "count") %>% dplyr::arrange(desc(count))
  expr[[the_file]] <- rbind.data.frame(paired_counts, single_counts) %>% group_by(acc) %>% summarize(count=sum(count))
  colnames(expr[[the_file]]) <- c("acc", the_file)
  class(expr[[the_file]]$acc) <- "character"
  
  all_files_list[[the_file]] <- reads
}

for (the_file in single_filelist) {
  reads <- try(read.table(paste0(the_file, reads_suffix), header=FALSE), silent = TRUE)
  if (inherits(reads, "try-error")) {
    reads <- data.frame(matrix(nrow = 0, ncol = 3)) 
    colnames(reads) <- c("read", "acc", "sample")
    print(paste0("file ", the_file, " not found or file is empty"))
  }
  else {
    reads$sample <- the_file
    colnames(reads) <- c("read", "acc", "sample")
  }
  single_reads <- reads
  single_ranks <- single_reads %>% group_by(acc) %>% tally() %>% dplyr::arrange(desc(n)) %>% rownames_to_column(var="order") %>% select(-n)
  single_ranks$order <- as.numeric(single_ranks$order)
  single_reads <- merge(single_reads, single_ranks, by="acc")
  single_counts <- single_reads %>% group_by(read) %>% slice_min(order_by = order, n=1, with_ties = FALSE) %>% group_by(acc) %>% tally(name = "count") %>% dplyr::arrange(desc(count))
  expr[[the_file]] <- single_counts %>% group_by(acc) %>% summarize(count=sum(count))
  colnames(expr[[the_file]]) <- c("acc", the_file)
  class(expr[[the_file]]$acc) <- "character"
  all_files_list[[the_file]] <- reads
}

#Merge into a big table with accessions as row names and samples as column names
mapping_table <- purrr::reduce(expr, full_join, by="acc") %>% replace(., is.na(.), 0) %>% column_to_rownames(var="acc")
#merge all read files into a big data frame
all_files <- as.data.frame(do.call(rbind, all_files_list))
colnames(all_files) <- c("read", "acc", "sample")

#Read in the number of reads per sample for normalization
reads <- read.table(readcountvector, header=FALSE, sep="\t", row.names=NULL)
colnames(reads) <- c("sample", "reads")
mapping_table <- as.data.frame(t(mapping_table)) %>% rownames_to_column(var="sample")
mapping_table <- left_join(reads, mapping_table, by="sample") %>% column_to_rownames(var="sample") %>% t() %>% as.data.frame()

#At this stage, filter out samples that you do not want to include in the analysis. For the Berlin case, all with non-standard processing (i.e. not Centricons)
mapping_table <- mapping_table %>% select(-contains(samples_to_remove)) %>% t() %>% as.data.frame()
all_files <- all_files %>% filter(!(str_detect(sample, paste(samples_to_remove, collapse="|"))))

#globalrank is the ranking of accession numbers by total number of raw counts in the entire dataset
globalrank <- all_files %>% group_by(acc) %>% tally() %>% dplyr::arrange(desc(n)) %>% rownames_to_column(var="order") %>% select(-n)
colnames(globalrank) <- c("order", "acc")

#Normalize counts with read counts
mapping_table_norm <- mapping_table/mapping_table$reads*1000000
mapping_table_norm$reads <- NULL

#Get number of samples that we have at this point
samples <- dim(mapping_table_norm)[[1]]

#List of accs to work with at this point
acc_list <- as.data.frame(t(mapping_table_norm)) %>%
  mutate(AT = rowSums(across(everything(), ~ . > min_sample))) %>%
  mutate(RS = rowSums(across(-contains("AT")))) %>%
  filter(AT>=min_nos | RS >= min_average*samples) %>% rownames()

#Now, for every acc, calculate the percentage of reads that also map to all other accs
#The values in the overlaps data frame mean: of all reads mapping to the accession in column acc, the value in the overlap column indicates how
#many percent of these also map to the accession in column acc2
all_files <- subset(all_files, acc %in% acc_list)
expr <- list()
acc_reads_list <- list()
for (the_acc in acc_list) {
  acc_reads_list[[the_acc]] <- subset(all_files, acc == the_acc) %>% pull(var="read")
}
for (the_acc in acc_list) {
  for (the_acc_2 in acc_list) {
    name_nor <- length(acc_reads_list[[the_acc]])
    overlap <- length(acc_reads_list[[the_acc]][acc_reads_list[[the_acc]] %in% acc_reads_list[[the_acc_2]]])
    overlap=overlap/name_nor*100
    expr[[paste(the_acc, the_acc_2, sep="_")]] <- data.frame(the_acc, the_acc_2, overlap)
    colnames(expr[[paste(the_acc, the_acc_2, sep="_")]]) <- c("acc", "acc2", "overlap")
  }
}
overlaps <- as.data.frame(do.call(rbind, expr))

#add the global ranking to the overlaps table
overlaps <- left_join(overlaps, globalrank, by="acc")
colnames(globalrank) <- c("order2", "acc2")
overlaps <- left_join(overlaps, globalrank, by="acc2")
overlaps$order <- as.numeric(overlaps$order)
overlaps$order2 <- as.numeric(overlaps$order2)

#After the initial filtering, we remove those accessions, that have more than overlap_threshold overlap with accessions that are higher in the globalrank
#The idea being - if the overlap with a more abundant accession number is high, then likely the reads rather belong to something closer to the more abundant accession
#Still, do not set overlap_threshold too low, in order not to loose to many accessions at this point
accs_to_remove <- overlaps %>% filter(acc != acc2 & overlap > overlap_threshold & order2<order) %>% pull(acc) %>% unique()
accs_filt <- setdiff(acc_list, accs_to_remove)

mapping_table_norm_all <- as.data.frame(t(mapping_table_norm)) %>% rownames_to_column(var="acc") %>%
  filter(acc %in% acc_list) %>% left_join(., acc_names, by="acc") %>% column_to_rownames(var="name") %>% mutate(acc=NULL)
mapping_table_norm_filt <- as.data.frame(t(mapping_table_norm)) %>% rownames_to_column(var="acc") %>%
  filter(acc %in% accs_filt) %>% left_join(., acc_names, by="acc") %>% column_to_rownames(var="name") %>% mutate(acc=NULL)

the_filename <- paste(prefix, "mapstat_samples_all.pdf", sep="_")
pheatmap::pheatmap(as.matrix(log10(mapping_table_norm_all+0.1)),  cluster_cols = FALSE, border_color = NA, treeheight_row = 4, treeheight_col = 4, cellwidth = 9, cellheight = 8, fontsize_row=7, fontsize_col = 6, filename = the_filename)

the_filename <- paste(prefix, "mapstat_samples_filt.pdf", sep="_")
pheatmap::pheatmap(as.matrix(log10(mapping_table_norm_filt+0.1)),  cluster_cols = FALSE, border_color = NA, treeheight_row = 4, treeheight_col = 4, cellwidth = 9, cellheight = 8, fontsize_row=7, fontsize_col = 6, filename = the_filename)

#Now with aggregation (dataset-specific code)
mapping_table_aggr <- mapping_table %>% mutate(month = gsub(pattern = "^WW_([0-9]{4}).*", replacement = "\\1", rownames(mapping_table), perl = TRUE))
mapping_table_aggr <- mapping_table_aggr %>% group_by(month) %>% summarize_all(sum) %>% column_to_rownames(var="month")
mapping_table_aggr_norm <- mapping_table_aggr/mapping_table_aggr$reads *1000000
mapping_table_aggr_norm$reads <- NULL

mapping_table_aggr_norm_all <- as.data.frame(t(mapping_table_aggr_norm)) %>% rownames_to_column(var="acc") %>%
  filter(acc %in% acc_list) %>% left_join(., acc_names, by="acc") %>% column_to_rownames(var="name") %>% mutate(acc=NULL) %>% select(-c("2107", "2202"))
mapping_table_aggr_norm_filt <- as.data.frame(t(mapping_table_aggr_norm)) %>% rownames_to_column(var="acc") %>%
  filter(acc %in% accs_filt) %>% left_join(., acc_names, by="acc") %>% column_to_rownames(var="name") %>% mutate(acc=NULL) %>% select(-c("2107", "2202"))

the_filename <- paste(prefix, "mapstat_aggr_all.pdf", sep="_")
pheatmap::pheatmap(as.matrix(log10(mapping_table_aggr_norm_all+0.1)),  cluster_cols = FALSE, border_color = NA, treeheight_row = 4, treeheight_col = 4, cellwidth = 9, cellheight = 8, fontsize_row=7, fontsize_col = 6, filename = the_filename)

the_filename <- paste(prefix, "mapstat_aggr_filt.pdf", sep="_")
pheatmap::pheatmap(as.matrix(log10(mapping_table_aggr_norm_filt+0.1)),  cluster_cols = FALSE, border_color = NA, treeheight_row = 4, treeheight_col = 4, cellwidth = 9, cellheight = 8, fontsize_row=7, fontsize_col = 6, filename = the_filename)

write.table(as.data.frame(accs_filt), filtered_list_file, row.names = FALSE, col.names = FALSE, quote=FALSE)

