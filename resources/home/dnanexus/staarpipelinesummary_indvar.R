args <- commandArgs(TRUE)
### mandatory
outfile <- args[1]
test.type <- args[2] # Gene_Centric_Coding, Gene_Centric_Noncoding, ncRNA, Sliding_Window
nullobj.file <- args[3]
chromosome <- as.numeric(args[4])
agds.file <- args[5]
### optional
annotation.name.catalog.file <- args[6]
known.varlist.4columns <- args[7]
gene.name <- args[8]
category <- args[9]
start.loc <- as.numeric(args[10])
end.loc <- as.numeric(args[11])
max.maf <- as.numeric(args[12])
QC_label <- args[13]
variant_type <- args[14]
geno_missing_imputation <- args[15]
Annotation_dir <- args[16]
Annotation_name <- args[17]

test.type.vals <- c("Gene_Centric_Coding", "Gene_Centric_Noncoding", "ncRNA", "Sliding_Window")
if(!test.type %in% test.type.vals) stop("Error: test.type must be Gene_Centric_Coding, Gene_Centric_Noncoding, ncRNA, or Sliding_Window")

cat("Output file prefix:", outfile, "\n")
cat("Type of test:", test.type, "\n")
cat("Null model object:", nullobj.file, "\n")
cat("Chromosome:", chromosome, "\n")
cat("Genotype and functional annotation all-in-one GDS (AGDS) file:", agds.file, "\n")
cat("Name and the corresponding channel name in the AGDS file:", annotation.name.catalog.file, "\n")
cat("Known variant list (4 columns):", known.varlist.4columns, "\n")
cat("Gene name of the gene category:", gene.name, "\n")
cat("Functional category of the gene category:", category, "\n")
cat("Start location of the genetic region:", start.loc, "\n")
cat("End location of the genetic region:", end.loc, "\n")
cat("Maximum minor allele frequency to be included for variant-set test:", max.maf, "\n")
cat("Channel name of the QC label in the AGDS file:", QC_label, "\n")
cat("Variants included in the analysis:", variant_type, "\n")
cat("Method of handling missing genotypes:", geno_missing_imputation, "\n")
cat("Channel name of the annotations in the AGDS file:", Annotation_dir, "\n")
cat("Annotations used in STAAR:", Annotation_name, "\n")

suppressMessages(library(gdsfmt))
suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(STAAR))
suppressMessages(library(SCANG))
suppressMessages(library(STAARpipeline))
suppressMessages(library(STAARpipelineSummary))

## Null Model
nullobj <- get(load(nullobj.file))
if(class(nullobj) == "GENESIS.nullMixedModel") {
  nullmod$sample.id <- row.names(nullmod$model.matrix)
  nullobj <- genesis2staar_nullmodel(nullmod)
  rm(nullmod); gc()
}

if(annotation.name.catalog.file == "NO_ANNOTATION_NAME_CATALOG_FILE") stop("Error: Annotation name catalog file cannot be missing")
Annotation_name_catalog <- read.csv(annotation.name.catalog.file, as.is=T)
Annotation_name <- unlist(strsplit(Annotation_name,","))

if(known.varlist.4columns != "NO_KNOWN_VARLIST_4COLUMNS") {
  known.varlist.4columns <- read.csv(known.varlist.4columns, as.is=T)
  known.varlist.4columns_info <- known.varlist.4columns[,c("CHR","POS","REF","ALT")]
  known.varlist.4columns_info <- known.varlist.4columns_info[known.varlist.4columns_info$CHR == chromosome,,drop=FALSE]
} else {
  known.varlist.4columns_info <- data.frame(CHR=character(),POS=character(),REF=character(),ALT=character(),stringsAsFactors = FALSE)
}

known_loci_genome <- known.varlist.4columns_info
known_loci_genome <- known_loci_genome[order(known_loci_genome$CHR, known_loci_genome$POS),]


## Main function for this applet
if(test.type == "Gene_Centric_Coding") {
  genofile <- seqOpen(agds.file)
  out <- Gene_Centric_Coding_Info(category=category,chr=chromosome,genofile=genofile,obj_nullmodel=nullobj,
                                  gene_name=gene.name,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                  Annotation_name=Annotation_name)
  seqClose(genofile)
}else if(test.type == "Gene_Centric_Noncoding") {
  genofile <- seqOpen(agds.file)
  out <- Gene_Centric_Noncoding_Info(category=category,chr=chromosome,genofile=genofile,obj_nullmodel=nullobj,
                                     gene_name=gene.name,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                     QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                     Annotation_name=Annotation_name)
  seqClose(genofile)
}else if(test.type == "ncRNA") {
  genofile <- seqOpen(agds.file)
  out <- Gene_Centric_Noncoding_Info(category="ncRNA",chr=chromosome,genofile=genofile,obj_nullmodel=nullobj,
                                     gene_name=gene.name,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                                     QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                     Annotation_name=Annotation_name)
  seqClose(genofile)
}else if(test.type == "Sliding_Window") {
  genofile <- seqOpen(agds.file)
  out <- Sliding_Window_Info(chr=chromosome,genofile=genofile,obj_nullmodel=nullobj,
                             start_loc=start.loc,end_loc=end.loc,known_loci=known_loci_genome,rare_maf_cutoff=max.maf,
                             QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                             Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                             Annotation_name=Annotation_name)
  seqClose(genofile)
}

rm(list=setdiff(ls(), c("out", "outfile"))); gc()
write.csv(out, file = paste0(outfile, ".csv"), row.names = FALSE)

