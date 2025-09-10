## header ---------------------------------------------------------------------

options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(data.table)
library(this.path)
library(scriptName)

## settings -------------------------------------------------------------------

### fixed settings
dirBase <- dirname(this.dir())
dirAnnoGff <- file.path(dirBase, "anno", "gff")
dirOut <- file.path(dirBase, "anno", "bed")
unlink(dirOut, recursive = T)
dir.create(dirOut, showWarnings = F, recursive = T)
idRef <- "SGDref"
hdGff <- c("Chr_id", "Strain_id", "Feat_type", "S_coord",
           "E_coord", "S_val", "Strand_id", "Frame_id", "Attribute_str")
vtClassSrt <- c("CDS")
# vtClassSrt <- c("gene",
#                 "pseudogene",
#                 "intron",
#                 "five_prime_UTR_intron",
#                 "noncoding_exon",
#                 "plus_1_translational_frameshift",
#                 "blocked_reading_frame",
#                 "external_transcribed_spacer_region",
#                 "internal_transcribed_spacer_region",
#                 "ARS",
#                 "ARS_consensus_sequence",
#                 "TY1",
#                 "TY1/TY2_soloLTR",
#                 "TY1_truncated",
#                 "TY2",
#                 "TY3_soloLTR",
#                 "TY4_soloLTR",
#                 "TSU4_soloLTR",
#                 "TY4_truncated",
#                 "TY5",
#                 "TY5_soloLTR",
#                 "LTR_retrotransposon",
#                 "long_terminal_repeat",
#                 "W_region",
#                 "Z1_region",
#                 "Z2_region",
#                 "centromere",
#                 "centromere_DNA_Element_I",
#                 "centromere_DNA_Element_II",
#                 "centromere_DNA_Element_III",
#                 "matrix_attachment_site",
#                 "X_element",
#                 "X_element_partial",
#                 "X_element_combinatorial_repeat",
#                 "X_region",
#                 "Y_prime_element",
#                 "Y_region",
#                 "recombination_enhancer",
#                 "silent_mating_type_cassette_array",
#                 "mating_type_region",
#                 "tRNA",
#                 "pseudogenic_transcript",
#                 "ncRNA",
#                 "ncRNA_gene",
#                 "non_transcribed_region",
#                 "rRNA",
#                 "rRNA_gene",
#                 "snRNA",
#                 "snRNA_gene",
#                 "snoRNA",
#                 "snoRNA_gene",
#                 "intein_encoding_region")
### just in case there's some redundancy
vtClassSrt <- unique(vtClassSrt)

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Transforming the gff to bed. ",
    "\n", sep = "")

### read strain-haplotypes ids from file
pathIds <- list.files(path = file.path(dirBase, "ids"), pattern = "ids-ps.txt",
                      full.names = T)
vtStrainHaplo <- as.character(fread(file = pathIds, header = F)[[1]])
vtRef <- grep(idRef, vtStrainHaplo, value = T)
### remove the reference
vtStrainHaplo <- grep(idRef, vtStrainHaplo, value = T, invert = T)

### dev indR <- vtRef[1]
### loop for reference annotations
for (indR in vtRef) {
  ### read the gff
  pathAnnoGff <- list.files(path = dirAnnoGff, pattern = indR,
                            full.names = T, recursive = T)
  dtGff <- fread(file = pathAnnoGff, sep = "\t", header = F, verbose = F)
  colnames(dtGff) <- hdGff
  ### id
  strIdPref <- paste0(sub(pattern = "-", replacement = "#", x = indR), "#")
  ### convert the start coordinate from 1-based (gff) to 0-based (bed)
  dtGff[, S_coord := S_coord - 1]
  ### check start and end coordinates: if start = end impg breaks
  nBad <- nrow(dtGff[S_coord >= E_coord])
  if (nBad > 0) {
    cat("[", myName, "] ",
        "Found ", nBad,
        " annotations with the start coordinate larger or equal to",
        " the end coordinate that will be removed.\n",
        sep = "")
    dtGff <- dtGff[S_coord < E_coord]
  }
  ### filter features (for CDS)
  dtGff <- dtGff[Feat_type %in% vtClassSrt, ]
  ### make the feature id column
  strFeatId <- sub("^.*Name=([^;]*).*$", "\\1", dtGff$Attribute_str)
  ### trim trailing redundancy, e.g. _CDS in YAL003W_CDS
  strFeatIdTrm <- sub("^([^_]*).*$", "\\1", strFeatId)
  ### if strFeatId = dtGff[, Feat_type] set strFeatId = "MN"
  indM <- which(dtGff[, Feat_type] == strFeatIdTrm)
  if (length(indM) != 0) {
    strFeatIdTrm[indM] <- "MN"
  }
  ### paste class and feature id (and the strand)
  strClassIdStrand <- paste0(dtGff[, Feat_type], ":",
                             strFeatIdTrm, "#", dtGff[, Strand_id])
  ### transform the gff into a bed file
  dtBed <- data.table(Chrom_id = paste0(strIdPref, dtGff[, Chr_id]),
                      Chrom_start = dtGff[, S_coord],
                      Chrom_end = dtGff[, E_coord],
                      Class_feat = strClassIdStrand)
  ### write the bed file
  nameOut <- sub(pattern = "-features.gff$", replacement = ".bed",
                 x = basename(pathAnnoGff))
  pathOutBed <- file.path(dirOut, nameOut)
  fwrite(file = pathOutBed, x = dtBed, sep = "\t",
         quote = F, row.names = F, col.names = F)
}

### dev
### indS <- vtStrainHaplo[1]
### indS <- "AIE-0"
### loop for phenovar annotations
for (indS in vtStrainHaplo) {
  ### read the gff
  pathAnnoGff <- list.files(path = dirAnnoGff, pattern = indS,
                            full.names = T, recursive = T)
  dtGff <- fread(file = pathAnnoGff, sep = "\t", header = F)
  colnames(dtGff) <- hdGff
  ### id
  strIdPref <- paste0(sub(pattern = "-", replacement = "#", x = indS), "#")
  ### convert the start coordinate from 1-based (gff) to 0-based (bed)
  dtGff[, S_coord := S_coord - 1]
  ### check start and end coordinates: impg breaks if start = end
  nBad <- nrow(dtGff[S_coord >= E_coord])
  if (nBad > 0) {
    cat("[", myName, "] ",
        "Found ", nBad,
        " annotations with the start coordinate larger or equal to",
        " the end coordinate that will be removed.\n",
        sep = "")
    dtGff <- dtGff[S_coord < E_coord]
  }
  
  ### choose valid genes: those genes that do not have a systematic name
  ### since systematic genes will be included via the reference
  dtGffGenes <- dtGff[Feat_type %in% "gene", ]
  ### make the feature id column
  strGeneSysName <- sub("^.*Name=([^;]*).*$", "\\1", dtGffGenes$Attribute_str)
  ### nuclear gene names without systematic names:
  ### e.g. this is not valid: ID=DBVPG6765_G0000030;Name=YAL067C
  strValidGenes <- grep(pattern = "^Y[A-P][L,R][0-9]{3}[W,C]",
                        x = strGeneSysName, value = T, invert = T)
  
  ### filter features (for CDS)
  dtGff <- dtGff[Feat_type %in% vtClassSrt, ]
  ### make the feature id column
  strFeatId <- sub("^.*Name=([^;]*).*$", "\\1", dtGff$Attribute_str)
  ### trim trailing redundancy, e.g. .mRNA.1.CDS.2 in YAL003W.mRNA.1.CDS.2
  strFeatIdTrm <- sub("^([^.]+)\\..*$", "\\1", strFeatId)
  ### if strFeatId = dtGff[, Feat_type] set strFeatId = "MN"
  indM <- which(dtGff[, Feat_type] == strFeatIdTrm)
  if (length(indM) != 0) {
    strFeatIdTrm[indM] <- "MN"
  }
  ### paste class and feature id (and the strand)
  strClassIdStrand <- paste0(dtGff[, Feat_type], ":",
                             strFeatIdTrm, "#", dtGff[, Strand_id])
  ### transform the gff into a bed file
  ### and filter out genes with a systematic name
  dtBed <- data.table(chrom = paste0(strIdPref, dtGff[, Chr_id]),
                      Chrom_start = dtGff[, S_coord],
                      Chrom_end = dtGff[, E_coord],
                      Class_feat = strClassIdStrand,
                      Name_id = strFeatIdTrm)
  ### remove invalid genes, e.g. ID=DBVPG6765_G0000030;Name=YAL067C
  dtBed <- dtBed[Name_id %in% strValidGenes, ]
  ### remove the the column to filter for valid genes
  dtBed <- dtBed[, 1:4]
  ### write the bed file
  nameOut <- sub(pattern = "-features.gff$", replacement = ".bed",
                 x = basename(pathAnnoGff))
  pathOutBed <- file.path(dirOut, nameOut)
  fwrite(file = pathOutBed, x = dtBed, sep = "\t",
         quote = F, row.names = F, col.names = F)
}
