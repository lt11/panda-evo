## header ---------------------------------------------------------------------

options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(data.table)
library(this.path)
library(scriptName)
library(GenomicRanges)
library(tictoc)

## function(s) ----------------------------------------------------------------

#' Split and Extract a Specific Sub-Column from a Character Column
#'
#' This function splits a character column using a specified delimiter 
#' and extracts the desired sub-column based on the given index.
#'
#' @param x A character vector, typically a column from a data table.
#' @param n An integer specifying the index of the sub-column 
#'          to extract after splitting.
#' @param s A character string representing the delimiter used for splitting.
#'
#' @return A character vector containing the extracted sub-column values.
#'
#' @examples
#' # example usage:
#' vec <- c("A|B|C", "D|E|F", "G|H|I")
#' # extract the second sub-column
#' SplitSubCol(vec, 2, "|")  
#'
#' @export
SplitSubCol <- function(x, n, s) {
  y <- sapply(strsplit(x, split = s, fixed = TRUE), "[[", n)
  return(y)
}

## settings -------------------------------------------------------------------

### size threshold: smaller overlaps are removed
sizeLim <- 25
### fixed settings
dirBase <- dirname(this.dir())
### dev dirBase <- "/Users/Lorenzo/dev/panda"
dirAnnoBed <- file.path(dirBase, "anno", "bed")
dirOut <- file.path(dirBase, "png")
unlink(dirOut, recursive = T)
dir.create(dirOut, showWarnings = F)
pathGen <- file.path(dirOut, "generators.txt")
cat("class:feature#strand", "\n", file = pathGen)
idRef <- "SGDref"
hdImpg <- c("Query_id", "Query_start", "Query_end",
            "Target_id", "Target_start", "Target_end",
            "Target_clsfeat", "N_score", "Query_alndir", "Target_alndir")

### the path to the paf file
dirAlnPaf <- file.path(dirBase, "aln")
pathAlnPaf <- list.files(path = dirAlnPaf, pattern = "paf$",
                         recursive = T, full.names = T)
if (length(pathAlnPaf) > 1) {
  cat("[", myName, "] ",
      "User error.\n",
      sep = "")
  ### stop prints "Error: " by default
  stop("multiple paf files found.")
} else if (length(pathAlnPaf) == 0) {
  cat("[", myName, "] ",
      "User error.\n",
      sep = "")
  ### stop prints "Error: " by default
  stop("no paf file found.")
}

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Making the pangenome. ",
    "\n", sep = "")

### read strain-haplotypes ids from file (the order is maintained in the output)
pathIds <- file.path(dirBase, "ids", "ids-ps.txt")
vtStrainHaplo <- as.character(fread(file = pathIds, header = F)[[1]])
vtStrainHaploHash <- gsub(pattern = "-", replacement = "#", x = vtStrainHaplo)
### initialise the output features table
nHaplos <- length(vtStrainHaploHash)
dtPanFeats <- data.table(matrix(character(length = 0), ncol = nHaplos))
setnames(dtPanFeats, vtStrainHaploHash)
### dev 
# vtStrainHaploHash <- c("SGDref#0", "AFI#0", "S288C#0",
#                        "XXX-h3", "stoca-zo", "DBVPG6765#0")

## impg run -------------------------------------------------------------------

### message
cat("[", myName, "] ",
    "Running impg. ",
    "\n", sep = "")
### dev indS <- vtStrainHaplo[1] 
indL <- 1
lsImpg <- list()
### for each strain-haplotypes (hyphen-separated)
for (indS in vtStrainHaplo) {
  
  ### pick the strain-haplotypes bed file (query)
  pathAnnoBed <- list.files(path = dirAnnoBed, pattern = indS, full.names = T)
  
  ### run "$ impg -p file.paf -b file.bed"  which is 
  ### (50x faster than going line-by-line)
  ### with direct load of the output
  
  ### string for bash (impg-0.2.0)
  strBashImpg <- paste0("impg -I -p ", pathAlnPaf, " -b ", pathAnnoBed)
  # ### string for bash (impg-0.2.3)
  # strBashImpg <- paste0("impg query -I -p ", pathAlnPaf, " -b ", pathAnnoBed)
  
  strOut <- system(strBashImpg, intern = T)
  ### make a single-file data-table
  lsImpgOne <- strsplit(strOut, split = "\t")
  dtImpgOne <- rbindlist(lapply(lsImpgOne,
                                function(x) as.data.table(as.list(x))))
  ### append the single-file data-table to a list
  lsImpg[[indL]] <- dtImpgOne
  sizeList <- format(object.size(lsImpg), units = "auto")
  cat("[", myName, "] ",
      "Size of impg's list: ", sizeList, " at iteration: ",
      indL, " out of ", length(vtStrainHaplo), ".",
      "\n", sep = "")
  indL <- indL + 1
}
cat("[", myName, "] ",
    "Running impg: done.",
    "\n", sep = "")

### get the redundant data-table
dtImpgAll <- rbindlist(lsImpg)
colnames(dtImpgAll) <- hdImpg
### add supplementary columns
dtImpgAll[, c("Query_start", "Query_end") := lapply(.SD, as.numeric),
          .SDcols = c("Query_start", "Query_end")]
dtImpgAll[, Target_cls := SplitSubCol(Target_clsfeat, 1, ":")]
dtImpgAll[, c("Target_start", "Target_end") := lapply(.SD, as.numeric),
          .SDcols = c("Target_start", "Target_end")]
dtImpgAll[, Target_len := Target_end - Target_start]
dtImpgAll[, Query_len := Query_end - Query_start]

## overlap size filter --------------------------------------------------------

### size filter
dtImpgAllSzFlt <- dtImpgAll[Target_len > sizeLim
                            & Query_len > sizeLim]
cat("[", myName, "] ",
    "Overlap size filter. ",
    "\n", sep = "")
nDtImpgAll <- nrow(dtImpgAll)
nDtImpgAllSzFlt <- nrow(dtImpgAllSzFlt)
nSmall <- nDtImpgAll - nDtImpgAllSzFlt
cat("[", myName, "] ",
    "Removing ", nSmall,
    " entries out of ", nDtImpgAll, ".",
    "\n", sep = "")

### garbage collection
rm(dtImpgAll)
invisible(gc())

## sorting and formatting -----------------------------------------------------

### sorting ascending or descending is not equivalent since the 
### overlap is done in two rounds and 
### is strand-aware (if two features have overlapping 
### coordinates but one is in the + strand and the other is in 
### the - strand they do not overlap), thus
### we sort in increasing (ascending) order to avoid spurious
### overlaps that may occur when multiple classes are analysed
### (a spurious overlap will produce a big block that cannot be easily
### decomposed, except for the different classes);
### the YFL066C locus (where we have the Y' region TEL06L) is a good
### example of a spurious overlap
setorder(dtImpgAllSzFlt, Target_len)

### transfer the strand data from Target_clsfeat to Query_id and Target_id,
### so that when we calculate the blocks these will be strand-aware
strStrand <- SplitSubCol(dtImpgAllSzFlt[, Target_clsfeat], 2, "#")
dtImpgAllSzFlt[, Query_id := paste0(Query_id, "#", strStrand)]
dtImpgAllSzFlt[, Target_id := paste0(Target_id, "#", strStrand)]

## blocks calculation ---------------------------------------------------------

cat("[", myName, "] ",
    "Blocks calculation. ",
    "\n", sep = "")
nBlocks <- 1
vtUnq <- unique(dtImpgAllSzFlt[, Target_clsfeat])
### dev
# indTarClsFeat <- "CDS:YAL003W#+"
# indTarClsFeat <- vtUnq[1]

### YAL005C overlaps with YAL004W, we do not want these two genes
### to generate a single sub-block, so we make the overlap strand-aware
# indTarClsFeat <- "gene:YAL004W#+"
# indTarClsFeat <- "gene:YAL005C#-"
for (indTarClsFeat in vtUnq) {
  
  ### dev print(nBlocks)
  
  ### get the target class-feature that will generate the block
  dtTarClsFeat <- dtImpgAllSzFlt[indTarClsFeat, on = "Target_clsfeat"]
  ### check if the corresponding interval is still in dtImpgAllSzFlt
  ### (if it was ovelapping with a former indTarClsFeat
  ### it has already been removed, but dtTarClsFeat will still have a line
  ### e.g. the following):
  ### Query_id Query_start Query_end Target_id Target_start Target_end
  ###   <char>       <num>     <num>    <char>        <num>      <num>
  ###     <NA>         NA        NA      <NA>           NA         NA
  ### Target_clsfeat N_score Query_alndir Target_alndir Target_cls Target_len
  ###         <char>  <char>       <char>        <char>     <char>      <num>
  ###    buciodeculo    <NA>         <NA>          <NA>       <NA>         NA
  
  if (any(is.na(dtTarClsFeat[, Query_id]))) {
    next
  }
  
  ### first round of overlap (on query coordinates)
  setkey(dtTarClsFeat, Query_id, Query_start, Query_end)
  dtIndOverOne <- foverlaps(dtImpgAllSzFlt, dtTarClsFeat,
                            nomatch = NULL, type = "any", which = T)
  dtBlock <- dtImpgAllSzFlt[dtIndOverOne[, xid], ]
  ### dev print(format(object.size(dtBlock), units = "auto"))
  
  ### second round of overlap (on target coordinates)
  setkey(dtBlock, Target_id, Target_start, Target_end)
  dtIndOverTwo <- foverlaps(dtImpgAllSzFlt, dtBlock, 
                            nomatch = NULL, type = "any", which = T)
  dtBlock <- rbindlist(list(dtBlock, dtImpgAllSzFlt[dtIndOverTwo[, xid], ]))
  ### dev print(format(object.size(dtBlock), units = "auto"))
  
  ### block classes
  allClassesInBlock <- unique(dtBlock[, Target_cls])
  
  ### dev
  # cat(indTarClsFeat, "\t",
  #     nBlocks, "\t", length(allClassesInBlock), "\n",
  #     file = "~/Desktop/check-nsb.txt", append = T)
  
  ## sub-blocks reduction (class-by-class, strand-aware) ----------------------
  
  ### dev
  # indSb <- "gene"
  # indSb <- "TY1/TY2_soloLTR"
  # indSb <- "Y_prime_element"
  for (indSb in allClassesInBlock) {
    ### reduce (same as bedtools merge) the rows for each class
    ### on the basis of the query coordinates
    grSblock <- GRanges(seqnames = dtBlock[Target_cls == indSb, Query_id],
                        ranges = IRanges(start = dtBlock[Target_cls == indSb,
                                                         Query_start],
                                         end = dtBlock[Target_cls == indSb,
                                                       Query_end]),
                        mcols = dtBlock[Target_cls == indSb, Target_clsfeat])
    grSblockRed <- reduce(grSblock)
    
    ### vector with the features names and class name
    ### cleaned of the "#+" or "#-" at the end of the string
    ### e.g. "gene:S288C_G0022800" and "gene:YGR296W" 
    vtClsFeat <- sub(pattern = "#[\\+\\-]$",
                     replacement = "",
                     x = unique(grSblock$mcols))
    ### features string
    strFeats <- paste(gsub(pattern = paste0(indSb, ":"), replacement = "",
                           x = vtClsFeat),
                      collapse = ",")
    
    dtSblockRed <- data.table(as.character(grSblockRed@seqnames),
                              as.character(grSblockRed@ranges))
    
    ## transformation of the reduced sub-block --------------------------------
    
    ### e.g. this dtSblockRed:
    ###                 V1              V2
    ###             <char>          <char>
    ###  SGDref#0#chrVII#+ 1084864-1090591
    ###   S288C#0#chrVII#+ 1085116-1090843
    ###
    ### becomes a two-column data-table with columns:
    ###                 SGDref#0                     S288C#0
    ###                   <char>                      <char>
    ### chrVII:1084864-1090591#+    chrVII:1085116-1090843#+
    
    ### transposition and collapsing explained with
    ### another example, this dtSblockRed:
    ###                      V1            V2
    ###                  <char>        <char>
    ###     DBVPG6765#0#chrIV#- 960288-966619
    ###         AFI#0#chrXIII#- 349539-349680
    ###         AFI#0#chrXIII#- 349713-349989
    ###   DBVPG6765#0#chrXIII#- 354903-355179
    ###           AFI#0#chrIV#- 960661-961002
    ###       S288C#0#chrXIII#- 384508-384784
    ###      SGDref#0#chrXIII#- 378732-379008
    ###
    ### must produce this column:
    ### DBVPG6765#0
    ### chrIV:960288-966619;chrXIII:354903-355179
    
    ### formatting dtSblockRed
    dtSblockRed[, Haplo_id := sub("#chr.*", "", dtSblockRed[, V1])]
    dtSblockRed[, Info_str := paste(SplitSubCol(x = dtSblockRed[, V1],
                                                n = 3, s = "#"),
                                    dtSblockRed[, V2],
                                    sep = ":")]
    ### collapsing with ";" all the Haplo_id elements of a Haplo_id
    dtSblockRedCo <- dtSblockRed[, paste(Info_str, collapse = ";"),
                                 by = Haplo_id]
    ### transpose, producing a data-table
    dtTra <- transpose(dtSblockRedCo)
    ### set column names and format
    setnames(dtTra, as.character(dtTra[1, ]))
    dtTra <- dtTra[-1]
    ### add Class_id and Features_id columns
    dtTra[, ':='(Class_id = rep(indSb, .N), Features_id = rep(strFeats, .N))]
    
    ### add the missing columns: not needed 
    ### since rbindlist makes it by default using fill = TRUE
    # newCols <- setdiff(vtStrainHaploHash, colnames(dtTra))
    # dtTra[, (newCols) := lapply(newCols, function(x) NA)]
    # setcolorder(dtTra, vtStrainHaploHash)
    
    ### append with rbindlist matching (default operation) column names
    dtPanFeats <- rbindlist(list(dtPanFeats, dtTra), fill = T)
    
    ### append the generator feature to the output
    cat(rep(indTarClsFeat, nrow(dtTra)), sep = "\n",
        append = T, file = pathGen)
    ### dev print(nrow(dtTra))
  }
  
  ### delete dtIndOverOne and dtIndOverTwo rows in dtImpgAllSzFlt,
  ### the (any(is.na(dtTarClsFeat[, Query_id]))) will avoid 
  ### the loop to break
  indOut <- unique(c(dtIndOverOne[, xid], dtIndOverTwo[, xid]))
  dtImpgAllSzFlt <- dtImpgAllSzFlt[-indOut]
  
  nBlocks <- nBlocks + 1
  
  ### stop prints "Error: " by default
  ### dev if (nBlocks == 20) stop("we did 20 iterations!")
}
cat("[", myName, "] ",
    "Blocks calculation left ",
    nrow(dtImpgAllSzFlt),
    " impg lines. ",
    "\n", sep = "")

cat("[", myName, "] ",
    "Blocks calculation: done. ",
    "\n", sep = "")

### format and write dtPanFeats
leftCols <- c("Class_id", "Features_id")
colOrder <- c(leftCols, setdiff(names(dtPanFeats), leftCols))
setcolorder(dtPanFeats, colOrder)
pathOutPanFeat <- file.path(dirOut, "pan-features.txt")
write.table(x = dtPanFeats, file = pathOutPanFeat, append = F, quote = F,
            sep = "\t", col.names = T, row.names = F, na = "MA")
save(dtPanFeats, file = file.path(dirOut, "pan-features.RData"))
