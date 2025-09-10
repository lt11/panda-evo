## header ---------------------------------------------------------------------

options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(data.table)
library(this.path)
library(scriptName)
library(stringr)
library(Biostrings)
library(BiocParallel)
library(gtools)

## function(s) ----------------------------------------------------------------

#' Prepend Column Names to Substrings in Selected Columns of a data.table
#'
#' This function modifies a `data.table` by prepending each value 
#' in the specified columns 
#' with the corresponding column name followed by `#`. 
#' Substrings within a column are assumed 
#' to be separated by semicolons (`;`). `NA` values remain unchanged.
#'
#' @param x A `data.table` containing the dataset to be modified.
#' @param y A numeric vector of column indices specifying 
#'          which columns should be modified.
#'
#' @return The input `data.table` is modified in-place, 
#' with the selected columns updated 
#' to have their column names prefixed to each substring.
#'
#' @examples
#' library(data.table)
#' 
#' # create example data.table
#' dt <- data.table(
#'   ID = 1:3,
#'   SGD#0 = c("chrIV:30-32", "chrXIV:386-541;chrXIV:260-415", NA),
#'   AAB#1 = c("", "chrV:35-80", "chrXII:99-111")
#' )
#' 
#' # define columns to modify 
#' # (assuming SGD#0 is in column 2 and AAB#1 in column 3)
#' indHapCols <- c(2, 3)
#'
#' # apply the function
#' AddColPref(dt, indHapCols)
#'
#' # expected output:
#' #     ID   SGD#0                                        AAB#1
#' # 1:  1    SGD#0#chrIV:30-32                            ""  (unchanged)
#' # 2:  2    SGD#0#chrXIV:386-541;SGD#0chrXIV:260-415     AAB#1#chrV:35-80
#' # 3:  3    NA (unchanged)                               AAB#1#chrXII:99-111
#'
#' @export
AddColPref <- function(x, y) {
  x[, (y) := lapply(y, function(indC) {
    col_name <- names(x)[indC]
    ifelse(!is.na(x[[indC]]),
           gsub("([^;]+)", paste0(col_name, "#\\1"), x[[indC]], perl = T),  
           x[[indC]])})]
  return(x)
}

#' Replace NA Values with Empty Strings in a data.table
#'
#' This function modifies a data.table replacing all NA values 
#' with empty strings ("").
#' You can choose whether to apply this only to character 
#' columns or to all columns
#' (which will coerce all types to character).
#'
#' @param x A `data.table` object to modify in-place.
#' @param y Logical. If `TRUE`, replaces NA in all columns 
#'          (and coerces to character).
#'          If `FALSE` (default), replaces only NA values
#'          in character columns.
#'
#' @return The modified `data.table` with NA values replaced by "".
#' @examples
#' library(data.table)
#' x <- data.table(A = c("apple", NA), B = c(NA, "orange"), C = c(1, NA))
#' ReplaceNAtoEmpty(x)         # only character columns
#' ReplaceNAtoEmpty(x, TRUE)   # all columns, coerced to character
ReplaceNAtoEmpty <- function(x, y = F) {
  if (!data.table::is.data.table(x)) {
    ### stop prints "Error: " by default
    stop("input must be a data.table.")
  }
  colsToMod <- if (y) names(x) else names(x)[sapply(x, is.character)]
  
  x[, (colsToMod) := lapply(.SD, function(x) {
    x[is.na(x)] <- ""
    x
  }), .SDcols = colsToMod]
  
  return(x)
}

#' Replace Empty Strings with NA Values in a data.table
#'
#' This function modifies a data.table by replacing all empty strings ("")
#' with NA values. You can choose whether to apply this only to character 
#' columns or to all columns (which will coerce all types to character).
#'
#' @param x A `data.table` object to modify in-place.
#' @param y Logical. If `TRUE`, replaces "" in all columns 
#'          (and coerces to character).
#'          If `FALSE` (default), replaces only "" values
#'          in character columns.
#'
#' @return The modified `data.table` with "" values replaced by NA.
#' @examples
#' library(data.table)
#' x <- data.table(A = c("apple", ""), B = c("", "orange"), C = c("1", ""))
#' ReplaceEmptyToNA(x)         # only character columns
#' ReplaceEmptyToNA(x, TRUE)   # all columns, coerced to character
ReplaceEmptyToNA <- function(x, y = F) {
  if (!data.table::is.data.table(x)) {
    ### stop prints "Error: " by default
    stop("input must be a data.table.")
  }
  colsToMod <- if (y) names(x) else names(x)[sapply(x, is.character)]
  
  x[, (colsToMod) := lapply(.SD, function(col) {
    col[col == ""] <- NA
    col
  }), .SDcols = colsToMod]
  
  return(x)
}

#' Read Multiple gzipped CDS FASTA Files into a data.table
#'
#' `ReadCds` reads all gzipped CDS FASTA files 
#' corresponding to the provided sample IDs,
#' extracting sequence names and sequences into a data.table fast lookup table.
#'
#' @param dir Character. Directory containing gzipped CDS FASTA files.
#' @param ids Character vector of sample identifiers 
#'            (e.g., "CMF#1", "SGDref#0").
#'
#' @return A data.table with columns: Source_id, Seq_id, 
#'         Sequence (single string).
#'         Keys: Source_id, Seq_id.
#' }
#'
#' @export
ReadCds <- function(dir, ids) {
  normIds <- gsub("#", "-", ids, fixed = T)
  dtAll <- data.table()
  
  for (genId in normIds) {
    oneFasta <- file.path(dir, paste0(genId, "-cds.fa.gz"))
    
    if (!file.exists(oneFasta)) {
      warning(sprintf("Missing: %s", oneFasta))
      next
    }
    
    ### wrap connection inside a local() block
    ### to ensure that it is closed once we have finished 
    ### working with that specific file (best practices)
    dtThis <- local({
      con <- gzfile(oneFasta, "r")
      on.exit(close(con))
      
      idSeq <- NULL
      seqLines <- character()
      records <- list()
      
      ### read one line at a time till EOF is found
      while (length(line <- readLines(con, n = 1, warn = F)) > 0) {
        if (startsWith(line, ">")) {
          if (!is.null(idSeq)) {
            records[[length(records) + 1]] <- list(
              Source_id = gsub("-", "#", genId, fixed = T),
              Seq_id = idSeq,
              Sequence = paste(seqLines, collapse = "")
            )
          }
          idSeq <- sub("^>", "", strsplit(line, "\\s+")[[1]][1])
          seqLines <- character()
        } else {
          seqLines <- c(seqLines, line)
        }
      }
      
      ### add the last entry of a fasta
      if (!is.null(idSeq)) {
        records[[length(records) + 1]] <- list(
          Source_id = gsub("-", "#", genId, fixed = T),
          Seq_id = idSeq,
          Sequence = paste(seqLines, collapse = "")
        )
      }
      ### append to dtThis
      rbindlist(records)
    })
    dtAll <- rbind(dtAll, dtThis)
  }
  return(dtAll)
}

#' Check if a DNA sequence is a valid protein-coding sequence
#'
#' This function verifies whether a given DNA sequence satisfies basic criteria
#' for being a valid open reading frame (ORF) in protein-coding genes.
#'
#' @param seq A character string representing a DNA sequence 
#'            (composed of A, T, G, and C).
#'
#' @return `TRUE` if the sequence:
#' \itemize{
#'   \item Has a length that is a multiple of 3 (i.e., full codons),
#'   \item Starts with a start codon (`ATG`),
#'   \item Ends with a stop codon (`TAA`, `TAG`, or `TGA`),
#'   \item Contains no internal stop codons (i.e., stop codons appear only 
#'         at the end).
#' }
#' Otherwise, it returns `FALSE`.
#'
#' @examples
#' IsValidCoding("ATGGACGACTAA")  # TRUE
#' IsValidCoding("ATGTAGGACTAA")  # FALSE due to internal stop codon
#'
#' @export
IsValidCoding <- function(seq) {
  strStartCod <- c("ATG")
  strStopCods <- c("TAA", "TAG", "TGA")
  if (nchar(seq) %% 3 != 0) return(F)
  codons <- str_sub(seq, seq(1, nchar(seq), 3), seq(3, nchar(seq), 3))
  if (!(codons[1] %in% strStartCod)) return(F)
  if (!(codons[length(codons)] %in% strStopCods)) return(F)
  if (any(codons[-length(codons)] %in% strStopCods)) return(F)
  return(T)
}

#' Reverse complement of a DNA sequence
#'
#' This function computes the reverse complement of a DNA string 
#' using the four standard bases.
#'
#' @param seq A character string representing a DNA sequence (e.g. "ATGC").
#' @param alphabet A named character vector defining the complement mapping.
#'                 Names are input bases, values are their complements.
#'                 Default is \code{c(A="T", T="A", G="C", C="G")}.
#'
#' @return A character string containing the reverse complement of \code{seq}.
#'         Or `NULL` if the input string is not valid.
#'
#' @examples
#' RevComp("CACCA")
#' # "TGGTG"
#'
#' RevComp("AATTGGCC")
#' # "GGCCAATT"
#'
#' # This will retutn NULL:
#' # RevComp("ATGN")
#'
#' @export
RevComp <- function(seq) {
  if (grepl("[^ATGC]", seq)) {
    message("Skipping invalid sequence.", "\n",
            seq, "\n")
    return(NULL)
  }
  alphabet = c(A="T", T="A", G="C", C="G")
  chars <- strsplit(seq, split = "")[[1]]
  comp <- alphabet[chars]
  rc <- paste(rev(comp), collapse = "")
  return(rc)
}

#' Append valid coding sequences to a FASTA file
#'
#' `AppendValidFa()` takes a set of sequences and their headers, checks whether
#' each sequence represents a valid coding sequence (using `IsValidCoding`),
#' and appends only the valid ones to a FASTA file. If a sequence is not valid
#' in the forward orientation, its reverse complement is checked instead.
#'
#' @param h A character vector of sequence headers (FASTA identifiers).
#'          Must be the same length as `s`. Must NOT contain the header
#'          starting character ">" (this is handled by Biostrings).
#' @param s A character vector of nucleotide sequences.
#' @param p A character string of length 1 giving the path to the output 
#'          FASTA file. The file is created if it does not exist, 
#'          otherwise valid sequences are appended.
#'
#' @details
#' This function requires the **Biostrings** package. It uses 
#' [Biostrings::DNAStringSet()], [Biostrings::reverseComplement()], and 
#' [Biostrings::writeXStringSet()] internally, but does not attach the package 
#' to the search path (via `library()`).
#'
#' The function checks validity of each sequence using `IsValidCoding`. Because
#' `IsValidCoding` is not vectorized, it is wrapped in a `vapply()` call for
#' safety. For sequences not valid in the forward orientation, their reverse
#' complement is tested. Only sequences that are valid in at least one 
#' orientation are written to the FASTA file. 
#'
#' Headers are supplied via `h` and written as standard FASTA identifiers 
#' (prefixed with ">" by Biostrings automatically).
#'
#' @return Invisibly returns the number of sequences written to the file.
#'         Returns `0L` (invisibly) if no sequence was valid.
#'
#' @examples
#' \dontrun{
#' headers <- c("seq1", "seq2", "seq3")
#' seqs <- c("ATGAAATGA", "TTTT", "ATGGGGGTAA")
#' out_file <- tempfile(fileext = ".fa")
#' AppendValidFa(headers, seqs, out_file)
#' }
#'
#' @seealso [Biostrings::DNAStringSet()], [Biostrings::writeXStringSet()],
#'          [IsValidCoding()]
#'
#' @importFrom Biostrings DNAStringSet reverseComplement writeXStringSet
#' @export
AppendValidFa <- function(h, s, p) {
  ### check for required package
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    ### stop prints "Error: " by default
    stop("package 'Biostrings' is required but not installed.")
  }
  
  ### check the output path is valid
  stopifnot(is.character(p), length(p) == 1L)
  
  if (!length(h)) return(invisible(0L))
  dssSeqs <- DNAStringSet(s)
  dssSeqsRc <- reverseComplement(dssSeqs)
  
  ### since IsValidCoding isn't vectorized we wrap it with vapply
  VcIsValidCoding <- function(x) vapply(x, IsValidCoding, logical(1))
  
  ### a named boolean with the sequence (name) and if it is valid
  nbooVaFw <- VcIsValidCoding(as.character(dssSeqs))
  ### initialise valid reverse-complement as FALSE
  nbooVaRc <- logical(length(dssSeqs))
  ### sequences that need to be the reverse-complement
  booNeedRc <- !nbooVaFw
  if (any(booNeedRc)) {
    nbooVaRc[booNeedRc] <- VcIsValidCoding(as.character(dssSeqsRc[booNeedRc]))
  }
  ### check if there is anything valid
  booKeep <- nbooVaFw | nbooVaRc
  if (!any(booKeep)) return(invisible(0L))
  
  ### choose forward by default, then replace where RC is needed
  out <- dssSeqs[booKeep]
  booUseRc <- !nbooVaFw[booKeep]
  out[booUseRc] <- dssSeqsRc[booKeep][booUseRc]
  ### header is formatted by writeXStringSet (it adds the ">")
  names(out) <- h[booKeep]
  
  writeXStringSet(out,
                  filepath = p,
                  format = "fasta",
                  append = T)
  
  return(invisible(length(out)))
}

#' Find Non-Low Outliers Based on Sequence Length (with Threshold)
#'
#' `FindNonLowOutliers()` identifies the indexes of sequences 
#' whose lengths are not considered
#' lower outliers, based on the interquartile range (IQR) method.
#'
#' This function removes `NA` values from the input, computes the length
#' of each valid sequence, and returns both the indexes of sequences whose
#' lengths are greater than or equal to the computed lower bound and the value
#' of that bound.
#'
#' @param seq A character vector of nucleotide or amino acid sequences.
#'            May contain `NA` values.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{indexes}{Integer vector of indexes of non-lower-outlier sequences, 
#'   relative to the original input vector.}
#'   \item{limit}{The numeric value of the lower bound 
#'   used to define outliers.}
#' }
#'
#' @examples
#' seqs <- c("ATTATATGTGATGAGC", "ATGCGTATGATCCATGAATGCGT", "AT", NA)
#' res <- FindNonLowOutliers(seqs)
#' res$indexes  # non-outlier indices
#' res$limit    # lower bound used
#'
#' @export
FindNonLowOutliers <- function(seq) {
  ### remove NAs, if any
  booValid <- !is.na(seq)
  seq <- seq[booValid]
  
  ### get sequence lengths
  lens <- nchar(seq)
  
  ### IQR-based lower bound
  q1 <- quantile(lens, 0.25)
  q3 <- quantile(lens, 0.75)
  iqr <- q3 - q1
  limit <- q1 - 1.5 * iqr
  
  ### return results as a list
  indBon <- which(lens >= limit)
  return(list(indexes = indBon, limit = limit))
}

#' Extract a genomic segment from a named genome object
#'
#' @description
#' Parses a coordinate string of the form `"chr:start-end"` and returns the
#' corresponding subsequence from a genome object.
#'
#' @details
#' This function expects `genome` to be a `DNAStringSet` object 
#' named by seqnames. Internally it:
#' 1) parses `coord` with a strict regex,
#' 2) subtracts 1 from the parsed start,
#' 3) calls [Biostrings::subseq()].
#'
#' **Coordinate convention note**: this implementation **adds 1**
#' to the parsed start position to convert the
#' “**BED-like 0-based**” input to “**GFF-like 1-based**”.
#'
#' @param coord A character scalar like `"chrVII:435624-435679"`.
#' @param genome A genome container with elements accessible 
#'               via `genome[[chr]]`, a named `DNAStringSet` (Biostrings).
#'
#' @return A `DNAString` (the extracted segment) or `NULL` if parsing fails,
#'   the chromosome is not found, or inputs are malformed.
#'
#' @section Dependencies:
#' - **Biostrings**: [Biostrings::subseq()] is used to slice sequences.
#'
#' @importFrom Biostrings subseq
#'
#' @examples
#' # genome <- Biostrings::DNAStringSet(c(chrI = "ACTG..."))
#' # ExtractExon("chrI:10-20", genome)
#'
#' @export
ExtractSequence <- function(coord, genome) {
  m <- regexec("^([^:]+):([0-9]+)-([0-9]+)$", coord)
  parts <- regmatches(coord, m)[[1]]
  if (length(parts) != 4) return(NULL)
  chr <- parts[2]
  start <- as.integer(parts[3]) + 1
  end <- as.integer(parts[4])
  seq <- genome[[chr]]
  if (is.null(seq)) return(NULL)
  bingo <- subseq(seq, start, end)
  return(bingo)
}

#' Build and validate coding sequences from one element (set of exons)
#'
#' @description
#' Given a semicolon-separated list of exon coordinates (e.g.,
#' `"chr:start-end;chr:start-end;..."`) and a genome identifier,
#' this function extracts each exon, concatenates them in **all permutations**
#' (each exon used at most once), evaluates validity on both strands using
#' `IsValidCoding()`, and returns the **longest valid** sequence (if any).
#' Imports:
#'   Biostrings,
#'   gtools
#'
#' @details
#' Steps:
#' 1) Split `coordinates` into exon strings 
#'    and extract each with [ExtractSequence()].
#' 2) Generate all permutations of exons (no repeats).
#' 3) Concatenate with [Biostrings::xscat()] to a candidate CDS.
#' 4) Evaluate forward (`+`) and reverse-complement (`-`) via `IsValidCoding()`
#'    and `RevComp()`.
#' 5) Return the longest valid sequence 
#'    and metadata (strand, exon order, length, coordinates).
#'
#' **Coordinate convention note**: `ExtractSequence()` adds 1 
#' to the parsed start, to convert from 0-based coordinates 
#' to 1-based coordinates.
#'
#' @param coordinates A character scalar with exon coordinates separated by `;`.
#' @param id Genome identifier used to select `genomes[[id]]`.
#' @param genomes A list-like container of genome objects; `genomes[[id]]` must
#'   yield a genome accessible by chromosome names (see [ExtractSequence()]).
#'
#' @return A named list with fields:
#' \itemize{
#'   \item `strand` (`"+"` or `"-"`)
#'   \item `ord` (comma-separated exon permutation)
#'   \item `seq` (the selected coding sequence as a character)
#'   \item `len` (sequence length, integer)
#'   \item `coo` (coordinates of the regions, character)
#' }
#' or `NULL` if no valid sequence is found or inputs are invalid.
#'
#' @section Dependencies:
#' - **Biostrings**: [Biostrings::xscat()] 
#'   to concatenate exons (XString objects).
#' - **gtools**: [gtools::permutations()] 
#'   to enumerate exon orders.
#' - **Custom**: `IsValidCoding()` and `RevComp()` 
#'  must be defined in scope.
#'
#' @importFrom Biostrings xscat
#' @importFrom gtools permutations
#'
#' @seealso [ExtractSequence()], `IsValidCoding()`, `RevComp()`
#'
#' @examples
#' # genomes <- list(SGDref = Biostrings::DNAStringSet(...))
#' # ProcessElement("chrI:10-20;chrI:30-50", "SGDref", genomes)
#'
#' @export
ProcessElement <- function(coordinates, id, genomes) {
  genome <- genomes[[id]]
  if (is.null(genome)) return(NULL)
  
  c <- unlist(strsplit(coordinates, ";", fixed = TRUE))
  nEx <- length(c)
  if (nEx > 3) return(NULL)
  
  exons <- lapply(c, ExtractSequence, genome = genome)
  if (any(sapply(exons, is.null))) return(NULL)
  
  ### assumption: any exon is picked only once
  mtPerms <- permutations(nEx, nEx, repeats.allowed = FALSE)
  
  lsValSeqs <- list()
  
  for (i in seq_len(nrow(mtPerms))) {
    ### pick one row of the matrix with all the permutations
    permExons <- exons[mtPerms[i, ]]
    cds <- do.call(xscat, permExons)
    cdsFw <- as.character(cds)
    cdsRc <- RevComp(cdsFw)
    
    ### forward strand
    if (!is.null(cdsFw) && IsValidCoding(cdsFw)) {
      lsValSeqs[[length(lsValSeqs) + 1]] <- list(
        strand = "+",
        ord = paste(mtPerms[i, ], collapse = ","),
        seq = cdsFw,
        len = nchar(cdsFw),
        coo = coordinates
      )
    }
    ### reverse strand
    if (!is.null(cdsRc) && IsValidCoding(cdsRc)) {
      lsValSeqs[[length(lsValSeqs) + 1]] <- list(
        strand = "-",
        ord = paste(mtPerms[i, ], collapse = ","),
        seq = cdsRc,
        len = nchar(cdsRc),
        coo = coordinates
      )
    }
  }
  
  # pick the longest valid sequence
  if (length(lsValSeqs) == 0) return(NULL)
  best <- lsValSeqs[[which.max(sapply(lsValSeqs, `[[`, "len"))]]
  return(best)
}

#' Process one table row across multiple genomes/columns
#'
#' @description
#' For a given input row (e.g., from a `data.table`), evaluate each coordinate
#' column with [ProcessElement()], collecting the best valid CDS (if any) per
#' genome/column.
#'
#' @details
#' - `row` must provide `"Features_id"` and the requested coordinate columns.
#' - Each `coords` string is passed to [ProcessElement()] with `id = col`.
#' - Results are returned as a list of small lists ready to be rbind-ed / bound
#'   into a data.table/data.frame.
#'
#' **Coordinate convention note**: the same **internal 1→0 index conversion**
#' applies via [ExtractExon()] (start is decremented by 1 before slicing).
#'
#' @param row A single-row list-like object supporting `[[` access 
#'            by column name (e.g., a `data.table` row obtained via `dt[i]`).
#' @param cols Character vector of column names that contain coordinate strings.
#' @param genomes A list-like container of genome objects; 
#'                see [ProcessElement()].
#'
#' @return A list of results; each element is a list with fields:
#'   `Feature_id`, `Genome_id`, `Strand_str`, `Exon_order`, `Seq_str`, 
#'   `Seq_length`.
#'   Returns an empty list if nothing valid is found.
#'
#' @section Dependencies:
#' - Relies on [ProcessElement()], which in turn uses 
#'   **Biostrings** and **gtools**,
#'   and requires `IsValidCoding()` and `RevComp()` to be defined.
#'
#' @seealso [ProcessElement()], [ExtractExon()]
#'
#' @examples
#' # res <- ProcessRow(dt[1], cols = c("SGDref#0","DBVPG6765#0"), 
#' #                   genomes = genomes)
#'
#' @export
ProcessRow <- function(row, cols, genomes) {
  id <- row[["Features_id"]]
  lsResults <- list()
  
  for (col in cols) {
    coords <- row[[col]]
    res <- ProcessElement(coords, col, genomes)
    if (is.null(res)) next
    lsResults[[length(lsResults) + 1]] <- list(
      Feature_id = id,
      Genome_id = col,
      Strand_str = res$strand,
      Exon_order = res$ord,
      Seq_str = res$seq,
      Seq_length = res$len,
      Coords_id = res$coo
    )
  }
  lsResults
}

#' Append sequences from list to a FASTA file (Biostrings-backed)
#'
#' This function appends a sequence to a FASTA file, including its 
#' header. The latter is built to be unique. The forward and the reverse 
#' strands are translated ("+" to "W", "-" to "C") to improve readability.
#' 
#' @param lst A list of elements each with `Genome_id`, 
#'            `Exon_order`, `Strand_str`, `Coords_id`, `Seq_str.`
#' @param outfile The path to FASTA file (created if missing).
#' @param width The wrap width (default 80).
#' @return invisibly, number of sequences written
AppendListToFasta <- function(lst, outfile, width = 80L) {
  stopifnot(is.list(lst), is.character(outfile), length(outfile) == 1L)
  if (!requireNamespace("Biostrings", quietly = T)) {
    stop("Biostrings is required.")
  }
  
  ### keep only well-formed elements
  need <- c("Genome_id", "Exon_order", "Strand_str", "Coords_id", "Seq_str")
  ok <- vapply(lst, function(el) all(need %in% names(el)), logical(1))
  if (!all(ok)) warning("Skipping ", sum(!ok), " malformed element(s).")
  lst <- lst[ok]
  if (!length(lst)) return(invisible(0L))
  
  headers <- vapply(lst, function(el) {
    strWC <- switch(el$Strand_str,
                    "+" = "W",
                    "-" = "C",
                    el$Strand_str)
    paste0(el$Genome_id, "-",
           el$Exon_order, "-",
           strWC, "-",
           el$Coords_id)
  }, character(1))
  
  seqs <- vapply(lst, function(el) el$Seq_str, character(1))
  
  dss <- Biostrings::DNAStringSet(seqs)
  names(dss) <- headers
  
  Biostrings::writeXStringSet(dss, filepath = outfile,
                              append = file.exists(outfile),
                              width = width)
  invisible(length(dss))
}

## settings -------------------------------------------------------------------

### reference genome
idRefHap <- "SGDref-0"
### pattern of systematic genes
ptnSys <- "Y[A-P][L,R][0-9]{3}[W,C](-[A-Z])?"
### pattern of random id genes
ptnRid <- "_G[0-9]{7}"
### the type of sub-blocks to process
strSblock <- "cds"
### minimal frequency of a gene to be labelled as core
freqCore <- 0.95
### if at least minPointsOutliers are found in the CDS fasta files 
### FindNonLowOutliers is used to calculate the minimum size to accept a CDS; 
### otherwise minCdsSize is used 
### (it is a number slightly lower than the shortest CDS known)
minPointsOutliers <- 30
minCdsSize <- 69

### counters
countN <- 1          # number of sub-blocks with one or more missing sequence
countSeqMissing <- 1 # number of missing sequences
countSeqFound <- 1   # number of valid sequences found in the genome fasta 

### fixed settings
dirBase <- dirname(this.dir())

### columns to skip when we make the fasta files
skippedCols <- c("Class_id", "Features_id", "Ν_pres", "F_pres", "N_feats",
                 "N_feats_sys", "N_feats_rid", "Sblock_type",
                 "N_total", "N_invalid", "F_pres_genome")

### input data
pathEsse <- file.path(Sys.getenv("HOME"), "data",
                      "SGD-essential-genes", "essentiality.txt")
pathInPan <- file.path(dirBase, "png", "pan-features.RData")
dirCds <- file.path(Sys.getenv("HOME"), "data",
                    "nano-assemblies-pansn-2024", "cds")
dirGen <- file.path(Sys.getenv("HOME"), "data",
                    "nano-assemblies-pansn-2024", "genomes")
### output data
dirOut <- file.path(dirBase, "seqs", strSblock)
unlink(dirOut, recursive = T)
dir.create(dirOut, recursive = T)

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Making the CDS files. ",
    "\n", sep = "")

### load the pangenome (dtPanFeats)
load(pathInPan)
### dev dtPanFeats <- dtPanFeats[1:20, 1:10, with = F]

### get the columns corresponding to genomes
bonCols <- setdiff(names(dtPanFeats), skippedCols)

### we need the first and last haplotype columns to correctly subset 
### the data-table when we calculate the statistics 
nHaplos <- ncol(dtPanFeats) - 2
indHapCols <- 3:ncol(dtPanFeats)

### load essential non-evolvable (esne) genes
dtEsse <- fread(pathEsse, sep = "\t", col.names = c("Gene_id", "Gene_type"))
### get the list of non-evolvable genes
vtEsne <- dtEsse[Gene_type == "non-evolvable", unique(Gene_id)]

### find all the CDS fasta files and their paths
pathAllCds <- list.files(dirCds, pattern = "cds.fa.gz", full.names = T)

### find all the genomes fasta files and their paths
pathGenomes <- list.files(path = dirGen, pattern = "-genome\\.fa\\.gz$",
                          recursive = T, full.names = T)
fileNames <- basename(pathGenomes)
patternGen <- paste0("^(", paste(sub("#" ,"-", bonCols), collapse = "|"), ")")
booMatches <- grepl(patternGen, fileNames)
matchedFiles <- pathGenomes[booMatches]
matchedFilesNames <- sub("-", "#",
                         sub("-genome.fa.gz", "", basename(matchedFiles)))

## collapse haplotypes into genomes -------------------------------------------

### keep dtPanFeats with e.g. chrIV:12-234 instead of S288C#0#chrIV:12-234
dtTmp <- copy(dtPanFeats)
AddColPref(dtTmp, indHapCols)
ReplaceNAtoEmpty(dtTmp)
### remove numeric suffix
vtGnmIds <- gsub("#[0-9]+$", "", names(dtTmp)[indHapCols])
### change column names and identify unique ids
names(dtTmp)[indHapCols] <- vtGnmIds
vtGnmIdsUniq <- unique(vtGnmIds)

### initialize the data.table with CDS presence based on genome
dtPanFeatsPresGnm <- dtTmp[, .(Class_id, Features_id)]

### process each unique prefix
for (indP in vtGnmIdsUniq) {
  vtColMatched <- which(names(dtTmp) == indP)
  
  if (length(vtColMatched) > 1) {
    colsBon <- dtTmp[, ..vtColMatched]
    dtPanFeatsPresGnm[, (indP) := do.call(paste, c(colsBon, sep = ";"))]
    dtPanFeatsPresGnm[, (indP) := gsub("^;|;$", "", get(indP))]
  } else if (length(vtColMatched) == 1) {
    dtPanFeatsPresGnm[, (indP) := dtTmp[[vtColMatched]]]
  }
}
ReplaceEmptyToNA(dtPanFeatsPresGnm)

### garbage collection
rm(dtTmp)
invisible(gc())

## presence (genome-based) ----------------------------------------------------

nGenomes <- ncol(dtPanFeatsPresGnm) - 2
indGnmCols <- 3:ncol(dtPanFeatsPresGnm)

### number of genomes with at least one region in the sub-block
dtPanFeatsPresGnm[, Ν_pres := rowSums(!is.na(as.matrix(.SD))),
                  .SDcols = indGnmCols]
### fraction of genomes with at least one region in the sub-block
dtPanFeatsPresGnm[, F_pres := c(Ν_pres / nGenomes) ]

### add F_pres in dtPanFeats
### by matching dtPanFeats[, "Features_id"]
### and dtPanFeatsPresGnm[, "Features_id"]
dtPanFeats[dtPanFeatsPresGnm,
           F_pres_genome := i.F_pres,
           on = .(Features_id = Features_id)]
rm(dtPanFeatsPresGnm)

## block annotation as non-evolvable, core, or dispensable --------------------

dtPanFeats[, Sblock_type := fifelse(
  grepl(paste(vtEsne, collapse = "|"), Features_id),
  "non-evolvable",
  fifelse(F_pres_genome > freqCore, "core", "dispensable")
)]

## read the genomes -----------------------------------------------------------

### a list of named genomes of class DNAStringSet,
### e.g. lsDssGenomes[["SGDref#0"]] gives the DNAStringSet
### with the SGDref#0 genome
lsDssGenomes <- setNames(lapply(matchedFiles, readDNAStringSet),
                         matchedFilesNames)

## read all the CDS files and build their table with the sequences ------------

### TODO: faster with readDNAStringSet (and requires less memory)
dtAllCds <- ReadCds(dirCds, bonCols)

## make fasta files (one per row) ---------------------------------------------

### dev indR <- 1002
for (indR in 1:nrow(dtPanFeats)) {
  ## find a meaningful unique id for the output file --------------------------
  
  ### get the first systematic name, if any
  strSysId <- regmatches(dtPanFeats[indR, Features_id],
                         regexpr(ptnSys, dtPanFeats[indR, Features_id]))
  if (length(strSysId) == 0) {
    ### get the first random name
    strSysId <- strsplit(dtPanFeats[indR, Features_id],
                         split = ",", fixed = T)[[1]][1]
  }
  pathFaOut <- file.path(dirOut,
                         paste0(indR, "-", strSysId,
                                "-", dtPanFeats[indR, Sblock_type], ".fa"))
  ### clean the output file if it already exists
  if (file.exists(pathFaOut)) {
    unlink(pathFaOut)
  }
  
  ## grab sequences from the CDS database -------------------------------------
  
  ### proteins of unknown function (e.g. YER189W)
  ### may be missing in the CDS fasta of SGDref
  ### even if they are reported in the gff
  
  ### we make a perl regexp with lookbehind and lookahead;
  ### e.g. S288C_G0016460,YER190W,YER189W is grepped with
  ### S288C_G0016460|YER190W|YER189W
  ### to avoid to get e.g. YEL075C/YER189W from YER189W;
  ### do not confuse the "|" in the patter (which acts as an OR)
  ### with "|" in x (dtAllCds[["Seq_id"]]) 
  ### which is just the default string separator in the input fasta files
  strPatterns <- paste0("(?<=^|\\|)(",
                        gsub(",", "|", dtPanFeats[indR, Features_id]),
                        ")$")
  ### Seq_id is e.g. SGDref_G0000040|SGDref_G0000040.mRNA.1|YAL067C
  indHit <- grep(strPatterns, dtAllCds[["Seq_id"]], value = F, perl = T)
  if (length(indHit) > 0) {
    strHds  <- paste0(dtAllCds[["Source_id"]][indHit],
                      " ",
                      dtAllCds[["Seq_id"]][indHit])
    strSeqs <- dtAllCds[["Sequence"]][indHit]
    ### filter out small sequences (e.g. wrong annotations like 
    ### CMF_HP1_G0006560|CMF_HP1_G0006560.mRNA.1|CMF_HP1_G0006560
    ### since otherwise it will end up with YNL233W)
    booEnoughCdsFasta <- length(indHit) > minPointsOutliers
    if (booEnoughCdsFasta) {
      lsNonLowOut <- FindNonLowOutliers(strSeqs)
      indLenVal <- lsNonLowOut[[1]]
      minLength <- lsNonLowOut[[2]]
      AppendValidFa(gsub(" ", "-", strHds[indLenVal]),
                    strSeqs[indLenVal],
                    pathFaOut)
    } else {
      minLength <- minCdsSize
      indLenVal <- which(nchar(strSeqs) > minLength)
      AppendValidFa(gsub(" ", "-", strHds[indLenVal]),
                    strSeqs[indLenVal],
                    pathFaOut)
    }
    
    doneCols <- sub("^([^ ]+).*", "\\1", strHds[indLenVal])
    naCols <- colnames(dtPanFeats)[is.na(dtPanFeats[indR, ])]
    caccaCols <- c(doneCols,
                   naCols,
                   skippedCols)
    ### columns of strains with regions reported in the pangenome data-table
    ### that can get a CDS from the genome fasta
    missCols <- setdiff(colnames(dtPanFeats), caccaCols)
    lsSeqsFound <- ProcessRow(dtPanFeats[indR, ], missCols, lsDssGenomes)
    if (booEnoughCdsFasta) {
      booBon <- vapply(lsSeqsFound,
                       function(lsElem)
                         lsElem$Seq_length > minLength,
                       logical(1))
      nSeqFound <- AppendListToFasta(lsSeqsFound[booBon], pathFaOut)
    } else {
      booBon <- vapply(lsSeqsFound,
                       function(lsElem)
                         lsElem$Seq_length > minCdsSize,
                       logical(1))
      nSeqFound <- AppendListToFasta(lsSeqsFound[booBon], pathFaOut)
    }
    ### as.integer() converts F to 0 and T to 1
    countN <- countN + as.integer(length(missCols) > 0)
    countSeqMissing <- countSeqMissing + length(missCols)
    countSeqFound <- countSeqFound + nSeqFound
  }
}
### here "missing" means that we found homology on the graph but 
### in the annotations the cds fasta is missing
### (thus we do not count the cases where no homology was detected 
### since these elements are NA)
cat("Number of sub-blocks with one (or more) missing or too short CDS: ",
    countN, "\n",
    "Number of missing or too short CDS fasta: ",
    countSeqMissing ,"\n",
    "Number of CDS built from their genome: ",
    countSeqFound, "\n",
    sep = "")
