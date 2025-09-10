## header ---------------------------------------------------------------------

### A script that runs an R script (provided as first argument) and
### performs its profiling. The output (*.Rprofvis) is written in 
### the same folder of the script being profiled.

library(this.path)
library(scriptName)
library(profvis)

## settings -------------------------------------------------------------------

myName <- current_filename()
cat("[", myName, "] ", "Running and profiling your script. ",
    "\n", sep = "")

argsVal <- commandArgs(trailingOnly = T)
if (length(argsVal) != 1) {
  cat("[", myName, "] ",
      "Usage: Rscript run-profvis.R <input.R>\n",
      sep = "")
  ### stop prints "Error: " by default
  stop("user error.", call. = F)
}
nameScript <- argsVal[1]

## clmnt ----------------------------------------------------------------------

dirWrk <- dirname(this.path())
pathIn <- file.path(dirWrk, nameScript)
pathOut <- paste0(pathIn, "profvis")
file.exists(pathOut) && file.remove(pathOut)

profvis({
  source(pathIn)
}, prof_output = pathOut)
