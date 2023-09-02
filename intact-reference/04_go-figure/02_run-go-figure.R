# Run GO-figure algorithm
# Brenda Pardo
# 2023-05-15

##This code has to be run in R but in the server, before opening R you have to activate conda environment
# cd /home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq/code/GOEnrichment/go-figure/
# conda activate gofigure
# R
library(here)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

##Run this in R
input_folder = here('intact-reference/04_go-figure/input')
output_folder = here('intact-reference/04_go-figure/output')


input_files = list.files(input_folder, pattern='input_file*')

for(input_file in input_files) {
  sample = gsub('input_file_', '', input_file)
  sample = gsub('.tsv', '', sample)
  
  out_folder = file.path(output_folder, sample)
  input_file.c = file.path(input_folder, input_file)
  if(!dir.exists(out_folder)) {
    dir.create(out_folder)
  }
  
  cmd = c(
    args[2], args[1],  '-i', input_file.c, 
    '-v 0.05 -a 20 -m 20 -si 0.4 -e 60 -q "pdf" -o', out_folder
  )
  cmd = paste(cmd, collapse=' ')
  system(cmd)
}
