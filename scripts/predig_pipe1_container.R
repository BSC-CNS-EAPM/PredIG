### PREDIG PIPELINE - HORUS
### INPUT: CSV with epitope, HLA_allele, uniprot_id

# libs --------------------------------------------------------------------
suppressPackageStartupMessages({
library(dplyr)
library(stringr)
library(argparse)
library(xgboost)
})

# command line arguments --------------------------------------------------
# Create a parser object
parser <- ArgumentParser(description = 'PREDIG PIPELINE: Input type > CSV with epitope, HLA_allele, uniprot ID')

# Add arguments
parser$add_argument('--input', type = 'character', required = TRUE, help = 'Path to the input file')
parser$add_argument('--out', type = 'character', required = TRUE, help = 'Path to the output file')
parser$add_argument('--model', type = 'character', required = TRUE, help = 'Specify which PredIG model to predict with: neoant, noncan or path')
parser$add_argument('--exp_name', type = 'character', required = TRUE, help = 'Experiment name')

# Parse the arguments
args <- parser$parse_args()

# Print the parsed arguments
print(paste("Input Path:", args$input))
print(paste("Output Path:", args$out))
print(paste("PredIG Model Choice:", args$model))
print(paste("Experiment Name:", args$exp_name))



# check if tmp directory exists and otherwise create it -------------------
if (!dir.exists("tmp")) {
  suppressWarnings(dir.create("tmp"))
}

# tmp file dir ------------------------------------------------------------
suppressWarnings(dir.create(paste("tmp/tmp_", args$exp_name, sep = "")))


# INPUT FORMAT CHECK -------------------------------------------------------------
# initial: epitope, HLA_allele, prot_seq
initial <- read.csv(args$input, stringsAsFactors = FALSE)

# try if initial_predig has columns epitope, HLA_allele, prot_seq
if (!all(c("epitope", "HLA_allele", "uniprot_id") %in% colnames(initial))) {
  stop("Input file must contain columns named: 'epitope', 'HLA_allele', 'uniprot_id'")
}

# create initial predig by selecting epitope, HLA_allele, prot_seq
initial_predig <- select(initial, c("epitope", "HLA_allele", "uniprot_id"))

# try if initial_predig has no NA values
if (any(is.na(initial_predig))) {
  stop("Input file contains NA values at one of these columns: epitope, HLA_allele, uniprot_id. Please, remove all NAs.")
}

# Check if initial contains other columns besides epitope, HLA_allele, uniprot_id
if (any(!colnames(initial) %in% c("epitope", "HLA_allele", "uniprot_id"))) {
  # Create initial_non_predig as a copy of initial
  initial_non_predig <- initial
  
  # Ensure epitope and HLA_allele columns exist before using them
  if ("epitope" %in% colnames(initial_non_predig) && "HLA_allele" %in% colnames(initial_non_predig)) {
    # Create phla_id at initial_non_predig
    initial_non_predig$phla_id <- paste(initial_non_predig$epitope, initial_non_predig$HLA_allele, sep = "_")
  }
  
  # Remove epitope, HLA_allele, and prot_seq columns
  initial_non_predig <- select(initial_non_predig, -c("epitope", "HLA_allele", "uniprot_id"))
}


# try if epitope is character
if (!is.character(initial_predig$epitope)) {
  stop("The 'epitope' column in your input file must be of class 'character'.")
}

# try if prot_seq contains any letter not in one letter amino acid code
# Identify rows with invalid characters
invalid_epitope_rows <- which(str_detect(initial_predig$epitope, "[^ACDEFGHIKLMNPQRSTVWY]"))

# Check if there are any invalid rows
if (length(invalid_epitope_rows) > 0) {
  # Print the row numbers of invalid rows
  # Stop execution with an error message
  stop(cat("The 'epitope' column in your input file contains non-amino acid characters in rows: ", paste(invalid_prot_rows, collapse = ", "),
           "\n Please correct these epitopes or remove their rows to prevent a PredIG crash.", sep = ""))
}

# Check if there are epitopes longer than 14 characters or shorter than 8
invalid_epitope_len <- which(str_length(initial_predig$epitope) < 8 | str_length(initial_predig$epitope) > 14)

# Check if there are any invalid rows
if (length(invalid_epitope_len) > 0) {
  # Print the row numbers of invalid rows
  # Stop execution with an error message
  stop(cat("The 'epitope' column in your input file contains epitopes with lengths outside the range of 8 to 14 amino acids in rows: ", 
           paste(invalid_epitope_len, collapse = ", "),
           "\n Please correct these epitopes or remove their rows to prevent a PredIG crash.", sep = ""))
}

# try if HLA_allele is character
if (!is.character(initial_predig$HLA_allele)) {
  stop("'HLA_allele' column in your input file must be of class character")
}


# try if HLA_allele is in the right 4-digit format for HLA-I alleles

# Define the regular expression for the valid HLA allele formats (HLA-A, HLA-B, HLA-C)
valid_hla_regex <- "^HLA-[ABC]\\*[0-9]{2,3}:[0-9]{2,3}$"

# Check if all HLA_allele entries match the valid format
if (!all(str_detect(initial_predig$HLA_allele, valid_hla_regex))) {
  # Identify rows with invalid HLA alleles
  invalid_hla_rows <- which(!str_detect(initial_predig$HLA_allele, valid_hla_regex))
  
  # # Print the row numbers of invalid rows
  # print(invalid_rows)
  
  # Construct the error message with invalid rows included
  error_message <- cat(
    "The 'HLA_allele' column must be any HLA-I allele in the right 4-digits HLA format, such as: HLA-A*01:01 or HLA-A*01:101. ",
    "\nInvalid entries found your input file at HLA_allele column. Precisely at rows: ", paste(invalid_hla_rows, collapse = ", "),
    "\nPlease correct these HLA alleles or remove their rows to prevent a PredIG crash.",
    sep = "")
  
  # Stop execution with the constructed error message
  stop(error_message)
}

# try if uniprot_id is character
if (!is.character(initial_predig$uniprot_id)) {
  stop("The 'uniprot_id' column in your input file must be of class 'character'")
}

print(initial_predig$uniprot_id)

# try if uniprot_id is in Uniprot ID format
valid_uniprot_regex <- "^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[O-Q][0-9][A-Z0-9]{3}[0-9]$|^[A-Z0-9]{3}[0-9][A-Z0-9]{6,}$"

# Check if all uniprot_id entries match the valid format
if (!all(str_detect(initial_predig$uniprot_id, valid_uniprot_regex))) {
  # Identify rows with invalid Uniprot IDs
  invalid_uniprot_rows <- which(!str_detect(initial_predig$uniprot_id, valid_uniprot_regex))
  
  # Construct the error message with invalid rows included
  error_message <- cat(
    "The 'uniprot_id' column must be any Uniprot ID in the right format, such as: P12345, Q9UJX6 or A0A12345678",
    "\nInvalid entries found your input file at uniprot_id column. Precisely at rows: ", paste(invalid_uniprot_rows, collapse = ", "),
    "\nPlease correct these Uniprot IDs or remove their rows to prevent a PredIG crash.",
    sep = "")
  
  # Stop execution with the constructed error message
  stop(error_message)
}


# INPUT PARSER ------------------------------------------------------------
# netcleave input parse ---------------------------------------------------------------
# generate input for NetCleave v2.0 - pred_input 2
netcleave_in <- select(initial, c("epitope", "uniprot_id"))

# export to specified input folder
# add project name to file name (later)
write.csv(netcleave_in, paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_netcleave_in.csv", sep = ""), row.names = FALSE)


# noah input parse --------------------------------------------------------------------
# generate input for NOAH: peptide, HLA
noah_in <- select(initial, c("epitope", "HLA_allele"))
noah_in <- rename(noah_in, peptide = epitope, HLA = HLA_allele)

# export to specified input folder
# add project name to file name (later)
write.table(noah_in, sep = ",", paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_noah_in.csv", sep = ""), col.names = FALSE, row.names = FALSE, quote = F)
#write.table(noah_in, paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_noah_in.csv", sep = ""), row.names = FALSE, quote = F)

# tap input parse ---------------------------------------------------------
# generate input for TAP: epitope
tap_in <- select(initial, c("epitope"))

# export to specified input folder
# add project name to file name (later)
# Write the column data to a text file without header
write.table(tap_in, file = paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_tap_in.csv", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)


# mhcflurry input parse ---------------------------------------------------------------
# generate input for MHCflurry: peptide, allele
mhcflurry_in <- select(initial, c("epitope", "HLA_allele"))

# rename
mhcflurry_in <- rename(mhcflurry_in, peptide = epitope, allele = HLA_allele)

# export to specified input folder
# add project name to file name (later)
write.csv(mhcflurry_in, file = paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_mhcflurry_in.csv", sep = ""), row.names = FALSE)

# pch input parse ---------------------------------------------------
# generate input for physicochemical properties: epitope
physicochem_in <- select(initial, c("epitope"))
physicochem_in <- rename(physicochem_in, peptide = "epitope")

# export to specified input folder
# add project name to file name (later)
write.csv(physicochem_in, paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_physicochem_in.csv", sep = ""), row.names = FALSE)



# RUNNERS -----------------------------------------------------------------
# Run the R scripts for each model


# run netcleave -----------------------------------------------------------
system(paste("cd predictors/netcleave &&",
             " python3 NetCleave.py --predict ", 
             "../../tmp/tmp_", args$exp_name, "/", args$exp_name, "_netcleave_in.csv",
             " --pred_input 2 --mhc_class I --mhc_allele HLA", " && ",
             "cd ../../", sep = ""))

file.rename(from = paste("predictors/netcleave/output/", args$exp_name, "_netcleave_in_NetCleave.csv", sep = ""), 
            to = paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_netcleave_out.csv", sep = ""))



# run NOAH ----------------------------------------------------------------
# export PYTHONPATH="/Users/rocbsc/bsc_phd/BSC/SCRIPTS/noah_run/NOAH:$PYTHONPATH"
system(paste('export PYTHONPATH="predictors/NOAH:$PYTHONPATH" && ',
             "cd predictors &&",
             "python3 NOAH/main_NOAH.py -i ", "../tmp/tmp_", args$exp_name, "/", args$exp_name, "_noah_in.csv",
             " -o ", "../tmp/tmp_", args$exp_name, "/", args$exp_name, "_noah_out.csv", 
             " -model NOAH/model.pkl" ,
             " && cd ..", sep = ""))


# run TAP -----------------------------------------------------------------
#system(paste("Rscript container_tap_runner.R --in tmp_", args$exp_name, "/", args$exp_name, "_tap_in.csv --out tmp_", args$exp_name, sep = ""))
# "/Users/rocbsc/bsc_phd/BSC/A_NEOANTIGENS/A_NeoA_DEVELOPMENT/20240801_predig_containers/a_immuno_plugins/Immunoinformatics-plugin-main/Immunoinformatics/Include/utils.py"

# generate length
tap_in$epitope <- as.character(tap_in$epitope)
tap_in$peptide_len <- str_length(tap_in$epitope)

# generate length dict
tap_len <- table(tap_in$peptide_len)[1]

# generate file per length
dict_sizes <- split(tap_in$epitope, tap_in$peptide_len)

# Function to write FASTA files from a list of sequences
write_fasta <- function(sequences, filename) {
  # Convert the sequences to a fasta format
  fasta_content <- sapply(names(sequences), function(name) {
    paste0(">", name, "\n", sequences[[name]])
  }, USE.NAMES = FALSE)
  
  # Write the fasta content to a file
  writeLines(fasta_content, filename)
  
  # Return the fasta content for printing
  return(fasta_content)
}

# Generate and write FASTA files for each peptide length
for (len in names(dict_sizes)) {
  sequences <- dict_sizes[[len]]
  names(sequences) <- paste0(sequences, "_", seq_along(sequences))  # Create sequence names
  filename <- paste0("tmp/tmp_", args$exp_name, "/",  args$exp_name, "_tap_in_length_", len, ".fasta")
  fasta_content <- write_fasta(sequences, filename)
  cat("Written", filename, "\n")
  # Print the fasta content
  print(fasta_content)
}

# create a list with full paths of fasta files
fasta_files <- list.files(path = paste("tmp/tmp_", args$exp_name, sep = ""), pattern = "*.fasta", full.names = TRUE)
print(fasta_files)

# Define the function to run the tapmap
run_tapmap <- function(fasta_files, tapmap_path, mat = NULL, alpha = NULL, precursor_len = NULL) {
  epitope <- c()
  tap <- c()
  
  # Iterate over each fasta file in the list
  for (file in fasta_files) {
    # Extract peptide length from the filename
    len <- gsub(".*_length_(\\d+)\\.fasta", "\\1", file)
    cat(sprintf("Running tapmap for peptides of size %s\n", len))
    
    command <- "predictors/netctlpan/tapmat_pred_fsa"
    #command <- tapmap_path
    args2 <- c()
    if (!is.null(mat)) {
      args2 <- c(args2, "-mat", mat)
    }
    if (!is.null(alpha)) {
      args2 <- c(args2, "-a", as.character(alpha))
    }
    args2 <- c(args2, "-l", len)
    if (!is.null(precursor_len)) {
      args2 <- c(args2, "-pl", as.character(precursor_len))
    }
    args2 <- c(args2, file)
    
    # Print the constructed command for debugging
    full_command <- paste(command, paste(args2, collapse = " "))
    cat("Executing command:", full_command, "\n")
    
    tryCatch({
      # Run the command and capture output
      # output <- system2(command, args = args2, stdout = TRUE, stderr = TRUE)
      # cat("Output:\n", output, "\n")
      # 
      # Parse the output
      cat("Parsing tapmap output\n")
      out_tab <- read.table(text = system2(command, args = args2, stdout = TRUE, stderr = TRUE),
                            sep = "", header = FALSE, comment.char = "#")
      print(glimpse(out_tab))
      if (nrow(out_tab) > 0) {
        epitope <- c(epitope, out_tab$V2)
        tap <- c(tap, out_tab$V3)
      } else {
        cat("No output from tapmap\n")
      }
    }, error = function(e) {
      cat(sprintf("An error occurred while running the tapmap size=%s: %s\n", len, e$message))
    })
  }
  
  # Create and save the dataframe
  df <- data.frame(epitope = epitope, TAP = tap)
  output_csv <- paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_output_tapmap.csv", sep = "")
  write.csv(df, output_csv, row.names = FALSE)
  cat("Saved TAPmap output to", output_csv, "\n")
  
  return(df)
}

# Example usage (replace with actual parameters)
# exp_name <- "example_experiment"  # Replace with your experiment name
# fasta_dir <- paste("tmp_", exp_name, sep = "")
# fasta_files <- list.files(path = fasta_dir, pattern = "*.fasta", full.names = TRUE)

#tapmap_path <- "/predictors/netctlpan/tapmat_pred_fsa"
mat_path <- "predictors/netctlpan/tap.logodds.mat"
tap_out <- run_tapmap(fasta_files, tapmap_path, mat = mat_path, alpha = 1.0, precursor_len = NULL)

print(tap_out)


# run mhcflurry -----------------------------------------------------------
system(paste("mhcflurry-predict ", 
             "tmp/tmp_", args$exp_name, "/", args$exp_name, "_mhcflurry_in.csv",
             " --out ", "tmp/tmp_", args$exp_name, "/", args$exp_name, "_mhcflurry_out.csv",
             " --no-throw --always-include-best-allele --no-flanking", sep = ""))


# run pch -----------------------------------------------------------------
system(paste("Rscript predictors/pch/predig_pch_calc.R --input tmp/tmp_", args$exp_name, "/", args$exp_name, 
             "_physicochem_in.csv", sep = ""))
             #" --out ","tmp/tmp_", args$exp_name, "/", args$exp_name, "_physicochem_out.csv", 
# test1_physicochem_in.csv

# OUTPUT PARSER -----------------------------------------------------------

# load results ------------------------------------------------------------

# NetCleave
netcleave_out <- read.csv(paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_netcleave_out.csv", sep = ""))

# NOAH
noah_out <- read.csv(paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_noah_out.csv", sep = ""))

# TAP
#tap_out <- read.table(comment.char = "#", sep = "\t", header = F, col.names = c("epitope", "TAP"), args$path_tap) 

# MHCflurry
print("here")
mhcflurry_out <- read.csv(paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_mhcflurry_out.csv", sep = ""))

# Physicochemical
physicochem_out <- read.csv(paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_physicochem_out.csv", sep = ""))


# parse netcleave out -----------------------------------------------------
# rename prediction to netcleave
netcleave_out <- rename(netcleave_out, netcleave = prediction)

# parse target cols
netcleave_pars <- select(netcleave_out, c("epitope", "netcleave"))

# to join
netcleave_tojoin <- netcleave_pars

# # export
# write.csv(netcleave_out, paste(args$out_path, "/", args$exp_name, "_netcleave_output_parsed.csv", sep = ""), row.names = FALSE)


# parse noah out ----------------------------------------------------------
# Define the function
noah_output_parser <- function(noah_results, df_name = "noah_pars") {
  library(dplyr, quietly = TRUE)
  library(argparser, quietly = TRUE)
  
  # Read input file
  noah_res <- read.table(file = noah_results, header = FALSE, sep = "\t", col.names = c("HLA_allele", "epitope", "NOAH"))
  
  # Store the dataframe in the global environment
  assign(df_name, noah_res, envir = .GlobalEnv)
  
  # Cheers message
  message("Cheers! Your NOAH output has been parsed. The dataframe should contain HLA_allele, epitope, and NOAH columns.")
}

# Parse NOAH and store output in the global environment
noah_output_parser(noah_results = noah_out)

# create phla id
noah_pars$phla_id <- paste(noah_pars$epitope, noah_pars$HLA_allele, sep = "_")

# to join
noah_tojoin <- select(noah_pars, c("phla_id", "NOAH"))


# parse tap out -----------------------------------------------------------
# tap results file is already parsed
# to join
tap_tojoin <- tap_out


# parse mhcflurry out -----------------------------------------------------
# select columns predig
mhcflurry_pars <- select(mhcflurry_out, c("peptide", "allele", "mhcflurry_affinity", "mhcflurry_affinity_percentile", "mhcflurry_processing_score", "mhcflurry_presentation_score"))

# rename columns
mhcflurry_pars <- rename(mhcflurry_pars, epitope = peptide, HLA_allele = allele)

# phla id
mhcflurry_pars <- mhcflurry_pars %>%
  mutate(phla_id = paste(epitope, HLA_allele, sep = "_"))

# to join
mhcflurry_tojoin <- select(mhcflurry_pars, -c("mhcflurry_presentation_percentile"))


# parse pch out -----------------------------------------------------------
# rename peptide to epitope
physicochem_pars <- rename(physicochem_out, epitope = peptide)

# to join
physicochem_tojoin <- select(physicochem_pars, -c("tcr_contact"))


# MERGE PARSED OUTPUTS ----------------------------------------------------

# merge all outputs
predig_joined <- left_join(netcleave_tojoin, mhcflurry_tojoin, by = "epitope") %>%
  left_join(tap_tojoin, by = "epitope") %>%
  left_join(noah_tojoin, by = "phla_id") %>%
  left_join(physicochem_tojoin, by = "epitope")


# reorder columns ---------------------------------------------------------
# netcleave	NOAH	TAP	mhcflurry_affinity	mhcflurry_affinity_percentile	mhcflurry_processing_score	mhcflurry_presentation_score	mw_peptide	hydroph_peptide	charge_peptide	stab_peptide	mw_tcr_contact	hydroph_tcr_contact	charge_tcr_contact
# reorder columns
predig_joined <- predig_joined %>% select(phla_id, epitope, HLA_allele, 
                                          netcleave, NOAH, TAP, 
                                          mhcflurry_affinity, mhcflurry_affinity_percentile, 
                                          mhcflurry_processing_score, mhcflurry_presentation_score,
                                          mw_peptide, hydroph_peptide, charge_peptide, stab_peptide,	
                                          mw_tcr_contact, hydroph_tcr_contact, charge_tcr_contact)
  


# PREDIG CALCULATION ------------------------------------------------------
# use if try loop to calculdate predig score
# depending on the model chosen at args$model

# prediction
if (args$model == "neoant") {
  predig_neoant <- xgb.load("models/predig_neoant.model")
  feat_space <- select(predig_joined, -c("epitope", "HLA_allele", "phla_id"))
  predig_joined <- predig_joined %>% mutate(predig_score_neoa = predict(predig_neoant, as.matrix(feat_space)))
} else if (args$model == "noncan") {
  predig_noncan <- xgb.load("models/predig_noncan.model")
  feat_space <- select(predig_joined, -c("epitope", "HLA_allele", "phla_id"))
  predig_joined <- predig_joined %>% mutate(predig_score_noncan = predict(predig_noncan, as.matrix(feat_space)))
} else if (args$model == "path") {
  predig_path <- xgb.load("models/predig_path.model")
  feat_space <- select(predig_joined, -c("epitope", "HLA_allele", "phla_id"))
  predig_joined <- predig_joined %>% mutate(predig_score_path = predict(predig_path, as.matrix(feat_space)))
} else {
  stop("Invalid model choice. Please choose from: neoant, noncan, path")
}

# remove non proprietary columns ------------------------------------------
predig_joined <- select(predig_joined, -c("TAP", "mhcflurry_affinity", "mhcflurry_affinity_percentile", "mhcflurry_processing_score", "mhcflurry_presentation_score"))

# left_join initial non predig columns if initial_non_predig exists
if (exists("initial_non_predig")) {
  predig_joined <- left_join(predig_joined, initial_non_predig, by = "phla_id")
}


# # export predig results -------------------------------------------------
write.csv(predig_joined, paste(args$out_path, "/", args$exp_name, "_predig_output.csv", sep = ""), row.names = FALSE)


# remove tmp files --------------------------------------------------------
unlink(paste("tmp/tmp_", args$exp_name, sep = ""), recursive = TRUE)



