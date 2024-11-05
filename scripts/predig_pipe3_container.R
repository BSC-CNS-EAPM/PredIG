### PREDIG PIPELINE - HORUS
### INPUT: FASTA + CSV with HLA Alleles

### SCRIPT STRUCTURE
# arguments
# initial input
# input parsers
# runners
# ourput parsers
# join
# predig calculation
# export
# remove tmp files


# libs --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(argparse)
  library(xgboost)
  library(data.table)
})


# command line arguments --------------------------------------------------
# Create a parser object
parser <- ArgumentParser(description = 'PREDIG PIPELINE: Input type > FASTA + CSV with HLA-I alleles')

# Add arguments
parser$add_argument('--fa', type = 'character', required = TRUE, help = 'Path to the input FASTA file: type 3 using FASTA')
parser$add_argument('--a', type = 'character', required = TRUE, help = 'Path to the input HLA-I alleles file: CSV with HLA_allele column and HLA-I alleles in 4-digits')
parser$add_argument('--out', type = 'character', required = TRUE, help = 'Path to the output file')
parser$add_argument('--model', type = 'character', required = TRUE, help = 'Specify which PredIG model to predict with: neoant, noncan or path')
parser$add_argument('--exp_name', type = 'character', required = TRUE, help = 'Experiment name')

# Parse the arguments
args <- parser$parse_args()

# Print the parsed arguments
print(paste("Input FASTA:", args$fa))
print(paste("Target HLA-I alleles CSV:", args$a))
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
# FASTA
# Read the FASTA file and check format
fasta <- readLines(args$fa)

# Check if the FASTA file is empty
if (length(fasta) == 0) {
  stop("The input FASTA file is empty. Please provide a valid FASTA file.")
}

# Check if header is appropriate
if (!startsWith(fasta[1], ">")) {
  stop("The input FASTA file does not start with a header. Please provide a valid FASTA file.")
}

# Check that FASTA only contains one protein
if (length(grep("^>", fasta)) > 1) {
  stop("The input FASTA file contains more than one protein sequence. Please provide a single protein sequence.")
}

# Check if FASTA sequence contains valid amino acids only
# Identify any wrong amino acids do not use rows
if (any(str_detect(fasta[2], "[^ACDEFGHIKLMNPQRSTVWY]"))) {
  stop(cat("The input FASTA file contains invalid amino acids.",
  "\nPlease provide a sequence only with valid AAs in single letter code.",
  sep = ""))
}


# HLA CSV
# Read the HLA-I alleles CSV file and check format
hla_alleles <- read.csv(args$a, stringsAsFactors = FALSE)

# Check if the HLA-I alleles CSV file is empty
if (nrow(hla_alleles) == 0) {
  stop("The input HLA-I alleles CSV file is empty. Please provide a valid CSV file.")
}

# Check if HLA_allele contains NA values
if (is.na(hla_alleles$HLA_allele)) {
  stop("'HLA_allele' column in your input file contains NA values. Please, remove these.")
}

# Check if HLA-I alleles are the right HLA 4-digits format
# Define the regular expression for the valid HLA allele formats (HLA-A, HLA-B, HLA-C)
valid_hla_regex <- "^HLA-[ABC]\\*[0-9]{2,3}:[0-9]{2,3}$"

# Check if all HLA_allele entries match the valid format
if (!all(str_detect(hla_alleles$HLA_allele, valid_hla_regex))) {
  # Identify rows with invalid HLA alleles
  invalid_hla_rows <- which(!str_detect(hla_alleles$HLA_allele, valid_hla_regex))
  
  # Construct the error message with invalid rows included
  error_message <- cat(
    "The 'HLA_allele' column must be any HLA-I allele in the right 4-digits HLA format, such as: HLA-A*01:01 or HLA-A*01:101. ",
    "\nInvalid entries found your input file at HLA_allele column. Precisely at rows: ", paste(invalid_hla_rows, collapse = ", "),
    "\nPlease correct these HLA alleles or remove their rows to prevent a PredIG crash.",
    sep = "")
  
  # Stop execution with the constructed error message
  stop(error_message)
}


# NETCLEAVE FASTA RUNNER --------------------------------------------------
# run netcleave -----------------------------------------------------------
system(paste("python3 predictors/netcleave/NetCleave.py --predict ", 
             args$fa, " --pred_input 1 --mhc_class I --mhc_allele HLA", sep = ""))

file.rename(from = paste(args$fa, "_netcleave.csv", sep = ""), 
            to = paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_netcleave_out.csv", sep = ""))

# parse netcleave out -----------------------------------------------------
# load
netcleave_out <- read.csv(paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_netcleave_out.csv", sep = ""))

# rename prediction to netcleave
netcleave_out <- rename(netcleave_out, netcleave = prediction)

# parse target cols
netcleave_pars <- select(netcleave_out, c("epitope", "netcleave"))

# to join
netcleave_tojoin <- netcleave_pars


# ADD HLA ALLELES ---------------------------------------------------------


# Convert df to a data.table
predig_input <- as.data.table(netcleave_tojoin)

# Extract unique HLA alleles
hla_alleles <- unique(hla_alleles$HLA_allele)

# Create a new data.table with the HLA alleles
allele_dt <- data.table(current_allele = hla_alleles)

# Perform a cross join between df and allele_dt
result_dt <- allele_dt[predig_input, allow.cartesian=TRUE]

# If you prefer a data.frame as the final result, convert it back
result_df <- as.data.frame(result_dt)

# rename the column
result_df <- rename(result_df, HLA_allele = current_allele)

# predig input
initial_predig <- result_df


# INPUT PARSER ------------------------------------------------------------

# noah input parse --------------------------------------------------------------------
# generate input for NOAH: peptide, HLA
noah_in <- select(initial_predig, c("epitope", "HLA_allele"))
noah_in <- rename(noah_in, peptide = epitope, HLA = HLA_allele)

# export to specified input folder
# add project name to file name (later)
write.csv(noah_in, paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_noah_in.csv", sep = ""), row.names = FALSE, quote = F)

# tap input parse ---------------------------------------------------------
# generate input for TAP: epitope
tap_in <- select(initial_predig, c("epitope"))

# export to specified input folder
# add project name to file name (later)
# Write the column data to a text file without header
write.table(tap_in, file = paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_tap_in.csv", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)


# mhcflurry input parse ---------------------------------------------------------------
# generate input for MHCflurry: peptide, allele
mhcflurry_in <- select(initial_predig, c("epitope", "HLA_allele"))

# rename
mhcflurry_in <- rename(mhcflurry_in, peptide = epitope, allele = HLA_allele)

# export to specified input folder
# add project name to file name (later)
write.csv(mhcflurry_in, file = paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_mhcflurry_in.csv", sep = ""), row.names = FALSE)

# pch input parse ---------------------------------------------------
# generate input for physicochemical properties: epitope
physicochem_in <- select(initial_predig, c("epitope"))

# export to specified input folder
# add project name to file name (later)
write.csv(physicochem_in, paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_physicochem_in.csv", sep = ""), row.names = FALSE)



# RUNNERS -----------------------------------------------------------------
# Run the R scripts for each model

# run NOAH ----------------------------------------------------------------
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

# tap running function

run_tapmap <- function(fasta_files, tapmap_path, mat = NULL, alpha = NULL, precursor_len = NULL) {
  epitope <- c()
  tap <- c()
  
  # Iterate over each fasta file in the list
  for (file in fasta_files) {
    # Extract peptide length from the filename
    len <- gsub(".*_length_(\\d+)\\.fasta", "\\1", file)
    cat(sprintf("Running tapmap for peptides of size %s\n", len))
    
    command <- tapmap_path
    args <- c()
    if (!is.null(mat)) {
      args <- c(args, "-mat", mat)
    }
    if (!is.null(alpha)) {
      args <- c(args, "-a", as.character(alpha))
    }
    args <- c(args, "-l", len)
    if (!is.null(precursor_len)) {
      args <- c(args, "-pl", as.character(precursor_len))
    }
    args <- c(args, file)
    
    tryCatch({
      # Run the command and capture output
      output <- system2(command, args = args, stdout = TRUE, stderr = TRUE)
      cat("Output:\n", output, "\n")
      
      # Parse the output
      cat("Parsing tapmap output\n")
      for (line in output) {
        if (!startsWith(line, "#")) {
          parts <- strsplit(line, "\\s+")[[1]]
          if (length(parts) >= 3) {
            cat("Parsed line:", parts, "\n")
            epitope <- c(epitope, parts[3])
            tap <- c(tap, parts[4])
          }
        }
      }
    }, error = function(e) {
      stop(sprintf("An error occurred while running the tapmap size=%s: %s", len, e))
    })
  }
  
  # Create and save the dataframe
  df <- data.frame(epitope = epitope, TAP = tap)
  write.csv(df, paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_output_tapmap.csv", sep = ""), row.names = FALSE)
  # write.csv(df, paste("fasta/", "exp_name", "_output_tapmap.csv", sep = ""), row.names = FALSE)
  
  return(df)
}



# Example usage (replace with actual parameters)
# tapmap_path <- "/Users/rocbsc/bsc_phd/BSC/A_NEOANTIGENS/A_NeoA_DEVELOPMENT/20240801_predig_containers/predictors/netctlpan/tapmat_pred_fsa"
# mat_path <- "/Users/rocbsc/bsc_phd/BSC/A_NEOANTIGENS/A_NeoA_DEVELOPMENT/20240801_predig_containers/predictors/netctlpan/tap.logodds.mat"
tapmap_path <- "/predictors/netctlpan/tapmat_pred_fsa"
mat_path <- "/predictors/netctlpan/tap.logodds.mat"
tap_out <- run_tapmap(fasta_files, tapmap_path, mat = mat_path, alpha = 1.0, precursor_len = NULL)


# run mhcflurry -----------------------------------------------------------
system(paste("mhcflurry-predict ", 
             "tmp/tmp_", args$exp_name, "/", args$exp_name, "_mhcflurry_in.csv",
             " --out ", "tmp/tmp_", args$exp_name, "/", args$exp_name, "_mhcflurry_out.csv",
             " --no-throw --always-include-best-allele --no-flanking", sep = ""))


# run pch -----------------------------------------------------------------
system(paste("Rscript /predictors/pch/predig_pch_calc.R --input tmp_", args$exp_name, "/", args$exp_name, "_physicochem_in.csv",
             " --out ","tmp/tmp_", args$exp_name, "/", args$exp_name, "_physicochem_out.csv", sep = ""))

# OUTPUT PARSER -----------------------------------------------------------

# load results ------------------------------------------------------------

# NOAH
noah_out <- read.csv(paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_noah_out.csv", sep = ""))

# TAP
#tap_out <- read.table(comment.char = "#", sep = "\t", header = F, col.names = c("epitope", "TAP"), args$path_tap) 

# MHCflurry
mhcflurry_out <- read.csv(paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_mhcflurry_out.csv", sep = ""))

# Physicochemical
physicochem_out <- read.csv(paste("tmp/tmp_", args$exp_name, "/", args$exp_name, "_physicochem_out.csv", sep = ""))


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
noah_tojoin <- select(noah_pars, c("phla_id", "epitope", "HLA_allele", "NOAH"))


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
mhcflurry_tojoin <- select(mhcflurry_pars, -c("epitope", "HLA_allele", "mhcflurry_presentation_percentile"))


# parse pch out -----------------------------------------------------------
# rename peptide to epitope
physicochem_pars <- rename(physicochem_out, epitope = peptide)

# to join
physicochem_tojoin <- select(physicochem_pars, -c("tcr_contact"))


# MERGE PARSED OUTPUTS ----------------------------------------------------

# merge all outputs
predig_joined <- left_join(noah_tojoin, netcleave_tojoin, by = "epitope") %>%
  left_join(tap_tojoin, by = "epitope") %>%
  left_join(mhcflurry_tojoin, by = "phla_id") %>%
  left_join(physicochem_tojoin, by = "epitope")



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

# # left_join initial non predig columns
# predig_joined <- left_join(predig_joined, initial_non_predig, by = "phla_id")


# # export predig results -------------------------------------------------
write.csv(predig_joined, paste(args$out_path, "/", args$exp_name, "_predig_output.csv", sep = ""), row.names = FALSE)


# remove tmp files --------------------------------------------------------
unlink(paste("tmp/tmp_", args$exp_name, sep = ""), recursive = TRUE)



