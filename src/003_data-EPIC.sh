#!/bin/bash
# Script: src/003_data-EPIC.sh
# Purpose: Prepare the subject ID linkage file (Subj_Id_2015.txt equivalent).

set -euo pipefail

# Robust environment sourcing
for env_file in ".env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../pipeline/.env"; do
  if [ -f "$env_file" ]; then
    set -a; source "$env_file"; set +a
    break
  fi
done

echo "=========================================="
echo " Generating EPIC Subject ID Linkage File"
echo "=========================================="

OUT_FILE="data/reference/Epic/Subj_Id_ALL.csv"
mkdir -p "$(dirname "$OUT_FILE")"

echo "Running R aggregation script..."

Rscript - <<EOF
library(dplyr)
library(readr)

# 1. Find all linkage files in data/genetics
link_files <- list.files("data/genetics", pattern = "Link_Ids|IDkey", recursive = TRUE, full.names = TRUE)

message("Found ", length(link_files), " linkage files.")

# 2. Function to read and harmonize
process_link <- function(f) {
  message("Processing: ", f)
  # Detect separator and read
  tryCatch({
    df <- if(grepl(".csv$", f)) read_csv(f, show_col_types = FALSE) else read_tsv(f, show_col_types = FALSE)
    
    # Harmonization logic: Map various names to common 'LabID' and 'EpicID'
    # This is a template; columns vary by study.
    # df <- df |> rename_with(~case_when(
    #   .x %in% c("Lab_ID", "Sample_ID", "IDBio") ~ "LabID",
    #   .x %in% c("EPIC_ID", "Project_ID") ~ "EpicID",
    #   TRUE ~ .x
    # ))
    
    return(df)
  }, error = function(e) {
    message("Error reading ", f, ": ", e$message)
    return(NULL)
  })
}

# 3. Merge all
# master_df <- lapply(link_files, process_link) |> bind_rows()

# 4. Save
# write_csv(master_df, "${OUT_FILE}")
message("Placeholder: Logic to merge study data into ${OUT_FILE} is ready for refinement.")
EOF

echo "Linkage file generation script template is ready."
