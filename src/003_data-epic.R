#!/usr/bin/env Rscript
#SBATCH --job-name=003_data_epic
#SBATCH --output=src/logs/003_data_epic.out
#SBATCH --error=src/logs/003_data_epic.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=low_p

options(stringsAsFactors = FALSE)

usage <- function() {
  cat(
    "Usage: Rscript src/003_data-epic.R [options]\n",
    "\n",
    "Build a tab-delimited EPIC study phenotype file from:\n",
    "  genetics_caco.sas7bdat, genetics_id.sas7bdat, genetics.sas7bdat\n",
    "\n",
    "Options:\n",
    "  --input-dir DIR            Directory containing the three SAS files.\n",
    "  --output FILE              Output txt path.\n",
    "  --id-column NAME           ID column to use; default: Idepic_Bio.\n",
    "  --missing-phenotype VALUE  PLINK missing phenotype code; default: -9.\n",
    "  --help                     Show this help.\n",
    "\n",
    sep = ""
  )
}

script_path <- function() {
  args <- commandArgs(FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  }
  normalizePath("src/003_data-epic.R", mustWork = FALSE)
}

parse_args <- function(args, defaults) {
  opts <- defaults
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (arg == "--help" || arg == "-h") {
      usage()
      quit(status = 0)
    }

    if (grepl("^--[^=]+=", arg)) {
      key <- sub("=.*$", "", arg)
      value <- sub("^[^=]+=", "", arg)
    } else if (grepl("^--", arg)) {
      key <- arg
      i <- i + 1
      if (i > length(args)) {
        stop("Missing value for ", key, call. = FALSE)
      }
      value <- args[[i]]
    } else {
      stop("Unexpected argument: ", arg, call. = FALSE)
    }

    if (key == "--input-dir") {
      opts$input_dir <- value
    } else if (key == "--output") {
      opts$output <- value
    } else if (key == "--id-column") {
      opts$id_column <- value
    } else if (key == "--missing-phenotype") {
      opts$missing_phenotype <- as.integer(value)
      if (is.na(opts$missing_phenotype)) {
        stop("--missing-phenotype must be an integer.", call. = FALSE)
      }
    } else {
      stop("Unknown option: ", key, call. = FALSE)
    }

    i <- i + 1
  }
  opts
}

require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Required R package is not installed: ", pkg,
      "\nInstall it before running this script.",
      call. = FALSE
    )
  }
}

read_sas_df <- function(path) {
  as.data.frame(haven::zap_labels(haven::read_sas(path)), stringsAsFactors = FALSE)
}

as_clean_character <- function(x) {
  x <- as.character(x)
  x[x == ""] <- NA_character_
  x
}

as_clean_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

collapse_character <- function(values) {
  values <- unique(stats::na.omit(as_clean_character(values)))
  if (length(values) == 0) {
    return(NA_character_)
  }
  values[[1]]
}

collapse_binary_any <- function(values) {
  values <- unique(stats::na.omit(as_clean_numeric(values)))
  if (length(values) == 0) {
    return(NA_integer_)
  }
  as.integer(max(values))
}

collapse_single_numeric <- function(values) {
  values <- unique(stats::na.omit(as_clean_numeric(values)))
  if (length(values) == 0) {
    return(NA_integer_)
  }
  if (length(values) == 1) {
    return(as.integer(values[[1]]))
  }
  NA_integer_
}

aggregate_by_id <- function(df, id_col, char_cols = character(), binary_cols = character(),
                            single_numeric_cols = character()) {
  df[[id_col]] <- as_clean_character(df[[id_col]])
  df <- df[!is.na(df[[id_col]]), , drop = FALSE]
  split_df <- split(df, df[[id_col]], drop = TRUE)

  out <- data.frame(ID = names(split_df), stringsAsFactors = FALSE)

  for (col in char_cols) {
    out[[col]] <- vapply(split_df, function(part) collapse_character(part[[col]]), character(1))
  }
  for (col in binary_cols) {
    out[[col]] <- vapply(split_df, function(part) collapse_binary_any(part[[col]]), integer(1))
  }
  for (col in single_numeric_cols) {
    out[[col]] <- vapply(split_df, function(part) collapse_single_numeric(part[[col]]), integer(1))
  }

  out
}

assert_columns <- function(df, cols, label) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(
      label, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

warn_unexpected_values <- function(df, cols, allowed, label) {
  for (col in cols) {
    values <- unique(stats::na.omit(as_clean_numeric(df[[col]])))
    unexpected <- setdiff(values, allowed)
    if (length(unexpected) > 0) {
      warning(
        label, " column ", col, " has unexpected value(s): ",
        paste(unexpected, collapse = ", "),
        call. = FALSE
      )
    }
  }
}

warn_collapsed_conflicts <- function(df, id_col, cols, label, numeric = FALSE) {
  df[[id_col]] <- as_clean_character(df[[id_col]])
  df <- df[!is.na(df[[id_col]]), , drop = FALSE]
  split_df <- split(df, df[[id_col]], drop = TRUE)

  for (col in cols) {
    conflict <- vapply(
      split_df,
      function(part) {
        values <- if (numeric) as_clean_numeric(part[[col]]) else as_clean_character(part[[col]])
        length(unique(stats::na.omit(values))) > 1
      },
      logical(1)
    )
    n_conflict <- sum(conflict)
    if (n_conflict > 0) {
      warning(
        label, " column ", col, " has conflicting non-missing values for ",
        n_conflict, " ID(s); those collapsed values are set to missing.",
        call. = FALSE
      )
    }
  }
}

collapse_numeric_by_id <- function(df, id_col, value_col) {
  df[[id_col]] <- as_clean_character(df[[id_col]])
  df <- df[!is.na(df[[id_col]]), , drop = FALSE]
  split_df <- split(df, df[[id_col]], drop = TRUE)

  data.frame(
    ID = names(split_df),
    value = vapply(split_df, function(part) collapse_single_numeric(part[[value_col]]), integer(1)),
    stringsAsFactors = FALSE
  )
}

plink_case_status <- function(in_study, case_control, missing_phenotype) {
  in_study <- as_clean_numeric(in_study)
  case_control <- as_clean_numeric(case_control)

  out <- rep.int(missing_phenotype, length(in_study))
  included <- !is.na(in_study) & in_study == 1

  out[included & !is.na(case_control) & case_control == 0] <- 1L
  out[included & !is.na(case_control) & case_control == 1] <- 2L
  as.integer(out)
}

study_case_control <- function(caco_raw, ids, included, study_id_col, caco_col, id_column) {
  caco_raw[[id_column]] <- as_clean_character(caco_raw[[id_column]])
  caco_raw$Proj_Acronym <- as_clean_character(caco_raw$Proj_Acronym)
  included <- as_clean_numeric(included)

  study_rows <- caco_raw[
    !is.na(caco_raw$Proj_Acronym) &
      caco_raw$Proj_Acronym == study_id_col &
      !is.na(caco_raw[[id_column]]),
    ,
    drop = FALSE
  ]
  out <- rep(NA_integer_, length(ids))

  if (nrow(study_rows) > 0) {
    warn_collapsed_conflicts(
      study_rows,
      id_column,
      caco_col,
      paste0("genetics_caco.sas7bdat Proj_Acronym=", study_id_col),
      numeric = TRUE
    )
    study_status <- collapse_numeric_by_id(study_rows, id_column, caco_col)
    out <- study_status$value[match(ids, study_status$ID)]
  } else {
    if (any(!is.na(included) & included == 1)) {
      warning(
        "genetics_caco.sas7bdat has no rows with Proj_Acronym=", study_id_col,
        "; included samples for this study will have missing phenotype.",
        call. = FALSE
      )
    }
  }

  out
}

study_map <- data.frame(
  study = c(
    "Brea_01_Erneg",
    "Brea_02_Onco",
    "Clrt_01_Gecco",
    "Ecvd_01",
    "Ecvd_02",
    "Ecvd_03",
    "Glbd_01",
    "Inte_01",
    "Inte_02",
    "Inte_03",
    "Kidn_01",
    "Kidn_02",
    "Lung_01",
    "Lymp_01",
    "Neuro_01",
    "Panc_01_PS1",
    "Panc_02_PS3",
    "Pros_01_Bpc3",
    "Pros_02_Icogs",
    "Pros_03_Onco",
    "Pros_04_P160555",
    "Ovar_01",
    "Stom_01",
    "Uadt_01"
  ),
  id_col = c(
    "GW_BREA_01",
    "GW_BREA_02",
    "GW_CLRT_01",
    "GW_ECVD_01",
    "GW_ECVD_02",
    "GW_ECVD_03",
    "GW_GLBD_01",
    "GW_INTE_01",
    "GW_INTE_02",
    "GW_INTE_03",
    "GW_KIDN_01",
    "GW_KIDN_02",
    "GW_LUNG_01",
    "GW_LYMP_01",
    "GW_NEUR_01",
    "GW_PANC_01",
    "GW_PANC_02",
    "GW_PROS_01",
    "GW_PROS_02",
    "GW_PROS_03",
    "GW_PROS_04",
    "GW_OVAR_01",
    "GW_STOM_01",
    "GW_UADT_01"
  ),
  caco_col = c(
    "Cncr_Caco_Brea_Orig",
    "Cncr_Caco_Brea_Orig",
    "Cncr_Caco_Clrt_Orig",
    "Ecvd_Case_Orig",
    "Ecvd_Case_Orig",
    "Ecvd_Case_Orig",
    "Cncr_Caco_Live_Orig",
    "InterAct_T2D_Status_Orig",
    "InterAct_T2D_Status_Orig",
    "InterAct_T2D_Status_Orig",
    "Cncr_Caco_Kidn_Orig",
    "Cncr_Caco_Kidn_Orig",
    "Cncr_Caco_Lung_Orig",
    "Cncr_Caco_Lymp_Orig",
    "Cncr_Caco_Neur_Orig",
    "Cncr_Caco_Panc_Orig",
    "Cncr_Caco_Panc_Orig",
    "Cncr_Caco_Pros_Orig",
    "Cncr_Caco_Pros_Orig",
    "Cncr_Caco_Pros_Orig",
    "Cncr_Caco_Pros_Orig",
    "Cncr_Caco_Ovar_Orig",
    "Cncr_Caco_Stom_Orig",
    "Cncr_Caco_Uadt_Orig"
  ),
  stringsAsFactors = FALSE
)

main <- function() {
  require_namespace("haven")

  root <- normalizePath(file.path(dirname(script_path()), ".."), mustWork = FALSE)
  default_input_dir <- file.path(root, "data", "reference", "Epic")
  defaults <- list(
    input_dir = Sys.getenv("EPIC_REF_DIR", default_input_dir),
    output = Sys.getenv(
      "EPIC_CASE_STATUS_FILE",
      file.path(default_input_dir, "EPIC_study_case_status.txt")
    ),
    id_column = Sys.getenv("EPIC_CASE_STATUS_ID_COLUMN", "Idepic_Bio"),
    missing_phenotype = as.integer(Sys.getenv("EPIC_MISSING_PHENOTYPE", "-9"))
  )
  opts <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
  if (is.na(opts$missing_phenotype)) {
    stop("Missing phenotype code must be an integer.", call. = FALSE)
  }

  input_dir <- normalizePath(opts$input_dir, mustWork = FALSE)
  output <- normalizePath(opts$output, mustWork = FALSE)

  caco_file <- file.path(input_dir, "genetics_caco.sas7bdat")
  id_file <- file.path(input_dir, "genetics_id.sas7bdat")
  sex_file <- file.path(input_dir, "genetics.sas7bdat")

  for (path in c(caco_file, id_file, sex_file)) {
    if (!file.exists(path)) {
      stop("Input file not found: ", path, call. = FALSE)
    }
  }

  message("==========================================")
  message(" Building EPIC Study Case Status File")
  message(" Input dir: ", input_dir)
  message(" Output:    ", output)
  message(" ID column: ", opts$id_column)
  message("==========================================")

  caco_raw <- read_sas_df(caco_file)
  id_raw <- read_sas_df(id_file)
  sex_raw <- read_sas_df(sex_file)

  assert_columns(id_raw, c("Country", "Center", "Idepic", opts$id_column, study_map$id_col), "genetics_id.sas7bdat")
  assert_columns(sex_raw, c(opts$id_column, "Sex"), "genetics.sas7bdat")
  assert_columns(caco_raw, c(opts$id_column, "Proj_Acronym", unique(study_map$caco_col)), "genetics_caco.sas7bdat")

  warn_unexpected_values(id_raw, study_map$id_col, c(0, 1), "genetics_id.sas7bdat")
  warn_unexpected_values(caco_raw, unique(study_map$caco_col), c(0, 1), "genetics_caco.sas7bdat")
  warn_unexpected_values(sex_raw, "Sex", c(1, 2), "genetics.sas7bdat")
  warn_collapsed_conflicts(id_raw, opts$id_column, c("Country", "Center", "Idepic"), "genetics_id.sas7bdat")
  warn_collapsed_conflicts(sex_raw, opts$id_column, "Sex", "genetics.sas7bdat", numeric = TRUE)

  id_agg <- aggregate_by_id(
    id_raw,
    opts$id_column,
    char_cols = c("Country", "Center", "Idepic"),
    binary_cols = study_map$id_col
  )

  sex_agg <- aggregate_by_id(
    sex_raw,
    opts$id_column,
    single_numeric_cols = "Sex"
  )

  sex_match <- match(id_agg$ID, sex_agg$ID)

  sex <- sex_agg$Sex[sex_match]
  sex[is.na(sex) | !(sex %in% c(1L, 2L))] <- 0L

  out <- data.frame(
    ID = id_agg$ID,
    country = id_agg$Country,
    centre = id_agg$Center,
    sex = as.integer(sex),
    stringsAsFactors = FALSE
  )

  summary_rows <- vector("list", nrow(study_map))
  for (i in seq_len(nrow(study_map))) {
    study <- study_map$study[[i]]
    id_col <- study_map$id_col[[i]]
    caco_col <- study_map$caco_col[[i]]

    included <- id_agg[[id_col]]
    case_control <- study_case_control(caco_raw, id_agg$ID, included, id_col, caco_col, opts$id_column)
    out[[study]] <- plink_case_status(included, case_control, opts$missing_phenotype)

    summary_rows[[i]] <- data.frame(
      study = study,
      included = sum(!is.na(included) & included == 1),
      controls = sum(out[[study]] == 1),
      cases = sum(out[[study]] == 2),
      missing_included = sum(!is.na(included) & included == 1 & out[[study]] == opts$missing_phenotype),
      stringsAsFactors = FALSE
    )
  }

  order_idx <- order(out$country, out$centre, out$ID, na.last = TRUE)
  out <- out[order_idx, , drop = FALSE]

  output_dir <- dirname(output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  write.table(out, output, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  summary_df <- do.call(rbind, summary_rows)
  message("Rows written: ", nrow(out))
  message("Study columns: ", nrow(study_map))
  message("Case/control summary:")
  print(summary_df, row.names = FALSE)
  message("Done.")
}

if (sys.nframe() == 0) {
  main()
}
