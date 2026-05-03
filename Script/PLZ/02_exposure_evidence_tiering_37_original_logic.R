#!/usr/bin/env Rscript
# Author: Zhaoxian Wang
# Purpose: Build FCC exposure-evidence tiering table for 37 plasticizers.

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(tibble)
})

get_arg_value <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name, "=")
  hit <- args[startsWith(args, key)]
  if (length(hit) == 0) return(default)
  sub(key, "", hit[[1]], fixed = TRUE)
}

`%coalesce_chr%` <- function(a, b) ifelse(is.na(a) | a == "", b, a)

resolve_project_root <- function() {
  from_arg <- get_arg_value("project-root", NA_character_)
  from_env <- Sys.getenv("CVD_MS1_PROJECT_ROOT", unset = "")
  cands <- unique(c(if (!is.na(from_arg)) from_arg else NULL, if (nzchar(from_env)) from_env else NULL, getwd()))
  for (s in cands) {
    p <- normalizePath(s, winslash = "/", mustWork = FALSE)
    for (i in 0:6) {
      r <- normalizePath(file.path(p, paste(rep("..", i), collapse = "/")), winslash = "/", mustWork = FALSE)
      if (dir.exists(file.path(r, "data"))) return(r)
    }
  }
  stop("Unable to resolve project root. Use --project-root=/path/to/repo")
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

project_root <- resolve_project_root()
out_dir <- get_arg_value("out-dir", file.path(project_root, "data", "RNA-Seq", "subtype", "reviewer_reply"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fcc_default <- first_existing(c(
  file.path(project_root, "data", "Toxicity", "FCC.xlsx"),
  file.path(project_root, "supp", "FCC.xlsx"),
  file.path(project_root, "data", "Zenodo", "Plasticizer", "FCC.xlsx")
))
fcc_file <- get_arg_value("fcc-file", fcc_default)
if (is.na(fcc_file) || !file.exists(fcc_file)) {
  stop("FCC file not found. Set --fcc-file=/path/to/FCC.xlsx")
}

ps <- read.xlsx(fcc_file, sheet = "Priority_Substances", check.names = TRUE)
tr <- read.xlsx(fcc_file, sheet = "Tier", check.names = TRUE)

ps$Abbr <- ps[[2]]
ps37 <- ps %>% filter(!is.na(Abbr), Abbr != "") %>% slice(1:37)

ps37$Abbr <- ifelse(ps37$Abbr == "E-PS", "EPS", ps37$Abbr)
tr$Abbr <- ifelse(tr[[1]] == "E-PS", "EPS", tr[[1]])

authorities_idx <- 17
vehicle_idx <- 21:30

tier_map <- tr %>%
  transmute(
    Abbr = Abbr,
    tier_info = tier_info,
    hazard_info = hazard_info
  )

df <- ps37 %>%
  transmute(
    Name = .[[1]],
    Abbr = Abbr,
    Pubchem_CID = .[[4]],
    Documenting_Authorities = suppressWarnings(as.numeric(.[[authorities_idx]])),
    tier_info_raw = NA_character_,
    Hazard_Info = NA_character_,
    exposure_fields_missing_n = rowSums(is.na(across(all_of(colnames(ps37)[vehicle_idx])))),
    exposure_vehicle_count = rowSums(across(all_of(colnames(ps37)[vehicle_idx]), ~ ifelse(is.na(.x), 0, as.numeric(.x))))
  ) %>%
  left_join(tier_map, by = "Abbr") %>%
  mutate(
    tier_info_raw = tier_info_raw %coalesce_chr% tier_info,
    Hazard_Info = Hazard_Info %coalesce_chr% hazard_info
  ) %>%
  select(-tier_info, -hazard_info)

df <- df %>%
  mutate(
    Exposure_Tier_Index = case_when(
      grepl("^Tier 1", ifelse(is.na(tier_info_raw), "", tier_info_raw)) ~ "top 1st",
      tier_info_raw == "Tier 2" ~ "top 2nd",
      TRUE ~ NA_character_
    ),
    Exposure_Vehicle = case_when(
      exposure_fields_missing_n == length(vehicle_idx) ~ "undocumented",
      exposure_vehicle_count > 0 ~ "documented",
      TRUE ~ "undocumented"
    ),
    Evidence_Data_Missing = ifelse(
      is.na(Documenting_Authorities) | is.na(Exposure_Tier_Index) | exposure_fields_missing_n > 0,
      "Yes", "No"
    ),
    Evidence_Tier = case_when(
      Exposure_Tier_Index == "top 1st" & Exposure_Vehicle == "documented" & !is.na(Documenting_Authorities) & Documenting_Authorities >= 20 ~ "High",
      Exposure_Tier_Index %in% c("top 1st", "top 2nd") & Exposure_Vehicle == "documented" & !is.na(Documenting_Authorities) & Documenting_Authorities >= 10 ~ "Moderate",
      TRUE ~ "Limited"
    ),
    Rule_Trace = case_when(
      Evidence_Tier == "High" ~ "top1st + documented vehicle + authorities>=20",
      Evidence_Tier == "Moderate" ~ "tier(top1st/top2nd) + documented vehicle + authorities>=10 and not High",
      TRUE ~ "undocumented vehicle OR authorities<10 OR missing/NA key exposure fields"
    ),
    Classification_Scope = "Exposure evidence strength only (not toxicity strength)"
  ) %>%
  arrange(factor(Evidence_Tier, levels = c("High", "Moderate", "Limited")), desc(Documenting_Authorities))

out_tbl <- df %>%
  select(
    Name, Abbr, Pubchem_CID,
    Exposure_Tier_Index,
    Documenting_Authorities,
    Exposure_Vehicle,
    exposure_vehicle_count,
    exposure_fields_missing_n,
    tier_info_raw,
    Evidence_Data_Missing,
    Evidence_Tier,
    Rule_Trace,
    Classification_Scope
  )

out_csv <- file.path(out_dir, "FCCdb_Exposure_Evidence_Tiering_37.csv")
write.csv(out_tbl, out_csv, row.names = FALSE)

summary_tbl <- out_tbl %>%
  count(Evidence_Tier, name = "n_chemicals") %>%
  arrange(factor(Evidence_Tier, levels = c("High", "Moderate", "Limited")))

rubric_tbl <- tibble::tribble(
  ~Rule_Level, ~Definition,
  "High", "Exposure Tier=top 1st AND Exposure vehicle=documented AND Documenting authorities>=20",
  "Moderate", "Not High, but Exposure Tier in {top 1st, top 2nd} AND Exposure vehicle=documented AND Documenting authorities>=10",
  "Limited", "Exposure vehicle=undocumented OR Documenting authorities<10 OR key exposure fields missing/NA"
)

note_tbl <- tibble::tribble(
  ~Note,
  "This tiering is based only on FCCdb exposure-evidence fields used in Fig.2 context (tier_info, documenting authorities, exposure vehicle fields).",
  "This tiering represents exposure evidence strength, not toxicity strength.",
  "ProTox/ADMETlab outputs are not used in this classification."
)

out_xlsx <- file.path(out_dir, "Exposure_Evidence_Tiering_37.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "FCCdb_Tiering_37")
writeData(wb, "FCCdb_Tiering_37", out_tbl)
addWorksheet(wb, "Tier_Summary")
writeData(wb, "Tier_Summary", summary_tbl)
addWorksheet(wb, "Rubric")
writeData(wb, "Rubric", rubric_tbl)
addWorksheet(wb, "Notes")
writeData(wb, "Notes", note_tbl)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("Generated:\n")
cat(" -", out_csv, "\n")
cat(" -", out_xlsx, "\n")
cat("\nTier counts:\n")
print(summary_tbl)
