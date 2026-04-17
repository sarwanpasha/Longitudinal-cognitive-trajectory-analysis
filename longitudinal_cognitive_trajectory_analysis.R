# =============================================================================
# LONGITUDINAL COGNITIVE TRAJECTORY ANALYSIS TOOLKIT
# =============================================================================
#
# Description:
#   A statistical framework for analyzing domain-specific cognitive decline
#   trajectories using linear mixed-effects models. Designed for longitudinal
#   neuropsychological data with repeated measures per participant.
#
# Key Features:
#   - Domain-specific cognitive trajectory modeling (memory, executive function,
#     language/visuospatial, processing speed, global composite)
#   - Time-anchored analysis relative to a clinical event or reference date
#   - Cognitive reserve moderation analysis (education, sex)
#   - Publication-ready figures and tables (Word .docx output)
#   - Battery/instrument stratification sensitivity analysis
#
# Input Requirements:
#   Two tab-delimited data files (see Section 1 for variable documentation):
#     (1) Diagnostic/clinical review file  — one row per participant per review
#     (2) Neuropsychological testing file  — one row per participant per exam
#
# Required R packages:
#   tidyverse, lme4, lmerTest, broom.mixed, emmeans, ggplot2, ggpubr,
#   patchwork, gt, gtsummary, flextable, splines, nlme, scales,
#   RColorBrewer, viridis, cowplot, officer, ggstance
#
# Usage:
#   1. Edit Section 0 (Configuration) to set your file paths and variable names
#   2. Source the script: source("longitudinal_cognitive_trajectory_analysis.R")
#   3. Outputs will be written to your working directory
#
# Outputs:
#   Tables (Word .docx): Table1 (sample characteristics), Table2 (LME results),
#                        Table3 (moderation), TableS1 (full fixed effects),
#                        TableS2 (random effects / model fit)
#   Figures (PDF + PNG): Figure1 (domain trajectories), Figure2 (divergence
#                        timeline), Figure3 (spaghetti plots), Figure4
#                        (education moderation), Figure5 (sex moderation),
#                        Figure6 (forest plot), FigureS1 (battery sensitivity),
#                        FigureS2 (individual test plots)
# =============================================================================


# =============================================================================
# SECTION 0: CONFIGURATION — EDIT THIS SECTION FOR YOUR DATA
# =============================================================================

# ---- File paths --------------------------------------------------------------
# Set to the paths of your two input data files.
# Files are expected to be tab-delimited with a header row.
# If your files have comment/metadata rows at the top, set SKIP_ROWS accordingly.

DIAG_FILE  <- "path/to/your_diagnostic_review_file.txt"   # clinical review file
NPD_FILE   <- "path/to/your_neuropsych_testing_file.txt"   # neuropsych exam file
SKIP_ROWS  <- 0   # number of header/comment rows to skip before the column names row

# ---- Variable name mapping ---------------------------------------------------
# Map the generic variable names used in this script to the actual column names
# in your data files. Edit the RIGHT-HAND side of each assignment.
#
# --- Diagnostic file variables (DIAG_FILE) ---
VAR_ID            <- "shareid"          # unique participant identifier
VAR_CASE_STATUS   <- "DEMRV103"         # case/control indicator
                                        #   0 = control, 1–4 = case, other = exclude
VAR_CDR           <- "DEMRV096"         # Clinical Dementia Rating score (optional)
VAR_DATE_IMPAIR   <- "IMPAIRMENT_DATE"  # SAS date: first clinical impairment
VAR_DATE_NORMAL   <- "NORMAL_DATE"      # SAS date: last known-normal review
VAR_DATE_MILD     <- "MILD_DATE"        # SAS date: mild stage (optional)
VAR_DATE_MODERATE <- "MODERATE_DATE"    # SAS date: moderate stage (optional)
VAR_DATE_SEVERE   <- "SEVERE_DATE"      # SAS date: severe stage (optional)
VAR_DATE_REVIEW   <- "REVIEW_DATE"      # SAS date: most recent review (controls anchor)

# --- Neuropsych file variables (NPD_FILE) ---
VAR_NPD_ID        <- "shareid"          # must match VAR_ID in diagnostic file
VAR_EXAM_DATE     <- "examdate"         # SAS date of neuropsychological exam
VAR_EDUCATION     <- "educg"            # education group: 0=<HS, 1=HS, 2=some college, 3=college+
VAR_SEX           <- "sex"             # sex: 1=male, 2=female
VAR_AGE           <- "age"             # age at exam (years)
VAR_BATTERY       <- "battery"         # test battery version (1 or 2); use 1 if only one battery

# ---- Neuropsychological test variable names ----------------------------------
# List the column names for each cognitive test in your NPD file.
# Assign NA (without quotes) for any test your data does not include.
# The script will exclude NA tests from domain composites automatically.

TEST_MEMORY_IMM   <- "LMi"      # immediate verbal memory (e.g., Logical Memory immediate)
TEST_MEMORY_DEL   <- "LMd"      # delayed verbal memory (e.g., Logical Memory delayed)
TEST_VISSP_IMM    <- "VRi"      # immediate visuospatial memory (e.g., Visual Reproduction imm.)
TEST_VISSP_DEL    <- "VRd"      # delayed visuospatial memory (e.g., Visual Reproduction del.)
TEST_LANGUAGE     <- "BNT36"    # confrontation naming (e.g., Boston Naming Test)
TEST_SPEED_A      <- "trailsA"  # processing speed (e.g., Trail Making Test Part A)
TEST_EXEC_B       <- "trailsB"  # executive function (e.g., Trail Making Test Part B)
TEST_WKMEM_FWD    <- "DSF"      # working memory forward (e.g., Digit Span Forward)
TEST_WKMEM_BWD    <- "DSB"      # working memory backward (e.g., Digit Span Backward)
TEST_FLUENCY      <- "FAS"      # verbal fluency — phonemic (e.g., FAS, COWAT)
TEST_VISORGAN     <- "HVOT"     # visual organization (e.g., Hooper Visual Organization Test)

# Tests where HIGHER scores = WORSE performance (will be sign-inverted for z-scoring)
# Edit this vector to match the tests in your battery where this applies.
TESTS_INVERT <- c("trailsA", "trailsB")   # Trail Making: higher time = worse

# ---- Analysis window --------------------------------------------------------
# Keep visits within this window (years) relative to the anchor date.
WINDOW_MIN <- -15   # earliest visit to include (e.g., 15 years before anchor)
WINDOW_MAX <-   5   # latest visit to include  (e.g., 5 years after anchor)

# ---- SAS date origin --------------------------------------------------------
# SAS stores dates as integers since this origin. Change only if your
# date variables use a different epoch (e.g., "1970-01-01" for Unix time).
SAS_ORIGIN <- "1960-01-01"


# =============================================================================
# SECTION 1: VARIABLE DOCUMENTATION
# =============================================================================
#
# This toolkit expects two longitudinal datasets with the following structure:
#
# DIAGNOSTIC / CLINICAL REVIEW FILE (one row per participant per review visit):
#   [VAR_ID]          — Participant ID (links to NPD file)
#   [VAR_CASE_STATUS] — Case status indicator
#                         0        = cognitively normal control
#                         1–4      = case (severity stages 1–4)
#                         5+/other = excluded (other diagnoses, ambiguous)
#   [VAR_CDR]         — CDR score at review (optional; used in Table 1 only)
#   [VAR_DATE_IMPAIR] — Date of first documented cognitive impairment (cases)
#   [VAR_DATE_NORMAL] — Date of last normal review (controls use as anchor)
#   [VAR_DATE_REVIEW] — Date of most recent review
#   [VAR_DATE_MILD, VAR_DATE_MODERATE, VAR_DATE_SEVERE] — Stage transition dates (optional)
#
# NEUROPSYCHOLOGICAL TESTING FILE (one row per participant per exam):
#   [VAR_NPD_ID]      — Participant ID (links to diagnostic file)
#   [VAR_EXAM_DATE]   — Date of neuropsychological exam
#   [VAR_EDUCATION]   — Education level (ordinal 0–3)
#   [VAR_SEX]         — Sex (1=male, 2=female)
#   [VAR_AGE]         — Age at exam (years)
#   [VAR_BATTERY]     — Battery version (1 or 2; set all to 1 if only one battery)
#   [TEST_* variables]— Raw scores for each neuropsychological test
#
# TIME ANCHOR:
#   Cases:    time = 0 at VAR_DATE_IMPAIR (estimated cognitive impairment onset)
#   Controls: time = 0 at VAR_DATE_REVIEW (most recent clinical review)
#   Visits prior to anchor have negative years_to_anchor (pre-diagnostic period).
# =============================================================================


# =============================================================================
# SECTION 2: SETUP — PACKAGES AND HELPERS
# =============================================================================

packages <- c(
  "tidyverse", "lme4", "lmerTest", "broom.mixed",
  "emmeans", "ggplot2", "ggpubr", "patchwork",
  "gt", "gtsummary", "flextable",
  "splines", "nlme", "scales", "RColorBrewer",
  "viridis", "cowplot", "officer", "ggstance"
)
to_install <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

# Helper: convert integer date (e.g., SAS) to R Date
int_to_date <- function(x, origin = SAS_ORIGIN) as.Date(x, origin = origin)

# Helper: read tab-delimited data file, skipping SKIP_ROWS comment rows
read_data <- function(path, skip = SKIP_ROWS) {
  read.delim(path, skip = skip, header = TRUE,
             sep = "\t", na.strings = c("", "NA"),
             stringsAsFactors = FALSE, check.names = FALSE)
}

# Resolve test variable vector, dropping any NAs (tests not in dataset)
resolve_tests <- function(...) {
  vals <- c(...)
  vals[!is.na(vals)]
}

# Build domain composite: row means of available z-scored tests
make_composite <- function(data, z_vars) {
  present <- intersect(paste0("z_", z_vars), names(data))
  if (length(present) == 0) return(rep(NA_real_, nrow(data)))
  rowMeans(data[present], na.rm = TRUE)
}


# =============================================================================
# SECTION 3: LOAD DATA
# =============================================================================
message("Loading data...")

diag_raw <- read_data(DIAG_FILE)
npd_raw  <- read_data(NPD_FILE)

message(sprintf("Diagnostic file: %d rows, %d columns", nrow(diag_raw), ncol(diag_raw)))
message(sprintf("NPD file:        %d rows, %d columns", nrow(npd_raw),  ncol(npd_raw)))


# =============================================================================
# SECTION 4: CLEAN AND DERIVE VARIABLES
# =============================================================================
message("Cleaning and deriving variables...")

## 4a. Diagnostic file: case status, CDR, event dates -----------------------
diag <- diag_raw %>%
  rename(
    shareid      = all_of(VAR_ID),
    case_raw     = all_of(VAR_CASE_STATUS),
    cdr          = all_of(VAR_CDR),
    d_impairment = all_of(VAR_DATE_IMPAIR),
    d_normal     = all_of(VAR_DATE_NORMAL),
    d_review     = all_of(VAR_DATE_REVIEW)
  ) %>%
  mutate(
    # case_status: 0 = control, 1 = case, NA = exclude
    case_status = case_when(
      case_raw == 0        ~ 0L,
      case_raw %in% 1:4   ~ 1L,
      TRUE                 ~ NA_integer_
    ),
    case_label  = factor(case_status, levels = c(0, 1),
                         labels = c("Control", "Case")),
    # Convert integer dates to R Dates
    d_impairment = int_to_date(d_impairment),
    d_normal     = int_to_date(d_normal),
    d_review     = int_to_date(d_review)
  )

# Optionally rename stage-transition date columns if they exist
if (VAR_DATE_MILD     %in% names(diag_raw)) diag$d_mild     <- int_to_date(diag_raw[[VAR_DATE_MILD]])
if (VAR_DATE_MODERATE %in% names(diag_raw)) diag$d_moderate <- int_to_date(diag_raw[[VAR_DATE_MODERATE]])
if (VAR_DATE_SEVERE   %in% names(diag_raw)) diag$d_severe   <- int_to_date(diag_raw[[VAR_DATE_SEVERE]])

## 4b. NPD file: demographics, test scores -----------------------------------
npd <- npd_raw %>%
  rename(
    shareid    = all_of(VAR_NPD_ID),
    exam_date  = all_of(VAR_EXAM_DATE),
    educg      = all_of(VAR_EDUCATION),
    sex        = all_of(VAR_SEX),
    age        = all_of(VAR_AGE),
    battery    = all_of(VAR_BATTERY)
  ) %>%
  mutate(
    exam_date  = int_to_date(exam_date),
    educ_group = factor(educg,
                        levels = 0:3,
                        labels = c("<HS", "HS grad", "Some college", "College+")),
    sex_label  = factor(sex, levels = c(1, 2), labels = c("Male", "Female"))
  )

## 4c. Merge files -----------------------------------------------------------
merged_raw <- npd %>%
  left_join(
    diag %>% select(shareid, case_status, case_label, cdr,
                    d_impairment, d_normal, d_review),
    by = "shareid"
  ) %>%
  filter(!is.na(case_status))   # retain only participants with definite status

message(sprintf("After merge: %d rows, %d unique participants",
                nrow(merged_raw), n_distinct(merged_raw$shareid)))

## 4d. Time-to-anchor variable -----------------------------------------------
# For cases:    anchor = estimated impairment date → captures pre-diagnostic decline
# For controls: anchor = most recent review date
# Positive values = after anchor; negative = before anchor
merged_raw <- merged_raw %>%
  mutate(
    anchor_date      = if_else(case_status == 1L, d_impairment, d_review),
    years_to_anchor  = as.numeric(exam_date - anchor_date) / 365.25,
    years_from_normal= as.numeric(exam_date - d_normal)    / 365.25  # sensitivity
  )

## 4e. Build test variable list and z-score all tests ------------------------
# Collect all test column names that are present in data
ALL_TESTS <- resolve_tests(
  TEST_MEMORY_IMM, TEST_MEMORY_DEL,
  TEST_VISSP_IMM,  TEST_VISSP_DEL,
  TEST_LANGUAGE,   TEST_SPEED_A,    TEST_EXEC_B,
  TEST_WKMEM_FWD,  TEST_WKMEM_BWD,
  TEST_FLUENCY,    TEST_VISORGAN
)
ALL_TESTS <- intersect(ALL_TESTS, names(merged_raw))  # keep only those present

# Normative stats from controls at their earliest available visit
controls_baseline <- merged_raw %>%
  filter(case_status == 0) %>%
  arrange(shareid, exam_date) %>%
  group_by(shareid) %>% slice(1) %>% ungroup()

norm_stats <- controls_baseline %>%
  summarise(across(
    all_of(ALL_TESTS),
    list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
    .names = "{.col}__{.fn}"
  ))

# Compute z-scores; invert sign for tests where higher = worse
for (tst in ALL_TESTS) {
  m     <- norm_stats[[paste0(tst, "__mean")]]
  s     <- norm_stats[[paste0(tst, "__sd")]]
  zname <- paste0("z_", tst)
  if (tst %in% TESTS_INVERT) {
    merged_raw[[zname]] <- -(merged_raw[[tst]] - m) / s
  } else {
    merged_raw[[zname]] <-  (merged_raw[[tst]] - m) / s
  }
}

## 4f. Domain composite scores -----------------------------------------------
# Each composite is the row mean of its constituent z-scored tests.
# Tests not present in the data are automatically excluded.
merged_raw <- merged_raw %>%
  mutate(
    dom_memory     = make_composite(., c(TEST_MEMORY_IMM, TEST_MEMORY_DEL,
                                         TEST_VISSP_IMM,  TEST_VISSP_DEL)),
    dom_executive  = make_composite(., c(TEST_EXEC_B, TEST_WKMEM_BWD, TEST_FLUENCY)),
    dom_language   = make_composite(., c(TEST_LANGUAGE, TEST_VISORGAN)),
    dom_processing = make_composite(., c(TEST_SPEED_A, TEST_WKMEM_FWD)),
    dom_global     = rowMeans(pick(dom_memory, dom_executive,
                                   dom_language, dom_processing), na.rm = TRUE)
  )

DOMAINS       <- c("dom_memory", "dom_executive", "dom_language",
                   "dom_processing", "dom_global")
DOMAIN_LABELS <- c("Memory", "Executive Function",
                   "Language/Visuospatial", "Processing Speed", "Global Composite")

## 4g. Restrict to analytical window -----------------------------------------
analytic <- merged_raw %>%
  filter(
    !is.na(years_to_anchor),
    !is.na(age),
    years_to_anchor >= WINDOW_MIN,
    years_to_anchor <= WINDOW_MAX
  )

message(sprintf(
  "Analytic sample: %d visits, %d participants (%d cases, %d controls)",
  nrow(analytic),
  n_distinct(analytic$shareid),
  n_distinct(analytic$shareid[analytic$case_status == 1]),
  n_distinct(analytic$shareid[analytic$case_status == 0])
))

## 4h. Baseline record per participant (for Table 1) -------------------------
baseline <- analytic %>%
  arrange(shareid, exam_date) %>%
  group_by(shareid) %>% slice(1) %>% ungroup()

## 4i. Education binary moderator --------------------------------------------
analytic <- analytic %>%
  mutate(
    educ_binary = case_when(
      educg %in% c(0, 1) ~ "Low (\u2264HS)",
      educg %in% c(2, 3) ~ "High (\u2265some college)",
      TRUE ~ NA_character_
    ),
    educ_binary = factor(educ_binary,
                         levels = c("Low (\u2264HS)", "High (\u2265some college)"))
  )

baseline <- baseline %>%
  mutate(
    educ_binary = case_when(
      educg %in% c(0, 1) ~ "Low (\u2264HS)",
      educg %in% c(2, 3) ~ "High (\u2265some college)",
      TRUE ~ NA_character_
    ),
    educ_binary = factor(educ_binary,
                         levels = c("Low (\u2264HS)", "High (\u2265some college)"))
  )


# =============================================================================
# SECTION 5: TABLE 1 — SAMPLE CHARACTERISTICS
# =============================================================================
message("Building Table 1...")

visit_summary <- analytic %>%
  group_by(shareid) %>%
  summarise(
    n_visits   = n(),
    followup_y = max(years_to_anchor) - min(years_to_anchor),
    .groups = "drop"
  )

tbl1_data <- baseline %>%
  left_join(visit_summary, by = "shareid") %>%
  mutate(
    `Age at baseline (years)` = age,
    `Sex (Female, %)`         = (sex == 2),
    `Education group`         = educ_group,
    `Visits per participant`  = n_visits,
    `Years of follow-up`      = followup_y,
    `CDR at baseline`         = cdr
  )

# Select available z-score columns for Table 1 display
z_cols_available <- intersect(
  c(paste0("z_", TEST_MEMORY_IMM), paste0("z_", TEST_VISSP_IMM),
    paste0("z_", TEST_LANGUAGE),   paste0("z_", TEST_SPEED_A),
    paste0("z_", TEST_EXEC_B),     paste0("z_", TEST_WKMEM_BWD)),
  names(tbl1_data)
)

tbl1 <- tbl1_data %>%
  select(case_label,
         `Age at baseline (years)`, `Sex (Female, %)`, `Education group`,
         CDR = cdr, `Visits per participant`, `Years of follow-up`,
         all_of(z_cols_available)) %>%
  tbl_summary(
    by = case_label,
    statistic = list(all_continuous()  ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p}%)"),
    digits  = all_continuous() ~ 2,
    missing = "no"
  ) %>%
  add_p() %>%
  add_overall() %>%
  bold_labels() %>%
  modify_caption("**Table 1. Baseline Characteristics by Case Status**") %>%
  modify_footnote(
    update = everything() ~
      "Continuous: mean (SD); Categorical: n (%). z-scores standardized to control baseline mean/SD. Inverted test z-scores: higher = better performance."
  )

tbl1_flex <- tbl1 %>%
  as_flex_table() %>%
  flextable::fontsize(size = 11, part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::autofit()

officer::read_docx() %>%
  flextable::body_add_flextable(tbl1_flex) %>%
  print(target = "Table1_Sample_Characteristics.docx")
message("Table 1 saved.")


# =============================================================================
# SECTION 6: PRIMARY ANALYSIS — DOMAIN-SPECIFIC TRAJECTORY MODELS
# =============================================================================
# Model:
#   domain_score ~ years_to_anchor * case_status + age + sex + battery
#                  + (1 + years_to_anchor | participant_id)
#
#   - years_to_anchor: continuous time variable (0 = clinical anchor date)
#   - case_status:     0 = control, 1 = case
#   - Interaction term (years × case): tests whether SLOPES differ between groups
#   - Random slope + intercept per participant: captures individual variability
#   - Battery covariate: adjusts for different test instruments across waves
# =============================================================================
message("Fitting primary LME trajectory models...")

fit_trajectory_model <- function(domain_var, data) {
  formula_str <- paste0(
    domain_var,
    " ~ years_to_anchor * case_status + age + sex + battery +",
    " (1 + years_to_anchor | shareid)"
  )
  df_clean <- data %>%
    filter(!is.na(.data[[domain_var]]),
           !is.na(years_to_anchor), !is.na(age),
           !is.na(sex), !is.na(battery)) %>%
    mutate(case_status = factor(case_status, levels = c(0, 1)))

  tryCatch(
    lmerTest::lmer(as.formula(formula_str), data = df_clean,
                   REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa",
                                         optCtrl  = list(maxfun = 2e5))),
    error = function(e) {
      message("  Model failed for ", domain_var, ": ", conditionMessage(e)); NULL
    }
  )
}

traj_models <- setNames(
  lapply(DOMAINS, fit_trajectory_model, data = analytic),
  DOMAINS
)

# Extract key fixed-effect terms from each model
extract_lme_table <- function(model, domain_label) {
  if (is.null(model)) return(NULL)
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(domain = domain_label) %>%
    filter(stringr::str_detect(term, "years_to_anchor|case_status"))
}

traj_results <- bind_rows(
  mapply(extract_lme_table, traj_models, DOMAIN_LABELS, SIMPLIFY = FALSE)
)

# Format for Table 2
primary_tbl <- traj_results %>%
  mutate(
    term_clean = case_when(
      term == "years_to_anchor"                ~ "Slope (years from anchor)",
      term == "case_status1"                   ~ "Case vs Control (intercept)",
      term == "years_to_anchor:case_status1"   ~ "Slope \u00d7 Case interaction",
      TRUE ~ term
    ),
    p_label = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    sig = case_when(
      p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
      p.value < 0.05  ~ "*",   p.value < 0.1  ~ "\u2020",
      TRUE ~ ""
    )
  ) %>%
  select(Domain = domain, Term = term_clean,
         Estimate = estimate, SE = std.error,
         `95% CI Low` = conf.low, `95% CI High` = conf.high,
         `p-value` = p_label, `Sig.` = sig)

ft_primary <- primary_tbl %>%
  flextable::flextable() %>%
  flextable::bold(part = "header") %>%
  flextable::colformat_double(j = c("Estimate","SE","95% CI Low","95% CI High"), digits = 3) %>%
  flextable::theme_booktabs() %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::autofit() %>%
  flextable::set_caption(
    "Table 2. Fixed effects from linear mixed-effects models of domain-specific cognitive trajectories"
  ) %>%
  flextable::footnote(
    i = 1, j = 1, part = "header",
    value = flextable::as_paragraph(
      "Models include random intercept and slope per participant. Time anchored to clinical event date (cases) or last review date (controls). Battery, age, and sex included as covariates. * p<0.05, ** p<0.01, *** p<0.001, \u2020 p<0.10."
    ),
    ref_symbols = "a"
  )

officer::read_docx() %>%
  flextable::body_add_flextable(ft_primary) %>%
  print(target = "Table2_LME_Primary_Results.docx")
message("Table 2 saved.")


# =============================================================================
# SECTION 7: FIGURE 1 — DOMAIN TRAJECTORIES (CASES vs CONTROLS)
# =============================================================================
message("Building Figure 1: Domain trajectories...")

time_grid <- seq(WINDOW_MIN + 3, WINDOW_MAX - 1, by = 0.5)

predict_trajectories <- function(model, domain_var, domain_label, data) {
  if (is.null(model)) return(NULL)
  df_clean <- data %>%
    filter(!is.na(.data[[domain_var]]), !is.na(years_to_anchor),
           !is.na(age), !is.na(sex), !is.na(battery)) %>%
    mutate(case_status = factor(case_status, levels = c(0, 1)))

  ref_age <- median(df_clean$age, na.rm = TRUE)
  ref_sex <- as.numeric(names(sort(table(df_clean$sex), decreasing = TRUE))[1])
  ref_bat <- as.numeric(names(sort(table(df_clean$battery), decreasing = TRUE))[1])

  pred_df <- expand.grid(
    years_to_anchor = time_grid,
    case_status     = factor(c(0, 1), levels = c(0, 1)),
    age = ref_age, sex = ref_sex, battery = ref_bat, shareid = NA
  )
  pred_df$fit       <- predict(model, newdata = pred_df, re.form = NA)
  pred_df$domain    <- domain_label
  pred_df$case_label <- factor(pred_df$case_status, levels = c(0, 1),
                                labels = c("Control", "Case"))
  pred_df
}

traj_preds <- bind_rows(
  mapply(predict_trajectories, traj_models, DOMAINS, DOMAIN_LABELS,
         MoreArgs = list(data = analytic), SIMPLIFY = FALSE)
) %>% filter(!is.na(fit))

# Observed yearly means (for overlaying on model lines)
obs_means <- analytic %>%
  mutate(year_bin  = round(years_to_anchor),
         case_label = factor(case_status, levels = c(0,1),
                             labels = c("Control","Case"))) %>%
  filter(year_bin >= WINDOW_MIN, year_bin <= WINDOW_MAX) %>%
  pivot_longer(cols = all_of(DOMAINS), names_to = "domain_var", values_to = "score") %>%
  mutate(domain = DOMAIN_LABELS[match(domain_var, DOMAINS)]) %>%
  group_by(domain, year_bin, case_label) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    se_score   = sd(score,   na.rm = TRUE) / sqrt(sum(!is.na(score))),
    n          = sum(!is.na(score)),
    .groups    = "drop"
  ) %>%
  filter(n >= 5)

GROUP_COLORS <- c("Control" = "#2166AC", "Case" = "#D6604D")

traj_preds$domain <- factor(traj_preds$domain, levels = DOMAIN_LABELS)
obs_means$domain  <- factor(obs_means$domain,  levels = DOMAIN_LABELS)

fig1 <- ggplot() +
  annotate("rect", xmin = WINDOW_MIN, xmax = 0,
           ymin = -Inf, ymax = Inf, alpha = 0.05, fill = "#D6604D") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_errorbar(data = obs_means,
                aes(x = year_bin,
                    ymin = mean_score - 1.96 * se_score,
                    ymax = mean_score + 1.96 * se_score,
                    color = case_label),
                width = 0.25, linewidth = 0.5, alpha = 0.6) +
  geom_point(data = obs_means,
             aes(x = year_bin, y = mean_score, color = case_label),
             size = 2.2, alpha = 0.85) +
  geom_line(data = traj_preds,
            aes(x = years_to_anchor, y = fit, color = case_label),
            linewidth = 1.1) +
  facet_wrap(~ domain, ncol = 3, scales = "free_y",
             labeller = label_wrap_gen(width = 22)) +
  scale_color_manual(values = GROUP_COLORS, name = NULL) +
  scale_x_continuous(
    breaks = seq(WINDOW_MIN, WINDOW_MAX, by = 4),
    name   = "Years relative to clinical anchor (cases: impairment onset; controls: last review)"
  ) +
  scale_y_continuous(name = "Cognitive score (z-score; higher = better)") +
  labs(
    title    = "Figure 1. Domain-specific cognitive trajectories by case status",
    subtitle = "LME model predictions (lines) \u00b1 observed group means (points \u00b1 95% CI). Shaded area = pre-event period.",
    caption  = "Time = 0: clinical anchor date. Scores standardized to control baseline mean/SD."
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    plot.caption  = element_text(size = 9,  color = "grey40"),
    strip.text    = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    axis.title    = element_text(size = 11),
    axis.text     = element_text(size = 10),
    legend.position  = "bottom",
    legend.text      = element_text(size = 12),
    legend.key.size  = unit(1.2, "lines"),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(1.2, "lines")
  )

ggsave("Figure1_Domain_Trajectories.pdf", fig1, width = 14, height = 10, device = cairo_pdf)
ggsave("Figure1_Domain_Trajectories.png", fig1, width = 14, height = 10, dpi = 300)
message("Figure 1 saved.")


# =============================================================================
# SECTION 8: FIGURE 2 — EFFECT SIZE DIVERGENCE TIMELINE (COHEN'S d)
# =============================================================================
message("Building Figure 2: Divergence timeline...")

# For each domain × year bin, compute Cohen's d (cases vs controls)
divergence_df <- analytic %>%
  mutate(year_bin   = floor(years_to_anchor),
         case_label = factor(case_status, levels = c(0,1),
                             labels = c("Control","Case"))) %>%
  filter(year_bin >= WINDOW_MIN, year_bin <= WINDOW_MAX) %>%
  pivot_longer(cols = all_of(DOMAINS[1:4]),
               names_to = "domain_var", values_to = "score") %>%
  mutate(domain = DOMAIN_LABELS[match(domain_var, DOMAINS)]) %>%
  filter(!is.na(score))

cohens_d <- divergence_df %>%
  group_by(domain, year_bin) %>%
  summarise(
    n_case    = sum(case_label == "Case"),
    n_ctrl    = sum(case_label == "Control"),
    mean_case = mean(score[case_label == "Case"],    na.rm = TRUE),
    mean_ctrl = mean(score[case_label == "Control"], na.rm = TRUE),
    sd_case   = sd(score[case_label == "Case"],      na.rm = TRUE),
    sd_ctrl   = sd(score[case_label == "Control"],   na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    pooled_sd = sqrt(((n_case - 1) * sd_case^2 + (n_ctrl - 1) * sd_ctrl^2) /
                       (n_case + n_ctrl - 2)),
    cohens_d  = (mean_ctrl - mean_case) / pooled_sd,   # positive = cases worse
    se_d      = sqrt((n_case + n_ctrl) / (n_case * n_ctrl) +
                       cohens_d^2 / (2 * (n_case + n_ctrl))),
    ci_low    = cohens_d - 1.96 * se_d,
    ci_high   = cohens_d + 1.96 * se_d,
    sufficient = (n_case >= 10 & n_ctrl >= 5)
  ) %>%
  filter(sufficient)

domain_colors <- c(
  "Memory"                = "#1B7837",
  "Executive Function"    = "#762A83",
  "Language/Visuospatial" = "#E08214",
  "Processing Speed"      = "#2166AC"
)

cohens_d$domain <- factor(cohens_d$domain,
                           levels = c("Memory","Executive Function",
                                      "Language/Visuospatial","Processing Speed"))

fig2 <- ggplot(cohens_d,
               aes(x = year_bin, y = cohens_d, color = domain, fill = domain)) +
  annotate("rect", xmin = WINDOW_MIN, xmax = 0,
           ymin = -Inf, ymax = Inf, alpha = 0.04, fill = "tomato") +
  geom_hline(yintercept = 0,             linetype = "solid",  color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = c(0.2, 0.5, 0.8),
             linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  annotate("text", x = WINDOW_MIN + 0.5, y = c(0.21, 0.51, 0.81),
           label = c("Small (d=0.2)","Medium (d=0.5)","Large (d=0.8)"),
           hjust = 0, size = 3.2, color = "grey50") +
  scale_color_manual(values = domain_colors, name = "Cognitive domain") +
  scale_fill_manual( values = domain_colors, name = "Cognitive domain") +
  scale_x_continuous(breaks = seq(WINDOW_MIN, WINDOW_MAX, by = 2),
                     name   = "Years relative to clinical anchor (cases)") +
  scale_y_continuous(name   = "Cohen's d (cases vs controls; positive = cases worse)",
                     limits = c(-0.2, NA)) +
  labs(
    title    = "Figure 2. Temporal sequence of cognitive domain divergence",
    subtitle = "Cohen's d per year bin. Shaded bands = 95% CI.",
    caption  = "Bins with \u226510 cases and \u22655 controls shown. Red shading = pre-event period."
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    plot.caption  = element_text(size = 9,  color = "grey40"),
    axis.title    = element_text(size = 12),
    axis.text     = element_text(size = 11),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(1.2, "lines"),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))

ggsave("Figure2_Divergence_Timeline.pdf", fig2, width = 12, height = 7, device = cairo_pdf)
ggsave("Figure2_Divergence_Timeline.png", fig2, width = 12, height = 7, dpi = 300)
message("Figure 2 saved.")


# =============================================================================
# SECTION 9: FIGURE 3 — INDIVIDUAL TRAJECTORY SPAGHETTI PLOTS
# =============================================================================
message("Building Figure 3: Spaghetti plots...")

set.seed(42)
case_ids    <- sample(unique(analytic$shareid[analytic$case_status == 1]),
                      min(60, sum(analytic$case_status == 1)))
control_ids <- sample(unique(analytic$shareid[analytic$case_status == 0]),
                      min(30, n_distinct(analytic$shareid[analytic$case_status == 0])))

spag_long <- analytic %>%
  filter(shareid %in% c(case_ids, control_ids)) %>%
  mutate(case_label = factor(case_status, levels = c(0,1),
                             labels = c("Control","Case"))) %>%
  select(shareid, case_label, years_to_anchor,
         dom_memory, dom_executive, dom_language, dom_processing) %>%
  pivot_longer(cols = starts_with("dom_"), names_to = "domain_var", values_to = "score") %>%
  mutate(domain = DOMAIN_LABELS[match(domain_var, DOMAINS)]) %>%
  filter(!is.na(score), !is.na(years_to_anchor))

spag_long$domain <- factor(spag_long$domain, levels = DOMAIN_LABELS[1:4])

fig3 <- ggplot(spag_long,
               aes(x = years_to_anchor, y = score, group = shareid)) +
  annotate("rect", xmin = WINDOW_MIN, xmax = 0,
           ymin = -Inf, ymax = Inf, alpha = 0.06, fill = "#D6604D") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "solid",  color = "grey50", linewidth = 0.4) +
  geom_line(aes(color = case_label), alpha = 0.30, linewidth = 0.4) +
  stat_smooth(aes(group = case_label, color = case_label),
              method = "lm", formula = y ~ x, se = TRUE,
              linewidth = 1.5, alpha = 0.15) +
  facet_grid(case_label ~ domain, labeller = label_wrap_gen(width = 20)) +
  scale_color_manual(values = GROUP_COLORS, name = NULL) +
  scale_x_continuous(breaks = seq(WINDOW_MIN, WINDOW_MAX, by = 4),
                     name   = "Years relative to anchor") +
  scale_y_continuous(name  = "Composite domain score (z)") +
  labs(
    title    = "Figure 3. Individual cognitive trajectories by domain and case status",
    subtitle = "Thin lines = individual participants (random sample); thick line = group linear trend \u00b1 95% CI.",
    caption  = "Red shading = pre-event period (before time = 0)."
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    plot.caption  = element_text(size = 9,  color = "grey40"),
    strip.text    = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    axis.title    = element_text(size = 11),
    axis.text     = element_text(size = 10),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(1.0, "lines")
  )

ggsave("Figure3_Spaghetti_Plots.pdf", fig3, width = 14, height = 9, device = cairo_pdf)
ggsave("Figure3_Spaghetti_Plots.png", fig3, width = 14, height = 9, dpi = 300)
message("Figure 3 saved.")


# =============================================================================
# SECTION 10: SECONDARY ANALYSIS — EDUCATION AND SEX MODERATION
# =============================================================================
# Model (for memory and global composite):
#   domain ~ years_to_anchor * case_status * educ_binary
#           + years_to_anchor * case_status * sex_label
#           + age + battery + (1 + years_to_anchor | participant_id)
#
#   Tests whether education or sex moderates the rate of case-related decline.
# =============================================================================
message("Running secondary moderation analyses...")

analytic_sec <- analytic %>% filter(!is.na(educ_binary), !is.na(sex_label))

fit_reserve_model <- function(domain_var, data) {
  formula_str <- paste0(
    domain_var,
    " ~ years_to_anchor * case_status * educ_binary +",
    " years_to_anchor * case_status * sex_label +",
    " age + battery + (1 + years_to_anchor | shareid)"
  )
  df_clean <- data %>%
    filter(!is.na(.data[[domain_var]]),
           !is.na(years_to_anchor), !is.na(age), !is.na(battery)) %>%
    mutate(case_status = factor(case_status, levels = c(0, 1)))

  tryCatch(
    lmerTest::lmer(as.formula(formula_str), data = df_clean,
                   REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa",
                                         optCtrl  = list(maxfun = 2e5))),
    error = function(e) {
      message("  Moderation model failed for ", domain_var, ": ", conditionMessage(e)); NULL
    }
  )
}

reserve_domains <- c("dom_memory", "dom_global")
reserve_labels  <- c("Memory", "Global Composite")
reserve_models  <- setNames(
  lapply(reserve_domains, fit_reserve_model, data = analytic_sec),
  reserve_domains
)

extract_reserve_table <- function(model, domain_label) {
  if (is.null(model)) return(NULL)
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(domain = domain_label) %>%
    filter(stringr::str_detect(term, "educ_binary|sex_label|years_to_anchor:case_status"))
}

reserve_results <- bind_rows(
  mapply(extract_reserve_table, reserve_models, reserve_labels, SIMPLIFY = FALSE)
)

reserve_tbl <- reserve_results %>%
  mutate(
    term_clean = term %>%
      str_replace("years_to_anchor", "Slope") %>%
      str_replace("case_status1", "Case") %>%
      str_replace("educ_binaryHigh.*", "High education") %>%
      str_replace("sex_labelFemale", "Female"),
    p_label = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    sig = case_when(
      p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
      p.value < 0.05  ~ "*",   p.value < 0.1  ~ "\u2020",
      TRUE ~ ""
    )
  ) %>%
  select(Domain = domain, Term = term_clean,
         Estimate = estimate, SE = std.error,
         `95% CI Low` = conf.low, `95% CI High` = conf.high,
         `p-value` = p_label, `Sig.` = sig)

ft_reserve <- reserve_tbl %>%
  flextable::flextable() %>%
  flextable::bold(part = "header") %>%
  flextable::colformat_double(j = c("Estimate","SE","95% CI Low","95% CI High"), digits = 3) %>%
  flextable::theme_booktabs() %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::autofit() %>%
  flextable::set_caption(
    "Table 3. Secondary analysis: education and sex as moderators of case-related cognitive decline"
  ) %>%
  flextable::footnote(
    i = 1, j = 1, part = "header",
    value = flextable::as_paragraph(
      "Three-way interaction model: Slope \u00d7 Case \u00d7 education and Slope \u00d7 Case \u00d7 sex. Interaction-relevant terms only. * p<0.05, ** p<0.01, *** p<0.001, \u2020 p<0.10."
    ),
    ref_symbols = "b"
  )

officer::read_docx() %>%
  flextable::body_add_flextable(ft_reserve) %>%
  print(target = "Table3_Moderation_Results.docx")
message("Table 3 saved.")


# =============================================================================
# SECTION 11: FIGURE 4 — EDUCATION × CASE TRAJECTORY PLOT
# =============================================================================
message("Building Figure 4: Education moderation plots...")

predict_reserve <- function(model, domain_var, domain_label, data) {
  if (is.null(model)) return(NULL)
  df_clean <- data %>%
    filter(!is.na(.data[[domain_var]]), !is.na(years_to_anchor),
           !is.na(age), !is.na(battery), !is.na(educ_binary), !is.na(sex_label)) %>%
    mutate(case_status = factor(case_status, levels = c(0, 1)))

  ref_age <- median(df_clean$age, na.rm = TRUE)
  ref_bat <- as.numeric(names(sort(table(df_clean$battery), decreasing = TRUE))[1])

  pred_df <- expand.grid(
    years_to_anchor = time_grid,
    case_status     = factor(c(0, 1), levels = c(0, 1)),
    educ_binary     = levels(df_clean$educ_binary),
    sex_label       = factor("Female", levels = c("Male","Female")),
    age = ref_age, battery = ref_bat, shareid = NA
  )
  pred_df$fit        <- predict(model, newdata = pred_df, re.form = NA)
  pred_df$domain     <- domain_label
  pred_df$case_label <- factor(pred_df$case_status, levels = c(0,1),
                                labels = c("Control","Case"))
  pred_df
}

reserve_preds <- bind_rows(
  mapply(predict_reserve, reserve_models, reserve_domains, reserve_labels,
         MoreArgs = list(data = analytic_sec), SIMPLIFY = FALSE)
) %>% filter(!is.na(fit))

educ_colors <- c("Low (\u2264HS)" = "#E08214", "High (\u2265some college)" = "#542788")

fig4 <- ggplot(reserve_preds,
               aes(x = years_to_anchor, y = fit,
                   color = educ_binary, linetype = case_label)) +
  annotate("rect", xmin = WINDOW_MIN, xmax = 0,
           ymin = -Inf, ymax = Inf, alpha = 0.05, fill = "tomato") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "solid",  color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 1.15) +
  facet_wrap(~ domain, ncol = 2, scales = "free_y") +
  scale_color_manual(values = educ_colors, name = "Education") +
  scale_linetype_manual(values = c("Control" = "dashed", "Case" = "solid"), name = "Group") +
  scale_x_continuous(breaks = seq(WINDOW_MIN, WINDOW_MAX, by = 4),
                     name   = "Years relative to anchor") +
  scale_y_continuous(name = "Composite domain score (z)") +
  labs(
    title    = "Figure 4. Education as a moderator of case-related cognitive decline",
    subtitle = "LME predictions at median age. Solid = cases; dashed = controls.",
    caption  = "Red shading = pre-event period. Education: Low = <HS or HS grad; High = some college or college+"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    plot.caption  = element_text(size = 9,  color = "grey40"),
    strip.text    = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    axis.title    = element_text(size = 12),
    axis.text     = element_text(size = 11),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(1.5, "lines"),
    panel.grid.minor = element_blank()
  )

ggsave("Figure4_Education_Moderation.pdf", fig4, width = 12, height = 7, device = cairo_pdf)
ggsave("Figure4_Education_Moderation.png", fig4, width = 12, height = 7, dpi = 300)
message("Figure 4 saved.")


# =============================================================================
# SECTION 12: FIGURE 5 — SEX × CASE TRAJECTORY PLOT
# =============================================================================
message("Building Figure 5: Sex moderation plots...")

predict_sex <- function(model, domain_var, domain_label, data) {
  if (is.null(model)) return(NULL)
  df_clean <- data %>%
    filter(!is.na(.data[[domain_var]]), !is.na(years_to_anchor),
           !is.na(age), !is.na(battery), !is.na(educ_binary), !is.na(sex_label)) %>%
    mutate(case_status = factor(case_status, levels = c(0, 1)))

  ref_age  <- median(df_clean$age, na.rm = TRUE)
  ref_bat  <- as.numeric(names(sort(table(df_clean$battery), decreasing = TRUE))[1])
  ref_educ <- names(sort(table(as.character(df_clean$educ_binary)), decreasing = TRUE))[1]

  pred_df <- expand.grid(
    years_to_anchor = time_grid,
    case_status     = factor(c(0, 1), levels = c(0, 1)),
    educ_binary     = factor(ref_educ, levels = levels(df_clean$educ_binary)),
    sex_label       = factor(c("Male","Female"), levels = c("Male","Female")),
    age = ref_age, battery = ref_bat, shareid = NA
  )
  pred_df$fit        <- predict(model, newdata = pred_df, re.form = NA)
  pred_df$domain     <- domain_label
  pred_df$case_label <- factor(pred_df$case_status, levels = c(0,1),
                                labels = c("Control","Case"))
  pred_df
}

sex_preds <- bind_rows(
  mapply(predict_sex, reserve_models, reserve_domains, reserve_labels,
         MoreArgs = list(data = analytic_sec), SIMPLIFY = FALSE)
) %>% filter(!is.na(fit))

sex_colors <- c("Male" = "#2166AC", "Female" = "#D6604D")

fig5 <- ggplot(sex_preds,
               aes(x = years_to_anchor, y = fit,
                   color = sex_label, linetype = case_label)) +
  annotate("rect", xmin = WINDOW_MIN, xmax = 0,
           ymin = -Inf, ymax = Inf, alpha = 0.05, fill = "tomato") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "solid",  color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 1.15) +
  facet_wrap(~ domain, ncol = 2, scales = "free_y") +
  scale_color_manual(values = sex_colors, name = "Sex") +
  scale_linetype_manual(values = c("Control" = "dashed", "Case" = "solid"), name = "Group") +
  scale_x_continuous(breaks = seq(WINDOW_MIN, WINDOW_MAX, by = 4),
                     name   = "Years relative to anchor") +
  scale_y_continuous(name = "Composite domain score (z)") +
  labs(
    title    = "Figure 5. Sex as a moderator of case-related cognitive decline",
    subtitle = "LME predictions at median age and modal education. Solid = cases; dashed = controls.",
    caption  = "Red shading = pre-event period."
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    plot.caption  = element_text(size = 9,  color = "grey40"),
    strip.text    = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    axis.title    = element_text(size = 12),
    axis.text     = element_text(size = 11),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(1.5, "lines"),
    panel.grid.minor = element_blank()
  )

ggsave("Figure5_Sex_Moderation.pdf", fig5, width = 12, height = 7, device = cairo_pdf)
ggsave("Figure5_Sex_Moderation.png", fig5, width = 12, height = 7, dpi = 300)
message("Figure 5 saved.")


# =============================================================================
# SECTION 13: FIGURE 6 — FOREST PLOT OF SLOPE × CASE INTERACTION EFFECTS
# =============================================================================
message("Building Figure 6: Forest plot...")

forest_data <- bind_rows(
  mapply(extract_lme_table, traj_models, DOMAIN_LABELS, SIMPLIFY = FALSE)
) %>%
  filter(term == "years_to_anchor:case_status1") %>%
  mutate(
    domain = factor(domain, levels = rev(DOMAIN_LABELS)),
    sig = case_when(
      p.value < 0.001 ~ "p < 0.001", p.value < 0.01 ~ "p < 0.01",
      p.value < 0.05  ~ "p < 0.05",  TRUE            ~ "p \u2265 0.05"
    ),
    sig = factor(sig, levels = c("p < 0.001","p < 0.01","p < 0.05","p \u2265 0.05"))
  )

sig_colors <- c(
  "p < 0.001" = "#B2182B", "p < 0.01"  = "#D6604D",
  "p < 0.05"  = "#F4A582", "p \u2265 0.05" = "#92C5DE"
)

fig6 <- ggplot(forest_data, aes(x = estimate, y = domain, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3, linewidth = 1.0) +
  geom_point(size = 4.5) +
  geom_text(
    aes(label = sprintf("\u03b2=%.3f\np=%s", estimate,
                        ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)))),
    hjust = -0.18, size = 3.8, color = "grey20", lineheight = 0.9
  ) +
  scale_color_manual(values = sig_colors, name = "Significance") +
  scale_x_continuous(
    name   = "\u03b2 coefficient (Slope \u00d7 Case interaction)\nNegative = steeper decline in cases",
    expand = expansion(mult = c(0.05, 0.35))
  ) +
  scale_y_discrete(name = NULL) +
  labs(
    title    = "Figure 6. Effect sizes for slope \u00d7 case interaction across cognitive domains",
    subtitle = "Negative estimates indicate steeper longitudinal decline in cases vs controls.",
    caption  = "Forest plot: fixed-effect estimate \u00b1 95% CI from LME models."
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title   = element_text(size = 15, face = "bold"),
    plot.subtitle= element_text(size = 11, color = "grey30"),
    plot.caption = element_text(size = 9,  color = "grey40"),
    axis.text.y  = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 11),
    axis.title.x = element_text(size = 11),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggsave("Figure6_Forest_Plot.pdf", fig6, width = 10, height = 6, device = cairo_pdf)
ggsave("Figure6_Forest_Plot.png", fig6, width = 10, height = 6, dpi = 300)
message("Figure 6 saved.")


# =============================================================================
# SECTION 14: SUPPLEMENTARY TABLE S1 — FULL FIXED EFFECTS
# =============================================================================
message("Building Table S1...")

full_fixed <- bind_rows(
  mapply(function(m, lbl) {
    if (is.null(m)) return(NULL)
    broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE) %>% mutate(domain = lbl)
  }, traj_models, DOMAIN_LABELS, SIMPLIFY = FALSE)
) %>%
  mutate(
    p_label = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
    sig = case_when(
      p.value < 0.001 ~ "***", p.value < 0.01 ~ "**",
      p.value < 0.05  ~ "*",   p.value < 0.1  ~ "\u2020",
      TRUE ~ ""
    )
  ) %>%
  select(Domain = domain, Term = term, Estimate = estimate,
         SE = std.error, `95% CI Low` = conf.low, `95% CI High` = conf.high,
         df, `t-value` = statistic, `p-value` = p_label, `Sig.` = sig)

officer::read_docx() %>%
  flextable::body_add_flextable(
    full_fixed %>%
      flextable::flextable() %>%
      flextable::bold(part = "header") %>%
      flextable::colformat_double(
        j = c("Estimate","SE","95% CI Low","95% CI High","t-value"), digits = 3) %>%
      flextable::theme_booktabs() %>%
      flextable::fontsize(size = 9, part = "all") %>%
      flextable::autofit() %>%
      flextable::set_caption("Table S1. Complete fixed effects from all primary trajectory models")
  ) %>%
  print(target = "TableS1_Full_Fixed_Effects.docx")
message("Table S1 saved.")


# =============================================================================
# SECTION 15: SUPPLEMENTARY TABLE S2 — RANDOM EFFECTS & MODEL FIT
# =============================================================================
message("Building Table S2...")

re_summary <- bind_rows(
  mapply(function(m, lbl) {
    if (is.null(m)) return(NULL)
    vcomp <- as.data.frame(VarCorr(m))
    tibble(
      Domain                = lbl,
      `Random intercept SD` = vcomp$sdcor[vcomp$var1 == "(Intercept)" & is.na(vcomp$var2)],
      `Random slope SD`     = vcomp$sdcor[vcomp$var1 == "years_to_anchor" & is.na(vcomp$var2)],
      `Intercept-slope r`   = vcomp$sdcor[!is.na(vcomp$var2)][1],
      `Residual SD`         = sigma(m),
      AIC                   = AIC(m),
      `N observations`      = nobs(m),
      `N participants`      = ngrps(m)["shareid"]
    )
  }, traj_models, DOMAIN_LABELS, SIMPLIFY = FALSE)
)

officer::read_docx() %>%
  flextable::body_add_flextable(
    re_summary %>%
      flextable::flextable() %>%
      flextable::bold(part = "header") %>%
      flextable::colformat_double(digits = 3) %>%
      flextable::theme_booktabs() %>%
      flextable::fontsize(size = 10, part = "all") %>%
      flextable::autofit() %>%
      flextable::set_caption("Table S2. Random effects and model fit statistics")
  ) %>%
  print(target = "TableS2_Random_Effects.docx")
message("Table S2 saved.")


# =============================================================================
# SECTION 16: SENSITIVITY ANALYSIS — BATTERY STRATIFICATION
# =============================================================================
# Verifies that the slope × case interaction is consistent within each
# test battery / instrument version used across different exam waves.
# Uses a simplified model (random intercept only) to ensure convergence
# in potentially smaller within-battery subsets.
# =============================================================================
message("Running sensitivity analysis: battery stratification...")

fit_sensitivity_simple <- function(domain_var, data) {
  formula_str <- paste0(
    domain_var,
    " ~ years_to_anchor * case_status + age + sex + (1 | shareid)"
  )
  df_clean <- data %>%
    filter(!is.na(.data[[domain_var]]), !is.na(years_to_anchor),
           !is.na(age), !is.na(sex)) %>%
    mutate(case_status = factor(case_status, levels = c(0, 1)))

  tryCatch(
    lmerTest::lmer(as.formula(formula_str), data = df_clean,
                   REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa",
                                         optCtrl  = list(maxfun = 2e5))),
    error = function(e) { message("  Failed: ", conditionMessage(e)); NULL }
  )
}

sens_results <- list()
for (bat in sort(unique(analytic$battery))) {
  dat_bat <- analytic %>% filter(battery == bat)
  message(sprintf("Battery %s: %d visits, %d participants",
                  bat, nrow(dat_bat), n_distinct(dat_bat$shareid)))

  for (i in seq_along(DOMAINS)) {
    m <- fit_sensitivity_simple(DOMAINS[i], dat_bat)
    if (!is.null(m)) {
      res <- broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE) %>%
        filter(term == "years_to_anchor:case_status1") %>%
        mutate(domain = DOMAIN_LABELS[i], battery = paste("Battery", bat))
      sens_results <- c(sens_results, list(res))
    }
  }
}

sens_df <- bind_rows(sens_results) %>%
  mutate(domain = factor(domain, levels = DOMAIN_LABELS))

figS1 <- ggplot(sens_df,
                aes(x = estimate, y = domain, color = battery, shape = battery)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.3, linewidth = 0.9,
                 position = position_dodgev(height = 0.5)) +
  geom_point(size = 3.5, position = position_dodgev(height = 0.5)) +
  scale_color_manual(values = c("Battery 1" = "#1F78B4", "Battery 2" = "#33A02C"),
                     name = NULL) +
  scale_shape_manual(values = c("Battery 1" = 16, "Battery 2" = 17), name = NULL) +
  scale_x_continuous(name = "\u03b2 coefficient (Slope \u00d7 Case interaction)") +
  scale_y_discrete(name = NULL) +
  labs(
    title    = "Figure S1. Sensitivity analysis: consistency across test batteries",
    subtitle = "Slope \u00d7 Case interaction estimated separately per battery.",
    caption  = "Points \u00b1 95% CI. Dodged vertically for clarity."
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title   = element_text(size = 14, face = "bold"),
    plot.subtitle= element_text(size = 11, color = "grey30"),
    axis.text.y  = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 11),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggsave("FigureS1_Sensitivity_Battery.pdf", figS1, width = 10, height = 6, device = cairo_pdf)
ggsave("FigureS1_Sensitivity_Battery.png", figS1, width = 10, height = 6, dpi = 300)
message("Figure S1 saved.")


# =============================================================================
# SECTION 17: SUPPLEMENTARY FIGURE S2 — INDIVIDUAL TEST Z-SCORE TRAJECTORIES
# =============================================================================
message("Building Figure S2: Individual test trajectories...")

indiv_tests  <- paste0("z_", intersect(ALL_TESTS, names(analytic)))
indiv_labels <- sub("^z_", "", indiv_tests)   # strip z_ prefix for display

obs_indiv <- analytic %>%
  mutate(year_bin   = round(years_to_anchor),
         case_label = factor(case_status, levels = c(0,1),
                             labels = c("Control","Case"))) %>%
  filter(year_bin >= WINDOW_MIN, year_bin <= WINDOW_MAX) %>%
  pivot_longer(cols = all_of(indiv_tests), names_to = "test_var", values_to = "z") %>%
  mutate(test_label = indiv_labels[match(test_var, indiv_tests)]) %>%
  group_by(test_label, year_bin, case_label) %>%
  summarise(
    mean_z = mean(z, na.rm = TRUE),
    se_z   = sd(z, na.rm = TRUE) / sqrt(sum(!is.na(z))),
    n      = sum(!is.na(z)),
    .groups = "drop"
  ) %>%
  filter(n >= 5)

obs_indiv$test_label <- factor(obs_indiv$test_label, levels = indiv_labels)

figS2 <- ggplot(obs_indiv, aes(x = year_bin, y = mean_z, color = case_label)) +
  annotate("rect", xmin = WINDOW_MIN, xmax = 0,
           ymin = -Inf, ymax = Inf, alpha = 0.06, fill = "#D6604D") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid",  color = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean_z - 1.96 * se_z, ymax = mean_z + 1.96 * se_z),
                width = 0.3, linewidth = 0.5, alpha = 0.6) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.0) +
  facet_wrap(~ test_label, ncol = 4, scales = "free_y") +
  scale_color_manual(values = GROUP_COLORS, name = NULL) +
  scale_x_continuous(breaks = c(WINDOW_MIN, WINDOW_MIN/2, 0, WINDOW_MAX),
                     name   = "Years relative to anchor") +
  scale_y_continuous(name = "Mean z-score (\u00b1 95% CI)") +
  labs(
    title    = "Figure S2. Trajectories of individual neuropsychological tests",
    subtitle = "Group means \u00b1 95% CI by year bin. Scores standardized to control baseline.",
    caption  = "Inverted tests: higher z = better performance."
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    strip.text    = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    axis.title    = element_text(size = 10),
    axis.text     = element_text(size = 9),
    legend.position  = "bottom",
    legend.text      = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

ggsave("FigureS2_Individual_Tests.pdf", figS2, width = 14, height = 11, device = cairo_pdf)
ggsave("FigureS2_Individual_Tests.png", figS2, width = 14, height = 11, dpi = 300)
message("Figure S2 saved.")


# =============================================================================
# SECTION 18: PARTICIPANT FLOW SUMMARY (CONSORT-STYLE)
# =============================================================================
message("Participant flow summary...")

cat("\n=== PARTICIPANT FLOW ===\n")
cat(sprintf("NPD file:             %d records, %d unique participants\n",
            nrow(npd_raw), n_distinct(npd_raw[[VAR_NPD_ID]])))
cat(sprintf("Diagnostic file:      %d unique participants\n",
            n_distinct(diag_raw[[VAR_ID]])))
cat(sprintf("Linked (both files):  %d participants\n",
            n_distinct(merged_raw$shareid)))
cat(sprintf("  Cases:              %d\n",
            n_distinct(merged_raw$shareid[merged_raw$case_status == 1])))
cat(sprintf("  Controls:           %d\n",
            n_distinct(merged_raw$shareid[merged_raw$case_status == 0])))
cat(sprintf("Analytic sample:      %d visits, %d participants (%d cases, %d controls)\n",
            nrow(analytic),
            n_distinct(analytic$shareid),
            n_distinct(analytic$shareid[analytic$case_status == 1]),
            n_distinct(analytic$shareid[analytic$case_status == 0])))
cat("========================\n\n")


# =============================================================================
# SECTION 19: SESSION INFO
# =============================================================================
message("\n=== ALL OUTPUTS SAVED ===")
message("Tables (Word .docx): Table1, Table2, Table3, TableS1, TableS2")
message("Figures (PDF + PNG): Figure1-6, FigureS1-S2")
message("\nRecommended manuscript order:")
message("  Table 1   — Sample characteristics")
message("  Figure 1  — Domain trajectory overview")
message("  Figure 2  — Divergence timeline (Cohen's d)")
message("  Figure 6  — Forest plot of interaction effects")
message("  Table 2   — LME fixed effects (primary)")
message("  Figure 4  — Education moderation")
message("  Figure 5  — Sex moderation")
message("  Table 3   — Moderation effects")
message("  Figure 3  — Individual spaghetti plots")
message("  Figure S1 — Battery sensitivity")
message("  Figure S2 — Individual test plots")
message("  Table S1  — Full fixed effects")
message("  Table S2  — Random effects / model fit")

writeLines(capture.output(sessionInfo()), "session_info.txt")
message("Session info saved to session_info.txt")
