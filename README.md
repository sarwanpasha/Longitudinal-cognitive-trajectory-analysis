# Longitudinal Cognitive Trajectory Analysis Toolkit

A statistical framework for analyzing **domain-specific cognitive decline trajectories** using linear mixed-effects (LME) models. Designed for longitudinal neuropsychological studies with repeated measures per participant, where the goal is to characterize cognitive change relative to a clinically meaningful anchor event (e.g., diagnosis, symptom onset, or last known-normal review).

---

## Overview

This toolkit provides a fully documented, configurable R script that takes two longitudinal data files — a clinical review/diagnostic file and a neuropsychological testing file — and produces a complete set of publication-ready figures and tables.

**Primary analysis:** Domain-specific cognitive trajectory modeling (LME with random intercept + slope per participant), testing whether slopes differ between cases and controls in the years surrounding a clinical anchor event.

**Secondary analysis:** Moderation of case-related decline by education level and sex (cognitive reserve framework).

**Sensitivity analysis:** Consistency of findings across test battery versions.

---

## Features

- ⚙️ **Fully configurable** — all variable names and file paths are set in a single `Section 0` block; no editing required elsewhere
- 📊 **5 cognitive domain composites** — memory, executive function, language/visuospatial, processing speed, global
- 📈 **8 publication-ready figures** (PDF + PNG): trajectory plots, divergence timelines, spaghetti plots, forest plots, moderation plots
- 📋 **5 formatted tables** (Word `.docx`): sample characteristics, LME results, moderation results, full fixed effects, random effects
- 🔁 **Battery/instrument sensitivity analysis** — verifies results hold within each test version
- 🧮 **Automatic z-scoring** against control baseline with sign-inversion for tests where higher = worse
- 🔒 **No data hardcoding** — the script works with any dataset that fits the described structure

---

## Requirements

### R version
R ≥ 4.1.0 recommended

### Required packages
The script will auto-install missing packages on first run.

```r
tidyverse, lme4, lmerTest, broom.mixed, emmeans,
ggplot2, ggpubr, patchwork, gt, gtsummary, flextable,
splines, nlme, scales, RColorBrewer, viridis,
cowplot, officer, ggstance
```

---

## Input Data Format

The script expects **two tab-delimited text files**. Each file should have one row per participant per visit/review, with a header row.

### File 1: Diagnostic / Clinical Review File
One row per participant per review visit.

| Variable (generic name) | Description |
|---|---|
| Participant ID | Unique participant identifier (links to NPD file) |
| Case status | 0 = control, 1–4 = case (stages), other = excluded |
| CDR score | Clinical Dementia Rating (optional; for Table 1) |
| Impairment date | Date of first documented impairment (cases' time anchor) |
| Last-normal date | Date of last normal cognitive review |
| Review date | Date of most recent review (controls' time anchor) |
| Stage transition dates | Mild / moderate / severe transition dates (optional) |

### File 2: Neuropsychological Testing File
One row per participant per exam.

| Variable (generic name) | Description |
|---|---|
| Participant ID | Links to diagnostic file |
| Exam date | Date of neuropsychological exam |
| Education | Ordinal: 0 = <HS, 1 = HS grad, 2 = some college, 3 = college+ |
| Sex | 1 = male, 2 = female |
| Age | Age at exam (years) |
| Battery | Test battery version (1 or 2; use 1 if only one version) |
| Test scores | Raw scores for each neuropsychological test (see below) |

### Supported Neuropsychological Tests
The script can accommodate the following cognitive domains. Any test not present in your data can be set to `NA` in the configuration block and will be excluded automatically.

| Domain | Tests supported |
|---|---|
| **Memory** | Immediate verbal memory, delayed verbal memory, immediate visuospatial memory, delayed visuospatial memory |
| **Executive Function** | Complex attention/switching (e.g., Trail Making B), working memory backward, phonemic verbal fluency |
| **Language / Visuospatial** | Confrontation naming, visual organization |
| **Processing Speed** | Simple psychomotor speed (e.g., Trail Making A), working memory forward |

> **Note on date formats:** By default the script assumes dates are stored as SAS integers (days since 1960-01-01). If your dates use a different epoch (e.g., Unix time: 1970-01-01), update `SAS_ORIGIN` in the configuration block.

---

## Usage

### 1. Clone the repository

```bash
git clone https://github.com/your-username/longitudinal-cognitive-trajectory-analysis.git
cd longitudinal-cognitive-trajectory-analysis
```

### 2. Place your data files

Copy your two data files into the project directory (or note their full paths).

### 3. Configure the script

Open `longitudinal_cognitive_trajectory_analysis.R` and edit **Section 0** only:

```r
# File paths
DIAG_FILE <- "path/to/your_diagnostic_review_file.txt"
NPD_FILE  <- "path/to/your_neuropsych_testing_file.txt"
SKIP_ROWS <- 0   # rows to skip before the header row

# Variable name mapping (right-hand side = your actual column names)
VAR_ID          <- "your_participant_id_column"
VAR_CASE_STATUS <- "your_case_status_column"
# ... (see script for full list)

# Test variable names
TEST_MEMORY_IMM <- "your_immediate_memory_test_column"
# Set to NA if a test is not in your data:
TEST_VISORGAN  <- NA
```

### 4. Run the script

```r
source("longitudinal_cognitive_trajectory_analysis.R")
```

All outputs will be written to the working directory.

---

## Output Files

### Tables (Word `.docx`)

| File | Contents |
|---|---|
| `Table1_Sample_Characteristics.docx` | Baseline demographics, cognitive scores by case status |
| `Table2_LME_Primary_Results.docx` | Fixed effects from primary trajectory models |
| `Table3_Moderation_Results.docx` | Education × case and sex × case interaction terms |
| `TableS1_Full_Fixed_Effects.docx` | Complete fixed-effect table for all models |
| `TableS2_Random_Effects.docx` | Random effects SDs, correlation, residual SD, AIC |

### Figures (PDF + PNG)

| File | Contents |
|---|---|
| `Figure1_Domain_Trajectories` | LME-predicted trajectories (cases vs controls) per domain |
| `Figure2_Divergence_Timeline` | Year-by-year Cohen's d showing when domains diverge |
| `Figure3_Spaghetti_Plots` | Individual participant trajectories by domain and group |
| `Figure4_Education_Moderation` | Education group × case trajectory plots |
| `Figure5_Sex_Moderation` | Sex × case trajectory plots |
| `Figure6_Forest_Plot` | Forest plot of slope × case interaction effect sizes |
| `FigureS1_Sensitivity_Battery` | Battery-stratified replication of primary interaction |
| `FigureS2_Individual_Tests` | Observed trajectories for each individual test |

---

## Statistical Methods

### Primary Model (per cognitive domain)

```
domain_score ~ years_to_anchor × case_status + age + sex + battery
               + (1 + years_to_anchor | participant_id)
```

- **Time variable (`years_to_anchor`):** continuous; 0 = clinical anchor date. Negative = pre-event; positive = post-event.
- **Anchor date:** For cases, the estimated date of clinical impairment onset. For controls, the most recent review date.
- **Key test:** The `years_to_anchor × case_status` interaction — tests whether longitudinal slopes differ between groups.
- **Random effects:** Random intercept and slope per participant to capture individual variability in baseline and rate of change.
- **Covariates:** Age, sex, battery version.
- **Estimation:** REML = FALSE (ML, for model comparison); optimizer = BOBYQA.

### Z-scoring
All test scores are standardized to the mean and SD of controls at their earliest (baseline) visit. Tests where higher scores indicate worse performance (e.g., timed tests) are sign-inverted so that higher z always means better performance.

### Domain Composites
Each domain composite is the row mean of its constituent z-scored tests, computed with `na.rm = TRUE` (participants with partial test completion are retained).

### Moderation Models (Secondary Analysis)
Three-way interaction models (`slope × case × moderator`) fit separately for each moderator (education, sex) on memory and global composite outcomes.

### Effect Size
Cohen's d is computed year-by-year as `(mean_control − mean_case) / pooled_SD`. Bins with fewer than 10 cases or 5 controls are excluded.

---

## Folder Structure

```
longitudinal-cognitive-trajectory-analysis/
├── longitudinal_cognitive_trajectory_analysis.R   # main script
├── README.md
└── session_info.txt                               # generated on run
```

---

## Adapting to Your Study

This toolkit is designed to be **domain-agnostic**. While the default domain groupings reflect a common neuropsychological battery structure, you can reassign any test to any domain by changing the `make_composite()` calls in Section 4f. You can also add domains beyond the four defaults by following the same pattern.

The time-anchor logic is equally flexible: any date variable in your diagnostic file can serve as the anchor by setting `VAR_DATE_IMPAIR` accordingly. For studies without a clinical event date, you can anchor all participants to a fixed visit number or calendar date.

---

## Citation

If you use this toolkit in your research, please cite it as:

```
Longitudinal Cognitive Trajectory Analysis Toolkit (2025).
GitHub: [https://github.com/your-username/longitudinal-cognitive-trajectory-analysis](https://github.com/sarwanpasha/Longitudinal-cognitive-trajectory-analysis)
```

---

## License

MIT License. See `LICENSE` for details.

---

## Contributing

Contributions, bug reports, and feature requests are welcome. Please open an issue or submit a pull request.
