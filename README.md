# HBC_SMB_PVA
Repo to support manuscript
Humpback Chub population viability model example R code

These five scripts demonstrate the population viability analysis (PVA) used in the Pine et al. paper. They contain the model assumptions and equations, but leave out the code used to format the tables and figures for NAJFM. They are examples and not an exact reproduction package for every published table and figure.

You need only R or RStudio. No extra R packages, input files, or directory setup are needed.
How to run
1. Download and unzip the folder. Keep the five .R files together.
2. Open RStudio and set the working directory to the folder containing these files (Session > Set Working Directory > Choose Directory).
3. In the R console, run these commands in order, in the same R session:
source("01_fixed_predation.R")
source("02_dynamic_predation.R")
source("03_sensitivity.R")
source("04_demographic_decomposition.R")
source("05_settled_abundance.R")


Do not restart R between scripts. Script 01 defines the model functions and population inputs used by the later scripts.

Each script prints progress in the console. The scripts automatically create a single results folder with CSV output files.

What each script does
01 — Fixed predation

Defines the model and runs fixed annual juvenile mortality rates for Dinosaur National Monument (DNM), Little Colorado River (LCR), and western Grand Canyon (WGC). Produces quasi-extinction risk and example adult-abundance trajectories.

02 — Dynamic predation

Uses a bootstrapped Smallmouth Bass recruitment series to generate temporally varying juvenile predation. Uses the paper's main Fmax values and population-specific packet sizes.

03 — Sensitivity

Changes the compensation ratio, size-dependent natural mortality, and juvenile predation packet size. Uses the paper's parameter ranges with a smaller subset of predation scenarios.

04 — Demographic decomposition

Holds LCR life history as a reference and changes starting abundance, compensation, or the survival schedule one factor at a time.

05 — Settled abundance

Summarizes average adult females over years 51–100, relative to each population's no-predation result, using the runs from 01 and 02.
Important model choices
- All abundance estimates represent females; adult females are age 4 and older.
- Recruitment is density dependent (Beverton–Holt) and includes lognormal annual variation.
- Fixed juvenile mortality occurs after density dependence and before fish enter age 1.
- Dynamic predation is related to a bootstrapped Smallmouth Bass abundance index. Packet mortality applies only to juveniles.
- Quasi-extinction is the event that adult female abundance falls to 5% or less of starting adult females at least once by a given year. It is not necessarily irreversible extinction.
- Each run has 100 annual transitions (year 0 through year 100). The random seed is 42.
- R0 is recalculated from the age-specific survival schedule, not entered as a fixed observed recruitment estimate.

Number of simulations
The default is 100 simulations per scenario so the examples are easier to run.

In 01_fixed_predation.R, change:
N_SIM <- 100L


to:
N_SIM <- 1000L


to use the paper's number of simulations.

The mortality levels used in script 03 are a subset of those in the manuscript. Simulation results will vary with Monte Carlo sample size and should not be taken as forecasts of actual Smallmouth Bass impacts.

Output
All CSV files appear in results/. Files 01–04 describe risk, trajectories, and sensitivities. File 05 describes long-term adult abundance.

The code does not reproduce the final publication styling, manuscript numbering, or map.

Note on the Smallmouth Bass input: The 18 recruitment values are included directly near the top of script 01, transcribed from the Middle Green River series used in the manuscript. This avoids any external input file.
