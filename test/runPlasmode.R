# add path to your data
# your data should have covariates columns, treatment and outcome.
# names of them do not matter as you can specify yourself later in the script
data_path <- "..."
final_cohort_cov <- "..." # read from data_path


# please provide column names of confounders you would like to use in the original data
# i.e. variables that affect both treatment and outcome
confColumns <- "..."

# please provide column names of instrumental variables you would like to use in the original data
# i.e., variables that only affect treatment
instrColumns <- "..."

# please provide column names of risk factors you would like to use in the original data
# i.e., variables that only affect outcome
rfColumns <- "..."

#------------------select confounders, risk factors and instrumental variables------------------
# specify the treatment/exposure column name
exp_cols <- "..."
# specify the outcome column name
out_cols <- "..."

# combine variables that affect treatment
exp_x_cols <- c(confColumns, instrColumns)

# combine variables that affect outcome
out_x_cols <- c(confColumns, rfColumns, exp_cols)

# specify number of simulations, sample size, effect size, id column name
nsim <- 100
nsize <- 10000
effectOR <- 1.5
idVar <- "id"

# create the exposure and outcome equations
formula_exp <- build_lm_formula(exp_x_cols, exp_cols)
formula_out <- build_lm_formula(out_x_cols, out_cols)

final_cohort_cov$id <- c(1:nrow(final_cohort_cov))

OriginalData <- final_cohort_cov %>% dplyr::mutate(id = row_number())
Bin_Form1 <- PlasmodeBin2(
  formulaOut = formula_out, objectOut = NULL, formulaExp = formula_exp,
  objectExp = NULL, data = final_cohort_cov, idVar = idVar, effectOR = effectOR,
  nsim = nsim, size = nsize, exposedPrev = NULL, eventRate = NULL
)

plasmodeData <- Bin_Form1$Sim_Data

#------------------save data------------------
for (i in 1:nsim)
{
  data <- NULL

  #specify outcome columns names
  idname <- "..." # please specify the ID column name eg paste("ID", i, sep = "")
  evtname <- "..." # please specify the event column name eg paste("EVENT", i, sep = "")
  expname <-"..." # please specify the outcome column name eg paste("EXPOSURE", i, sep = "")

  data <- plasmodeData %>%
    select(evtname, expname, idname) %>%
    rename(id = idname)

  data <- OriginalData %>% left_join(data, by = "id")

  # data file names saved each time with different i
  form <- sprintf("data_%s.csv", i)
  write.csv(data, file = form)
}
