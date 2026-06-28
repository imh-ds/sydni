
# -- EXPERIMENT 1 -- #

library(tidyverse)
# library(sydni)
devtools::load_all(".")
library(cscore)
library(qol)

# save.image("3_Experiments/Experiment_01/1_Data/Exp01_Working.rdata")
# load("3_Experiments/Experiment_01/1_Data/Exp01_Working.rdata")

# DATA PREPARATION --------------------------------------------------------

# Read in data
data <- readr::read_csv(
  "tests/3_Experiments/Experiment_01/1_Data/cc_data_s2_data.csv"
) |> 
  
  # Filter to a priori attentive data
  filter(
    attention == 1,
    Late == 0,
    gen != 3
  ) |> 
  
  # Select relevant variables
  dplyr::select(
    ResponseId, age, gen, DAD, MOM,
    starts_with("mss"),
    starts_with("dcs"),
    starts_with("hssb"),
    -ends_with("DO")
  ) %>%
  
  # Create new ID
  dplyr::mutate(
    ResponseId = sprintf("id_%03d", 
                         seq(nrow(.)))
  )

# Make variable names lower-case
names(data) <- tolower(names(data))


# METADATA PARAMETERS -----------------------------------------------------

# Composite list
composite_list <- cscore::composite_list(
  
  # Socioeconomic Status
  ses = c("mom", "dad"),
  
  # Consumer attitudes
  material = c("mss_1", "mss_2", "mss_3", "mss_4"),
  dissatis = c("dcs_1", "dcs_2", "dcs_3", "dcs_4"),
  
  # Merit
  hardwork = c("hssb_1", "hssb_2"),
  ability = c("hssb_3", "hssb_4"),
  merit = c("hardwork", "ability"),
  
  # Opportunity
  socialconnections = c("hssb_5", "hssb_6"),
  privilege = c("hssb_7", "hssb_8"),
  opportunity = c("socialconnections", "privilege"),
  
  # Chance
  luck = c("hssb_9", "hssb_10"),
  fate = c("hssb_11", "hssb_12"),
  chance = c("luck", "fate")
)

# Variables
vars <- c(sprintf("mss_%01d", seq(4)),
          sprintf("dcs_%01d", seq(4)),
          sprintf("hssb_%01d", seq(12)))
covars <- c("age", "gen", "dad", "mom")

# Synthetic data parameters
id <- "responseid"
min_class <- "gen"
min_class_value <- 2
n <- 142
dup_size <- 6

# Impute data with knn imputation
kdata <- knn_impute(
  data = data
)

# Create copy for SYDNI augmentation
sdata <- kdata %>% 
  dplyr::select(-dplyr::all_of(id)) %>% 
  dplyr::filter(
    !!sym(min_class) == min_class_value
  )

# Create copy for observed data
odata <- kdata %>% 
  dplyr::filter(
    !!sym(min_class) == min_class_value
  ) %>% 
  cscore::composite_score(
    composite_list
  )

# Create copy with composite variables
cdata <- kdata %>% 
  cscore::composite_score(
    composite_list
  )


# -- BENCHMARK PARAMETERS -- #

# Correlations performance benchmark parameters
bench_vars <- c("material", "dissatis", "merit", "opportunity", "chance")
bench_var_names <- c("Materialism", "Dissatisfaction", "Merit", "Opportunity", "Chance")
synthetic_order <- c("sydni", "synthpop", "gans", "smote", "adasyn")


# GENERATE SYNTHETIC DATA -------------------------------------------------

# # -- SYDNI -- #
# sydni_list <- sydni(
#   data = sdata,
#   vars = vars,
#   covariates = covars,
#   aug_n = 142,
#   tree_n = 100,
#   iter_n = 20,
#   mean_index_weight = 0.3,
#   sd_index_weight = 0.2,
#   cor_index_weight = 0.4,
#   skew_index_weight = 0.05,
#   kurtosis_index_weight = 0.05,
#   return_all_mods = TRUE,
#   replace = TRUE
# )

# one_var <- vars[1]
# two_var <- vars[1:2]
families <- stats::setNames(
  rep(list(gaussian()), length(vars)),
  vars
)

results <- sydni(
    data = sdata,
    variables = vars,
    covariates = covars,
    engine = "copula",
    # measurement_blocks = composite_list,
    n_aug = nrow(cdata)-nrow(sdata),
    families = families,
    # predictor_k = 2,
    # iter = 500,
    # warmup = 25,
    chains = 1,
    cores = min(4, max(1, parallel::detectCores() - 1)),
    algorithm = "fullrank",
    backend = "cmdstanr",
    seed = 123
  )


kdata %>% 
  cscore::composite_score(
    composite_list=composite_list,
    return_metrics=TRUE
  ) |> 
  purrr::pluck('validity')

results$data %>% 
  cscore::composite_score(
    composite_list,
    return_metrics=TRUE
  ) |> 
  purrr::pluck('validity')


obs <- kdata |> 
  cscore::composite_score(
    composite_list = composite_list
  )
syd <- results$data |> 
  cscore::composite_score(
    composite_list = composite_list
  )

qol::cor_matrix(
  obs,
  vars = c("ses", "material", "dissatis", "merit", "opportunity", "chance")
)
qol::cor_matrix(
  syd,
  vars = c("ses", "material", "dissatis", "merit", "opportunity", "chance")
)

m <- results$models$all

m$structural_correlation_target
m$structural_correlation_generated
m$structural_correlation_error_raw
m$structural_correlation_error_calibrated
m$within_family_correlation_max_change


# # See which model is best
# head(sydni_list$accuracy_table)

# # Extract SYDNI data
# data_sydni <- sydni_list$augmented_data$model_55$aug_data


# Generate synthetic SMOTE, ADASYN, GANs, and SynthPop data
synthetic_data <- synthetic_data_generation(
  data = kdata,
  id = id,
  min_class = min_class,
  min_class_value = min_class_value,
  n = n,
  dup_size = dup_size
)

# Merge SYDNI data into synthetic data list
synthetic_data <- synthetic_data_merge(
  synthetic_data,
  data_sydni,
  id = id,
  min_class = min_class,
  min_class_value = min_class_value,
  n = n
)


# SYNTHETIC DATA METRICS --------------------------------------------------

# Calculate benchmarks
benchmark_metrics <- benchmark(
  observed_data = odata,
  synthetic_data = synthetic_data,
  composite_list = composite_list,
  vars = bench_vars,
  var_names = bench_var_names,
  synthetic_order = synthetic_order
)
