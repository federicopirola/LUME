# Set Working Directory to Source File Location
#setwd("/path/to/your/project")

# Load Functions
source('distributions.R')
source('data_model_generation.R')
source('model_parameter_update.R')
source('missing_imputation.R')
source('group_update.R')
source('LUME.R')

# Install and Load Libraries
packages <- c(
  "ggplot2",
  "emdbook",
  "geomtextpath",
  "reshape2",
  "MASS",
  "dplyr",
  "coda",
  "zoo"
)
missing_packages <- packages[!(packages %in% rownames(installed.packages()))]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}
invisible(lapply(packages, library, character.only = TRUE))

# Initialize Data Generation Hyperarameters
N<- 100 # number of samples
G <- 5 # number of groups
k <- 10 # number of nodes
time_frames <- 5 # number of temporal frames
missingness_rate <- 0.2 # missingness probability

# Define Data 
base <- rep(1, G)
remaining <- N - G
alpha <- rep(10, G)
p <- rgamma(G, shape = alpha, rate = 1)
p <- p / sum(p)
extra <- as.vector(rmultinom(1, remaining, p))
sampled_sizes <- base + extra

# Generate GH-DBN 
DBNs <- vector("list", G)
datasets <- vector("list", G)
DBN_star <- DBNGeneration(k = k, a1 = 0.4, a2 = 0.8)
for (g in 1:G) {
  DBNs[[g]] <- DBN_star
  for (var in sample(names(DBN_star$TransitionNetwork))) {
    n_par <- length(DBNs[[g]]$TransitionNetwork[[var]]$parents)
    DBNs[[g]]$TransitionNetwork[[var]]$params <-
      MultivariateNormalSample(N = 1, mu = DBNs[[g]]$TransitionNetwork[[var]]$params, Sigma = diag(1,n_par))
  }
  datasets[[g]] <- DataGenerationFromDBN(sampled_sizes[g], time_frames, DBNs[[g]])
  datasets[[g]]$Sample_Id <- paste0(LETTERS[[g]],datasets[[g]]$Sample_Id)
  datasets[[g]] <- datasets[[g]] %>%
    mutate(across(-c(1, 2), ~ ifelse(runif(n()) < missingness_rate, NA, .)))
}
data_all <- do.call(rbind, datasets)

#LEARN HOMOGENEOUS DBN FROM DATA

# Initialize group allocation vector: all in one group (homogeneous)
init_grouping_list <- rep(1,length(unique(data_all$Sample_Id))) 
names(init_grouping_list) <- as.character(unique(data_all$Sample_Id))

# LUME-DBN Algorithm
posterior_DBN <- LUME(
  data_all, # dataset
  sample_id = 'Sample_Id', # sample id column
  time_id = 'Time', # time id column
  group_segmentation_var = 'Sample_Id', # group id column (sample id)
  alpha = 2, # psi2 Inverse Gamma scale parameter
  beta = 0.2, # psi2 Inverse Gamma rate parameter
  a = 0.01, # sigma2 Inverse Gamma scale parameter
  b = 0.01, # sigma2 Inverse Gamma scale parameter
  epochs = 10000, # number of epochs
  parameters_moves = 'Gibbs', # Gibbs Step for parameter updates
  model_moves = 'MH-Marginal', # Metropolis-Hasting Step based on Marginal Likelihoods for structural updates
  type_model_prior = 'Poisson', # Poisson Prior for parent sets 
  lambda = 1, # Parent Set Cardinality Poisson prior parameter
  max_parents = 5, # Maximum number of parents per node
  standardized = FALSE, # If data standardization should be applied
  segmentation = FALSE, # If temporal segmentation should be enabled
  group_segmentation = FALSE, # If group segmentation should be enabled
  group_segmentation_list = init_grouping_list # Initial group allocation vector
)

#LEARN GROUP-HETEROGENEOUS DBN FROM DATA
# Initialize group allocation vector: every sample in a distinct group
init_grouping_list <- 1:length(unique(data_all$Sample_Id)) 
names(init_grouping_list) <- as.character(unique(data_all$Sample_Id))

# LUME-GHDBN Algorithm
posterior_GHDBN <- LUME(
  data_all, # dataset
  sample_id = 'Sample_Id', # sample id column
  time_id = 'Time', # time id column
  group_segmentation_var = 'Sample_Id', # group id column (sample id)
  alpha = 2, # psi2 Inverse Gamma scale parameter
  beta = 0.2, # psi2 Inverse Gamma rate parameter
  a = 0.01, # sigma2 Inverse Gamma scale parameter
  b = 0.01, # sigma2 Inverse Gamma scale parameter
  epochs = 10000, # number of epochs
  parameters_moves = 'Gibbs', # Gibbs Step for parameter updates
  model_moves = 'MH-Marginal', # Metropolis-Hasting Step based on Marginal Likelihoods for structural updates
  type_model_prior = 'Poisson', # Poisson Prior for parent sets 
  lambda = 1, # Parent Set Cardinality Poisson prior parameter
  max_parents = 5, # Maximum number of parents per node
  standardized = FALSE, # If data standardization should be applied
  segmentation = FALSE, # If temporal segmentation should be enabled
  group_segmentation = TRUE, # If group segmentation should be enabled
  group_segmentation_list = init_grouping_list, # Initial group allocation vector
  coupling = TRUE, # If parameter coupling (hierarchical model) should be imposed
  lambda_group = 0.1 # Group Cardinality Poisson prior parameter
)
