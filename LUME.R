LUME <- function(df,
                               sample_id = 'Sample_Id',
                               time_id = 'Time',
                               parameters_moves = 'Gibbs',
                               model_moves = 'MH-Marginal',
                               type_model_prior = 'Uniform',
                               type_segment_prior = 'Uniform_card', # type of prior for the number of segments ('Uniform_card', 'Geometric', 'Poisson_card' or 'Geometric+Poisson')
                               type_coupling = 'global', # type of parameter coupling ('piece-wise' or 'global')
                               coupling = FALSE, # if parameter coupling should be applied
                               standardized = FALSE,
                               time_lags = 1,
                               self_loop = FALSE,
                               p_MCMC = 1, # probability of doing an MCMC move
                               p_RJMCMC = 1,
                               p_ChangePoint = 1, # probability of doing a Change Point Detection move
                               missing_update_period = 10, #how many epochs between missing updates
                               group_update_epochs = 10, 
                               alpha = 0.01,
                               beta = 0.01,
                               alpha_xi = 0.01, #xi2 scale parameter
                               beta_xi = 0.01, #xi2 rate parameter
                               a = 0.01,
                               b = 0.01,
                               v = 0.1,
                               r = 1,
                               lambda = 1,
                               lambda_group = 1, 
                               rho = 0.1, #geometric probability parameter for segment distance prior
                               iota = 1, # expected value of number of successes of Poisson segments' granularity prior
                               mu_cross = 'zeros', # global mean vector hyperprior mean (only in globally-coupled settings) - default is a 0 vector for each outcome variable
                               Sigma_cross = 'ones', # global mean vector hyperprior covariance matrix (only in globally-coupled settings) - default is an identity matrix for each outcome variable
                               intercept_move = FALSE,
                               max_parents = ncol(df) - 3 + self_loop,
                               max_changepoints = NaN, #max number of changepoints
                               epochs = 1000,
                               init_sigma2 = 1, 
                               init_delta2 = 1,
                               init_B = 'zeros',
                               init_M = 'random',
                               init_A = 'no changepoints',
                               init_MU = 'zeros',
                               init_group_segmentation = 'all_separated', #could be 'all_separated', 'all_in_one', 'random'
                               init_xi2 = 1, # initial xi2
                               init_allocation_vector = NaN, # initial segment set
                               segmentation = FALSE, # if data segmentation should be applied
                               group_segmentation = FALSE, # if sample grouping is applied
                               segmentation_var = time_id, # variable on which data is segmented
                               group_segmentation_var = sample_id,
                               group_segmentation_list = NaN, # Fixed group segmentation vector (if desired)
                               n_segments = NaN, # Fixed number of segments (if desired)
                               segmentation_vector = NaN, # Fixed segmentation vector (if desired)
                               log = TRUE,
                               sigma2_in_betas_prior = FALSE, # Inclusion of sigma^2 in the Prior - and FCD - of beta
                               max_missing_input = +Inf
) {
  # Get variable names (assumes all columns except Sample_Id and Time are variables)
  tot_time <- list(params = 0, structure = 0, missing = 0, data_alloc = 0)
  posterior_samples <- list()
  missing_values <- list()
  group_segmentation_lists <- list()
  log_marginals <- list()
  variables <- sort(setdiff(names(df), c(sample_id, time_id)))
  k <- length(variables)
  T_max <- max(df[,segmentation_var])
  # Define covariates as all variables at t-1
  if (!self_loop) {
    covars <- lapply(variables, function(x) {c('1',unlist(sapply(1:time_lags, function(lag) paste0(setdiff(variables, x), "_t-", lag))))})
  }
  else{
    covars <- lapply(variables, function(x) {c('1',unlist(sapply(1:time_lags, function(lag) paste0(variables, "_t-", lag))))})
  }
  var_names <- setdiff(names(df), c(sample_id, time_id))
  
  #mean_0 <- sapply(var_names, function(x) mean(df[df[[time_id]] == 0, x], na.rm = TRUE))
  df <- cbind(df[,c(sample_id, time_id)], '1' = 1, df[, sort(setdiff(names(df), c(sample_id, time_id, "1")))])
  missing_mask <- df %>%
    mutate(across(-c(sample_id, time_id), ~ ifelse(is.na(.), TRUE, FALSE)))
  if (standardized){
    df[,variables] <- StandardizeData(df[,variables])
  }
  mean_0 <- sapply(var_names, function(x) {
    time0_vals <- df[df[[time_id]] == 0, x]
    if (any(!is.na(time0_vals))) {
      mean(time0_vals, na.rm = TRUE)
    } else {
      time1_vals <- df[df[[time_id]] == 1, x]
      if (any(!is.na(time1_vals))) {
        mean(time1_vals, na.rm = TRUE)
      } else {
        mean(df[[x]], na.rm = TRUE)
      }
    }
  })
  names(mean_0) <- var_names  # ensures names are retained
  sd_0 <- sapply(var_names, function(x) {1}) #sd(df[df[[time_id]] == 0, x], na.rm = TRUE)
  names(sd_0) <- var_names
  
  df <-InitialMissingImputation(df, s = 5)
  vars_ <- colnames(df)[!colnames(df) %in% c(sample_id,time_id)]
  # Create lagged dataset for DBN
  M <- M_Init(df, init = init_M, sample_id = sample_id, time_id = time_id, intercept_move = intercept_move, max_parents = max_parents - intercept_move)
  B <- B_Init_group(df, init = init_B, sample_id = sample_id, time_id = time_id)
  MU <- B_Init_group(df, init = init_MU, sample_id = sample_id, time_id = time_id)
  A <- A_Init(df, init = init_A, sample_id = sample_id, time_id = time_id) #matrix(0, nrow = k, ncol = max(df[,segmentation_var]) - 1)
  segmentation_vector
  Sigma <- rep(init_sigma2, k)
  Delta <- rep(init_delta2, k)
  for (j in 1:k){
    posterior_samples[[variables[j]]] <- list('betas' = list(), 'mus' = list(), 'sigma2s' = list(),'delta2s' = list(),'models' = list(),'allocation_vectors' = list(),'covariates'= gsub("_t-1", "", covars[[j]]))
  }
  lagged_data <- LaggedDatasetGeneration(df[,!(names(df) %in% c('1'))], lag = time_lags, sample_id = sample_id, time_id = time_id)
  lagged_data <- cbind(lagged_data[, 1:2], '1' = 1, lagged_data[, 3:ncol(lagged_data)])
  # Fit one model for each variable at time t
  for (epoch in 1:round(epochs/missing_update_period)) {
    group_segmentation_vector <- sapply(lagged_data[[sample_id]], function(x) {group_segmentation_list[[as.character(x)]]})
    for (j in 1:k) {
      outcome_var <- paste0(variables[j], "_t")
      allocation_vector <- c(1, unname(which(A[variables[j], ] == 1)), max(df[,time_id]))
      res <- RJMCMC(
        data = lagged_data[,c(sample_id, time_id, covars[[j]], outcome_var)],
        epochs = missing_update_period,
        v = v,
        type_segment_prior = type_segment_prior, # type of prior for the number of segments ('Uniform_card', 'Geometric', 'Poisson_card' or 'Geometric+Poisson')
        type_coupling = type_coupling, # type of parameter coupling ('piece-wise' or 'global')
        coupling = coupling, # if parameter coupling should be applied
        p_MCMC = p_MCMC, # probability of doing a parameter move
        p_RJMCMC = p_RJMCMC, # probability of doing a model move
        p_ChangePoint = p_ChangePoint, # probability of doing a Change Point Detection move
        alpha = alpha,
        beta =  beta,
        alpha_xi = alpha_xi, #xi2 scale parameter
        beta_xi = beta_xi, #xi2 rate parameter
        a = a,
        b = b,
        r = r,
        lambda = lambda,
        rho = rho, #geometric probability parameter for segment distance prior
        iota = iota, # expected value of number of successes of Poisson segments' granularity prior
        mu_cross = rep(0, ifelse(self_loop,k+1,k)),
        Sigma_cross = diag(1, nrow = ifelse(self_loop,k+1,k)),
        max_parents = max_parents,
        max_changepoints = max_changepoints, #max number of changepoints
        intercept_move = intercept_move,
        outcome = outcome_var,
        covariates = covars[[j]],
        init_model = M[j,ifelse(self_loop,all(),-(j+1))],
        init_betas = lapply(1:max(group_segmentation_vector), function(x) {B[j,ifelse(self_loop,all(),-(j+1)),x]}),
        init_allocation_vector = allocation_vector,
        betas_mean = lapply(1:max(group_segmentation_vector), function(x) {MU[j,ifelse(self_loop,all(),-(j+1)),x]}),
        # this parameters depend on the type of coupling adopted - no coupling renders just to the initial hyperprior parameters, global coupling takes the values from MU and sequential coupling update them for each time slice based on the current parameter values
        init_sigma2 = Sigma[j],
        init_delta2 = Delta[j],
        parameters_moves = parameters_moves,
        model_moves = model_moves,
        type_model_prior = type_model_prior,
        segmentation = segmentation, # if data segmentation should be applied
        group_segmentation = group_segmentation,
        segmentation_var = segmentation_var, # variable on which data is segmented
        group_segmentation_var = group_segmentation_var,
        n_segments = n_segments, # Fixed number of segments (if desired)
        segmentation_vector = segmentation_vector, # Fixed segmentation vector
        group_segmentation_vector = group_segmentation_vector,
        log = log,
        sigma2_in_betas_prior = sigma2_in_betas_prior
      )
      for (x in c('betas', 'mus', 'sigma2s', 'delta2s', 'models', 'allocation_vectors')) {
        posterior_samples[[variables[j]]][[x]] <-  append(posterior_samples[[variables[j]]][[x]], res[[x]])
      }
    }
    posterior_samples_mat <- PosteriorSamplesMatricesG(posterior_samples, T = T_max, group_segmentation_list = group_segmentation_list)
    M <- posterior_samples_mat$M
    B <- posterior_samples_mat$B
    A <- posterior_samples_mat$A
    MU <- posterior_samples_mat$MU
    Sigma <- posterior_samples_mat$Sigma
    Delta <- posterior_samples_mat$Delta
    missing_structure <- MissingFCDComponentsG(B, M, Sigma, mean_0, sd_0)
    df <- MissingValuesUpdateG(data = df, missing_mask = missing_mask, structure = missing_structure, group_segmentation_list = group_segmentation_list, mean_0 = mean_0, sd_0 = sd_0, sample_id = sample_id, time_id = time_id)
    df[vars_][abs(df[vars_]) >= max_missing_input] <- 0
    if (epoch %% 100 == 0) {
      cat("Epoch:", (epoch * missing_update_period), "\n")
    }
    missing_values[[epoch]] <- as.matrix(df[,-c(1, 2)])[which(as.matrix(missing_mask[,-c(1, 2)]), arr.ind = TRUE)]
    lagged_data <- LaggedDatasetGeneration(df[,!(names(df) %in% c('1'))], lag = time_lags, sample_id = sample_id, time_id = time_id)
    lagged_data <- cbind(lagged_data[, 1:2], '1' = 1, lagged_data[, 3:ncol(lagged_data)])
    if (group_segmentation){
      for (s in 1:group_update_epochs){
        posterior_group <- GroupUpdate(grouping_list = group_segmentation_list, lagged_data = lagged_data, M = M, MU = MU, Sigma = Sigma, Delta = Delta, mu_cross = rep(0, (k+1)), Sigma_cross = diag(1, (k+1)), group_var = group_segmentation_var, lambda_group = lambda_group, coupling = coupling)
        group_segmentation_list <- posterior_group$grouping_list
        MU <- posterior_group$MU
        B <- posterior_group$B
        log_marginal <- posterior_group$log_marginal
        group_segmentation_lists[[(epoch-1)*group_update_epochs+s]] <- group_segmentation_list
        log_marginals[[(epoch-1)*group_update_epochs+s]] <- log_marginal
      }
    }
  }
  list(posterior_samples = posterior_samples, missings = missing_values, group_segmentation_lists = group_segmentation_lists, log_marginals = log_marginals)
}
