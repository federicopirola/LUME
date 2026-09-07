# In case of NH-DBNs with sequentially coupled parameters these vectors is not directly required because the current value of the parameters are considered instead of the prior mean
BetasMeanVector <- function(Y, X, S,
                            model = rep(1, ncol(X)),
                            allocation_vector = c(1, max(S)),
                            init_betas = rep(0, ncol(X)),
                            betas_mean = rep(0, ncol(X)),
                            sigma2 = 1, delta2 = 1, xi2 = 1, coupled = FALSE, type_coupling = 'global',
                            group_segmentation = FALSE, segment_indexes = rep(1, length(Y))) {
  
  model_logical <- model == 1
  k <- sum(model_logical)
  if (group_segmentation){
    H <- length(segment_indexes) # Number of groups
    return(betas_mean)
  }
  else{
    H <- length(allocation_vector) - 1  # Number of segments
    
  }
  
  if (coupled & type_coupling == 'sequential'){
    betas_mean <- vector("list", H - 1)
    X_included <- X[, model_logical, drop = FALSE]
    init_betas_included <- init_betas[model_logical]
    
    # h = 1 → initial beta vector
    betas_mean[[1]] <- setNames(init_betas, colnames(X))  # already in full space
    
    # h = 2
    if (H >= 2) {
      X_curr <- as.matrix(X_included[S >= allocation_vector[2] & S < allocation_vector[3], , drop = FALSE])
      Y_curr <- Y[S >= allocation_vector[2] & S < allocation_vector[3]]
      
      beta_h2 <- solve((1 / delta2) * diag(k) + t(X_curr) %*% X_curr) %*%
        ((1 / delta2) * init_betas_included + t(X_curr) %*% Y_curr)
      
      beta_full <- rep(0, ncol(X))
      beta_full[model_logical] <- beta_h2[, 1]
      betas_mean[[2]] <- setNames(beta_full, colnames(X))
    }
    
    # h >= 3
    if (H >= 3) {
      for (h in 3:H) {
        X_curr <- as.matrix(X_included[S >= allocation_vector[h] & S < allocation_vector[h + 1], , drop = FALSE])
        Y_curr <- Y[S >= allocation_vector[h] & S < allocation_vector[h + 1]]
        
        beta_prev <- betas_mean[[h - 1]][model_logical]  # extract relevant subset
        
        beta_h <- solve((1 / xi2) * diag(k) + t(X_curr) %*% X_curr) %*%
          ((1 / xi2) * beta_prev + t(X_curr) %*% Y_curr)
        
        beta_full <- rep(0, ncol(X))
        beta_full[model_logical] <- beta_h[, 1]
        betas_mean[[h]] <- setNames(beta_full, colnames(X))
      }
    }
  }
  else if (coupled & type_coupling == 'global'){
    betas_mean <- betas_mean
  }
  else if (any(allocation_vector != c(1, max(S)))){
    betas_mean <- betas_mean
  }
  else{
    betas_mean <- vector("list", H - 1)
    named_vector_template <- setNames(rep(0, ncol(X)), names(model))
    if (H == 1){
      betas_mean[[1]] <- named_vector_template
    }
    else{
      betas_mean <- replicate(H-1, named_vector_template, simplify = FALSE)
    }
  }
  return(betas_mean)
}

Gibbs_move = function(Y, X, S = rep(1, length(Y)), model = rep(1,ncol(X)), betas = rep(0,ncol(X)), sigma2 = 1, delta2 = 1, alpha = 0.01, beta = 0.01, a = 0.01, b = 0.01, betas_mean = rep(0,ncol(X)), sigma2_in_betas_prior = TRUE, collapsed = FALSE, allocation_vector = c(1, max(S)), coupling = FALSE, group_segmentation = FALSE, segment_indexes = rep(1, length(Y))){
  betas_included <- lapply(if (is.list(betas)) betas else list(betas), function(b) b[model == 1])
  X <- as.matrix(X[,model == 1])
  betas_mean  <- lapply(if (is.list(betas_mean)) betas_mean else list(betas_mean), function(b) b[model == 1])
  if (!group_segmentation){
    segment_indexes <- if(length(allocation_vector) > 2){
      cut(S, breaks = allocation_vector, labels = FALSE, right = FALSE, include.lowest = TRUE)
    } else{
      rep(1, length(Y))
    }
  }
  # Sample sigma2 | betas, Y or sigma2 | Y
  sigma2 <- Sigma2FCDSample(N = 1, Y = Y, X = X, betas = betas_included, delta2 = delta2, a = a, b = b, mu_betas = betas_mean, sigma2_in_betas_prior = sigma2_in_betas_prior, collapsed = collapsed, segment_indexes = segment_indexes)
  
  # Sample betas | sigma2, delta2, Y
  betas_included <- BetasFCDSample(N = 1, Y = Y, X = X, sigma2 = sigma2, delta2 = delta2, mu_betas = betas_mean, sigma2_in_betas_prior = sigma2_in_betas_prior, segment_indexes = segment_indexes)
  
  # Sample delta2 | betas
  delta2 <- Delta2FCDSample(N = 1, betas = betas_included, sigma2 = sigma2, alpha = alpha, beta = beta, mu_betas = betas_mean, sigma2_in_betas_prior = sigma2_in_betas_prior, segment_indexes = segment_indexes)
  
  if (coupling){
    xi2 <- 1 # Add also xi2 update for delta2
  }
  
  #betas <- replace(numeric(length(betas)), model == 1, betas_included)
  if (is.list(betas_included)) {
    betas <- lapply(betas_included, function(b_in) replace(numeric(length(model)), model == 1, b_in))
  } else {
    betas <- replace(numeric(length(model)), model == 1, betas_included)
  }
  
  
  list(betas = betas, sigma2 = sigma2, delta2 = delta2)
}

RJ_MH_Marginal_move = function(Y,
                               X,
                               S = rep(1, length(Y)),
                               model = rep(1, ncol(X)),
                               betas = rep(0, ncol(X)),
                               sigma2 = 1,
                               delta2 = 1,
                               alpha = 0.01,
                               beta = 0.01,
                               a = 0.01,
                               b = 0.01,
                               r = 1,
                               lambda = 1,
                               mu_cross = rep(0, ncol(X)), 
                               Sigma_cross = diag(1, ncol(X)),
                               betas_mean = rep(0, ncol(X)),
                               sigma2_in_betas_prior = TRUE,
                               log = TRUE,
                               intercept_move = FALSE,
                               type_model_prior = 'Uniform',
                               max_parents = ncol(X),
                               allocation_vector = c(1, max(S)),
                               coupling = FALSE,
                               type_coupling = 'global',
                               group_segmentation = FALSE, 
                               segment_indexes = rep(1, length(Y))) {
  # Update Model
  updated_model <-
    ModelUpdate(
      k = ncol(X),
      model = model,
      max_parents = max_parents,
      intercept_move = intercept_move
    )
  model_star <- updated_model$model
  move_type <- updated_model$move
  
  if (!group_segmentation){
    segment_indexes <- if (length(allocation_vector) > 2) {
      cut(
        S,
        breaks = allocation_vector,
        labels = FALSE,
        right = FALSE,
        include.lowest = TRUE
      )
    } else{
      rep(1, length(Y))
    }
  }
  
  
  # Construct submatrices for current and proposed models
  X_included <- as.matrix(X[, model == 1])
  X_star <- as.matrix(X[, model_star == 1])
  
  betas_mean_included <-
    lapply(if (is.list(betas_mean))
      betas_mean
      else
        list(betas_mean), function(b)
          b[model == 1])
  # in case of global coupling here betas_mean_star is sampled from its FCD
  if (coupling & type_coupling == 'global'){
    betas_mean_star <- MuFCDSample(N = 1, Y = Y, X = X_star, segment_indexes = segment_indexes, sigma2 = sigma2, delta2 = delta2, mu_cross = mu_cross[model_star == 1], Sigma_cross = Sigma_cross[model_star == 1,model_star == 1])
    betas_mean_star <-
      lapply(1:max(segment_indexes), function(i) {
        betas_mean_star
      })
  }
  else if(!coupling){
    betas_mean_star <-
      lapply(if (is.list(betas_mean))
        betas_mean
        else
          list(betas_mean), function(b)
            b[model_star == 1])
  }
  A = AcceptanceRatio(
    Y = Y,
    X = X_included,
    X_star = X_star,
    model = model,
    model_star = model_star,
    MCMC_type = 'RJ-Marginal',
    delta2 = delta2,
    delta2_star = delta2,
    sigma2 = sigma2,
    sigma2_star = sigma2,
    betas_0 = betas_mean_included,
    betas_0_star = betas_mean_star,
    likelihood = FALSE,
    marginal = TRUE,
    model_prior = ifelse(type_model_prior == 'Uniform', FALSE, TRUE),
    beta_prior = FALSE,
    mu_prior = ifelse(coupling & type_coupling == 'global', TRUE, FALSE),
    mu_posterior = ifelse(coupling & type_coupling == 'global', TRUE, FALSE),
    a = a,
    b = b,
    r = r,
    lambda = lambda,
    mu_cross = mu_cross,
    Sigma_cross = Sigma_cross,
    log = log,
    sigma2_in_betas_prior = sigma2_in_betas_prior,
    move_type = move_type,
    intercept_move = intercept_move,
    type_model_prior = type_model_prior,
    max_parents = max_parents,
    segment_indexes = segment_indexes,
    segment_indexes_star = segment_indexes
  )
  # Accept/reject move
  if (min(A, 1) >= runif(1)) {
    model <- model_star
    
    # Resample betas under accepted model
    betas_star_included <-
      BetasFCDSample(
        N = 1,
        Y = Y,
        X = X_star,
        sigma2 = sigma2,
        delta2 = delta2,
        mu_betas = betas_mean_star,
        sigma2_in_betas_prior = TRUE,
        segment_indexes = segment_indexes
      )
    
    # Form full betas_star vector (with zeros for excluded covariates)
    if (is.list(betas_star_included)) {
      betas <-
        lapply(betas_star_included, function(b_in)
          replace(numeric(length(model_star)), model_star == 1, b_in))
    } else {
      betas <-
        replace(numeric(length(model_star)), model_star == 1, betas_star_included)
    }
    if (is.list(betas_mean_star)) {
      betas_mean <-
        lapply(betas_mean_star, function(b_in)
          replace(numeric(length(model_star)), model_star == 1, b_in))
    } else {
      betas_mean <-
        replace(numeric(length(model_star)), model_star == 1, betas_mean_star)
    }
  }
  list(betas = betas, model = model, betas_mean = betas_mean)
}


RJMCMC = function(data,
                  # data.frame object
                  segmentation = FALSE,
                  # if data segmentation should be applied
                  group_segmentation = FALSE,
                  # if grouping should be applied
                  segmentation_var = 'Time',
                  # variable on which data is segmented
                  group_segmentation_var = 'Sample_Id',
                  # variable on which data are grouped
                  n_segments = NaN,
                  # Fixed number of segments (if desired)
                  segmentation_vector = NaN,
                  # Fixed segmentation vector
                  group_segmentation_vector = NaN,
                  # Fixed grouping vector
                  epochs = 1000,
                  # number of epochs of RJMCMC sampling
                  v = 0.2,
                  # parameter for MH moves
                  p_MCMC = 1,
                  # probability of doing an MCMC move
                  p_RJMCMC = 0.1,
                  # probability of doing an RJMCMC move
                  p_ChangePoint = 0.1,
                  # probability of doing a Change Point Detection move
                  alpha = 0.01,
                  # delta2 scale parameter
                  beta =  0.01,
                  # delta2 rate parameter
                  alpha_xi = 0.01,
                  #xi2 scale parameter (only in sequentially-coupled settings)
                  beta_xi = 0.01,
                  #xi2 rate parameter (only in sequentially-coupled settings)
                  a = 0.01,
                  # sigma2 scale parameter
                  b = 0.01,
                  # sigma2 rate parameter
                  r = 1,
                  # number of failures of Neg-Bin model prior
                  lambda = 1,
                  # expected value of number of successes of Poisson/Neg-Bin model prior
                  rho = 0.1,
                  #geometric probability parameter for segment distance prior
                  iota = 1,
                  # expected value of number of successes of Poisson segments' granularity prior
                  intercept_move = FALSE,
                  # if intercept could be deleted too
                  outcome = "Y",
                  # outcome variable
                  covariates = colnames(data)[!colnames(data) %in% outcome],
                  # set of covariates (default: all the variables except outcome)
                  max_parents = length(covariates),
                  # max number of parents
                  max_changepoints = NaN,
                  #max number of changepoints
                  init_model = character(0),
                  # initial model vector
                  init_betas = rep(0, (length(covariates))),
                  # initial parameters betas
                  betas_mean = rep(0, (length(covariates))),
                  # prior mean of the betas
                  mu_cross = rep(0, (length(covariates))),
                  # global mean vector hyperprior mean (only in globally-coupled settings)
                  Sigma_cross = diag(1, nrow = length(covariates)),
                  # global mean vector hyperprior covariance matrix (only in globally-coupled settings)
                  init_sigma2 = 1,
                  # initial sigma2
                  init_delta2 = 1,
                  # initial delta2
                  init_xi2 = 1,
                  # initial xi2 (only in sequentially-coupled setting)
                  # init_psi2 = 1
                  # initial psi2 (only if grouping is applied)
                  init_allocation_vector = NaN,
                  # initial segment set
                  parameters_moves = 'Gibbs',
                  # type of move on the parameter domain ('MH' or 'Gibbs')
                  model_moves = 'MH-Marginal',
                  # type of move on the model domain ('MH-Prior', 'MH-FCD' or 'MH-Marginal')
                  type_model_prior = 'Poisson',
                  # type of model prior ('Uniform', 'Uniform_trunc', 'Uniform_card', 'Poisson' or 'Neg-Bin')
                  type_segment_prior = 'Uniform_card',
                  # type of prior for the number of segments ('Uniform_card', 'Geometric', 'Poisson_card' or 'Geometric+Poisson')
                  type_coupling = 'global',
                  # type of parameter coupling ('sequential' or 'global')
                  log = TRUE,
                  # if logarithm should be applied in acceptance ratio computation
                  coupling = FALSE,
                  # if parameter coupling should be applied
                  #group_coupling = FALSE
                  # if parameter coupling over groups is applied
                  sigma2_in_betas_prior = TRUE) {
  # if sigma2 is present in betas prior
  segment_indexes = group_segmentation_vector
  # Experimental Conditions
  X <- as.matrix(data[, (names(data) %in% covariates)])
  Y <- data[, outcome]
  k <- ncol(X)
  
  # Temporal Vector Initialization
  if (segmentation) {
    S <- data[, segmentation_var]
  }
  else{
    S <- rep(1, length(Y))
  }
  
  if (group_segmentation) {
    S <- data[, group_segmentation_var]
  }
  else{
    S <- rep(1, length(Y))
  }
  # Changepoint Maximum Number
  max_changepoints <-
    ifelse(is.nan(max_changepoints), max(unique(data[, segmentation_var])), max_changepoints)
  # Model Initialization
  if (length(init_model) == 0) {
    init_model <- rep(0, k)
    if (max_parents > 0 && max_parents <= k) {
      selected_indices <- sample(k, max_parents)
      init_model[selected_indices] <- 1
    }
  }
  if (intercept_move == FALSE) {
    init_model[1] <- 1
    max_parents <- max_parents + 1
  }
  
  # Allocation Vector initialization
  if (length(init_allocation_vector) == 0){
    if (is.nan(n_segments) &
        any(is.nan(segmentation_vector))) {
      # unknown changepoints number and changepoints locations
      init_allocation_vector <- c(1, max(S)) # DBN
    }
    else if (!any(is.nan(segmentation_vector))) {
      # known changepoints number and changepoints locations
      init_allocation_vector <-
        c(1, segmentation_vector, max(S)) # NH-DBN with fixed changepoints
      n_segments <- length(segmentation_vector) + 1
    }
    else{
      # known changepoints number and unknown changepoints locations
      init_allocation_vector <-
        round(seq(1, max(S), length.out = (n_segments + 1))) # NH-DBN with equi-spaced changepoints
    }
  }
  # Betas Mean Initialization
  betas_mean <-
    BetasMeanVector(
      Y = Y,
      X = X,
      S = S,
      model = init_model,
      allocation_vector = init_allocation_vector,
      init_betas = init_betas,
      betas_mean = betas_mean,
      sigma2 = init_sigma2,
      delta2 = init_delta2,
      xi2 = init_xi2,
      coupled = coupling,
      type_coupling = type_coupling,
      group_segmentation = group_segmentation,
      segment_indexes = segment_indexes
    )
  # Initial Parameters
  model = init_model
  sigma2 = init_sigma2
  betas = init_betas
  #betas = list(init_betas)
  #betas = append(betas, replicate(length(init_allocation_vector) - 2, rep(0, ncol(X)), simplify = FALSE))
  delta2 = init_delta2
  xi2 = init_xi2
  allocation_vector <- init_allocation_vector
  betas_list = list() #list(betas)
  sigma2s <- list() #list(sigma2)
  delta2s <- list() #list(delta2)
  models <- list() #list(model)
  mus <- list() #list(betas_mean)
  allocation_vectors <- list() #list(allocation_vector)
  # Sampling Beta parameters
  for (epoch in 1:epochs) {
    # Randomly select between an RJMCMC move and a standard MH move
    if (p_MCMC >= runif(1, min = 0, max = 1)) {
      if (parameters_moves == 'Gibbs') {
        param_set <-
          Gibbs_move(
            Y = Y,
            X = X,
            S = S,
            model = model,
            betas = betas,
            sigma2 = sigma2,
            delta2 = delta2,
            alpha = alpha,
            beta = beta,
            a = a,
            b = b,
            betas_mean = betas_mean,
            sigma2_in_betas_prior = sigma2_in_betas_prior,
            collapsed = ifelse(model_moves == 'MH-Marginal', TRUE, FALSE),
            allocation_vector = allocation_vector,
            coupling = coupling,
            group_segmentation = group_segmentation,
            segment_indexes = segment_indexes
          )
      }
      else if (parameters_moves == 'MH') {
        param_set <-
          MH_move(
            Y = Y,
            X = X,
            S = S,
            model = model,
            betas = betas,
            sigma2 = sigma2,
            delta2 = delta2,
            alpha = alpha,
            beta = beta,
            a = a,
            b = b,
            betas_mean = betas_mean,
            sigma2_in_betas_prior = sigma2_in_betas_prior,
            log = log,
            v = v,
            allocation_vector = allocation_vector,
            coupling = coupling
          )
      }
      sigma2 <- param_set$sigma2
      delta2 <- param_set$delta2
      betas <- param_set$betas
    }
    if (p_RJMCMC >= runif(1, min = 0, max = 1)) {
      if (model_moves == 'MH-Prior') {
        param_set <-
          RJ_MH_move(
            Y = Y,
            X = X,
            S = S,
            model = model,
            betas = betas,
            sigma2 = sigma2,
            delta2 = delta2,
            alpha = alpha,
            beta = beta,
            a = a,
            b = b,
            r = r,
            lambda = lambda,
            betas_mean = betas_mean,
            sigma2_in_betas_prior = sigma2_in_betas_prior,
            log = log,
            intercept_move = intercept_move,
            type_model_prior = type_model_prior,
            max_parents = max_parents,
            allocation_vector = allocation_vector,
            coupling = coupling,
            type_coupling = type_coupling,
            mu_cross = rep(0, ncol(X)),
            Sigma_cross = diag(1, ncol(X))
          )
      }
      else if (model_moves == 'MH-FCD') {
        param_set <-
          RJ_MH_FCD_move(
            Y = Y,
            X = X,
            S = S,
            model = model,
            betas = betas,
            sigma2 = sigma2,
            delta2 = delta2,
            alpha = alpha,
            beta = beta,
            a = a,
            b = b,
            r = r,
            lambda = lambda,
            betas_mean = betas_mean,
            sigma2_in_betas_prior = sigma2_in_betas_prior,
            log = log,
            intercept_move = intercept_move,
            type_model_prior = type_model_prior,
            max_parents = max_parents,
            allocation_vector = allocation_vector,
            coupling = coupling,
            type_coupling = type_coupling,
            mu_cross = rep(0, ncol(X)),
            Sigma_cross = diag(1, ncol(X))
          )
      }
      else if (model_moves == 'MH-Marginal') {
        param_set <-
          RJ_MH_Marginal_move(
            Y = Y,
            X = X,
            S = S,
            model = model,
            betas = betas,
            sigma2 = sigma2,
            delta2 = delta2,
            alpha = alpha,
            beta = beta,
            a = a,
            b = b,
            r = r,
            lambda = lambda,
            betas_mean = betas_mean,
            sigma2_in_betas_prior = sigma2_in_betas_prior,
            log = log,
            intercept_move = intercept_move,
            type_model_prior = type_model_prior,
            max_parents = max_parents,
            allocation_vector = allocation_vector,
            coupling = coupling,
            type_coupling = type_coupling,
            mu_cross = mu_cross,
            Sigma_cross = Sigma_cross,
            group_segmentation = group_segmentation,
            segment_indexes = segment_indexes
          )
      }
      model <- param_set$model
      betas <- param_set$betas
      betas_mean <- param_set$betas_mean
      # betas_mean should be updated here, it represents the value of mu in case of global coupling
    }
    # if segmentation is not selected or segmentation vector is fixed than the changepoint move is avoided
    if (((segmentation) &
         (any(is.nan(
           segmentation_vector
         )))) & (p_ChangePoint >= runif(1, min = 0, max = 1))) {
      param_set <-
        ChangePoint_move(
          Y = Y,
          X = X,
          S = S,
          model = model,
          betas = betas,
          sigma2 = sigma2,
          delta2 = delta2,
          a = a,
          b = b,
          rho = rho,
          iota = iota,
          betas_mean = betas_mean,
          sigma2_in_betas_prior = sigma2_in_betas_prior,
          log = log,
          max_changepoints = max_changepoints,
          allocation_vector = allocation_vector,
          coupling = coupling,
          type_segment_prior = type_segment_prior,
          n_segments = n_segments,
          mu_cross = rep(0, ncol(X)),
          Sigma_cross = diag(1, ncol(X))
        )
      allocation_vector <- param_set$allocation_vector
      betas <- param_set$betas
      betas_mean <- param_set$betas_mean
      # betas_mean should be updated here, not only due to segment dependency but due to the joint update with the changepoint set
    }
    # Updating Markov Chain
    betas_list[[length(betas_list) + 1]] <- betas
    mus[[length(mus) + 1]] <- betas_mean #substitute with betas_mean vectors
    sigma2s[[length(sigma2s) + 1]] <- sigma2
    delta2s[[length(delta2s) + 1]] <- delta2
    models[[length(models) + 1]] <- model
    allocation_vectors[[length(allocation_vectors) + 1]] <-
      allocation_vector
    if (epoch %% 1000 == 0) {
      cat("Epoch:", epoch, "\n")
    }
  }
  #cat('\n\n')
  # Return thinned samples post-burn-in
  list(
    betas = betas_list,
    mus = mus,
    sigma2s = sigma2s,
    delta2s = delta2s,
    models = models,
    allocation_vectors = allocation_vectors,
    covariates = covariates
  )
}

B_Init_group <- function(df, init = 'zeros', sample_id = 'Sample_Id', time_id = 'Time') {
  # Extract variable names (excluding sample_id, time_id, intercept column)
  variable_names <- setdiff(colnames(df), c(sample_id, time_id, "1"))
  
  k <- length(variable_names)
  
  # Get unique time points and their count
  N <- length(unique(df[[sample_id]]))
  
  # Create 3D array: [target variables, intercept + predictors, time]
  if (init == 'zeros') {
    B <- array(0, dim = c(k, k + 1, N),
               dimnames = list(variable_names, c("1", variable_names), 1:N))
  } 
  else {
    B <- init
  }
  
  B
}

PosteriorSamplesMatricesG <- function(posterior_samples, T, group_segmentation_list) {
  
  vars <- names(posterior_samples)
  k <- length(vars)
  G <- max(group_segmentation_list)
  
  all_covariates <- sort(unique(unlist(lapply(posterior_samples, function(x) x$covariates))))
  p <- length(all_covariates)
  
  B <- array(0, dim = c(k, p, G),
             dimnames = list(vars, all_covariates, as.character(1:G)))
  
  A <- matrix(0, nrow = k, ncol = T - 1,
              dimnames = list(vars, as.character(1:(T - 1))))
  
  M <- matrix(0, nrow = k, ncol = p,
              dimnames = list(vars, all_covariates))
  
  MU <- array(0, dim = c(k, p, G),
              dimnames = list(vars, all_covariates, as.character(1:G)))
  
  Sigma <- setNames(rep(0, k), vars)
  Delta <- setNames(rep(0, k), vars)
  
  for (i in seq_along(vars)) {
    
    v <- vars[i]
    post <- posterior_samples[[v]]
    covars <- post$covariates
    
    ## ---------------------------
    ## Change points (A matrix)
    ## ---------------------------
    allocation <- tail(post$allocation_vectors, 1)[[1]]
    
    if (length(allocation) > 2) {
      cp_positions <- allocation[-c(1, length(allocation))]
      A[i, cp_positions] <- 1
    }
    
    ## ---------------------------
    ## Betas + Mus (FIXED SECTION)
    ## ---------------------------
    last_betas_list <- tail(post$betas, 1)[[1]]
    last_mus_list   <- tail(post$mus, 1)[[1]]
    
    for (s in 1:G) {
      
      beta_segment <- last_betas_list[[s]]
      mu_segment   <- last_mus_list[[s]]
      
      # SAFETY CHECK: lengths must match
      if (length(beta_segment) != length(covars)) {
        stop(paste("Length mismatch in beta for variable:", v, "segment:", s))
      }
      
      if (length(mu_segment) != length(covars)) {
        stop(paste("Length mismatch in mu for variable:", v, "segment:", s))
      }
      
      # POSITIONAL mapping (SAFE)
      for (j in seq_along(covars)) {
        cov <- covars[j]
        
        B[i, cov, s]  <- beta_segment[j]
        MU[i, cov, s] <- mu_segment[j]
      }
    }
    
    ## ---------------------------
    ## Model matrix M
    ## ---------------------------
    last_model <- tail(post$models, 1)[[1]]
    included_idx <- which(last_model == 1)
    
    if (length(included_idx) > 0) {
      M[i, covars[included_idx]] <- 1
    }
    
    ## ---------------------------
    ## Variance parameters
    ## ---------------------------
    Sigma[v] <- tail(post$sigma2s, 1)[[1]]
    Delta[v] <- tail(post$delta2s, 1)[[1]]
  }
  
  list(B = B, MU = MU, M = M, A = A, Sigma = Sigma, Delta = Delta)
}

MuFCDSample = function(N = 1, Y, X, segment_indexes, sigma2, delta2, mu_cross, Sigma_cross){
  # Posterior covariance
  Precision_post = Reduce('+', lapply(1:max(segment_indexes), function(h){
    X_current <- as.matrix(X[segment_indexes == h, , drop = FALSE])
    Y_current <- Y[segment_indexes == h]
    nh <- nrow(X_current)
    
    #Vh <- sigma2 * (diag(nh) + delta2 * X_current %*% t(X_current))
    
    # Solve Vh %*% Z = X_current via Cholesky
    #chol_Vh <- chol(Vh)
    #Vinv_X <- backsolve(chol_Vh, forwardsolve(t(chol_Vh), X_current)) # Vh^{-1} X
    
    #t(X_current) %*% Vinv_X
    
    ## PROVA
    XtX <- crossprod(X_current)
    A   <- diag(ncol(X_current)) + delta2 * XtX
    Ainv <- chol2inv(chol(A))
    
    (1 / sigma2) * (XtX - delta2 * XtX %*% Ainv %*% XtX)
  })) + solve(Sigma_cross)
  
  Sigma_cross2 <- solve(Precision_post)
  
  # Posterior mean
  Linear_post = Reduce('+', lapply(1:max(segment_indexes), function(h){
    X_current <- as.matrix(X[segment_indexes == h, , drop = FALSE])
    Y_current <- Y[segment_indexes == h]
    nh <- nrow(X_current)
    
    XtX <- crossprod(X_current)
    XtY <- crossprod(X_current, Y_current)
    
    A   <- diag(ncol(X_current)) + delta2 * XtX
    Ainv <- chol2inv(chol(A))
    
    (1 / sigma2) * (XtY - delta2 * XtX %*% Ainv %*% XtY)
  })) + solve(Sigma_cross) %*% mu_cross
  
  mu_cross2 <- Sigma_cross2 %*% Linear_post
  MultivariateNormalSample(N = N, mu = mu_cross2, Sigma = Sigma_cross2)
}

MuFCDDensity = function(mu_current, Y, X, segment_indexes, sigma2, delta2, mu_cross, Sigma_cross, log = TRUE){
  # Posterior covariance
  Precision_post = Reduce('+', lapply(1:max(segment_indexes), function(h){
    X_current <- as.matrix(X[segment_indexes == h, , drop = FALSE])
    Y_current <- Y[segment_indexes == h]
    nh <- nrow(X_current)
  
    XtX <- crossprod(X_current)
    A   <- diag(ncol(X_current)) + delta2 * XtX
    Ainv <- chol2inv(chol(A))
    
    (1 / sigma2) * (XtX - delta2 * XtX %*% Ainv %*% XtX)
  })) + solve(Sigma_cross)
  
  Sigma_cross2 <- solve(Precision_post)
  
  # Posterior mean
  Linear_post = Reduce('+', lapply(1:max(segment_indexes), function(h){
    X_current <- as.matrix(X[segment_indexes == h, , drop = FALSE])
    Y_current <- Y[segment_indexes == h]
    nh <- nrow(X_current)
    
    XtX <- crossprod(X_current)
    XtY <- crossprod(X_current, Y_current)
    
    A   <- diag(ncol(X_current)) + delta2 * XtX
    Ainv <- chol2inv(chol(A))
    
    (1 / sigma2) * (XtY - delta2 * XtX %*% Ainv %*% XtY)
  })) + solve(Sigma_cross) %*% mu_cross
  
  mu_cross2 <- Sigma_cross2 %*% Linear_post
  MultivariateNormalDensity(X = as.vector(mu_current), mu = as.vector(mu_cross2), Sigma = Sigma_cross2, log = log)
}

PosteriorRatio = function(Y,
                          X,
                          X_star = X,
                          model = rep(1, ncol(X)),
                          model_star = model,
                          betas = rep(0, ncol(X)),
                          betas_star = betas,
                          sigma2 = 1,
                          sigma2_star = sigma2,
                          delta2 = 1,
                          delta2_star = delta2,
                          betas_0 = rep(0, ncol(X)),
                          betas_0_star = betas_0,
                          a = 0.01,
                          b = 0.01,
                          alpha = 0.01,
                          beta = 0.01,
                          lambda = 1,
                          r = 1,
                          rho = 0.01,
                          iota = 1,
                          mu_cross = rep(0, (ncol(X))),
                          Sigma_cross = diag(1, nrow = ncol(X)),
                          type_model_prior = 'Uniform',
                          type_segment_prior = 'Uniform_card',
                          intercept_move = intercept_move,
                          max_parents = length(model) - 1 + intercept_move,
                          max_changepoints = 0,
                          likelihood = TRUE,
                          marginal = FALSE,
                          model_prior = FALSE,
                          segment_prior = FALSE,
                          beta_prior = TRUE,
                          sigma_prior = FALSE,
                          delta_prior = FALSE,
                          mu_posterior = FALSE,
                          mu_prior = FALSE,
                          log = TRUE,
                          sigma2_in_betas_prior = TRUE,
                          segment_indexes = rep(1, length(Y)),
                          segment_indexes_star = segment_indexes,
                          allocation_vector = c(1, max(1)),
                          allocation_vector_star = allocation_vector) {
  if (log) {
    MarginalRate = ifelse(
      marginal,
      LogMarginal(
        Y = Y,
        X = X_star,
        delta2 = delta2_star,
        mu_beta = betas_0_star,
        a = a,
        b = b,
        segment_indexes = segment_indexes_star
      ) - LogMarginal(
        Y = Y,
        X = X,
        delta2 = delta2,
        mu_beta = betas_0,
        a = a,
        b = b,
        segment_indexes = segment_indexes
      ),
      0
    )
    LikelihoodRate = ifelse(
      likelihood,
      Likelihood(
        Y = Y,
        X = X_star,
        betas = betas_star,
        sigma2 = sigma2_star,
        log = TRUE,
        segment_indexes = segment_indexes_star
      ) - Likelihood(
        Y = Y,
        X = X,
        betas = betas,
        sigma2 = sigma2,
        log = TRUE,
        segment_indexes = segment_indexes
      ),
      0
    )
    BetaPriorRate = ifelse(
      beta_prior,
      BetaPrior(
        betas = betas_star,
        betas_0 = betas_0_star,
        sigma2 = sigma2_star,
        delta2 = delta2_star,
        sigma2_in_betas_prior = sigma2_in_betas_prior,
        log = TRUE,
        segment_indexes = segment_indexes_star
      ) - BetaPrior(
        betas = betas,
        betas_0 = betas_0,
        sigma2 = sigma2,
        delta2 = delta2,
        sigma2_in_betas_prior = sigma2_in_betas_prior,
        log = TRUE,
        segment_indexes = segment_indexes
      ),
      0
    )
    SigmaPriorRate = ifelse(
      sigma_prior,
      InverseGammaDensity(
        sigma2_star,
        a = a,
        b = b,
        log = TRUE
      ) - InverseGammaDensity(
        sigma2,
        a = a,
        b = b,
        log = TRUE
      ),
      0
    )
    DeltaPriorRate = ifelse(
      delta_prior,
      InverseGammaDensity(
        delta2_star,
        a = alpha,
        b = beta,
        log = TRUE
      ) - InverseGammaDensity(
        delta2,
        a = alpha,
        b = beta,
        log = TRUE
      ),
      0
    )
    ModelPriorRate = ifelse(
      model_prior,
      ModelPrior(
        model_star,
        type = type_model_prior,
        lambda = lambda,
        r = r,
        max_parents = max_parents,
        intercept_move = intercept_move,
        log = TRUE
      ) - ModelPrior(
        model,
        type = type_model_prior,
        lambda = lambda,
        r = r,
        max_parents = max_parents,
        intercept_move = intercept_move,
        log = TRUE
      ),
      0
    )
    SegmentPriorRate = ifelse(
      segment_prior,
      SegmentPrior(
        allocation_vector = allocation_vector_star,
        type_segment_prior = type_segment_prior,
        rho = rho,
        iota = iota,
        max_changepoints = max_changepoints,
        log = TRUE
      ) - SegmentPrior(
        allocation_vector = allocation_vector,
        type_segment_prior = type_segment_prior,
        rho = rho,
        iota = iota,
        max_changepoints = max_changepoints,
        log = TRUE
      ),
      0
    )
    #print(Sigma_cross[model== 1, model == 1])
    #print(Sigma_cross[model_star == 1, model_star == 1])
    MuPosteriorRate = ifelse(
      mu_posterior,
      MuFCDDensity(
        mu_current = betas_0[[1]],
        Y = Y,
        X = X,
        segment_indexes = segment_indexes,
        sigma2 = sigma2,
        delta2 = delta2,
        mu_cross = mu_cross[model == 1],
        Sigma_cross = Sigma_cross[model == 1, model == 1],
        log = TRUE
      ) - MuFCDDensity(
        mu_current = betas_0_star[[1]],
        Y = Y,
        X = X_star,
        segment_indexes = segment_indexes_star,
        sigma2 = sigma2_star,
        delta2 = delta2_star,
        mu_cross = mu_cross[model_star == 1],
        Sigma_cross = Sigma_cross[model_star == 1, model_star == 1],
        log = TRUE
      ),
      0
    )
    #print(Sigma_cross[model== 1, model == 1])
    #print(Sigma_cross[model_star == 1, model_star == 1])
    MuPriorRate = ifelse(
      mu_prior,
      MultivariateNormalDensity(
        betas_0_star[[1]],
        mu = mu_cross[model_star == 1],
        Sigma = Sigma_cross[model_star == 1, model_star == 1],
        log = TRUE
      ) - MultivariateNormalDensity(
        betas_0[[1]],
        mu = mu_cross[model == 1],
        Sigma = Sigma_cross[model == 1, model == 1],
        log = TRUE
      ),
      0
    )
    #if(all(allocation_vector_star == c(1,10))){
    #  cat(sprintf(
    #    "Marginal: %.4f | Likelihood: %.4f | BetaPrior: %.4f | SigmaPrior: %.4f | DeltaPrior: %.4f | ModelPrior: %.4f | SegmentPrior: %.4f | MuPosterior: %.4f | MuPrior: %.4f | TOTAL: %.4f\n",
    #    MarginalRate, LikelihoodRate, BetaPriorRate, SigmaPriorRate, DeltaPriorRate,
    #    ModelPriorRate, SegmentPriorRate, MuPosteriorRate, MuPriorRate,
    #    MarginalRate + LikelihoodRate + BetaPriorRate + SigmaPriorRate +
    #      DeltaPriorRate + ModelPriorRate + SegmentPriorRate +
    #      MuPosteriorRate + MuPriorRate
    #  ))
    #}
    PosteriorRate = MarginalRate + LikelihoodRate + BetaPriorRate + SigmaPriorRate + DeltaPriorRate + ModelPriorRate + SegmentPriorRate + MuPosteriorRate + MuPriorRate
  }
  else{
    MarginalRate = ifelse(
      marginal,
      Marginal(
        Y = Y,
        X = X_star,
        delta2 = delta2,
        mu_beta = betas_0_star,
        a = a,
        b = b
      ) / Marginal(
        Y = Y,
        X = X,
        delta2 = delta2,
        mu_beta = betas_0,
        a = a,
        b = b
      ),
      1
    )
    LikelihoodRate = ifelse(
      likelihood,
      MultivariateNormalDensity(
        Y,
        mu = X_star %*% betas_star,
        Sigma = sigma2_star,
        log = FALSE
      ) / MultivariateNormalDensity(
        Y,
        mu = X %*% betas,
        Sigma = sigma2,
        log = FALSE
      ),
      1
    )
    BetaPriorRate = ifelse(
      beta_prior,
      MultivariateNormalDensity(
        betas_star,
        mu = betas_0,
        Sigma = BetaNoise_star,
        log = FALSE
      ) / MultivariateNormalDensity(
        betas,
        mu = betas_0,
        Sigma = BetaNoise,
        log = FALSE
      ),
      1
    )
    SigmaPriorRate = ifelse(
      sigma_prior,
      InverseGammaDensity(
        sigma2_star,
        a = a,
        b = b,
        log = FALSE
      ) / InverseGammaDensity(
        sigma2,
        a = a,
        b = b,
        log = FALSE
      ),
      1
    )
    DeltaPriorRate = ifelse(
      delta_prior,
      InverseGammaDensity(
        delta2_star,
        a = alpha,
        b = beta,
        log = FALSE
      ) / InverseGammaDensity(
        delta2,
        a = alpha,
        b = beta,
        log = FALSE
      ),
      1
    )
    ModelPriorRate = ifelse(
      model_prior,
      ModelPrior(
        model_star,
        type = type_model_prior,
        lambda = lambda,
        r = r,
        max_parents = max_parents,
        intercept_move = intercept_move,
        log = FALSE
      ) / ModelPrior(
        model,
        type = type_model_prior,
        lambda = lambda,
        r = r,
        max_parents = max_parents,
        intercept_move = intercept_move,
        log = FALSE
      ),
      1
    )
    SegmentPriorRate = ifelse(
      segment_prior,
      SegmentPrior(
        allocation_vector = allocation_vector_star,
        type_segment_prior = type_segment_prior,
        rho = rho,
        iota = iota,
        log = FALSE
      ) / SegmentPrior(
        allocation_vector = allocation_vector,
        type_segment_prior = type_segment_prior,
        rho = rho,
        iota = iota,
        log = FALSE
      ),
      1
    )
    MuPosteriorRate = ifelse(
      mu_posterior,
      MuFCDDensity(
        mu_current = betas_0[[1]],
        Y = Y,
        X = X,
        segment_indexes = segment_indexes,
        sigma2 = sigma2,
        delta2 = delta2,
        mu_cross = mu_cross[model == 1],
        Sigma_cross = Sigma_cross[model == 1,model == 1],
        log = FALSE
      ) / MuFCDDensity(
        mu_current = betas_0_star[[1]],
        Y = Y,
        X = X_star,
        segment_indexes = segment_indexes_star,
        sigma2 = sigma2_star,
        delta2 = delta2_star,
        mu_cross = mu_cross[model_star == 1],
        Sigma_cross = Sigma_cross[model_star == 1,model_star == 1],
        log = FALSE
      ),
      1
    )
    MuPriorRate = ifelse(
      mu_prior,
      MultivariateNormalDensity(
        betas_0_star[[1]],
        mu = mu_cross[model_star == 1],
        Sigma = Sigma_cross[model_star == 1,model_star == 1],
        log = FALSE
      ) / MultivariateNormalDensity(
        betas_0[[1]],
        mu = mu_cross[model == 1],
        Sigma = Sigma_cross[model = 1,model == 1],
        log = FALSE
      ),
      1
    )
    PosteriorRate = LikelihoodRate * BetaPriorRate * SigmaPriorRate * DeltaPriorRate * ModelPriorRate * SegmentsPriorRate * MuPosteriorRate * MuPriorRate
  }
  PosteriorRate
}

AcceptanceRatio = function(Y,
                           X,
                           X_star = X,
                           model = rep(1, ncol(X)),
                           model_star = model,
                           MCMC_type = 'Gibbs',
                           betas = rep(0, ncol(X)),
                           betas_star = rep(0, ncol(X_star)),
                           sigma2 = 1,
                           sigma2_star = 1,
                           delta2 = 1,
                           delta2_star = 1,
                           betas_0 = rep(0, ncol(X)),
                           betas_0_star = betas_0,
                           a = 0.01,
                           b = 0.01,
                           alpha = 0.01,
                           beta = 0.01,
                           lambda = 1,
                           r = 1,
                           rho = 0.01,
                           iota = 1,
                           mu_cross = rep(0, (ncol(X))),
                           Sigma_cross = diag(1, nrow = ncol(X)),
                           likelihood = TRUE,
                           marginal = FALSE,
                           model_prior = FALSE,
                           segment_prior = FALSE,
                           beta_prior = TRUE,
                           sigma_prior = FALSE,
                           delta_prior = FALSE,
                           mu_prior = FALSE,
                           mu_posterior = FALSE,
                           log = TRUE,
                           sigma2_in_betas_prior = TRUE,
                           move_type = 1,
                           type_model_prior = 'Uniform',
                           type_segment_prior = 'Uniform_card',
                           intercept_move = FALSE,
                           max_parents = length(model) - 1 + intercept_move,
                           max_changepoints = 0,
                           segment_indexes = rep(1, length(Y)),
                           segment_indexes_star = segment_indexes,
                           allocation_vector = c(1, max(1)),
                           allocation_vector_star = allocation_vector) {
  if (log) {
    exp(
      PosteriorRatio(
        Y = Y,
        X = X,
        X_star = X_star,
        model = model,
        model_star = model_star,
        betas = betas,
        betas_star = betas_star,
        sigma2 = sigma2,
        sigma2_star = sigma2_star,
        delta2 = delta2,
        delta2_star = delta2_star,
        betas_0 = betas_0,
        betas_0_star = betas_0_star,
        a = a,
        b = b,
        alpha = alpha,
        beta = beta,
        lambda = lambda,
        r = r,
        mu_cross = mu_cross,
        Sigma_cross = Sigma_cross,
        likelihood = likelihood,
        marginal = marginal,
        model_prior = model_prior,
        beta_prior = beta_prior,
        sigma_prior = sigma_prior,
        delta_prior = delta_prior,
        log = TRUE,
        sigma2_in_betas_prior = sigma2_in_betas_prior,
        type_model_prior = type_model_prior,
        intercept_move = intercept_move,
        max_parents = max_parents,
        segment_indexes = segment_indexes,
        segment_indexes_star = segment_indexes_star,
        allocation_vector = allocation_vector,
        allocation_vector_star = allocation_vector_star,
        segment_prior = segment_prior,
        mu_prior = mu_prior,
        mu_posterior = mu_posterior,
        rho = rho,
        iota = iota,
        type_segment_prior = type_segment_prior,
        max_changepoints = max_changepoints
      ) + ProposalRatio(
        Y = Y,
        X = X,
        X_star = X_star,
        model = model,
        betas = betas,
        betas_star = betas_star,
        sigma2 = sigma2,
        delta2 = delta2,
        betas_0 = betas_0,
        betas_0_star = betas_0_star,
        type = MCMC_type,
        move_type = move_type,
        intercept_move = intercept_move,
        max_parents = max_parents,
        log = TRUE,
        segment_indexes = segment_indexes,
        segment_indexes_star = segment_indexes_star,
        allocation_vector = allocation_vector
      )
    )
  }
  else{
    PosteriorRatio(
      Y = Y,
      X = X,
      X_star = X_star,
      model = model,
      model_star = model_star,
      betas = betas,
      betas_star = betas_star,
      sigma2 = sigma2,
      sigma2_star = sigma2_star,
      delta2 = delta2,
      delta2_star = delta2_star,
      betas_0 = betas_0,
      betas_0_star = betas_0_star,
      a = a,
      b = b,
      alpha = alpha,
      beta = beta,
      lambda = lambda,
      r = r,
      likelihood = likelihood,
      marginal = marginal,
      model_prior = model_prior,
      beta_prior = beta_prior,
      sigma_prior = sigma_prior,
      delta_prior = delta_prior,
      log = FALSE,
      sigma2_in_betas_prior = sigma2_in_betas_prior,
      type_model_prior = type_model_prior,
      intercept_move = intercept_move,
      max_parents = max_parents,
      segment_indexes = segment_indexes,
      segment_indexes_star = segment_indexes_star,
      allocation_vector = allocation_vector,
      allocation_vector_star = allocation_vector_star,
      segment_prior = segment_prior,
      rho = rho,
      iota = iota,
      type_segment_prior = type_segment_prior,
      max_changepoints = max_changepoints
    ) * ProposalRatio(
      Y = Y,
      X = X,
      X_star = X_star,
      model = model,
      betas = betas,
      betas_star = betas_star,
      sigma2 = sigma2,
      delta2 = delta2,
      betas_0 = betas_0,
      betas_0_star = betas_0_star,
      type = MCMC_type,
      move_type = move_type,
      intercept_move = intercept_move,
      max_parents = max_parents,
      log = FALSE,
      segment_indexes = segment_indexes,
      segment_indexes_star = segment_indexes_star,
      allocation_vector = allocation_vector
    )
  }
}

ChangePoint_move <-
  function(X,
           Y,
           S = rep(1, length(Y)),
           model = rep(1, ncol(X)),
           betas = rep(0, ncol(X)),
           sigma2 = 1,
           delta2 = 1,
           a = 0.01,
           b = 0.01,
           rho = 0.01,
           iota = 1,
           betas_mean = rep(0, ncol(X)),
           type_segment_prior = 'Uniform_card',
           allocation_vector = c(1, max(S)),
           max_changepoints = max(S),
           n_segments = NaN,
           sigma2_in_betas_prior = TRUE,
           log = TRUE,
           coupling = FALSE,
           type_coupling = 'global',
           mu_cross = rep(0, ncol(X)), 
           Sigma_cross = diag(1, ncol(X))) {
    # Update Changepoints
    params <-
      UpdateChangepoints(
        allocation_vector = allocation_vector,
        max_changepoints = max_changepoints,
        n_segments = n_segments
      )
    move_type <- params$move
    allocation_vector_star <- params$allocation_vector
    segment_indexes <- if (length(allocation_vector) > 2) {
      cut(
        S,
        breaks = allocation_vector,
        labels = FALSE,
        right = FALSE,
        include.lowest = TRUE
      )
    } else{
      rep(1, length(Y))
    }
    segment_indexes_star <- if (length(allocation_vector_star) > 2) {
      cut(
        S,
        breaks = allocation_vector_star,
        labels = FALSE,
        right = FALSE,
        include.lowest = TRUE
      )
    } else{
      rep(1, length(Y))
    }
    X_included <- X[, model == 1, drop = FALSE]
    betas_included <-
      lapply(if (is.list(betas))
        betas
        else
          list(betas), function(b)
            b[model == 1])
    betas_mean_included <-
      lapply(if (is.list(betas_mean))
        betas_mean
        else
          list(betas_mean), function(b)
            b[model == 1])
    
    if (coupling & type_coupling == 'global'){
      betas_mean_star <- MuFCDSample(N = 1, Y = Y, X = X_included, segment_indexes = segment_indexes_star, sigma2 = sigma2, delta2 = delta2, mu_cross = mu_cross[model == 1], Sigma_cross = Sigma_cross[model == 1,model == 1])
      betas_mean_star <-
        lapply(1:max(segment_indexes_star), function(i) {
          betas_mean_star
        })
    }
    else if(!coupling){
      betas_mean_star <-
        lapply(1:max(segment_indexes_star), function(i) {
          betas_mean_included[[1]]
        })
    }
    #print(allocation_vector_star)
    # betas mean star in case of global coupling should be sampled and reformatted here
    A = AcceptanceRatio(
      Y = Y,
      X = X_included,
      X_star = X_included,
      model = model,
      model_star = model,
      MCMC_type = 'ChangePoint',
      delta2 = delta2,
      delta2_star = delta2,
      sigma2 = sigma2,
      sigma2_star = sigma2,
      betas_0 = betas_mean_included,
      betas_0_star = betas_mean_star,
      likelihood = FALSE,
      marginal = TRUE,
      beta_prior = FALSE,
      segment_prior = TRUE,
      mu_prior = ifelse(coupling & type_coupling == 'global', TRUE, FALSE),
      mu_posterior = ifelse(coupling & type_coupling == 'global', TRUE, FALSE),
      a = a,
      b = b,
      rho = rho,
      iota = iota,
      log = log,
      mu_cross = mu_cross,
      Sigma_cross = Sigma_cross,
      sigma2_in_betas_prior = sigma2_in_betas_prior,
      move_type = move_type,
      type_segment_prior = type_segment_prior,
      max_changepoints = max_changepoints,
      segment_indexes = segment_indexes,
      segment_indexes_star = segment_indexes_star,
      allocation_vector = allocation_vector,
      allocation_vector_star = allocation_vector_star
    )
    #if(all(allocation_vector_star == c(1,10))){
    #  print(A)
    #}
    # Accept/reject move
    if (min(A, 1) >= runif(1)) {
      allocation_vector <- allocation_vector_star
      # Resample betas under accepted model
      betas_included <-
        BetasFCDSample(
          N = 1,
          Y = Y,
          X = X_included,
          sigma2 = sigma2,
          delta2 = delta2,
          mu_betas = betas_mean_star,
          sigma2_in_betas_prior = TRUE,
          segment_indexes = segment_indexes_star
        )
      
      # Form full betas_star vector (with zeros for excluded covariates)
      if (is.list(betas_included)) {
        betas <-
          lapply(betas_included, function(b_in)
            replace(numeric(length(model)), model == 1, b_in))
      } else {
        betas <- replace(numeric(length(model)), model == 1, betas_included)
      }
      if (is.list(betas_mean)) {
        betas_mean <-
          lapply(betas_mean_star, function(b_in)
            replace(numeric(length(model)), model == 1, b_in))
      } else {
        betas_mean <-
          replace(numeric(length(model)), model == 1, betas_mean_star)
      }
    }
    list(
      allocation_vector = allocation_vector,
      betas = betas,
      betas_mean = betas_mean
    )
  }

B_Init <- function(df, init = 'zeros', sample_id = 'Sample_Id', time_id = 'Time') {
  # Extract variable names (excluding sample_id, time_id, intercept column)
  variable_names <- setdiff(colnames(df), c(sample_id, time_id, "1"))
  
  k <- length(variable_names)
  
  # Get unique time points and their count
  time_points <- sort(unique(df[[time_id]]))
  T <- length(time_points)
  
  # Create 3D array: [target variables, intercept + predictors, time]
  if (init == 'zeros') {
    B <- array(0, dim = c(k, k + 1, T),
               dimnames = list(variable_names, c("1", variable_names), time_points))
  } 
  else {
    B <- init
  }
  
  B
}


A_Init <- function(df, init = 'no changepoints', sample_id = 'Sample_Id', time_id = 'Time', max_changepoints = NULL) {
  # Extract variable names
  variable_names <- setdiff(colnames(df), c(sample_id, time_id, "1"))
  k <- length(variable_names)
  
  # Extract unique sorted time points
  time_points <- sort(unique(df[[time_id]]))
  T <- length(time_points)
  
  # Number of transitions (T - 1)
  n_transitions <- T - 1
  
  # Default max_changepoints
  if (is.null(max_changepoints)) {
    max_changepoints <- n_transitions
  }
  
  # Initialize matrix A (k × (T - 1)) for changepoint indicators
  if (init == 'no changepoints') {
    A <- matrix(0, nrow = k, ncol = n_transitions,
                dimnames = list(variable_names, as.character(1:n_transitions)))
  } else if (init == 'random') {
    A <- matrix(0, nrow = k, ncol = n_transitions)
    for (i in 1:k) {
      n_cp <- sample(0:min(max_changepoints, n_transitions), 1)
      if (n_cp > 0) {
        cp_locs <- sample(1:n_transitions, n_cp)
        A[i, cp_locs] <- 1
      }
    }
    dimnames(A) <- list(variable_names, as.character(1:n_transitions))
  } else {
    A <- init
  }
  
  A
}

ModelUpdate <- function(k, model = rep(1, k), max_parents = length(model), intercept_move = FALSE) {
  model_star <- model
  idx_start <- 2 - intercept_move
  model_subset <- model[idx_start:k]
  model_sum <- sum(model)
  
  available_add <- which(model_subset == 0) + (1 - intercept_move)
  available_del <- which(model_subset == 1) + (1 - intercept_move)
  
  if (model_sum == 1 - intercept_move) {
    move_type <- 1  # ADD
    sampled <- sample(available_add, 1)
    model_star[sampled] <- 1
  } else if (model_sum == k - intercept_move) {
    move_type <- 2  # DELETE
    sampled <- sample(available_del, 1)
    model_star[sampled] <- 0
  } else {
    move_type <- ifelse(model_sum == max_parents, sample(2:3, 1), sample(1:3, 1))
    if (move_type == 1) {
      sampled <- if (length(available_add) == 1) available_add else sample(available_add, 1)
      model_star[sampled] <- 1
    } else if (move_type == 2) {
      sampled <- if (length(available_del) == 1) available_del else sample(available_del, 1)
      model_star[sampled] <- 0
    } else {  # EXCHANGE
      add_sample <- if (length(available_add) == 1) available_add else sample(available_add, 1)
      del_sample <- if (length(available_del) == 1) available_del else sample(available_del, 1)
      model_star[add_sample] <- 1
      model_star[del_sample] <- 0
    }
  }
  
  list(model = model_star, move = move_type)
}

BinaryChangepoints <- function(changepoints, T = max(changepoints)) {
  internal_cp <- changepoints[changepoints > 1 & changepoints < T]
  binary_vec <- integer(T - 2)
  binary_vec[internal_cp - 1] <- 1
  binary_vec
}

ChangepointsList <- function(binary_vec) {
  internal_cps <- which(binary_vec == 1) + 1
  c(1, internal_cps, length(binary_vec) + 2)
}

ReallocationIndex <- function(x, i) {
  # Find left changepoint (exclusive)
  l <- i - 1
  while (l > 0 && x[l] == 0) l <- l - 1
  left_bound <- l + 1  # first non-changepoint to the right of previous CP
  
  # Find right changepoint (exclusive)
  r <- i + 1
  while (r <= length(x) && x[r] == 0) r <- r + 1
  right_bound <- r - 1  # last non-changepoint to the left of next CP
  
  candidates <- setdiff(left_bound:right_bound, i)
  if (length(candidates) == 0) {
    NaN
  }
  else if (length(candidates) == 1){
    candidates
  }  
  else{
    sample(candidates, 1)
  }
}

UpdateChangepoints <- function(allocation_vector, max_changepoints = (max(allocation_vector) - 2), n_segments = NaN) {
  allocation_vector_star <- BinaryChangepoints(allocation_vector)
  #print(allocation_vector_star)
  k <- length(allocation_vector_star)
  current_changepoints <- sum(allocation_vector_star)
  
  # Determine move type: 1 = BIRTH, 2 = DEATH, 3 = REALLOCATION
  if (!is.nan(n_segments)) {
    move_type <- 3  # Only reallocation allowed when n_segments is fixed
  } else if (current_changepoints == 0) {
    move_type <- 1  # Only BIRTH possible
  } else if (current_changepoints >= max_changepoints) {
    move_type <- sample(2:3, 1)  # Only DEATH or REALLOCATION
  } else {
    move_type <- sample(1:3, 1)  # Any move
  }
  
  # Perform the move
  if (move_type == 1) {  # BIRTH
    zero_indices <- which(allocation_vector_star == 0)
    if (length(zero_indices) > 0) {
      sampled_changepoint <- if (current_changepoints == (k - 1)) {
        zero_indices
      } else {
        sample(zero_indices, 1)
      }
      allocation_vector_star[sampled_changepoint] <- 1
    }
  } else if (move_type == 2) {  # DEATH
    one_indices <- which(allocation_vector_star == 1)
    if (length(one_indices) > 0) {
      sampled_changepoint <- if (current_changepoints == 1) {
        one_indices
      } else {
        sample(one_indices, 1)
      }
      allocation_vector_star[sampled_changepoint] <- 0
    }
  } else if (move_type == 3) {  # REALLOCATION
    one_indices <- which(allocation_vector_star == 1)
    if (length(one_indices) > 0) {
      sampled_changepoint <- if (length(one_indices) == 1) {
        one_indices
      } else {
        sample(one_indices, 1)
      }
      sampled_changepoint_2 <- ReallocationIndex(allocation_vector_star, sampled_changepoint)
      if (!is.nan(sampled_changepoint_2)) {
        allocation_vector_star[sampled_changepoint] <- 0
        allocation_vector_star[sampled_changepoint_2] <- 1
      }
    }
  }
  
  allocation_vector <- ChangepointsList(allocation_vector_star)
  list(allocation_vector = allocation_vector, move = move_type)
}

ProposalRatio <- function(Y, X, X_star = X, model = rep(1, ncol(X)), betas = rep(0, ncol(X)), betas_star = betas, sigma2 = 1, delta2 = 1, betas_0 = rep(0, ncol(X)), betas_0_star = betas_0, type = 'Gibbs', move_type = 1, max_parents = length(model) - 1 + intercept_move, intercept_move = FALSE, log = TRUE, segment_indexes = rep(1, length(Y)), segment_indexes_star = segment_indexes, allocation_vector = c(1, 1)){
  p <- sum(model)
  k <- length(model)
  TT <- max(allocation_vector)
  tau_card <- length(allocation_vector) - 2
  if (type %in% c('MH','Gibbs')){
    if (log){
      0
    }
    else{
      1
    }
  }
  else if (type %in% c('RJ-FCD')){
    if (move_type == 1){ # addition
      cov_proposal <- ifelse(log, log(k - p) - log(p - intercept_move), (k - p) / (p - intercept_move))
    }
    else if (move_type == 2){ # deletion
      cov_proposal <- ifelse(log, log(p - 1 + intercept_move) - log(k - p + 1), (p - 1 + intercept_move) / (k - p + 1))
    }
    else{ # exchange
      cov_proposal <- 0
    }
    # Log-posterior ratio (likelihood + prior + proposal ratio)
    (- BetasFCDDensity(Y = Y, X = X_star, betas = betas_star, sigma2 = sigma2, delta2 = delta2, mu_betas = betas_0_star, sigma2_in_betas_prior = FALSE, log = log, segment_indexes = segment_indexes_star) +
        BetasFCDDensity(Y = Y, X = X, betas = betas, sigma2 = sigma2, delta2 = delta2, mu_betas = betas_0, sigma2_in_betas_prior = FALSE, log = log, segment_indexes = segment_indexes) + cov_proposal)
  }
  else if (type %in% c('RJ-Prior','RJ-Marginal')){
    if (move_type == 1){ # addition
      ifelse(log, log(k - p) - log(p - intercept_move), (k - p) / (p - intercept_move))
    }
    else if (move_type == 2){ # deletion
      ifelse(log, log(p - 1 + intercept_move) - log(k - p + 1), (p - 1 + intercept_move) / (k - p + 1))
    }
    else{ # reversal
      ifelse(log,0,1)
    }
  }
  else if (type %in% c('ChangePoint')){
    if (move_type == 1){ #birth
      ifelse(log, log(TT - 2 - tau_card) - log(tau_card + 1), (TT - tau_card) / (tau_card + 1))
    }
    else if (move_type == 2){ #death
      ifelse(log, log(tau_card) - log(TT - 1 - tau_card), (tau_card) / (TT - 1 - tau_card))
    }
    else{ # reallocation
      ifelse(log,0,1)
    }
  }
}

