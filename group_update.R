GroupUpdate <- function(grouping_list, lagged_data, M, MU, Delta, Sigma, mu_cross, Sigma_cross, lambda_group = 1, group_var = 'Sample_Id', coupling = TRUE){
  N <- length(grouping_list)
  G <- max(grouping_list)
  outcomes <- names(M[,1])
  sampled_group <- sample(1:G, 1)
  sampled_idx <- sample(names(which(grouping_list == sampled_group)),1, replace = TRUE)
  only_member = length(which(grouping_list == sampled_group)) == 1
  if (only_member) {
    grouping_list <- grouping_list - as.integer(grouping_list > sampled_group)
    grouping_list[sampled_idx] = 1
    G = G - 1
  }
  reallocations <- sapply(1:(G + 1), function(x) {
    grouping_list_star <- grouping_list
    grouping_list_star[sampled_idx] = x
    last <- x == (G + 1)
    segmentation_vec <- sapply(lagged_data[[group_var]], function(x) {grouping_list_star[[as.character(x)]]})
    logMarginal <- sapply(outcomes, function(outcome){
      outcome_t <- paste0(outcome, '_t')
      covariates <- names(which(M[outcome,] == 1))
      covariates_t <- c('1', colnames(lagged_data)[colnames(lagged_data) %in% paste0(covariates, '_t-1')])
      logM <- LogMarginal(
        Y = lagged_data[[outcome_t]],
        X = as.matrix(lagged_data[, covariates_t]),
        delta2 = Delta[[outcome]],
        mu_beta = lapply(1:(G + last), function(x) MU[outcome,covariates,1]),
        segment_indexes = segmentation_vec
      )
      logM
    })
    ifelse(last, (log(G)), (log(N - G) - log(G) + log(sum(grouping_list_star == x) + 1))) + PoissonDensity(G + last, lambda = lambda_group) + sum(logMarginal)
  })
  reallocation_scores <- exp(reallocations - max(reallocations)) 
  x = sample(1:(G + 1),1,replace=TRUE,prob=reallocation_scores/sum(reallocation_scores))
  grouping_list[sampled_idx] = x
  G = max(grouping_list)
  MU <- array(
    0,
    dim = c(length(outcomes), (length(outcomes) + 1), G),
    dimnames = list(
      outcomes,
      c(1, outcomes),
      1:G
    )
  )
  if (coupling){
    for (outcome in outcomes){
      outcome_t <- paste0(outcome, '_t')
      covariates <- names(which(M[outcome,] == 1))
      covariates_t <- c('1', colnames(lagged_data)[colnames(lagged_data) %in% paste0(covariates, '_t-1')])
      model_inc <- M[outcome,] == 1
      mu_s <- MuFCDSample(N = 1,
                          Y = lagged_data[[outcome_t]],
                          X = as.matrix(lagged_data[, covariates_t]),
                          delta2 = Delta[[outcome]],
                          sigma2 = Sigma[[outcome]],
                          segment_indexes = sapply(lagged_data[[group_var]], function(x) {grouping_list[[as.character(x)]]}), #init_grouping_list
                          mu_cross = mu_cross[model_inc],
                          Sigma_cross = Sigma_cross[model_inc,model_inc]
      )
      MU[outcome, , ] <- sapply(1:G, function(i) {
        replace(
          numeric(length(M[outcome, ])),
          model_inc,
          mu_s
        )
      })
    }
  }
  list(grouping_list = grouping_list,
       MU = MU,
       B = MU,
       log_marginal = reallocations[x])
}