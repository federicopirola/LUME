## R Files

* **`distributions.R`** — Contains the probability distributions used in LUME-DBN and its extensions, including functions for sampling and evaluating their densities.

* **`data_model_generation.R`** — Contains functions for generating DBN, NH-DBN, and GH-DBN models, as well as simulating data from homogeneous and non-homogeneous models.

* **`model_parameter_update.R`** — Contains functions for parameter updates using Gibbs sampling and structural updates using Metropolis–Hastings moves.

* **`missing_imputation.R`** — Contains functions for the initialization and sampling of missing values using LUME-DBN and its non-homogeneous extensions.

* **`group_update.R`** — Contains functions for group and patient selection and for sampling group allocation vectors.

* **`LUME.R`** — Contains the main functions for generating posterior samples of DBNs, NH-DBNs, and GH-DBNs using LUME-DBN and its extensions.

* **`postprocessing.R`** — Contains functions for post-processing MCMC chains, including burn-in removal, thinning, visualization, convergence diagnostics, and structural accuracy assessment.

* **`main.R`** — Provides examples of data generation from a GH-DBN and model learning using the LUME-DBN algorithm and the extended MCMC method for GH-DBNs.
