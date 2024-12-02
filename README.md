# R package deeptrafo

[![R build status](https://github.com/neural-structured-additive-learning/deeptransformation/workflows/R-CMD-check/badge.svg)](https://github.com/neural-structured-additive-learning/deeptransformation/actions) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This package provides function to fit

- Deep Conditional Transformation Models (Baumann et al., 2021)
- Ordinal neural network transformation models (Kook et al., 2022)
- Autoregressive Transformation Models (Ruegamer et al., 2022)
- Deep Interpretable Ensembles (Kook et al., 2022)

using `deeptrafo`.

# Installation

To install `deeptrafo` run
```r
remotes::install_github("neural-structured-additive-learning/deeptrafo")
```
for the development version or
```r
install.packages("deeptrafo")
```
for the stable version from CRAN. If you set up a Python environment
for the first time, install `reticulate` and run the `check_and_install()`
function from the `deepregression` package. This tries to install miniconda, TF
2.10.0, TFP 0.16 and keras 2.8.0.

# Troubleshooting

See
[here](https://github.com/neural-structured-additive-learning/deepregression/blob/main/README.md#troubleshooting)
for troubleshooting your Python/R installation.

# Mini example

We fit a small regression model with smooth `s` and neural network predictor 
`nn` for the ordinal outcome `rating` in the `wine` data from package 
`ordinal`. After fitting the model for a couple of epochs, we illustrate the 
methods for `'deeptrafo'` objects, cross-validation and ensembles.

```r
library("deeptrafo")

# Prepare the data
data("wine", package = "ordinal")
wine$z <- rnorm(nrow(wine))
wine$x <- rnorm(nrow(wine))

# Set up neural network architecture
nn <- \(x) x |>
    layer_dense(input_shape = 1L, units = 2L, activation = "relu") |>
    layer_dense(1L)

# Model formula and definition
fml <- rating ~ 0 + temp + contact + s(z, df = 3) + nn(x)
m <- deeptrafo(fml, wine, latent_distr = "logistic", monitor_metric = NULL,
    return_data = TRUE, list_of_deep_models = list(nn = nn))

# Overview
print(m)

# Fit
m %>% fit(epochs = 10, batch_size = nrow(wine))

# Coefficients for structured predictors
coef(m, which_param = "interacting")
coef(m, which_param = "shifting")

# Predictions in-sample
predict(m, type = "pdf")

# Prediction over the whole range of the response
predict(m, type = "pdf", newdata = wine[, -2])

# In- and out-of-sample log-likelihood
logLik(m)
logLik(m, newdata = wine[1:10, ])

# Plot the smooth predictor
plot(m)

# Run a cross-validation
mcv <- cv(m, cv_folds = 3)

# Fit an ensemble of `m` with coefficients 
ens <- ensemble(m, n_ensemble = 3)
coef(ens)

# Transformation ensemble predictions
unlist(logLik(ens))
```

# How to cite `deeptrafo`

When using the software, please cite

```
@Article{kook2024deeptrafo,
  title = {Estimating Conditional Distributions with Neural Networks Using {R} Package {deeptrafo}},
  author = {Lucas Kook and Philipp F. M. Baumann and Oliver D\"urr and Beate Sick and David R\"ugamer},
  journal = {Journal of Statistical Software},
  year = {2024},
  volume = {111},
  number = {10},
  pages = {1--36},
  doi = {10.18637/jss.v111.i10},
}
```

For methodology please cite an appropriate selection of:

```
@InProceedings{baumann2020deep,
  doi = {10.1007/978-3-030-86523-8\_1},
  year = {2021},
  publisher = {Springer International Publishing},
  pages = {3--18},
  author = {Philipp F. M. Baumann and Torsten Hothorn and David R\"{u}gamer},
  title = {Deep Conditional Transformation Models},
  booktitle = {Machine Learning and Knowledge Discovery in Databases. Research Track}
}
```

2. Probabilistic time series forecasts with autoregressive transformation models

```
@article{rugamer2021timeseries,
  doi = {10.48550/arXiv.2110.08248},
  year = {2021},
  journal = {arXiv preprint},
  note = {To appear in \emph{Statistics \& Computing}},
  author = {David R\"ugamer and Philipp FM Baumann and Thomas Kneib and Torsten Hothorn},
  title = {Probabilistic Time Series Forecasts with Autoregressive Transformation Models}
}
```

3. Ordinal neural network transformation models

```
@article{kook2020ordinal,
  doi = {10.1016/j.patcog.2021.108263},
  year = {2022},
  publisher = {Elsevier {BV}},
  volume = {122},
  pages = {108263},
  author = {Lucas Kook and Lisa Herzog and Torsten Hothorn and Oliver D\"{u}rr and Beate Sick},
  title = {Deep and interpretable regression models for ordinal outcomes},
  journal = {Pattern Recognition}
}
```

4. Transformation ensembles

```
@article{kook2022deep,
  title={Deep interpretable ensembles},
  author={Kook, Lucas and G{\"o}tschi, Andrea and Baumann, Philipp FM and Hothorn, Torsten and Sick, Beate},
  journal={arXiv preprint arXiv:2205.12729},
  year={{2022}},
  doi={10.48550/arXiv.2205.12729}
}
```

