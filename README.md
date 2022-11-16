# R package deeptransformation

[![R build status](https://github.com/neural-structured-additive-learning/deeptransformation/workflows/R-CMD-check/badge.svg)](https://github.com/neural-structured-additive-learning/deeptransformation/actions) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This package provides function to fit

-   Deep Conditional Transformation Models (Baumann et al., 2021)
-   Autoregressive Transformation Models (Ruegamer et al., 2021)
-   Deep Interpretable Ensembles (Kook et al., 2022)

using `deepregression`.

# Installation

To install the package, use the following command:

``` r
remotes::install_github("neural-structured-additive-learning/deeptransformation")
```

Note that the installation requires additional packages (see below) including `deepregression`.

# How to cite this?

For deep conditional transformation models, please cite the following:

    @InProceedings{baumann2021dctm,
      author={Baumann, Philipp F. M. and Hothorn, Torsten, and R{\"u}gamer, David},
      title="Deep Conditional Transformation Models",
      booktitle="Machine Learning and Knowledge Discovery in Databases (ECML-PKDD)",
      year="2021",
      publisher="Springer International Publishing",
      pages="3--18"
    }

For autoregressive transformation models, please cite:

    @article{ruegamer2021atm,
      title={Transforming autoregression for expressive forecasts with parametric uncertainty quantification},
      author={R{\"u}gamer, David and Baumann, Philipp and Hothorn, Torsten},
      year={2021},
      eprint={2110.08248},
      archivePrefix={arXiv},
      journal={arXiv preprint arXiv:2110.08248}
    }

For transformation ensembles, please cite:

    @article{kook2022,
      title={Deep interpretable ensembles},
      author={Kook, Lucas and G{\"o}tschi, Andrea and Baumann, Philipp and Hothorn, Torsten and Sick, Beate},
      year={2022},
      eprint={2205.12729},
      archivePrefix={arXiv},
      journal={arXiv preprint arXiv:2205.12729}
    }
