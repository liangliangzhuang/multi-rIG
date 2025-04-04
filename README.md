# Multivariate Reparameterized Inverse Gaussian Processes with Common Effects for Degradation-based Reliability Prediction



The repository related to this paper is structured into three primary directories, each designed for specific functions as outlined below:

-   **case**: This directory includes:

    -   The "Fatigue-crack-size.xlsx" dataset, originally introduced in [Meeker et al. (2022)](https://www.wiley.com/en-us/Statistical+Methods+for+Reliability+Data%2C+2nd+Edition-p-9781118115459) and further processed as described in Appendix H of [Fang et al. (2022)](https://www.sciencedirect.com/science/article/abs/pii/S0377221721008985).
    -   The "crack.R" script, which performs parameter estimation across various models and calculates their respective AIC.
    -   The "results" folder, where the final analytical outputs are stored.

-   **simulation**: "Integral_appr.R": the primary R script for replicating Figure 2. This script evaluates the efficacy of various numerical integration techniques in approximating the CDF of failure time.

-   **utility**: This directory contains a collection of essential functions for computational analysis:

    -   "appr.R": Functions dedicated to numerical integration approximations.

    -   "em.R": Functions related to the EM algorithm.

    -   "fct.R": A collection of auxiliary functions regularly employed throughout the analyses.

> For optimal interaction with these codes, it is recommended to open "multi-rIG.Rproj" using [RStudio](https://posit.co/download/rstudio-desktop/), install all necessary packages as initially specified, and proceed to execute the code sequentially, section by section.

## Note 

This work has been accepted for publication by Journal of Quality Technology.  If you use the provided code in this project, please remember to cite the paper accordingly. Detailed citation information is

```bibtex
@article{zhuang2025multivariate,
  title={Multivariate reparameterized inverse Gaussian processes with common effects for degradation-based reliability prediction},
  author={Zhuang, Liangliang and Xu, Ancha and Fang, Guanqi and Tang, Yincai},
  journal={Journal of Quality Technology},
  volume={57},
  number={1},
  pages={51--67},
  year={2025},
  publisher={Taylor \& Francis}
}
```

If you have any questions or need help with the code, please submit them in the [issue](https://github.com/liangliangzhuang/multi-rIG/issues).


