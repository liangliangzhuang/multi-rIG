# Multivariate inverse Gaussian process with common effects

The online code repository associated with this paper is organized into three main directories, each tailored for specific purposes as detailed below:

-   **case**: This directory houses

    -   The "Fatigue-crack-size.xlsx" dataset, originally introduced in [Meeker et al. (2022)](https://www.wiley.com/en-us/Statistical+Methods+for+Reliability+Data%2C+2nd+Edition-p-9781118115459) and further processed as described in Appendix H of [Fang et al. (2022)](https://www.sciencedirect.com/science/article/abs/pii/S0377221721008985).
    -   The "crack.R" script, which performs parameter estimation across various models and calculates their respective AIC.
    -   The "results" folder, where the final analytical outputs are systematically stored.

-   **simulation**: "Integral_appr.R": the primary R script for replicating Figure 2 in the paper. This script evaluates the efficacy of various numerical integration techniques in approximating the CDF of failure times.

-   **utility**: Houses a suite of functions crucial for the computational analyses:

    -   "appr.R": Functions dedicated to numerical integration approximations.

    -   "em.R": Functions related to the EM algorithm.

    -   "fct.R": A collection of auxiliary functions regularly employed throughout the analyses.

For optimal interaction with these resources, it is recommended to open "multi-rIG.Rproj" using [RStudio](https://posit.co/download/rstudio-desktop/), install all necessary packages as initially specified, and proceed to execute the code sequentially, section by section.
