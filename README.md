# Multivariate inverse Gaussian process with common effects
 
The online code repository associated with this paper is organized into three main directories, each tailored for specific purposes as detailed below:

- \textbf{case}: This directory houses 

 		- The ``Fatigue-crack-size.xlsx" dataset, originally introduced in \cite{meeker2022statistical} and further processed as described in Appendix H of \cite{Fang2020}.
   - The ``crack.R" script, which performs parameter estimation across various models and calculates their respective AIC.
   - The ``results" folder, where the final analytical outputs are systematically stored.
 
- \textbf{simulation}: ``Integral\_appr.R": the primary R script for replicating Figure 2 in the paper. This script evaluates the efficacy of various numerical integration techniques in approximating the CDF of failure times.

  
- \textbf{utility}: Houses a suite of functions crucial for the computational analyses:
  
 		- ``appr.R": Functions dedicated to numerical integration approximations.
  
 		- ``em.R": Functions related to the EM algorithm.
  
 		- ``fct.R": A collection of auxiliary functions regularly employed throughout the analyses.

For optimal interaction with these resources, it is recommended to open ``multi-rIG.Rproj" using RStudio\footnote{https://posit.co/download/rstudio-desktop/}, install all necessary packages as initially specified, and proceed to execute the code sequentially, section by section. 
