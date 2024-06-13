# Multivariate inverse Gaussian process with common effects
 
The online code repository associated with this paper is organized into three main directories, each tailored for specific purposes as detailed below:
 \begin{itemize}
 	\item \textbf{case}: This directory houses 
 	\begin{itemize}
 		\item The ``Fatigue-crack-size.xlsx" dataset, originally introduced in \cite{meeker2022statistical} and further processed as described in Appendix H of \cite{Fang2020}.
 		\item The ``crack.R" script, which performs parameter estimation across various models and calculates their respective AIC.
 		\item The ``results" folder, where the final analytical outputs are systematically stored.
  	\end{itemize}
 	\item \textbf{simulation}: ``Integral\_appr.R": the primary R script for replicating Figure 2 in the paper. This script evaluates the efficacy of various numerical integration techniques in approximating the CDF of failure times.
 	\item \textbf{utility}: Houses a suite of functions crucial for the computational analyses:
 	\begin{itemize}
 		\item ``appr.R": Functions dedicated to numerical integration approximations.
 		\item ``em.R": Functions related to the EM algorithm.
 		\item ``fct.R": A collection of auxiliary functions regularly employed throughout the analyses.
 	\end{itemize}
 \end{itemize}
For optimal interaction with these resources, it is recommended to open ``multi-rIG.Rproj" using RStudio\footnote{https://posit.co/download/rstudio-desktop/}, install all necessary packages as initially specified, and proceed to execute the code sequentially, section by section. 
