# Robust Estimation and Inference for High-Dimensional IV Cointegration Models 

### Description 

The R package [‘Robust_high_dimensional_IV_Cointegration’](https://github.com/christiskatsouris/Robust_high_dimensional_IV_Cointegration) (under development by Christis G. Katsouris) implements robust econometric estimation and inference methodologies for high-Dimensional IV Cointegration Models with autoregressive roots across the spectrum of stationarity and nonstationarity, with persistence types as defined by Magdalinos and Phillips (2020). In particular, this R package builds on the [‘ivxPredictive’](https://github.com/christiskatsouris/ivxPredictive) R package prepared by Christis G. Katsouris. 

<p align="center">
  
<img src="https://github.com/christiskatsouris/ivxPredictive/blob/main/data/persistence.jpg" width="460"/>

</p>  

We consider the following persistence classes:
   
#### (P1). nearly stable processes: 

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \zeta = - \infty \ \ \text{and it holds that} \ \ \theta_n \to | \theta | < 1.$$
    
#### (P2). nearly unstable processes:   

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \ \zeta \equiv c \in \mathbb{R} \ \ \text{and it holds that} \ \ \theta_n \to \theta = 1.$$

    
#### (P3). nearly explosive processes:   

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \ \zeta = + \infty \ \ \text{and it holds that} \ \ \theta_n \to | \theta | > 1.$$

  
### Methodology  
  
This R package implements a novel endogenous instrumentation approach based on the IVX estimator examined by [Phillips and Magdalinos (2009)](https://ideas.repec.org/p/skb/wpaper/cofie-06-2009.html) and [Kostakis, Magdalinos and Stamatogiannis (2015)](https://academic.oup.com/rfs/article/28/5/1506/1867633?login=true). The current procedure has a similar construction to the IV instrumentation proposed in the recent working paper of Magdalinos and Petrova (2022), with the aim to provide uniform and robust to the nuisance parameter of persistence inference across the spectrum of stationary and nonstationary roots, specifically for quantile autoregressive processes. We call this variant of the original IVX estimator, IVX-P, which can be employed to both conditional mean and conditional quantile functional forms when the model includes either univariate or multivariate regressors. The novelty of the IVX-P estimator is that is a 'hybrid estimator' which in contrast to the classical least squares estimator has desirable asymptotic theory properties and is constructed based on the underline nonstationary stochastic processes using information both within the admissible parameter space as well as outside the usual parameter space.  

Furthermore, this R package implements robust estimation and testing for high-Dimensional IV Cointegration Models with either a conditional mean or conditional quantile specification function. 
  
## Installation (under development) 

The R package ‘Robust_high_dimensional_IV_Cointegration’ will be able to be installed from Github.

## Usage 

```R

# After development the package will be able to be installed using
install.packages("Robust_high_dimensional_IV_Cointegration")
library("Robust_high_dimensional_IV_Cointegration")

```

## Key References:

- Katsouris, C. (2022b). "Asymptotic Theory for Moderate Deviations from the Unit Boundary in Quantile Autoregression and Predictive Regression Models". University of Southampton, Working paper.  
- Katsouris, C. (2022a). "Asymptotic Theory for Moderate Deviations from the Unit Boundary in Quantile Autoregressive Time Series". arXiv preprint [arXiv:2204.02073](https://arxiv.org/abs/2204.02073).
- Katsouris, C. (2021d). "Testing for Structural Breaks in Predictive Regression Models". University of Southampton, Working paper.  
- Katsouris, C. (2021e). "Bootstrapping Nonstationary Autoregressive Processes in Predictive Regression". University of Southampton, Working paper.   
- Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2015). "Robust econometric inference for stock return predictability". The Review of Financial Studies, 28(5), 1506-1553.
- Lee, J. H. (2016). "Predictive quantile regression with persistent covariates: IVX-QR approach". Journal of Econometrics, 192(1), 105-118.
- Fan, R., & Lee, J. H. (2019). Predictive quantile regressions under persistence and conditional heteroskedasticity. Journal of Econometrics, 213(1), 261-280.
- Magdalinos, T. (2016). Least squares and IVX limit theory in systems of predictive regressions with GARCH innovations. Econometric Theory, 1-38.
- Magdalinos, T., & Petrova, K. (2022). "Uniform and distribution-free inference with general autoregressive processes". University of Southampton, Working paper.
- Magdalinos, T., & Phillips, P. C. B. (2009). "Limit theory for cointegrated systems with moderately integrated and moderately explosive regressors". Econometric Theory, 25(2), 482-526.
- Phillips, P. C. B., & Magdalinos, T. (2009). "Econometric inference in the vicinity of unity". Singapore Management University, CoFie Working paper, 7.
- Phillips, P. C. B., & Magdalinos, T. (2008). Limit theory for explosively cointegrated systems. Econometric Theory, 24(4), 865-887.
- Phillips, P. C. B., & Magdalinos, T. (2007). Limit theory for moderate deviations from a unit root. Journal of Econometrics, 136(1), 115-130.
