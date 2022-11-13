# Robust Estimation and Inference for High-Dimensional IV Cointegration Models 

### Description 

The R package [‘Robust_high_dimensional_IV_Cointegration’](https://github.com/christiskatsouris/Robust_high_dimensional_IV_Cointegration) (under development by Christis G. Katsouris) implements robust econometric estimation and inference methodologies for high-Dimensional IV Cointegration Models with autoregressive roots under both the regimes of stationarity and nonstationarity and persistence types as defined by Phillips and Magdalinos (2009): "Econometric inference in the vicinity of unity" (see, also Magdalinos and Phillips (2020)). The current package builds on the [‘ivxPredictive’](https://github.com/christiskatsouris/ivxPredictive) package prepared by Christis G. Katsouris. 

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
  
The [‘ivxPredictive’](https://github.com/christiskatsouris/ivxPredictive) R package implements a novel endogenous instrumentation approach based on the IVX estimator proposed by [Phillips and Magdalinos (2009)](https://ideas.repec.org/p/skb/wpaper/cofie-06-2009.html) and [Kostakis, Magdalinos and Stamatogiannis (2015)](https://academic.oup.com/rfs/article/28/5/1506/1867633?login=true). Notice that the current procedure has a similar construction to the IV instrumentation proposed in the recent working paper of Magdalinos and Petrova (2022), with the aim to provide uniform and robust to the nuisance parameter of persistence inference across the spectrum of stationary and nonstationary roots, specifically for quantile autoregressive processes. We call this variant of the original IVX estimator, IVX-P, which can be employed to both conditional mean and conditional quantile functional forms when the model includes either univariate or multivariate regressors. The novelty of the IVX-P estimator is that is a 'hybrid estimator' which in contrast to the classical least squares estimator has desirable asymptotic theory properties and is constructed based on the underline nonstationary stochastic processes using information both within the admissible parameter space as well as outside the usual parameter space of autoregressive prcosesses.  

Furthermore, this R package implements robust estimation and testing for high-Dimensional IV Cointegration Models with either a conditional mean or conditional quantile specification function. Specifically, the user can choose the type of the functional form of the cointegration model. 
  
## Installation (under development) 

The R package [‘Robust_high_dimensional_IV_Cointegration’](https://github.com/christiskatsouris/Robust_high_dimensional_IV_Cointegration) will be able to be installed from Github.

## Usage: 

```R

# After development the package will be able to be installed using
install.packages("Robust_high_dimensional_IV_Cointegration")
library("Robust_high_dimensional_IV_Cointegration")

```

## Estimation and Inference Examples:

```R



```

## Monte Carlo Simulation Studies:

Here we present the R coding and some preliminary results based on a Monte Carlo Simulation study of the finite-sample properties of the inferential procedure that is under development for this research project.  

Consider the following Toeplitz structure of the covariance matrix. 

$\textbf{\Sigma} = \rho^{| j - k |}.$ 

```R



```

## Main References:

- Katsouris, C. (2022c). "Estimation and Inference in Quantile Predictive Regression Systems" (Chapter 4, PhD thesis, School of Economic, Social and Political Sciences, University of Southampton.
- Katsouris, C. (2022b). "Asymptotic Theory for Moderate Deviations from the Unit Boundary in Quantile Autoregression and Predictive Regression Models". University of Southampton, Working paper.  
- Katsouris, C. (2022a). "Asymptotic Theory for Moderate Deviations from the Unit Boundary in Quantile Autoregressive Time Series". arXiv preprint [arXiv:2204.02073](https://arxiv.org/abs/2204.02073).
- Katsouris, C. (2021e). "Bootstrapping Nonstationary Autoregressive Processes in Predictive Regression". University of Southampton, Working paper.   
- Katsouris, C. (2021d). "Testing for Structural Breaks in Predictive Regression Models". University of Southampton, Working paper.  
- Ke-Li Xu & Junjie Guo (2022). "A New Test for Multiple Predictive Regression". Journal of Financial Econometrics.
- Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2015). "Robust econometric inference for stock return predictability". The Review of Financial Studies, 28(5), 1506-1553.
- Koo, B., Anderson, H. M., Seo, M. H., & Yao, W. (2020). "High-dimensional predictive regression in the presence of cointegration". Journal of Econometrics, 219(2), 456-477.
- Lee, J. H. (2016). "Predictive quantile regression with persistent covariates: IVX-QR approach". Journal of Econometrics, 192(1), 105-118.
- Lee, J. H., Shi, Z., & Gao, Z. (2022). "On LASSO for predictive regression". Journal of Econometrics, 229(2), 322-349.
- Fan, R., & Lee, J. H. (2019). "Predictive quantile regressions under persistence and conditional heteroskedasticity". Journal of Econometrics, 213(1), 261-280.
- Chen, W. W., Deo, R. S., & Yi, Y. (2013). "Uniform inference in predictive regression models". Journal of Business & Economic Statistics, 31(4), 525-533.
- Magdalinos, T. (2016). "Least squares and IVX limit theory in systems of predictive regressions with GARCH innovations". Econometric Theory, 1-38.
- Magdalinos, T., and Petrova, K. (2022). "Uniform and distribution-free inference with general autoregressive processes". University of Southampton, Working paper.
- Magdalinos, T., and Phillips, P. C. B. (2009). "Limit theory for cointegrated systems with moderately integrated and moderately explosive regressors". Econometric Theory, 25(2), 482-526.
- Phillips, P. C. B., and Magdalinos, T. (2009). "Econometric inference in the vicinity of unity". Singapore Management University, CoFie Working paper, 7.
- Phillips, P. C. B., and Magdalinos, T. (2008). "Limit theory for explosively cointegrated systems". Econometric Theory, 24(4), 865-887.
- Phillips, P. C. B., and Magdalinos, T. (2007). "Limit theory for moderate deviations from a unit root". Journal of Econometrics, 136(1), 115-130.
- Wagner, M., Grabarczyk, P., & Hong, S. H. (2020). "Fully modified OLS estimation and inference for seemingly unrelated cointegrating polynomial regressions and the environmental Kuznets curve for carbon dioxide emissions". Journal of Econometrics, 214(1), 216-255.
- Yousuf, K., & Ng, S. (2021). "Boosting high dimensional predictive regressions with time varying parameters". Journal of Econometrics, 224(1), 60-87.
- Zhu, F., Cai, Z., & Peng, L. (2014). "Predictive regressions for macroeconomic data". The Annals of Applied Statistics, 8(1), 577-594.

## Literature on high-dimensional IV regression:

- Belloni, A., Hansen, C., & Newey, W. (2017). Simultaneous confidence intervals for high-dimensional linear models with many endogenous variables. arXiv preprint arXiv:1712.08102.
- Gold, D., Lederer, J., & Tao, J. (2020). "Inference for high-dimensional instrumental variables regression". Journal of Econometrics, 217(1), 79-111.
- Ning, Y., & Liu, H. (2017). A general theory of hypothesis tests and confidence regions for sparse high dimensional models. The Annals of Statistics, 45(1), 158-195.
- Van de Geer, S., Bühlmann, P., Ritov, Y. A., & Dezeure, R. (2014). On asymptotically optimal confidence regions and tests for high-dimensional models. The Annals of Statistics, 42(3), 1166-1202.

# Acknowledgments

The author greatfully acknowledges financial support from the [Department of Economics](http://business-school.exeter.ac.uk/about/departments/economics/) of the [Faculty of Environment, Science and Economy](https://www.exeter.ac.uk/departments/ese/) at the University of Exeter, United Kingdom. 

Christis G. Katsouris is a Lecturer in Economics at the [University of Exeter Business School](http://business-school.exeter.ac.uk/). He is also a member of the [Time Series and Machine Learning Group](https://www.personal.soton.ac.uk/cz1y20/Reading_Group/mlts-group-2022.html) at the [School of Mathematical Sciences](https://www.southampton.ac.uk/about/faculties-schools-departments/school-of-mathematical-sciences) (Statistics Division) of the University of Southampton. 

# Declarations

The author (Christis G. Katsouris) declares no conflicts of interest.

Notice that the academic research presented here is considered to be as open access to the academic and non-academic community. Therefore, we would appreciate it if appropriate acknolwedgement is given to statistical methodologies and econometric procedures developed by academic researchers and made available to the wider applied data scientist community.   

# Historical Background

#### Harald Cramér 

Harald Cramér  (25 September 1893 – 5 October 1985) was a Swedish mathematician, actuary, and statistician, specializing in mathematical statistics and probabilistic number theory. John Kingman described him as "one of the giants of statistical theory". A large portion of Cramér's work concerned the field of actuarial science and insurance mathematics. In 1929, Cramér was appointed to a newly created chair in Stockholm University, becoming the first Swedish professor of Actuarial Mathematics and Mathematical Statistics. Cramér retained this position up until 1958. During his tenure at Stockholm University, Cramér was a PhD advisor for 10 students, most notably Herman Wold and Kai Lai Chung. In 1950 he was elected as a Fellow of the American Statistical Association. Starting in 1950, Cramér took on the additional responsibility of becoming the President of Stockholm University. In 1958, he was also appointed to be Chancellor of the entire Swedish university system (Source: Wikepedia). 

