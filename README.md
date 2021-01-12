# **realPSD**: **R**obust and **E**fficient C**al**ibration of Parametric Power Spectral Density (PSD) Models 

*Feiyu Zhu, Martin Lysy*

---

### Description
**realPSD** provides **R**obust and **E**fficient c**AL**ibration of any parametric power spectral density (**PSD**) models, especially for dynamical processes with extremely large number (at least millions) of observations. Three built-in estimation methods are provided: maximum likelihood estimation (MLE) based on Whittle-type likelihood, log-periodogram estimation and nonlinear least squares estimation (NLS). The most effective of these combines the simplicity of nonlinear least-squares with the statistical efficiency of maximum likelihood. The technical details about these three estimation methods for high-throughput data (at extremely high sampling frequency and for extended durations) are extensively explained in this [preprint](). **realPSD** also provides a routine for removing most electronic noise in an automated pre-processing step, which is also discussed in this [preprint](). This feature would be very useful for spectral analysis based on recordings contaminated by various sources of instrumental noise in real-world applications.

**realPSD** relies upon [**TMB**](https://github.com/kaskr/adcomp.git) internally for [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) (AD) to enable efficient and accurate numerical gradient calculations. But users are not required to be familiar with **TMB**. **realPSD** itself can be used as a user-friendly platform that handles the **TMB** interface to R for users. All you need to do is to focus on writing your model and estimating it. That being said, there are some basic conventions required by **realPSD** to write your own model, see the [quick tutorial](http://htmlpreview.github.io/?https://github.com/mlysy/realPSD/blob/devel-ferris-prerelease/doc/realPSD-quicktut.html). It would be easier for users to quickly grasp these steps if users are familiar with any object-oriented programming (OOP) language like C++ or Java. Understanding of the basic data types provided by the C++ [**Eigen**](http://eigen.tuxfamily.org/index.php?title=Main_Page) library is preferred but not required.

### Installation

**realPSD** relies on **TMB**. Please make sure that **TMB** is correctly installed on your machine. The installation of **TMB** can be found [here](https://github.com/kaskr/adcomp/wiki/Download).

Once you have installed **TMB** successfully, **realPSD** can be installed from GitHub via the command:

```r
devtools::install_github("mlysy/realPSD")
```

### Test the installation

There are many built-in unit tests in **realPSD**. Once you have installed the package, run the following command to confirm that it can work properly on your machine.

```r
require(realPSD)
require(testthat)
test_pacakge("realPSD")
```

### Tutorials

- [`realPSD-quicktut`](http://htmlpreview.github.io/?https://github.com/mlysy/realPSD/blob/master/doc/realPSD-quicktut.html): A quick tutorial teaches you how to use **realPSD** by showing you how to build and estimate an Ornstein-Uhlenbeck (OU) model step by step.
- [`realPSD-advancedtut`](http://htmlpreview.github.io/?https://github.com/mlysy/realPSD/blob/master/doc/realPSD-advancedtut.html): An Advanced tutorial explains the theoretical results behind estimation methods supplied by the package and also the design pattern of **realPSD** that handles the **TMB** interface to R effectively, which allows users to build any model without restrictions.

