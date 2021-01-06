# **realPSD**: **R**obust and **E**fficient C**al**ibration of Parametric PSD Models 

*Martin Lysy, Feiyu Zhu*

---

### Description
**realPSD** allows users to easily build and estimate any parametric power spectral density (PSD) model in a robust and efficient way, especially for dynamical processes with extremely large (at least millions) observations. Three built-in estimation methods are provided, maximum likelihood estimation (MLE) based on Whittle-type likelihood, log-periodogram estimation and nonlinear least squares estimation (NLS). The most effective of these combines the simplicity of nonlinear least-squares with the statistical efficiency of maximum likelihood. The technical details about these three estimation methods for parametric PSD calibration based on high-throughput data (at high sampling frequency and for extended durations) are extensively explained in this [preprint](). **realPSD** also provides a routine for removing most electronic noise in an automated pre-processing step, which is also discussed in this [preprint]().

**realPSD** relies upon [**TMB**](https://github.com/kaskr/adcomp.git) internally for efficient and accurate numerical automatic differentiation. But users are not required to be familiar with **TMB**. **realPSD** itself can be used as a user-friendly platform that handles the **TMB** interface to R for users. All you need to do is to focus on writing your model and estimating it. That being said, there are some basic convetions required by **realPSD** to write your own model, see the [quick tutorial](http://htmlpreview.github.com/?https://github.com/mlysy/realPSD/blob/devel-ferris-prerelease/doc/realPSD-quicktut.Rmd). It would be easier for users to quickly grasp these steps if users are familiar with any object-oriented programming (OOP) language like C++ or Java. Understanding of the basic data types provided by the C++ [**Eigen**](http://eigen.tuxfamily.org/index.php?title=Main_Page) library is preferred but not required.

### Installation

**TMB** uses the C++ library **Eigen**. Thus, **realPSD** also relies on it. Please make sure that **Eigen** is installed on your machine. The installation of **Eigen** can be found [here](https://eigen.tuxfamily.org/dox/GettingStarted.html).

Once you have installed **Eigen** correctly, **realPSD** can be installed from GitHub via the command:

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

- [`realPSD-quicktut`](http://htmlpreview.github.com/?https://github.com/mlysy/realPSD/blob/devel-ferris-prerelease/doc/realPSD-quicktut.Rmd): A quick tutorial teaches you how to use **realPSD** by showing you how to build and estimate an Ornstein-Uhlenbeck (OU) model step by step.
- [`realPSD-advacnedtut`](http://htmlpreview.github.com/?https://github.com/mlysy/realPSD/blob/devel-ferris-prerelease/doc/realPSD-advancedtut.Rmd): An Advanced tutorial explains the theoretical results behind estimation methods supplied by the package and also the design pattern of **realPSD** that handles the **TMB** interface to R effectively, which allows users to build any model without restrictions.

