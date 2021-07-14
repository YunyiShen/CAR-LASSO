# Contributing to CAR-LASSO R package

The following guidelines are designed for contributors to the CAR-LASSO R package. 

## Reporting Issues

For reporting a bug or a failed function or requesting a new feature, you can simply open an issue in the [issue tracker](https://github.com/YunyiShen/CAR-LASSO/issues). If you are reporting a bug, please also include a minimal code example or all relevant information for us to replicate the issue.

## Contributing Code

To make contributions to CAR-LASSO R package, you need to set up your [GitHub](https://github.com) 
account if you do not have and sign in, and request your change(s) or contribution(s) via 
a pull request against the ``dev``
branch of the [CAR-LASSO repository](https://github.com/YunyiShen/CAR-LASSO). 

Please use the following steps:

1. Open a new issue for new feature or failed function in the [Issue tracker](https://github.com/solislemuslab/bioklustering/issues)
2. Fork the [CAR-LASSO repository](https://github.com/YunyiShen/CAR-LASSO) to your GitHub account
3. Clone your fork locally:
```
$ git clone https://github.com/your-username/CAR-LASSO.git
```   
4. Make your change(s) in the `master` (or `development`) branch of your cloned fork
5. Make sure that all tests pass without any errors
6. Push your change(s) to your fork in your GitHub account
7. [Submit a pull request](https://github.com/YunyiShen/CAR-LASSO/pulls) describing what problem has been solved and linking to the issue you had opened in step 1

Your contribution will be checked and merged into the original repository. You will be contacted if there is any problem in your contribution.

Make sure to include the following information in your pull request:

* **Code** which you are contributing to this package

* **Documentation** of this code if it provides new functionality. This should be a
  description of new functionality added to the `vignettes` folder and documentations directly in the `R` folder written in `Roxygen`.

- **Tests** of this code to make sure that the previously failed function or the new functionality now works properly

---

_These Contributing Guidelines have been adapted from the [Contributing Guidelines](https://github.com/atomneb/AtomNeb-py/blob/master/CONTRIBUTING.md) of [The Turing Way](https://github.com/atomneb/AtomNeb-py)! (License: MIT)_