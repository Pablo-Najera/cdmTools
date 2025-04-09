# cdmTools 1.0.6
* Added - `FCGDINA` function to estimate the G-DINA model with forced-choice items
* Added - `simFCGDINA` function to simulate forced-choice data from the G-DINA model
* Added - `is.CDMid` function to check model identifiability via post-hoc simulations
* Changed - `RDINA` function to enable empirical Bayes estimation of structural parameters

# cdmTools 1.0.5
* Added - `RDINA2GDINA` function to increase compatibility with the `GDINA` package
* Changed - `CA.MI` and `personFit` have been updated to be compatible with `RDINA` objects
* Kevin Santos included as a contributor

# cdmTools 1.0.4
* Changed - `is.Qid` has been updated to include new identifiability criteria
* Added - `personfit` function to identify aberrant response patterns

# cdmTools 1.0.3
* Fixed - `CA.MI` had a small issue that has been corrected
* Added - `RDINA` function to estimate the R-DINA model for small sample size settings

# cdmTools 1.0.2
* Added - `CA.MI` function to calculate corrected classification accuracy with multiple imputation
* Rodrigo S. Kreitchmann included as a contributor

# cdmTools 1.0.1
* Added - `GNPC` function to estimate attribute profiles
* Added - `print` methods for classes

# cdmTools 1.0.0
* CRAN release

# cdmTools 0.1.1
* Fixed - bug in `modelcompK` when a stop criterion is used and exploreK does not include K = 1
* Fixed - bug in `paK` when missing values are present in correlation matrix
* Fixed - bug in `modelcompK` when M2 statistic cannot be computed for a model

# cdmTools 0.1.0
* GitHub release
