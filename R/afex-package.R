#' Analysis of Factorial Experiments.
#'
#' \tabular{ll}{
#' Package: \tab afex\cr
#' Type: \tab Package\cr
#' Version: \tab 0.2-26\cr
#' Date: \tab 2012-08-09\cr
#' Depends: \tab R (>= 2.14.0), car, lme4, pbkrtest, reshape2, stringr\cr
#' Encoding: \tab UTF-8\cr
#' License: \tab GPL (>=3)\cr
#' URL: \tab http://www.psychologie.uni-freiburg.de/Members/singmann/R/afex\cr
#' }
#'
#' Provides convenience functions for analyzing factorial experiments using ANOVA or mixed-models. Functions ez.glm() and aov.car() allow convenient specification of pure-between, pure-within (i.e., repeated-measures), and mixed between-within (i.e., split-plot) ANOVAs and ANCOVAs with type 2 and type 3 sums of squares for data in the long format (i.e., one observation per row) wrapping car::Anova() (aggregating more then one observation per individual and cell of the design). Function nice.anova() produces publication ready ANOVA tables. Function mixed() fits a mixed model using lme4::lmer() and computes p-values for all effects in the model (using the Kenward-Rogers approximation of degrees of freedom). afex uses type 3 sums of squares as default (imitating commercial statistical software) and sets the default contrasts to contr.sum.
#'
#' @aliases afex-package afex
#' @name afex-package
#' @docType package
#' @title Analysis of Factorial Experiments.
#' @author Henrik Singmann \email{henrik.singmann@@psychologie.uni-freiburg.de}
#' @keywords package
NULL
