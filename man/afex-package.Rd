\docType{package}
\name{afex-package}
\alias{afex}
\alias{afex-package}
\title{The afex Package}
\description{
  Analysis of Factorial Experiments.
}
\details{
  \tabular{ll}{ Package: \tab afex\cr Type: \tab Package\cr
  Version: \tab 0.7-90\cr Date: \tab 2013-12-04\cr Depends:
  \tab R (>= 2.14.0), coin, car, lme4 (>= 1.0.4), pbkrtest
  (>= 0.3-6)\cr Encoding: \tab UTF-8\cr License: \tab GPL
  (>=3)\cr URL: \tab
  http://www.psychologie.uni-freiburg.de/Members/singmann/R/afex\cr
  }

  Provides convenience functions for analyzing factorial
  experiments using ANOVA or mixed-models. ez.glm() and
  aov.car() allow convenient calculation of between, within
  (i.e., repeated-measures), or mixed between-within (i.e.,
  split-plot) ANOVAs for data in the long format (i.e., one
  observation per row) wrapping car::Anova() (aggregating
  more then one observation per individual and cell of the
  design), per default returning a print ready ANOVA table.
  Function mixed() fits a mixed model using lme4::lmer()
  and computes p-values for all effects in the model using
  either the Kenward-Rogers approximation of degrees of
  freedom (LMM only), parametric bootstrap (LMMs and GLMMs)
  or likelihood ratio tests (LMMs and GLMMs). afex uses
  type 3 sums of squares as default (imitating commercial
  statistical software) and sets the default contrasts to
  contr.sum. Furthermore, compare.2.vectors() conveniently
  compares two vectors using a variety of tests.
}
\author{
  Henrik Singmann
  \email{henrik.singmann@psychologie.uni-freiburg.de}
}
\keyword{package}

