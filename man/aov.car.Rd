\name{aov.car}
\alias{aov.car}
\alias{ez.glm}
\alias{univariate}
\title{Convenience wrappers for car::Anova using either a formula or factor based interface.}
\usage{
  aov.car(formula, data, fun.aggregate = NULL, type = 3,
    ...)

  ez.glm(id, dv, data, between = NULL, within = NULL,
    covariate = NULL, fun.aggregate = NULL, type = 3, ...,
    print.formula = FALSE)

  univariate(object)
}
\arguments{
  \item{formula}{A formula specifying the ANOVA model
  similar to \code{\link{aov}}. Should include an error
  term (i.e., \code{Error( / )}). Note that the
  within-subject factors do not need to be outside the
  Error term (this contrasts with \code{aov}). See
  Details.}

  \item{id}{\code{character} vector (of length 1)
  indicating the subject identifier column in \code{data}.}

  \item{dv}{\code{character} vector (of length 1)
  indicating the column containing the \strong{dependent
  variable} in \code{data}.}

  \item{between}{\code{character} vector indicating the
  \strong{between}-subject(s) factor(s)/column(s) in
  \code{data}. Default is \code{NULL} indicating no
  between-subjects factors.}

  \item{within}{\code{character} vector indicating the
  \strong{within}-subject(s) factor(s)/column(s) in
  \code{data}.  Default is \code{NULL} indicating no
  within-subjects factors.}

  \item{covariate}{\code{character} vector indicating the
  between-subject(s) covariate(s) (i.e., column(s)) in
  \code{data}. Default is \code{NULL} indicating no
  covariates.}

  \item{data}{A \code{data.frame} containing the data.
  Mandatory.}

  \item{fun.aggregate}{The function for aggregating the
  data before running the ANOVA if there is more than one
  obervation per individuum and cell of the design. The
  default \code{NULL} issues a warning if aggregation is
  necessary and uses \code{\link{mean}}.}

  \item{type}{The type of sums of squares for the ANOVA.
  \strong{Defaults to 3}. Passed to
  \code{\link[car]{Anova}}. Possible values are
  \code{"II"}, \code{"III"}, \code{2}, or \code{3}.}

  \item{print.formula}{\code{ez.glm} is a wrapper for
  \code{aov.car}. This boolean argument indicates whether
  the formula in the call to \code{car.aov} should be
  printed.}

  \item{...}{Further arguments passed to
  \code{fun.aggregate}.}

  \item{object}{An object of class \code{Anova.mlm} as
  returned by \code{aov.car}, \code{ez.glm}, or
  \code{\link[car]{Anova}}.}
}
\value{
  \code{aov.car} and \code{ez.glm} are wrappers and
  therfore return the same as \code{\link[car]{Anova}}.
  Usually an object of class \code{"Anova.mlm"} (with
  within-subjects factors) or of class \code{c("anova",
  "data.frame")}.

  \code{univariate} returns a \code{list} of
  \code{data.frame}s containing the univariate results
  (i.e., the classical ANOVA results) from an object of
  class \code{"Anova.mlm"}. This is essentially the output
  from \code{summary.Anova.mlm} with \code{multivariate =
  FALSE}, e.g. \code{summary(aov.car(...), multivriate =
  FALSE)}, as a list instead of printed to the console.\cr
  For objects of class \code{"anova"} (i.e., the object
  returned by \code{car::Anova} for a purely
  between-subjects ANOVA) the object is returned unaltered.

  The elements of the list returned by \code{univariate}
  are: \code{anova}, \code{mauchly}, and
  \code{spehricity.correction} (containing both,
  Greenhouse-Geisser and Hyundt-Feldt correction).
}
\description{
  These functions allow convenient access to
  \code{\link[car]{Anova}} (from the \pkg{car} package) for
  data in the \strong{long} format (i.e., one observation
  per row), possibly aggregating the data if there is more
  than one obersvation per individuum and cell. Hence,
  mixed between-within ANOVAs can be calculated
  conveniently without using the rather unhandy format of
  \code{car::Anova}. \code{aov.car} can be called using a
  formula similar to \code{\link{aov}} specifying an error
  strata for the within-subject factor(s). \code{ez.glm} is
  called specifying the factors as character vectors.
}
\details{
  \strong{Type 3 sums of squares are default in
  \pkg{afex}.} Note that type 3 sums of squares are said to
  be dangerous and/or problematic. On the other side they
  are the default in in SPSS and SAS and recommended by
  e.g. Maxwell and Delaney (2004). For a brief discussion
  see
  \href{http://stats.stackexchange.com/q/6208/442}{here}.

  However, note that lower order effects (e.g., main
  effects) in type 3 ANOVAs are only meaningful with
  \href{http://www.ats.ucla.edu/stat/mult_pkg/faq/general/effect.htm}{effects
  coding}. That is, contrasts should be set to
  \code{\link{contr.sum}} via
  \code{options(contrasts=c('contr.sum','contr.poly'))}.
  This should be done automatically when loading \pkg{afex}
  and \pkg{afex} will issue a warning when running type 3
  SS and
  \href{http://www.ats.ucla.edu/stat/r/library/contrast_coding.htm}{other
  coding schemes}. You can check the coding with
  \code{options("contrasts")}.

  The \code{formula} for \code{aov.car} must contain a
  single \code{Error} term specyfying the \code{ID} column
  and potential within-subject factors (you may use
  \code{\link{mixed}} with multiple error terms). Factors
  outside the \code{Error} term are treated as
  between-subject factors (the within-subject factors
  specified in the \code{Error} term are ignored outside
  the \code{Error} term, i.e., it is not necessary to
  specify them outside the \code{Error} term, see
  Examples).

  Suppressing the intercept (i.e, via \code{0 +} or \code{-
  1}) is ignored. Specific specifications of effects (e.g.,
  excluding terms with \code{-} or using \code{^}) could be
  okay but is not tested. Using the \code{\link{I}} or
  \code{\link{poly}} function within the formula is not
  tested and not supported!

  For \code{ez.glm} either \code{between} or \code{within}
  must not be \code{NULL}.

  \code{ez.glm} will concatante all between-subject factors
  using \code{*} (i.e., producing all main effects and
  interactions) and all covariates by \code{+} (i.e.,
  adding only the main effects to the existing
  between-subject factors). The within-subject factors do
  fully interact with all between-subject factors and
  covariates. This is essentially identical to the behavior
  of SPSS's \code{glm} function.
}
\note{
  Variables entered as within-subjects (i.e., repeated
  measures) factors are silently converted to factors and
  unused levels dropped.

  Contrasts attached to a factor as an attribute are
  probably not preserved and not supported.
}
\examples{
# exampel using obk.long (see ?obk.long), a long version of the OBrienKaiser dataset from car.
data(obk.long)

# run univariate mixed ANCOVA for the full design:
univariate(aov.car(value ~ treatment * gender + age + Error(id/phase*hour), data = obk.long))
univariate(ez.glm("id", "value", obk.long, c("treatment", "gender"), c("phase", "hour"), "age"))

# both calls return the same:

## $anova
##                                      SS num Df  Error SS den Df           F       Pr(>F)
## (Intercept)                 6454.236987      1 215.65658      9 269.3547893 5.152317e-08
## treatment                    171.399953      2 215.65658      9   3.5765187 7.193619e-02
## gender                        94.598340      1 215.65658      9   3.9478742 7.818280e-02
## age                           12.398975      1 215.65658      9   0.5174466 4.901885e-01
## treatment:gender              61.531858      2 215.65658      9   1.2839551 3.231798e-01
## phase                        134.586005      2  59.72439     18  20.2810632 2.448505e-05
## treatment:phase               80.604542      4  59.72439     18   6.0732385 2.826803e-03
## gender:phase                   1.634246      2  59.72439     18   0.2462681 7.843036e-01
## age:phase                     20.553392      2  59.72439     18   3.0972362 6.982439e-02
## treatment:gender:phase        21.254421      4  59.72439     18   1.6014379 2.170946e-01
## hour                         108.513510      4  47.59543     36  20.5192290 7.001584e-09
## treatment:hour                 7.547869      8  47.59543     36   0.7136275 6.779072e-01
## gender:hour                    3.746135      4  47.59543     36   0.7083708 5.915285e-01
## age:hour                      14.904567      4  47.59543     36   2.8183608 3.926421e-02
## treatment:gender:hour          6.235198      8  47.59543     36   0.5895186 7.798264e-01
## phase:hour                     9.762579      8  88.62706     72   0.9913814 4.501348e-01
## treatment:phase:hour           6.579092     16  88.62706     72   0.3340505 9.915014e-01
## gender:phase:hour              8.851396      8  88.62706     72   0.8988515 5.222336e-01
## age:phase:hour                 7.539611      8  88.62706     72   0.7656409 6.339004e-01
## treatment:gender:phase:hour   12.822199     16  88.62706     72   0.6510416 8.307936e-01
##
## $mauchly
##                             Test statistic    p-value
## phase                         0.8217571566 0.45600959
## treatment:phase               0.8217571566 0.45600959
## gender:phase                  0.8217571566 0.45600959
## age:phase                     0.8217571566 0.45600959
## treatment:gender:phase        0.8217571566 0.45600959
## hour                          0.0966749877 0.04923980
## treatment:hour                0.0966749877 0.04923980
## gender:hour                   0.0966749877 0.04923980
## age:hour                      0.0966749877 0.04923980
## treatment:gender:hour         0.0966749877 0.04923980
## phase:hour                    0.0002379741 0.08651564
## treatment:phase:hour          0.0002379741 0.08651564
## gender:phase:hour             0.0002379741 0.08651564
## age:phase:hour                0.0002379741 0.08651564
## treatment:gender:phase:hour   0.0002379741 0.08651564
##
## $sphericity.correction
##                                GG eps   Pr(>F[GG])    HF eps   Pr(>F[HF])
## phase                       0.8487215 8.383485e-05 1.0252867 2.448505e-05
## treatment:phase             0.8487215 5.159591e-03 1.0252867 2.826803e-03
## gender:phase                0.8487215 7.493990e-01 1.0252867 7.843036e-01
## age:phase                   0.8487215 8.073373e-02 1.0252867 6.982439e-02
## treatment:gender:phase      0.8487215 2.279698e-01 1.0252867 2.170946e-01
## hour                        0.5341747 1.302016e-05 0.7054545 8.046331e-07
## treatment:hour              0.5341747 6.010781e-01 0.7054545 6.342676e-01
## gender:hour                 0.5341747 5.137213e-01 0.7054545 5.478398e-01
## age:hour                    0.5341747 8.155027e-02 0.7054545 6.211130e-02
## treatment:gender:hour       0.5341747 6.843526e-01 0.7054545 7.263729e-01
## phase:hour                  0.4355822 4.186799e-01 0.7444364 4.402119e-01
## treatment:phase:hour        0.4355822 9.317848e-01 0.7444364 9.787985e-01
## gender:phase:hour           0.4355822 4.651930e-01 0.7444364 5.020890e-01
## age:phase:hour              0.4355822 5.395151e-01 0.7444364 5.992844e-01
## treatment:gender:phase:hour 0.4355822 7.100921e-01 0.7444364 7.878433e-01
##
## Warning message:
## In univariate(aov.car(value ~ treatment * gender + age + Error(id/phase *  :
##   HF eps > 1 treated as 1

# To get a nicer ANOVA table use function nice.anova (see ?noce.anova):
nice.anova(ez.glm("id", "value", obk.long, c("treatment", "gender"), c("phase", "hour"), "age"))

##                         Effect          df   MSE         F     p
## 1                    treatment        2, 9 23.96    3.58 +   .07
## 2                       gender        1, 9 23.96    3.95 +   .08
## 3                          age        1, 9 23.96      0.52   .49
## 4             treatment:gender        2, 9 23.96      1.28   .32
## 5                        phase  1.7, 15.28  3.91 20.28 *** <.001
## 6              treatment:phase 3.39, 15.28  3.91   6.07 **  .005
## 7                 gender:phase  1.7, 15.28  3.91      0.25   .75
## 8                    age:phase  1.7, 15.28  3.91    3.10 +   .08
## 9       treatment:gender:phase 3.39, 15.28  3.91      1.60   .23
## 10                        hour 2.14, 19.23  2.48 20.52 *** <.001
## 11              treatment:hour 4.27, 19.23  2.48      0.71   .60
## 12                 gender:hour 2.14, 19.23  2.48      0.71   .51
## 13                    age:hour 2.14, 19.23  2.48    2.82 +   .08
## 14       treatment:gender:hour 4.27, 19.23  2.48      0.59   .68
## 15                  phase:hour 3.48, 31.36  2.83      0.99   .42
## 16        treatment:phase:hour 6.97, 31.36  2.83      0.33   .93
## 17           gender:phase:hour 3.48, 31.36  2.83      0.90   .47
## 18              age:phase:hour 3.48, 31.36  2.83      0.77   .54
## 19 treatment:gender:phase:hour 6.97, 31.36  2.83      0.65   .71

# replicating ?Anova using aov.car:
aov.car(value ~ treatment * gender + Error(id/phase*hour), data = obk.long, type = 2)
# in contrast to aov you do not need the within-subject factors outside Error()

# replicating ?Anova using ez.glm:
ez.glm("id", "value", obk.long, c("treatment", "gender"), c("phase", "hour"), type = 2)

#both return:
## Type II Repeated Measures MANOVA Tests: Pillai test statistic
##                             Df test stat approx F num Df den Df       Pr(>F)
## (Intercept)                  1     0.970      318      1     10 0.0000000065 ***
## treatment                    2     0.481        5      2     10      0.03769 *
## gender                       1     0.204        3      1     10      0.14097
## treatment:gender             2     0.364        3      2     10      0.10447
## phase                        1     0.851       26      2      9      0.00019 ***
## treatment:phase              2     0.685        3      4     20      0.06674 .
## gender:phase                 1     0.043        0      2      9      0.82000
## treatment:gender:phase       2     0.311        1      4     20      0.47215
## hour                         1     0.935       25      4      7      0.00030 ***
## treatment:hour               2     0.301        0      8     16      0.92952
## gender:hour                  1     0.293        1      4      7      0.60237
## treatment:gender:hour        2     0.570        1      8     16      0.61319
## phase:hour                   1     0.550        0      8      3      0.83245
## treatment:phase:hour         2     0.664        0     16      8      0.99144
## gender:phase:hour            1     0.695        1      8      3      0.62021
## treatment:gender:phase:hour  2     0.793        0     16      8      0.97237
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# aggregating over one within-subjects factor (phase) with warning:

aov.car(value ~ treatment * gender + age + Error(id/hour), data = obk.long)

ez.glm("id", "value", obk.long, c("treatment", "gender"), "hour", "age")


# runs with "numeric" factors
obk.long$hour2 <- as.numeric(as.character(obk.long$hour))

aov.car(value ~ treatment * gender + Error(id/hour2), data = obk.long, type = 2)

# only between
aov.car(value ~ treatment * gender + age + Error(id), data = obk.long, type = 2)
aov.car(value ~ treatment * gender + Error(id), data = obk.long, type = 2)

ez.glm("id", "value", obk.long, c("treatment", "gender"), within = NULL, covariate = "age", type = 2, print.formula = TRUE)

ez.glm("id", "value", obk.long, c("treatment", "gender"), within = NULL, type = 2, print.formula = TRUE)

# only within

univariate(aov.car(value ~ Error(id/phase*hour), data = obk.long, type = 2))

univariate(ez.glm("id", "value", obk.long,  NULL, c("phase", "hour"), type = 2, print.formula = TRUE))
}
\author{
  \code{univariate} is basically a copy of
  \code{\link[car]{summary.Anova.mlm}} written by John
  Fox.\cr The other functions were written by Henrik
  Singmann.

  The design of these functions is heavily influenced by
  \code{\link[ez]{ezANOVA}} from package \pkg{ez}.
}
\references{
  Maxwell, S. E., & Delaney, H. D. (2004). \emph{Designing
  Experiments and Analyzing Data: A Model-Comparisons
  Perspective}. Mahwah, N.J.: Lawrence Erlbaum Associates.
}
\seealso{
  \code{\link{nice.anova}} is a function for creating nice
  ANOVA tables (including sphercitiy corrections) from
  objects returned by \code{ez.glm} and \code{aov.car}.

  \code{\link{mixed}} provides a (formula) interface for
  obtaining p-values for mixed-models via \pkg{lme4}.

  \code{\link{obk.long}} describes the long version of the
  \code{OBrienKaiser} dataset used in the examples.
}

