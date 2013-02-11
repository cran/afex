\name{aov.car}
\alias{aov.car}
\alias{ez.glm}
\alias{univ}
\title{Convenience wrappers for car::Anova using either a formula or factor based interface.}
\usage{
  aov.car(formula, data, fun.aggregate = NULL, type = 3,
    return = "nice", observed = NULL, args.return = list(),
    ...)

  ez.glm(id, dv, data, between = NULL, within = NULL,
    covariate = NULL, observed = NULL, fun.aggregate =
    NULL, type = 3, return = "nice", args.return = list(),
    ..., print.formula = FALSE)

  univ(object)
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

  \item{observed}{\code{character} vector indicating which
  of the variables are observed (i.e, measured) as compared
  to experimentally manipulated. The default behavior is to
  return a ANOVA table with generalized eta squared effect
  size for which this information is necessary (see
  \code{\link{nice.anova}}).}

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

  \item{return}{What should be returned? If \code{"nice"}
  (the default) will return a nice ANOVA table (produced by
  \code{\link{nice.anova}}. Possible values are
  \code{c("Anova", "lm", "data", "nice", "full", "all",
  "univariate")} (possibly abbreviated).}

  \item{args.return}{\code{list} of further arguments
  passed to the function which produces the return value.
  Currently only supports \code{return = "nice"} (the
  default) which then passes arguments to
  \code{\link{nice.anova}} (see examples).}

  \item{...}{Further arguments passed to
  \code{fun.aggregate}.}

  \item{object}{An object of class \code{Anova.mlm} as
  returned by \code{aov.car}, \code{ez.glm}, or
  \code{\link[car]{Anova}}.}
}
\value{
  \code{aov.car} and \code{ez.glm} are wrappers to
  \code{\link[car]{Anova}}, the return value is dependent
  on the \code{return} argument. When argument
  \code{return} is \code{"nice"} (the default) a nice ANOVA
  table is returnd (\code{\link{nice.anova}}). If
  \code{return = "full"} or \code{return = "all"} a list
  \code{list} with the following elements:

  \describe{ \item{"Anova"}{the same as
  \code{\link[car]{Anova}}. Usually an object of class
  \code{"Anova.mlm"} (with within-subjects factors) or of
  class \code{c("anova", "data.frame")}. Also returned if
  \code{return = "Anova"}.} \item{"lm"}{the object fitted
  with \code{lm} and passed to \code{Anova} (i.e., an
  object of class \code{"lm"} or \code{"mlm"}). Also
  returned if \code{return = "lm"}.} \item{"data"}{the data
  used to fit the \code{lm} object. Also returned if
  \code{return = "data"}.} \item{"idata"}{if within-subject
  factors are present, the \code{idata} argument passed to
  \code{Anova}.} }

  If \code{return = "univariate"} the object returned from
  \code{univ}.

  \code{univ} returns a \code{list} of \code{data.frame}s
  containing the univariate results (i.e., the classical
  ANOVA results) from an object of class
  \code{"Anova.mlm"}. This is essentially the output from
  \code{summary.Anova.mlm} with \code{multivariate =
  FALSE}, e.g. \code{summary(aov.car(...), multivriate =
  FALSE)}, as a list instead of printed to the console.\cr
  For objects of class \code{"anova"} (i.e., the object
  returned by \code{car::Anova} for a purely
  between-subjects ANOVA) the object is returned unaltered.

  The elements of the list returned by \code{univ} are:
  \code{anova}, \code{mauchly}, and
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

  Note that the default behavior is to return a
  \code{\link{nice.anova}} \code{data.frame}. This includes
  calculation of generalized eta squared for which
  \strong{all non manipluated (i.e., observed)} variables
  need to be specified via the \code{observed} argument.
  Changing the effect size to \code{"pes"} (partial
  eta-squared) via \code{args.return} or the return value
  via \code{return} removes this necessity.
}
\note{
  Variables entered as within-subjects (i.e., repeated
  measures) factors are silently converted to factors and
  unused levels dropped.

  Contrasts attached to a factor as an attribute are
  probably not preserved and not supported.

  Function \code{univ} was called \code{univariate} in
  prior versions, but there was a function with similar
  name in package \pkg{multcomp} leading to bugs and
  unexpected behavior.
}
\examples{

# exampel using obk.long (see ?obk.long), a long version of the OBrienKaiser dataset from car.

data(obk.long, package = "afex")

# run univariate mixed ANCOVA for the full design:
aov.car(value ~ treatment * gender + age + Error(id/phase*hour), data = obk.long, observed = c("gender", "age"))

ez.glm("id", "value", obk.long, between = c("treatment", "gender"), within = c("phase", "hour"), covariate = "age", observed = c("gender", "age"))

# both calls return the same:
##                         Effect          df   MSE         F  ges     p
## 1                    treatment        2, 9 23.96    3.58 +  .20   .07
## 2                       gender        1, 9 23.96    3.95 +  .14   .08
## 3                          age        1, 9 23.96      0.52  .02   .49
## 4             treatment:gender        2, 9 23.96      1.28  .09   .32
## 5                        phase  1.7, 15.28  3.91 20.28 ***  .17 <.001
## 6              treatment:phase 3.39, 15.28  3.91   6.07 **  .11  .005
## 7                 gender:phase  1.7, 15.28  3.91      0.25 .002   .75
## 8                    age:phase  1.7, 15.28  3.91    3.10 +  .03   .08
## 9       treatment:gender:phase 3.39, 15.28  3.91      1.60  .03   .23
## 10                        hour 2.14, 19.23  2.48 20.52 ***  .14 <.001
## 11              treatment:hour 4.27, 19.23  2.48      0.71  .01   .60
## 12                 gender:hour 2.14, 19.23  2.48      0.71 .006   .51
## 13                    age:hour 2.14, 19.23  2.48    2.82 +  .02   .08
## 14       treatment:gender:hour 4.27, 19.23  2.48      0.59 .009   .68
## 15                  phase:hour 3.48, 31.36  2.83      0.99  .01   .42
## 16        treatment:phase:hour 6.97, 31.36  2.83      0.33 .010   .93
## 17           gender:phase:hour 3.48, 31.36  2.83      0.90  .01   .47
## 18              age:phase:hour 3.48, 31.36  2.83      0.77  .01   .54
## 19 treatment:gender:phase:hour 6.97, 31.36  2.83      0.65  .02   .71

# replicating ?Anova using aov.car:
aov.car(value ~ treatment * gender + Error(id/phase*hour), data = obk.long, type = 2, return = "Anova")
# in contrast to aov you do not need the within-subject factors outside Error()

# replicating ?Anova using ez.glm:
ez.glm("id", "value", obk.long, c("treatment", "gender"), c("phase", "hour"), type = 2, return = "Anova")

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
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 


# aggregating over one within-subjects factor (phase) with warning:

aov.car(value ~ treatment * gender + age + Error(id/hour), data = obk.long, observed = c("gender", "age"))

ez.glm("id", "value", obk.long, c("treatment", "gender"), "hour", "age", observed = c("gender", "age"))


# runs with "numeric" factors
obk.long$hour2 <- as.numeric(as.character(obk.long$hour))

aov.car(value ~ treatment * gender + Error(id/hour2), data = obk.long, type = 2,observed = c("gender"))

# only between
aov.car(value ~ treatment * gender + age + Error(id), data = obk.long, type = 2,observed = c("gender", "age"))
aov.car(value ~ treatment * gender + Error(id), data = obk.long, type = 2, observed = c("gender"))

ez.glm("id", "value", obk.long, c("treatment", "gender"), within = NULL, covariate = "age", type = 2, print.formula = TRUE, observed = c("gender", "age"))

ez.glm("id", "value", obk.long, c("treatment", "gender"), within = NULL, type = 2, print.formula = TRUE, observed = c("gender"))

# only within

aov.car(value ~ Error(id/phase*hour), data = obk.long, type = 2)

ez.glm("id", "value", obk.long,  NULL, c("phase", "hour"), type = 2, print.formula = TRUE)

# using return = "full":

str(aov.car(value ~ Error(id/phase*hour), data = obk.long, return = "full"), 1)

## List of 4
##  $ Anova:List of 14
##   ..- attr(*, "class")= chr "Anova.mlm"
##  $ lm   :List of 11
##   ..- attr(*, "class")= chr [1:2] "mlm" "lm"
##  $ data :'data.frame':  16 obs. of  16 variables:
##  $ idata:'data.frame':  15 obs. of  2 variables:

# use args.return arguments:
aov.car(value ~ treatment * gender + age + Error(id/phase*hour), data = obk.long, args.return = list(correction = "none", es = "pes"))

aov.car(value ~ treatment * gender + age + Error(id/phase*hour), data = obk.long,observed = c("gender", "age"), args.return = list(correction = "none", MSE = FALSE))

}
\author{
  \code{univ} is basically a copy of
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
  \code{\link{nice.anova}} creats the nice ANOVA tables
  which are by default returned. See also there for a
  slightly longer discussion of effect sizes.

  \code{\link{mixed}} provides a (formula) interface for
  obtaining p-values for mixed-models via \pkg{lme4}.

  \code{\link{obk.long}} describes the long version of the
  \code{OBrienKaiser} dataset used in the examples.
}

