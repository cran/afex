\name{mixed}
\alias{mixed}
\title{Obtain p-values for a mixed-model from lmer().}
\usage{
  mixed(formula, data, type = 3, method = c("KR", "PB"),
    ...)
}
\arguments{
  \item{formula}{a formula describing the full mixed-model
  to be fitted. As this formula is passed to \code{lmer},
  it needs at least one random term.}

  \item{data}{data.frame containing the data. Should have
  all the variables present in \code{fixed}, \code{random},
  and \code{dv} as columns.}

  \item{type}{type of test on which effects are based.
  Currently only type 3 tests (\code{3} or \code{"III"})
  are correctly implemented (see Details).}

  \item{method}{character vector indicating which methods
  for obtaining p-values should be used. Currently only
  \code{"KR"} is implemented corresponding to the
  Kenward-Rogers approximation for degrees of freedom.}

  \item{...}{further arguments (such as \code{weights})
  passed to \code{\link{lmer}}.}
}
\value{
  An object of class \code{"mixed"} (i.e., a list) with the
  following elements:

  \enumerate{ \item \code{anova.table} a data.frame
  containing the statistics returned from
  \code{\link[pbkrtest]{KRmodcomp}}. \item
  \code{full.model} the \code{"mer"} object returned from
  fitting the full mixed model. \item
  \code{restricted.models} a list of \code{"mer"} objects
  from fitting the restricted models (i.e., each model
  lacks the corresponding effect) \item \code{tests} a list
  of objects returned by the function for obtaining the
  p-values (objects are of class \code{"KRmodcomp"} when
  \code{method = "KR"}). \item \code{type} The \code{type}
  argument used when calling this function. \item
  \code{method} The \code{method} argument used when
  calling this function. }

  The following methods exist for objects of class
  \code{"mixed"}: \code{print} (which uses rounding and
  only returns the results wiuth \code{F.scaling = 1}),
  \code{summary}, and \code{anova} (the latter two return
  the same data.frame).
}
\description{
  Fits and calculates p-values for all effects in a mixed
  model fitted with \code{\link[lme4]{lmer}}. The default
  behavior calculates type 3 like p-values using the
  Kenward-Rogers approximation for degrees-of-freedom
  implemented in \code{\link[pbkrtest]{KRmodcomp}}.
  \code{print}, \code{summary}, and \code{anova} methods
  for the returned object of class \code{"mixed"} are
  available (all return the same data.frame).
}
\details{
  Type 3 tests are obtained by comparing a model in which
  only the corresponding effect is missing with the full
  model (containing all effects). This corresponds to the
  (type 3) Wald tests given by \code{car::Anova} for
  \code{"mer"} models (from version 2.0-13).

  Type 2 tests are obtained by comparing a model in which
  the corresponding effect and all higher oder effect
  (e.g., all three-way interactions for a two-way
  interaction) are missing with a model in which all
  effects of the relevant order are present and all higher
  order effects absent. Consequently, the results for lower
  order effects are identical of wether or not higher order
  effects are part of the model or not, which is rather
  dubious (but
  \href{https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q3/018992.html}{I
  didn't find a better way} of implementing the Type 2
  tests). This \strong{does not} correspond to the (type 2)
  Wald Test reported by \code{car::Anova}. If you want type
  2 tests, use \code{car::Anova} with \code{test = "F"}
  (from version 2.0-13) instead of this function.

  For an introduction to mixed-modeling for experimental
  designs using p-values see Judd, Westfall, and Kenny
  (2012). Further introductions to mixed-modeling for
  experimental designs are given by Baayen and colleagues
  (Baayen, 2008; Baayen, Davidson & Bates, 2008; Baayen &
  Milin, 2010).
}
\note{
  This function is not thoroughly tested so please report
  all bugs to henrik.singmann (at)
  psychologie.uni-freiburg.de \cr There might be problems
  with rather big models when constructing the model matrix
  to fit the \code{lmer} models (potentially problematic
  with Type 2 sums of squares). If you find any such bug,
  please send an example including code and data!

  This functions needs a lot of RAM and rather long time
  especially with complex random structures. The RAM demand
  is a problem especially on 32 bit Windows which only
  supports up to 2 or 3GB RAM (see
  \href{http://cran.r-project.org/bin/windows/base/rw-FAQ.html}{R
  Windows FAQ}).

  This function calls \code{lme4:::nobars} for dealing with
  the formula. So any significant changes to \pkg{lme4} or
  \code{lme4:::nobars} may disrupt its functionality.
}
\examples{
\dontrun{
# example data from package languageR:
# Lexical decision latencies elicited from 21 subjects for 79 English concrete nouns, with variables linked to subject or word.
data(lexdec, package = "languageR")

# using the simplest model
m1 <- mixed(RT ~ Correct + Trial + PrevType * meanWeight + Frequency + NativeLanguage * Length + (1|Subject) + (1|Word), data = lexdec)

m1
# gives:
##                   Effect df1       df2      Fstat p.value
## 1            (Intercept)   1   96.6379 13573.1410  0.0000
## 2                Correct   1 1627.7303     8.1452  0.0044
## 3                  Trial   1 1592.4301     7.5738  0.0060
## 4               PrevType   1 1605.3939     0.1700  0.6802
## 5             meanWeight   1   75.3919    14.8545  0.0002
## 6              Frequency   1   76.0821    56.5348  0.0000
## 7         NativeLanguage   1   27.1213     0.6953  0.4117
## 8                 Length   1   75.8259     8.6959  0.0042
## 9    PrevType:meanWeight   1 1601.1850     6.1823  0.0130
## 10 NativeLanguage:Length   1 1555.4858    14.2445  0.0002
}
}
\author{
  Henrik Singmann with contributions from
  \href{http://stackoverflow.com/q/11335923/289572}{Ben
  Bolker and Joshua Wiley}.
}
\references{
  Baayen, R. H. (2008). \emph{Analyzing linguistic data: a
  practical introduction to statistics using R}. Cambridge,
  UK; New York: Cambridge University Press.

  Baayen, R. H., Davidson, D. J., & Bates, D. M. (2008).
  Mixed-effects modeling with crossed random effects for
  subjects and items. \emph{Journal of Memory and
  Language}, 59(4), 390-412. doi:10.1016/j.jml.2007.12.005

  Baayen, R. H., & Milin, P. (2010). Analyzing Reaction
  Times. \emph{International Journal of Psychological
  Research}, 3(2), 12-28.

  Judd, C. M., Westfall, J., & Kenny, D. A. (2012).
  Treating stimuli as a random factor in social psychology:
  A new and comprehensive solution to a pervasive but
  largely ignored problem. \emph{Journal of Personality and
  Social Psychology}, 103(1), 54–69. doi:10.1037/a0028347
}
\seealso{
  \code{\link{ez.glm}} and \code{\link{aov.car}} for
  convenience functions to analyze experimental deisgns
  with classical ANOVA or ANCOVA wrapping
  \code{\link[car]{Anova}}.
}

