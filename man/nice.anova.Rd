\name{nice.anova}
\alias{nice.anova}
\title{Make nice ANOVA table for printing.}
\usage{
  nice.anova(object, es = NULL, correction = c("GG", "HF",
    "none"), sig.symbols = c(" +", " *", " **", " ***"),
    MSE = TRUE)
}
\arguments{
  \item{object}{An object of class \code{"Anova.mlm"} or
  \code{"anova"} as returned from \code{\link[car]{Anova}},
  \code{\link{ez.glm}}, or \code{\link{aov.car}}.}

  \item{es}{Effect Size to be reported. Currently none
  implemented.}

  \item{correction}{Character. Which sphericity correction
  on the degrees of freedom should be reported for the
  within-subject factors. The default \code{c("GG", "HF",
  "none")} corresponds to the Greenhouse-Geisser
  correction.}

  \item{sig.symbols}{Character. What should be the symbols
  designating significance? When entering an vector with
  \code{length(sig.symbol) < 4} only those elements of the
  default (\code{c(" +", " *", " **", " ***")}) will be
  replaced. \code{sig.symbols = ""} will display the stars
  but not the \code{+}, \code{sig.symbols = rep("", 4)}
  will display no symbols.}

  \item{MSE}{logical. Should the column containing the Mean
  Sqaured Error (MSE) be displayed? Default is
  \code{TRUE}.}
}
\value{
  A \code{data.frame} with the ANOVA table consisting of
  characters. The columns that are always present are:
  \code{Effect}, \code{df} (degrees of freedom), \code{F},
  and \code{p}.
}
\description{
  These functions produce a nice ANOVA table best for
  prointing. \code{nice.anova} takes an object from
  \code{\link[car]{Anova}} possible created by the
  convenience functions \code{\link{ez.glm}} or
  \code{\link{aov.car}}. When within-subject factors are
  present, either sphericity corrected or uncorrected
  degrees of freedom can be reported.
}
\details{
  The returned \code{data.frame} is print-ready when adding
  to a document with proper methods. I recommend
  \pkg{ascii} and \pkg{xtable}. \pkg{ascii} provides
  conversion to
  \href{http://www.methods.co.nz/asciidoc/}{AsciiDoc} but
  most notably to \href{http://orgmode.org/}{org-mode} (see
  \code{\link[ascii]{ascii}} and
  \code{\link[ascii]{print-ascii}}). \pkg{xtable} converts
  a \code{data.frame} into LaTeX code with many possible
  options (e.g., allowing for \code{"longtable"} or
  \code{"sidewaystable"}), see \code{\link[xtable]{xtable}}
  and \code{\link[xtable]{print.xtable}}. See Examples.

  Conversion functions to other formats (such as HTML, ODF,
  or Word) can be found at the
  \href{http://cran.r-project.org/web/views/ReproducibleResearch.html}{Reproducible
  Research Task View}.

  Effect sizes are the next thing to implement!
}
\examples{
# exampel using obk.long (see ?obk.long), a long version of the OBrienKaiser dataset from car.

data(obk.long)

# run univariate mixed ANCOVA for the full design:
nice.anova(aov.car(value ~ treatment * gender + age + Error(id/phase*hour), data = obk.long))

nice.anova(ez.glm("id", "value", obk.long, c("treatment", "gender"), c("phase", "hour"), "age"))

# no between
nice.anova(ez.glm("id", "value", obk.long, NULL, c("phase", "hour")))

# no within
nice.anova(ez.glm("id", "value", obk.long, c("treatment", "gender")))

nice.anova(ez.glm("id", "value", obk.long, c("treatment", "gender")), sig.symbol = rep("", 4))

\dontrun{
# use package ascii or xtable for formatting of tables ready for printing.

full <- nice.anova(ez.glm("id", "value", obk.long, c("treatment", "gender"), c("phase", "hour"), "age"))

require(ascii)
print(ascii(full, include.rownames = FALSE, caption = "ANOVA 1"), type = "org")

require(xtable)
print.xtable(xtable(full, caption = "ANOVA 2"), include.rownames = FALSE)
}
}
\seealso{
  \code{\link{ez.glm}} and \code{\link{aov.car}} are the
  convenience functions to create the object appropriate
  for \code{nice.anova}.
}

