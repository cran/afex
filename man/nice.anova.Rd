\encoding{UTF-8}
\name{nice.anova}
\alias{nice.anova}
\title{Make nice ANOVA table for printing.}
\usage{
nice.anova(object, es = "ges", observed = NULL, correction = c("GG", "HF", "none"),
     MSE = TRUE, intercept = FALSE, sig.symbols = c(" +", " *", " **", " ***"))
}
\arguments{
  \item{object}{An object of class \code{"Anova.mlm"} or
  \code{"anova"} as returned from \code{\link[car]{Anova}},
  \code{\link{ez.glm}}, or \code{\link{aov.car}}.}

  \item{es}{Effect Size to be reported. Default is
  \code{"ges"}, which reports generalized eta-squared (see
  details). Also supported is partial eta-squared
  (\code{"pes"}) or \code{"none"}.}

  \item{observed}{character vector referring to the
  observed (i.e., non manipulated) variables/effects in the
  design. Important for calculation of generalized
  eta-squared (ignored if \code{es} is not \code{"ges"}),
  see details.}

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

  \item{intercept}{logical. Should intercept (if present)
  be printed (default is \code{FALSE} which suppresses
  printing of the intercept)}
}
\value{
A \code{data.frame} with the ANOVA table consisting of
characters. The columns that are always present are:
\code{Effect}, \code{df} (degrees of freedom), \code{F},
and \code{p}.

\code{ges} contains the generalized eta-squared effect size
measure (Bakeman, 2005), \code{pes} contains partial
eta-squared (if requested).
}
\description{
These functions produce a nice ANOVA table best for
prointing. \code{nice.anova} takes an object from
\code{\link[car]{Anova}} possible created by the
convenience functions \code{\link{ez.glm}} or
\code{\link{aov.car}}. When within-subject factors are
present, either sphericity corrected or uncorrected degrees
of freedom can be reported.
}
\details{
The returned \code{data.frame} is print-ready when adding
to a document with proper methods. I recommend \pkg{ascii}
and \pkg{xtable}. \pkg{ascii} provides conversion to
\href{http://www.methods.co.nz/asciidoc/}{AsciiDoc} but
most notably to \href{http://orgmode.org/}{org-mode} (see
\code{\link[ascii]{ascii}} and
\code{\link[ascii]{print-ascii}}). \pkg{xtable} converts a
\code{data.frame} into LaTeX code with many possible
options (e.g., allowing for \code{"longtable"} or
\code{"sidewaystable"}), see \code{\link[xtable]{xtable}}
and \code{\link[xtable]{print.xtable}}. See Examples.

Conversion functions to other formats (such as HTML, ODF,
or Word) can be found at the
\href{http://cran.r-project.org/web/views/ReproducibleResearch.html}{Reproducible
Research Task View}.

The default reports generalized eta squared (Olejnik &
Algina, 2003), the "recommended effect size for repeated
measured designs" (Bakeman, 2005). Note that it is
important that all measured variables (as opposed to
experimentally manipulated variables), such as e.g., age,
gender, weight, ..., must be declared via \code{observed}
to obtain the correct effect size estimate. Partial eta
squared (\code{"pes"}) does not require this.
}
\examples{
## example from Olejnik & Algina (2003)
# "Repeated Measures Design" (pp. 439):
data(md_12.1)
# create object of class Anova:
rmd <- ez.glm("id", "rt", md_12.1, within = c("angle", "noise"),
              return = "Anova")
# use different es:
nice.anova(rmd, es = "pes") # noise: .82
nice.anova(rmd, es = "ges") # noise: .39

# exampel using obk.long (see ?obk.long), a long version of the OBrienKaiser dataset from car.
data(obk.long)
# create object of class Anova:
tmp.aov <- aov.car(value ~ treatment * gender + Error(id/phase*hour),
             data = obk.long, return = "Anova")

nice.anova(tmp.aov, observed = "gender")

nice.anova(tmp.aov, observed = "gender", sig.symbol = rep("", 4))

\dontrun{
# use package ascii or xtable for formatting of tables ready for printing.

full <- nice.anova(tmp.aov, observed = "gender")

require(ascii)
print(ascii(full, include.rownames = FALSE, caption = "ANOVA 1"), type = "org")

require(xtable)
print.xtable(xtable(full, caption = "ANOVA 2"), include.rownames = FALSE)
}
}
\author{
The code for calculating generalized eta-squared was
written by Mike Lawrence.\cr Everything else was written by
Henrik Singmann.
}
\references{
Bakeman, R. (2005). Recommended effect size statistics for
repeated measures designs. \emph{Behavior Research
Methods}, 37(3), 379-384. doi:10.3758/BF03192707

Olejnik, S., & Algina, J. (2003). Generalized Eta and Omega
Squared Statistics: Measures of Effect Size for Some Common
Research Designs. \emph{Psychological Methods}, 8(4),
434-447. doi:10.1037/1082-989X.8.4.434
}
\seealso{
\code{\link{ez.glm}} and \code{\link{aov.car}} are the
convenience functions to create the object appropriate for
\code{nice.anova}.
}

