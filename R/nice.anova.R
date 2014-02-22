#' Make nice ANOVA table for printing.
#'
#' These functions produce a nice ANOVA table best for prointing. \code{nice.anova} takes an object from \code{\link[car]{Anova}} possible created by the convenience functions \code{\link{ez.glm}} or \code{\link{aov.car}}. When within-subject factors are present, either sphericity corrected or uncorrected degrees of freedom can be reported.
#' 
#' @usage nice.anova(object, es = "ges", observed = NULL, correction = c("GG", "HF", "none"), 
#'      MSE = TRUE, intercept = FALSE, sig.symbols = c(" +", " *", " **", " ***"))
#' 
#'
#' @param object An object of class \code{"Anova.mlm"} or \code{"anova"} as returned from \code{\link[car]{Anova}},  \code{\link{ez.glm}}, or \code{\link{aov.car}}.
#' @param es Effect Size to be reported. Default is \code{"ges"}, which reports generalized eta-squared (see details). Also supported is partial eta-squared (\code{"pes"}) or \code{"none"}.
#' @param observed character vector referring to the observed (i.e., non manipulated) variables/effects in the design. Important for calculation of generalized eta-squared (ignored if \code{es} is not \code{"ges"}), see details.
#' @param correction Character. Which sphericity correction on the degrees of freedom should be reported for the within-subject factors. The default \code{c("GG", "HF", "none")} corresponds to the Greenhouse-Geisser correction.
#' @param sig.symbols Character. What should be the symbols designating significance? When entering an vector with \code{length(sig.symbol) < 4} only those elements of the default (\code{c(" +", " *", " **", " ***")}) will be replaced. \code{sig.symbols = ""} will display the stars but not the \code{+}, \code{sig.symbols = rep("", 4)} will display no symbols.
#' @param MSE logical. Should the column containing the Mean Sqaured Error (MSE) be displayed? Default is \code{TRUE}.
#' @param intercept logical. Should intercept (if present) be printed (default is \code{FALSE} which suppresses printing of the intercept)
#'
#' @return A \code{data.frame} with the ANOVA table consisting of characters. The columns that are always present are: \code{Effect}, \code{df} (degrees of freedom), \code{F}, and \code{p}.
#'
#' \code{ges} contains the generalized eta-squared effect size measure (Bakeman, 2005), \code{pes} contains partial eta-squared (if requested).
#'
#' @details The returned \code{data.frame} is print-ready when adding to a document with proper methods. I recommend \pkg{ascii} and \pkg{xtable}. \pkg{ascii} provides conversion to \href{http://www.methods.co.nz/asciidoc/}{AsciiDoc} but most notably to \href{http://orgmode.org/}{org-mode} (see \code{\link[ascii]{ascii}} and \code{\link[ascii]{print-ascii}}). \pkg{xtable} converts a \code{data.frame} into LaTeX code with many possible options (e.g., allowing for \code{"longtable"} or \code{"sidewaystable"}), see \code{\link[xtable]{xtable}} and \code{\link[xtable]{print.xtable}}. See Examples.
#'
#' Conversion functions to other formats (such as HTML, ODF, or Word) can be found at the \href{http://cran.r-project.org/web/views/ReproducibleResearch.html}{Reproducible Research Task View}.
#'
#' The default reports generalized eta squared (Olejnik & Algina, 2003), the "recommended effect size for repeated measured designs" (Bakeman, 2005). Note that it is important that all measured variables (as opposed to experimentally manipulated variables), such as e.g., age, gender, weight, ..., must be declared via \code{observed} to obtain the correct effect size estimate. Partial eta squared (\code{"pes"}) does not require this.
#'
#' @seealso \code{\link{ez.glm}} and \code{\link{aov.car}} are the convenience functions to create the object appropriate for \code{nice.anova}.
#'
#' @author The code for calculating generalized eta-squared was written by Mike Lawrence.\cr Everything else was written by Henrik Singmann.
#'
#' @references Bakeman, R. (2005). Recommended effect size statistics for repeated measures designs. \emph{Behavior Research Methods}, 37(3), 379-384. doi:10.3758/BF03192707

#'
#' Olejnik, S., & Algina, J. (2003). Generalized Eta and Omega Squared Statistics: Measures of Effect Size for Some Common Research Designs. \emph{Psychological Methods}, 8(4), 434-447. doi:10.1037/1082-989X.8.4.434
#' 
#' @name nice.anova
#' @export nice.anova
#'
#' @examples
#'
#' ## example from Olejnik & Algina (2003)
#' # "Repeated Measures Design" (pp. 439):
#' data(md_12.1)
#' # create object of class Anova:
#' rmd <- ez.glm("id", "rt", md_12.1, within = c("angle", "noise"),
#'               return = "Anova")
#' # use different es:
#' nice.anova(rmd, es = "pes") # noise: .82
#' nice.anova(rmd, es = "ges") # noise: .39

#'
#' # exampel using obk.long (see ?obk.long), a long version of the OBrienKaiser dataset from car.
#' data(obk.long)
#' # create object of class Anova:
#' tmp.aov <- aov.car(value ~ treatment * gender + Error(id/phase*hour), 
#'              data = obk.long, return = "Anova")
#' 
#' nice.anova(tmp.aov, observed = "gender")
#' 
#' nice.anova(tmp.aov, observed = "gender", sig.symbol = rep("", 4))
#' 
#' \dontrun{
#' # use package ascii or xtable for formatting of tables ready for printing.
#' 
#' full <- nice.anova(tmp.aov, observed = "gender")
#' 
#' require(ascii)
#' print(ascii(full, include.rownames = FALSE, caption = "ANOVA 1"), type = "org")
#' 
#' require(xtable)
#' print.xtable(xtable(full, caption = "ANOVA 2"), include.rownames = FALSE)
#' }
#' 
#' 

nice.anova <- function(object, es = "ges", observed = NULL, correction = c("GG", "HF", "none"), MSE = TRUE, intercept = FALSE, sig.symbols = c(" +", " *", " **", " ***")) {
  # internal functions:
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  # 	round.ps <- function(x) {
  # 		as.character(ifelse(x < 0.001, "<.001", substr(ifelse(x < 0.01, formatC(x, digits = 3, format = "f"), ifelse(round(x, 2) == 1, " >.99", formatC(x, digits = 2, format = "f"))), 2, 5)))
  # 	}
  make.fs <- function(anova, symbols) {
    ifelse(anova[["Pr(>F)"]] < 0.001, str_c(formatC(anova[["F"]], digits = 2, format = "f"), symbols[4]), 
           ifelse(anova[["Pr(>F)"]] < 0.01, str_c(formatC(anova[["F"]], digits = 2, format = "f"), symbols[3]), 
                  ifelse(anova[["Pr(>F)"]] < 0.05, str_c(formatC(anova[["F"]], digits = 2, format = "f"), symbols[2]), 
                         ifelse(anova[["Pr(>F)"]] < 0.1, str_c(formatC(anova[["F"]], digits = 2, format = "f"), symbols[1]), formatC(anova[["F"]], digits = 2, format = "f")))))
  }
  # check arguments
  es <- match.arg(es, c("none", "ges", "pes"), several.ok = TRUE)
  #browser()
  if (class(object)[1] == "Anova.mlm") {
    tmp <- suppressWarnings(univ(object))
    t.out <- tmp[["anova"]]
    if (correction[1] == "GG") {
      t.out[row.names(tmp[["sphericity.correction"]]), "num Df"] <- t.out[row.names(tmp[["sphericity.correction"]]), "num Df"] * tmp[["sphericity.correction"]][,"GG eps"]
      t.out[row.names(tmp[["sphericity.correction"]]), "den Df"] <- t.out[row.names(tmp[["sphericity.correction"]]), "den Df"] * tmp[["sphericity.correction"]][,"GG eps"]
      t.out[row.names(tmp[["sphericity.correction"]]), "Pr(>F)"] <- tmp[["sphericity.correction"]][,"Pr(>F[GG])"]
    } else {
      if (correction[1] == "HF") {
        if (any(tmp[["sphericity.correction"]][,"HF eps"] > 1)) warning("HF eps > 1 treated as 1")
        t.out[row.names(tmp[["sphericity.correction"]]), "num Df"] <- t.out[row.names(tmp[["sphericity.correction"]]), "num Df"] * pmin(1, tmp[["sphericity.correction"]][,"HF eps"])
        t.out[row.names(tmp[["sphericity.correction"]]), "den Df"] <- t.out[row.names(tmp[["sphericity.correction"]]), "den Df"] * pmin(1, tmp[["sphericity.correction"]][,"HF eps"])
        t.out[row.names(tmp[["sphericity.correction"]]), "Pr(>F)"] <- tmp[["sphericity.correction"]][,"Pr(>F[HF])"]
      } else {
        if (correction[1] == "none") {
          TRUE
        } else stop("None supported argument to correction.")
      }
    }
    tmp.df <- t.out		
  } else {
    if (class(object)[1] == "anova") {
      #browser()
      #class(object) <- "data.frame"
      tmp.df <- cbind(object[-nrow(object),], data.frame("Error SS" = object[nrow(object), "Sum Sq"], "den Df" = object[nrow(object), "Df"], check.names = FALSE))
      colnames(tmp.df)[1:3] <- c("SS", "num Df", "F")
    } else stop("Non-supported object passed. Object must be of class 'Anova.mlm' or 'anova'.")
  }
  tmp2 <- as.data.frame(tmp.df)
  tmp2[,"df"] <- paste(ifelse(is.wholenumber(tmp2[,"num Df"]), tmp2[,"num Df"], formatC(tmp2[,"num Df"], digits = 2, format = "f")),  ifelse(is.wholenumber(tmp2[,"den Df"]),tmp2[,"den Df"], formatC(tmp2[,"den Df"], digits = 2, format = "f")), sep = ", ")
  tmp2[,"MSE"] <- tmp2[,"Error SS"]/tmp2[,"den Df"]
  symbols.use <-  c(" +", " *", " **", " ***")
  symbols.use[seq_along(sig.symbols)] <- sig.symbols
  df.out <- data.frame(Effect = row.names(tmp2), df = tmp2[,"df"], stringsAsFactors = FALSE)
  if (MSE) df.out <- cbind(df.out, data.frame(MSE = formatC(tmp2[,"MSE"], digits = 2, format = "f"), stringsAsFactors = FALSE))
  df.out <- cbind(df.out, data.frame(F = make.fs(tmp2, symbols.use), stringsAsFactors = FALSE))
  # calculate es
  if ("pes" %in% es) {
    df.out <- cbind(df.out, pes = round_ps(tmp2$SS/(tmp2$SS + tmp2[,"Error SS"])), stringsAsFactors = FALSE)
  }
  if ("ges" %in% es) {
    # This code is basically a copy from ezANOVA by Mike Lawrence!
    if(!is.null(observed)){
      obs <- rep(FALSE,nrow(tmp2))
      for(i in observed){
        if (!any(str_detect(rownames(tmp2),str_c("\\<",i,"\\>")))) stop(str_c("Observed variable not in data: ", i))
        obs <- obs | str_detect(rownames(tmp2),str_c("\\<",i,"\\>"))
      }
      obs_SSn1 <- sum(tmp2$SS*obs)
      obs_SSn2 <- tmp2$SS*obs
    }else{
      obs_SSn1 <- 0
      obs_SSn2 <- 0
    }
    df.out$ges <- round_ps(tmp2$SS/(tmp2$SS+sum(unique(tmp2[,"Error SS"]))+obs_SSn1-obs_SSn2))
  }
  df.out$p  <-  round_ps(tmp2[,"Pr(>F)"])
  if (!intercept) if (df.out[1,1] == "(Intercept)")  df.out <- df.out[-1,, drop = FALSE]
  rownames(df.out) <- NULL
  df.out
}

