#' Convenience wrappers for car::Anova using either a formula or factor based interface.
#'
#' These functions allow convenient access to \code{\link[car]{Anova}} (from the \pkg{car} package) for data in the \strong{long} format (i.e., one observation per row), possibly aggregating the data if there is more than one obersvation per individuum and cell. Hence, mixed between-within ANOVAs can be calculated conveniently without using the rather unhandy format of \code{car::Anova}. \code{aov.car} can be called using a formula similar to \code{\link{aov}} specifying an error strata for the within-subject factor(s). \code{ez.glm} is called specifying the factors as character vectors.
#'
#' @usage aov.car(formula, data, fun.aggregate = NULL, type = 3, 
#'      factorize = TRUE, check.contrasts = TRUE, 
#'      return = "nice", observed = NULL, args.return = list(), ...)
#'
#' ez.glm(id, dv, data, between = NULL, within = NULL, covariate = NULL, 
#'      observed = NULL, fun.aggregate = NULL, type = 3, 
#'      factorize = TRUE, check.contrasts = TRUE, 
#'      return = "nice", args.return = list(), ..., print.formula = FALSE)
#' 
#' univ(object)
#'
#' @param formula A formula specifying the ANOVA model similar to \code{\link{aov}}. Should include an error term (i.e., \code{Error( / )}). Note that the within-subject factors do not need to be outside the Error term (this contrasts with \code{aov}). See Details.
#' @param id \code{character} vector (of length 1) indicating the subject identifier column in \code{data}.
#' @param dv \code{character} vector (of length 1) indicating the column containing the \strong{dependent variable} in \code{data}.
#' @param between \code{character} vector indicating the \strong{between}-subject(s) factor(s)/column(s) in \code{data}. Default is \code{NULL} indicating no between-subjects factors.
#' @param within \code{character} vector indicating the \strong{within}-subject(s) factor(s)/column(s) in \code{data}.  Default is \code{NULL} indicating no within-subjects factors.
#' @param covariate \code{character} vector indicating the between-subject(s) covariate(s) (i.e., column(s)) in \code{data}. Default is \code{NULL} indicating no covariates.
#' @param observed \code{character} vector indicating which of the variables are observed (i.e, measured) as compared to experimentally manipulated. The default behavior is to return a ANOVA table with generalized eta squared effect size for which this information is necessary (see \code{\link{nice.anova}}).
#' @param data A \code{data.frame} containing the data. Mandatory.
#' @param fun.aggregate The function for aggregating the data before running the ANOVA if there is more than one obervation per individuum and cell of the design. The default \code{NULL} issues a warning if aggregation is necessary and uses \code{\link{mean}}.
#' @param type The type of sums of squares for the ANOVA. \strong{Defaults to 3}. Passed to \code{\link[car]{Anova}}. Possible values are \code{"II"}, \code{"III"}, \code{2}, or \code{3}.
#' @param factorize logical. Should between subject factors be factorized (with note) before running the analysis. Default is \code{TRUE}. If one wants to run an ANCOVA, needs to be set to \code{FALSE} (in which case centering on 0 is checked on numeric variables).
#' @param check.contrasts \code{logical}. Should contrasts for between-subject factors be checked and (if necessary) changed to be \code{"contr.sum"}. See details.
#' @param print.formula \code{ez.glm} is a wrapper for \code{aov.car}. This boolean argument indicates whether the formula in the call to \code{car.aov} should be printed. 
#' @param return What should be returned? If \code{"nice"} (the default) will return a nice ANOVA table (produced by \code{\link{nice.anova}}. Possible values are \code{c("Anova", "lm", "data", "nice", "full", "all", "univariate", "lme4")} (possibly abbreviated).
#' @param args.return \code{list} of further arguments passed to the function which produces the return value. Currently only supports \code{return = "nice"} (the default) which then passes arguments to \code{\link{nice.anova}} (see examples).
#' @param ... Further arguments passed to \code{fun.aggregate}.
#' @param object An object of class \code{Anova.mlm} as returned by \code{aov.car}, \code{ez.glm}, or \code{\link[car]{Anova}}.
#'
#' @return \code{aov.car} and \code{ez.glm} are wrappers to \code{\link[car]{Anova}}, the return value is dependent on the \code{return} argument. When argument \code{return} is \code{"nice"} (the default) a nice ANOVA table is returnd (\code{\link{nice.anova}}) with the following columns: \code{Effect}, \code{df}, \code{MSE} (mean-squared errors), \code{F} (potentially with significant symbols), \code{ges} (generalized eta-squared), \code{p}.
#'
#' If \code{return = "full"} or \code{return = "all"} a list \code{list} with the following elements:
#'
#' \describe{
#'   \item{"Anova"}{the same as \code{\link[car]{Anova}}. Usually an object of class \code{"Anova.mlm"} (with within-subjects factors) or of class \code{c("anova", "data.frame")}. Also returned if \code{return = "Anova"}.}
#'   \item{"lm"}{the object fitted with \code{lm} and passed to \code{Anova} (i.e., an object of class \code{"lm"} or \code{"mlm"}). Also returned if \code{return = "lm"}.}
#'   \item{"data"}{the data used to fit the \code{lm} object. Also returned if \code{return = "data"}.}
#'   \item{"idata"}{if within-subject factors are present, the \code{idata} argument passed to \code{Anova}.}
#' }
#'
#' If \code{return = "univariate"} the object returned from \code{univ}.
#' 
#' \code{univ} returns a \code{list} of \code{data.frame}s containing the univariate results (i.e., the classical ANOVA results) from an object of class \code{"Anova.mlm"}. This is essentially the output from \code{summary.Anova.mlm} with \code{multivariate = FALSE}, e.g. \code{summary(aov.car(...), multivriate = FALSE)}, as a list instead of printed to the console.\cr
#' For objects of class \code{"anova"} (i.e., the object returned by \code{car::Anova} for a purely between-subjects ANOVA) the object is returned unaltered.
#'
#' The elements of the list returned by \code{univ} are: \code{anova}, \code{mauchly}, and \code{spehricity.correction} (containing both, Greenhouse-Geisser and Hyundt-Feldt correction).
#' 
#' If \code{return = "lme4"} The data (possibly aggregated to have one observation per cell) is fitted with \code{lme4::lmer} using all within-subjects factors as fully crossed random slopes and the obtained \code{mer} object returned for further analysis. (This behavhior is rather experimental!)
#' 
#' @details \strong{Type 3 sums of squares are default in \pkg{afex}.} Note that type 3 sums of squares are said to be dangerous and/or problematic. On the other side they are the default in in SPSS and SAS and recommended by e.g. Maxwell and Delaney (2004). For a brief discussion see \href{http://stats.stackexchange.com/q/6208/442}{here}. 
#'
#' However, note that lower order effects (e.g., main effects) in type 3 ANOVAs are only meaningful with \href{http://www.ats.ucla.edu/stat/mult_pkg/faq/general/effect.htm}{effects coding}. That is, contrasts should be set to \code{\link{contr.sum}} via \code{options(contrasts=c('contr.sum','contr.poly'))}. This should be done automatically when loading \pkg{afex} and \pkg{afex} will issue a warning when running type 3 SS and \href{http://www.ats.ucla.edu/stat/r/library/contrast_coding.htm}{other coding schemes}. You can check the coding with \code{options("contrasts")}. 
#' 
#' The \code{formula} for \code{aov.car} must contain a single \code{Error} term specyfying the \code{ID} column and potential within-subject factors (you may use \code{\link{mixed}} with multiple error terms). Factors outside the \code{Error} term are treated as between-subject factors (the within-subject factors specified in the \code{Error} term are ignored outside the \code{Error} term, i.e., it is not necessary to specify them outside the \code{Error} term, see Examples).
#'
#' Suppressing the intercept (i.e, via \code{0 +} or \code{- 1}) is ignored. Specific specifications of effects (e.g., excluding terms with \code{-} or using \code{^}) could be okay but is not tested. Using the \code{\link{I}} or \code{\link{poly}} function within the formula is not tested and not supported!
#'
#' For \code{ez.glm} either \code{between} or \code{within} must not be \code{NULL}.
#'
#' \code{ez.glm} will concatante all between-subject factors using \code{*} (i.e., producing all main effects and interactions) and all covariates by \code{+} (i.e., adding only the main effects to the existing between-subject factors). The within-subject factors do fully interact with all between-subject factors and covariates. This is essentially identical to the behavior of SPSS's \code{glm} function.
#'
#' To run an ANCOVA you need to set \code{factorize = FALSE} and make sure that all variables have the correct type (i.e., factors are factors and numeric variables are numeric and centered).
#'
#' Note that the default behavior is to return a \code{\link{nice.anova}} \code{data.frame}. This includes calculation of generalized eta squared for which \strong{all non manipluated (i.e., observed)} variables need to be specified via the \code{observed} argument. Changing the effect size to \code{"pes"} (partial eta-squared) via \code{args.return} or the return value via \code{return} removes this necessity.
#' 
#' If \code{check.contrasts = TRUE}, contrasts will be set to \code{"contr.sum"} for all between-subject factors if default contrasts are not equal to \code{"contr.sum"} or \code{attrib(factor, "contrasts") != "contr.sum"}. (within-subject factors are hard-coded \code{"contr.sum"}.)
#'
#' @author \code{univ} is basically a copy of \code{\link[car]{summary.Anova.mlm}} written by John Fox.\cr The other functions were written by Henrik Singmann.
#'
#' The design of these functions is heavily influenced by \code{\link[ez]{ezANOVA}} from package \pkg{ez}.
#'
#' @note Variables entered as within-subjects (i.e., repeated measures) factors are silently converted to factors and unused levels dropped.
#'
#' Contrasts attached to a factor as an attribute are probably not preserved and not supported.
#'
#' Function \code{univ} was called \code{univariate} in prior versions, but there was a function with similar name in package \pkg{multcomp} leading to bugs and unexpected behavior.
#'
#' @seealso \code{\link{nice.anova}} creats the nice ANOVA tables which are by default returned. See also there for a slightly longer discussion of effect sizes.
#'
#' \code{\link{mixed}} provides a (formula) interface for obtaining p-values for mixed-models via \pkg{lme4}.
#'
#' \code{\link{obk.long}} describes the long version of the \code{OBrienKaiser} dataset used in the examples.
#'
#' @references Maxwell, S. E., & Delaney, H. D. (2004). \emph{Designing Experiments and Analyzing Data: A Model-Comparisons Perspective}. Mahwah, N.J.: Lawrence Erlbaum Associates.
#'
#' @name aov.car
#' @aliases aov.car ez.glm univ
#' @export aov.car ez.glm univ
#' @import car
#' @importFrom stringr str_c str_detect str_replace_all
#' @importFrom reshape2 dcast
#' @example examples/examples.aov.car.R
#'

aov.car <- function(formula, data, fun.aggregate = NULL, type = 3, factorize = TRUE, check.contrasts = TRUE, return = "nice", observed = NULL, args.return = list(), ...) {
  #browser()
  return <- match.arg(return, c("Anova", "lm", "data", "nice", "full", "all", "univariate", "lme4"))
  # stuff copied from aov:
  Terms <- terms(formula, "Error", data = data)
  indError <- attr(Terms, "specials")$Error
  if (length(indError) > 1L) 
    stop(sprintf(ngettext(length(indError), "there are %d Error terms: only 1 is allowed", 
                          "there are %d Error terms: only 1 is allowed"), length(indError)), 
         domain = NA)
  # from here, code by Henrik Singmann:
  vars <- all.vars(formula)
  dv <- vars[1]
  vars <- vars[-1]
  parts <- attr(terms(formula, "Error", data = data), "term.labels")
  error.term <- parts[str_detect(parts, "^Error\\(")]
  id <- all.vars(parse(text = error.term))[1]
  within <- all.vars(parse(text = error.term))[-1]
  between <- vars[!(vars %in% c(id, within))]
  effect.parts <- parts[!str_detect(parts, "^Error\\(")]
  effect.parts.no.within <- effect.parts[!str_detect(effect.parts, str_c("\\<",within,"\\>", collapse = "|"))]
  # factorize if necessary
  if (factorize) {
    if (any(!vapply(data[, between, drop = FALSE], is.factor, TRUE))) {
      to.factor <- between[!vapply(data[,between, drop = FALSE], is.factor, TRUE)]
      message(str_c("Converting to factor: ", str_c(to.factor, collapse = ", ")))
      for (tmp.c in to.factor) {
        data[,tmp.c] <- factor(data[,tmp.c])
      }
    }
  } else {
    # check if numeric variables are centered.
    c.ns <- between[vapply(data[, between, drop = FALSE], is.numeric, TRUE)]
    if (length(c.ns) > 0) {
      non.null <- c.ns[!abs(vapply(data[, c.ns, drop = FALSE], mean, 0)) < .Machine$double.eps ^ 0.5]
      if (length(non.null) > 0) warning(str_c("Numerical variables NOT centered on 0 (i.e., likely bogus results): ", str_c(non.null, collapse = ", ")))
    }
  }  
  # make formulas
  rh2 <- if (length(between) > 0) str_c(effect.parts.no.within, collapse = "+") else "1"
  lh1 <- str_c(id, if (length(between) > 0) str_c(between, collapse = "+") else NULL, sep = "+")
  rh1 <- str_c(within, collapse = "+")
  rh3 <- str_c(within, collapse = "*")
  # converting all within subject factors to factors and adding a leading charcter (x) if starting with a digit.
  new.factor.levels <- c(letters, LETTERS)
  for (within.factor in within) {
    data[,within.factor] <- factor(make.names(as.character(data[,within.factor])))
    #if (length(levels(data[,within.factor])) <= length(new.factor.levels)) levels(data[,within.factor]) <- new.factor.levels[seq_along(levels(data[,within.factor]))]
  }
  # Check if each id is in only one between subjects cell.
  between.factors <- between[vapply(data[, between, drop = FALSE], is.factor, TRUE)]
  if (length(between.factors) > 0) {
    split.data <- split(data, lapply(between.factors, function(x) data[,x]))
    ids.per.condition <- lapply(split.data, function(x) unique(as.character(x[,id])))
    ids.in.more.condition <- unique(unlist(lapply(seq_along(ids.per.condition), function(x) unique(unlist(lapply(ids.per.condition[-x], function(y, z = ids.per.condition[[x]]) intersect(z, y)))))))
    if (length(ids.in.more.condition) > 0) stop(str_c("Following ids are in more than one between subjects condition:\n", str_c(ids.in.more.condition, collapse = ", ")))
  }
  # Is fun.aggregate NULL and aggregation necessary?
  if (is.null(fun.aggregate)) {
    if (any(xtabs(as.formula(str_c("~", id, if (length(within) > 0) "+", rh1)), data = data) > 1)) {
      warning("More than one observation per cell, aggregating the data using mean (i.e, fun.aggregate = mean)!")
      fun.aggregate <- mean
    }
  }
  # Is Type == 3 and contrasts != contr.sum and check.contrasts == FALSE?
  if ((type == 3 | type == "III") & options("contrasts")[[1]][1] != "contr.sum" & !check.contrasts) warning(str_c("Calculating Type 3 sums with contrasts = ", options("contrasts")[[1]][1], "\n  Results likely bogus or not interpretable!\n  You should use check.contrasts = TRUE or options(contrasts=c('contr.sum','contr.poly'))"))
  # if return = "lme4" return the (aggregated) data fitted with lmer!
  if (return == "lme4") {
    warning("lme4 return is experimental!\nAlso: Missing values and contrasts not checked for return = 'lme4'!")
    n.dat <- dcast(data, formula = as.formula(str_c(lh1, if (length(within) > 0) paste0("+", rh1) else "", "~ .", sep = "")), fun.aggregate = fun.aggregate, ..., value.var = dv)
    colnames(n.dat)[length(colnames(n.dat))] <- "value"
    f.within.new <- str_replace_all(rh1, pattern="\\+", replacement="*")
    return(lmer(as.formula(str_c("value~", rh2, if (length(within) > 0) paste0("*", f.within.new) else "", "+ (1", if (length(within) > 0) paste0("+", f.within.new) else "", "|", id, ")" , sep = "")), data = n.dat))
  }
  # prepare the data:
  tmp.dat <- dcast(data, formula = as.formula(str_c(lh1, if (length(within) > 0) rh1 else ".", sep = "~")), fun.aggregate = fun.aggregate, ..., value.var = dv)
  #browser()
  # check for missing values:
  if (any(is.na(tmp.dat))) {
    missing.values <- apply(tmp.dat, 1, function(x) any(is.na(x)))
    warning(str_c("Missing values for following ID(s):\n", str_c(tmp.dat[missing.values,1], collapse = ", "), "\nRemoving those cases from the analysis."))        
  }
  if (length(between) > 0) {
    if (check.contrasts) {
      resetted <- NULL
      for (i in between) {
        if (is.factor(tmp.dat[,i])) {
          if (is.null(attr(tmp.dat[,i], "contrasts")) & (options("contrasts")[[1]][1] != "contr.sum")) {
            contrasts(tmp.dat[,i]) <- "contr.sum"
            resetted  <- c(resetted, i)
          }
          else if (!is.null(attr(tmp.dat[,i], "contrasts")) && attr(tmp.dat[,i], "contrasts") != "contr.sum") {
            contrasts(tmp.dat[,i]) <- "contr.sum"
            resetted  <- c(resetted, i)
          }
        }
      }
    if (!is.null(resetted)) message(str_c("Contrasts set to contr.sum for the following variables: ", str_c(resetted, collapse=", ")))
    }
  }
  data.l <- list(data = tmp.dat)
  if (return == "data") return(tmp.dat)
  # branching based on type of ANOVA
  if (length(within) > 0) {  # if within-subject factors are present:
    # make idata argument
    if (length(within) > 1) {
      within.levels <- lapply(lapply(data[,within], levels), factor)
      idata <- rev(expand.grid(rev(within.levels)))
    } else {
      idata <- data.frame(levels(data[,within]))
      colnames(idata) <- within
    }
    # print(as.formula(str_c("cbind(",str_c(colnames(tmp.dat[-(seq_along(c(id, between)))]), collapse = ", "), ") ~ ", rh2)))
    # browser()
    tmp.lm <- do.call("lm", list(formula = as.formula(str_c("cbind(",str_c(colnames(tmp.dat[-(seq_along(c(id, between)))]), collapse = ", "), ") ~ ", rh2)), data = tmp.dat))
    if (return == "lm") return(tmp.lm)
    Anova.out <- Anova(tmp.lm, idata = idata, idesign = as.formula(str_c("~", rh3)), type = type)
    data.l <- c(data.l, idata = list(idata))
  } else { # if NO within-subjetc factors are present (i.e., purley between ANOVA):
    colnames(tmp.dat)[ncol(tmp.dat)] <- "dv"
    tmp.lm <- do.call("lm", list(formula = as.formula(str_c("dv ~ ", rh2)), data = tmp.dat))
    if (return == "lm") return(tmp.lm)
    Anova.out <- Anova(tmp.lm, type = type)
  }
  if (return == "Anova") return(Anova.out)
  else if ((return == "full")  | (return == "all")) return(c("Anova" = list(Anova.out), "lm" = list(tmp.lm), data.l))
  else if (return == "univariate") return(univ(Anova.out))
  else if (return == "nice") return(do.call("nice.anova", args = c(object = list(Anova.out), observed = list(observed), args.return)))
}



ez.glm <- function(id, dv, data, between = NULL, within = NULL, covariate = NULL, observed = NULL, fun.aggregate = NULL, type = 3, factorize = TRUE, check.contrasts = TRUE, return = "nice", args.return = list(), ..., print.formula = FALSE) {
  if (is.null(between) & is.null(within)) stop("Either between or within need to be non-NULL!")
  if (!is.null(covariate)) covariate <- str_c(covariate, collapse = "+")
  #browser()
  rh <- if (!is.null(between) || !is.null(covariate)) str_c(if (!is.null(between)) str_c(between, collapse = " * ") else NULL, covariate, sep = " + ") else "1"
  error <- str_c(" + Error(", id, if (!is.null(within)) "/" else "", str_c(within, collapse = " * "), ")")
  formula <- str_c(dv, " ~ ", rh, error)
  if (print.formula) message(str_c("Formula send to aov.car: ", formula))
  aov.car(formula = as.formula(formula), data = data, fun.aggregate = fun.aggregate, type = type, return = return, factorize = factorize, check.contrasts = check.contrasts, observed = observed, args.return = args.return, ...)
}


univ <- function(object) { 
	if (all(class(object) == c("anova", "data.frame"))) return(object)
	# This function is basically a cropped copy of car::summary.Anova.mlm written by John Fox returning the output as a list (instead of printing it).
	GG <- function(SSPE, P){ # Greenhouse-Geisser correction
		p <- nrow(SSPE)
		if (p < 2) return(NA) 
		lambda <- eigen(SSPE %*% solve(t(P) %*% P))$values
		lambda <- lambda[lambda > 0]
		((sum(lambda)/p)^2)/(sum(lambda^2)/p)
	}
	HF <- function(gg, error.df, p){ # Huynh-Feldt correction
		((error.df + 1)*p*gg - 2)/(p*(error.df - p*gg))
	}
	mauchly <- function (SSD, P, df) {
		# most of this function borrowed from stats:::mauchly.test.SSD
		if (nrow(SSD) < 2) return(c(NA, NA))
		Tr <- function (X) sum(diag(X))
		p <- nrow(P)
		I <- diag(p)
		Psi <- t(P) %*% I %*% P 
		B <- SSD 
		pp <- nrow(SSD) 
		U <- solve(Psi, B)
		n <- df 
		logW <- log(det(U)) - pp * log(Tr(U/pp))
		rho <- 1 - (2 * pp^2 + pp + 2)/(6 * pp * n)
		w2 <- (pp + 2) * (pp - 1) * (pp - 2) * (2 * pp^3 + 6 * pp^2 + 
					3 * p + 2)/(288 * (n * pp * rho)^2)
		z <- -n * rho * logW
		f <- pp * (pp + 1)/2 - 1
		Pr1 <- pchisq(z, f, lower.tail = FALSE)
		Pr2 <- pchisq(z, f + 4, lower.tail = FALSE)
		pval <- Pr1 + w2 * (Pr2 - Pr1)
		c(statistic = c(W = exp(logW)), p.value = pval)
	}
	nterms <- length(object$terms)
	error.df <- object$error.df
	table <- matrix(0, nterms, 6)
	table2 <- matrix(0, nterms, 4)
	table3 <- matrix(0, nterms, 2)
	rownames(table3) <- rownames(table2) <- rownames(table) <- object$terms
	colnames(table) <- c("SS", "num Df", "Error SS", "den Df", "F", "Pr(>F)")
	colnames(table2) <- c("GG eps", "Pr(>F[GG])",  "HF eps", "Pr(>F[HF])")
	colnames(table3) <- c("Test statistic", "p-value")
	#browser()
	for (term in 1:nterms){
		SSP <- object$SSP[[term]]
		SSPE <- object$SSPE[[term]]
		P <- object$P[[term]]
		p <- ncol(P)
		PtPinv <- solve(t(P) %*% P)
		gg <- GG(SSPE, P)
		table[term, "SS"] <- sum(diag(SSP %*% PtPinv))
		table[term, "Error SS"] <- sum(diag(SSPE %*% PtPinv))
		table[term, "num Df"] <- object$df[term] * p
		table[term, "den Df"] <- error.df * p
		table[term, "F"] <-  (table[term, "SS"]/table[term, "num Df"])/
				(table[term, "Error SS"]/table[term, "den Df"])
		table[term, "Pr(>F)"] <- pf(table[term, "F"], table[term, "num Df"],
				table[term, "den Df"], lower.tail=FALSE)
		table2[term, "GG eps"] <- gg
		table2[term, "HF eps"] <- HF(gg, error.df, p)
		table3[term,] <- mauchly(SSPE, P, object$error.df)
	}
	results <- list(anova = table)
	table3 <- na.omit(table3)
	if (nrow(table3) > 0){
		table2[,"Pr(>F[GG])"] <- pf(table[,"F"], table2[,"GG eps"]*table[,"num Df"],
				table2[,"GG eps"]*table[,"den Df"], lower.tail=FALSE)
		table2[,"Pr(>F[HF])"] <- pf(table[,"F"], 
				pmin(1, table2[,"HF eps"])*table[,"num Df"],
				pmin(1, table2[,"HF eps"])*table[,"den Df"], lower.tail=FALSE)
		table2 <- na.omit(table2)
		if (any(table2[,"HF eps"] > 1)) 
			warning("HF eps > 1 treated as 1")
		attributes(table2)[["na.action"]] <- NULL
		attributes(table3)[["na.action"]] <- NULL
		results <- c(results, mauchly = list(table3), sphericity.correction = list(table2))
	}
	results
}



