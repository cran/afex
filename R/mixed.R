#' Obtain p-values for a mixed-model from lmer().
#'
#' Fits and calculates p-values for all effects in a mixed model fitted with \code{\link[lme4]{lmer}}. The default behavior calculates type 3 like p-values using the Kenward-Rogers approximation for degrees-of-freedom implemented in \code{\link[pbkrtest]{KRmodcomp}} (for LMMs only), but also allows for parametric bootstrap (\code{method = "PB"}) (for LMMs and GLMMs). \code{print}, \code{summary}, and \code{anova} methods for the returned object of class \code{"mixed"} are available (all return the same data.frame).
#'
#' @usage mixed(formula, data, type = 3, method = c("KR", "PB"), args.test = list(), ...)
#'
#' @param formula a formula describing the full mixed-model to be fitted. As this formula is passed to \code{lmer}, it needs at least one random term.
#' @param data data.frame containing the data. Should have all the variables present in \code{fixed}, \code{random}, and \code{dv} as columns.
#' @param type type of test on which effects are based. Only type 3 tests (\code{3} or \code{"III"}) are correctly implemented (see Details).
#' @param method character vector indicating which methods for obtaining p-values should be used. \code{"KR"} (the default) corresponds to the Kenward-Rogers approximation for degrees of freedom (only working with linear mixed models). \code{"PB"} calculates p-values based on parametric bootstrap.
#' @param args.test \code{list} of arguments passed to the function calculating the p-values. See details.
#' @param ... further arguments (such as \code{weights}) passed to \code{\link{lmer}}.
#'
#' @return An object of class \code{"mixed"} (i.e., a list) with the following elements:
#'
#' \enumerate{
#' \item \code{anova.table} a data.frame containing the statistics returned from \code{\link[pbkrtest]{KRmodcomp}}.
#' \item \code{full.model} the \code{"mer"} object returned from fitting the full mixed model.
#' \item \code{restricted.models} a list of \code{"mer"} objects from fitting the restricted models (i.e., each model lacks the corresponding effect)
#' \item \code{tests} a list of objects returned by the function for obtaining the p-values.
#' \item \code{type} The \code{type} argument used when calling this function.
#' \item \code{method} The \code{method} argument used when calling this function.
#' }
#'
#' The following methods exist for objects of class \code{"mixed"}: \code{print} (which uses rounding and only returns the results with \code{F.scaling = 1}), \code{summary}, and \code{anova} (the latter two return the same data.frame).
#'
#' @details Type 3 tests are obtained by comparing a model in which only the corresponding effect is missing with the full model (containing all effects). This corresponds to the (type 3) Wald tests given by \code{car::Anova} for \code{"mer"} models (from version 2.0-13).
#'
#' Type 2 tests are obtained by comparing a model in which the corresponding effect and all higher oder effect (e.g., all three-way interactions for a two-way interaction) are missing with a model in which all effects of the relevant order are present and all higher order effects absent. Consequently, the results for lower order effects are identical of wether or not higher order effects are part of the model or not, which is rather dubious (but \href{https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q3/018992.html}{I didn't find a better way} of implementing the Type 2 tests). This \strong{does not} correspond to the (type 2) Wald Test reported by \code{car::Anova}. If you want type 2 tests, use \code{car::Anova} with \code{test = "F"} (from version 2.0-13) instead of this function.
#'
#' For an introduction to mixed-modeling for experimental designs using p-values see Judd, Westfall, and Kenny (2012). Further introductions to mixed-modeling for experimental designs are given by Baayen and colleagues (Baayen, 2008; Baayen, Davidson & Bates, 2008; Baayen & Milin, 2010). 
#'
#' p-values are calculated via methods from \pkg{pbkrtest}. When \code{method = "KR"}, the Kenward-Rogers approximation for degrees-of-freedom is calculated using \code{\link[pbkrtest]{KRmodcomp}}, which is only applicable to linear-mixed models. The possible argument to \code{pbkrtest} via \code{args.test} is \code{details}.
#'
#' \code{method = "PB"} calculates p-values using parametric bootstrap using \code{\link[pbkrtest]{PBmodcomp}}. This can be used for linear and also generalized linear mixed models (GLMM) by specifiying a \code{\link[stats]{family}} argument to \code{mixed}. Note that you should specify further arguments to \code{PBmodcomp} via \code{args.test}, especially \code{nsim} (the number of simulations to form the reference distribution. For other arguments see \code{\link[pbkrtest]{PBmodcomp}}. Note that \code{REML} (argument to \code{lmer}) will be set to \code{FALSE} if method is \code{PB}.
#'
#' @note This function is not thoroughly tested so please report all bugs to henrik.singmann (at) psychologie.uni-freiburg.de \cr
#' There might be problems with rather big models when constructing the model matrix to fit the \code{lmer} models (potentially problematic with Type 2 tests). If you find any such bug, please send an example including code and data!
#'
#' This functions needs a lot of RAM and rather long time especially with complex random structures. The RAM demand is a problem especially on 32 bit Windows which only supports up to 2 or 3GB RAM (see \href{http://cran.r-project.org/bin/windows/base/rw-FAQ.html}{R Windows FAQ}).
#'
#' This function calls \code{lme4:::nobars} for dealing with the formula. So any significant changes to \pkg{lme4} or \code{lme4:::nobars} may disrupt its functionality.
#'
#' @author Henrik Singmann with contributions from \href{http://stackoverflow.com/q/11335923/289572}{Ben Bolker and Joshua Wiley}.
#'
#' @seealso \code{\link{ez.glm}} and \code{\link{aov.car}} for convenience functions to analyze experimental deisgns with classical ANOVA or ANCOVA wrapping \code{\link[car]{Anova}}. 
#'
#' @references Baayen, R. H. (2008). \emph{Analyzing linguistic data: a practical introduction to statistics using R}. Cambridge, UK; New York: Cambridge University Press.
#'
#' Baayen, R. H., Davidson, D. J., & Bates, D. M. (2008). Mixed-effects modeling with crossed random effects for subjects and items. \emph{Journal of Memory and Language}, 59(4), 390-412. doi:10.1016/j.jml.2007.12.005
#' 
#' Baayen, R. H., & Milin, P. (2010). Analyzing Reaction Times. \emph{International Journal of Psychological Research}, 3(2), 12-28.
#'
#' Judd, C. M., Westfall, J., & Kenny, D. A. (2012). Treating stimuli as a random factor in social psychology: A new and comprehensive solution to a pervasive but largely ignored problem. \emph{Journal of Personality and Social Psychology}, 103(1), 54-69. doi:10.1037/a0028347
#'
#' @export mixed
#' @S3method print mixed
#' @S3method summary mixed
#' @S3method anova mixed
#'
#' @examples
#' \dontrun{
#' # example data from package languageR:
#' # Lexical decision latencies elicited from 21 subjects for 79 English concrete nouns, with variables linked to subject or word. 
#' data(lexdec, package = "languageR")
#' 
#' # using the simplest model
#' m1 <- mixed(RT ~ Correct + Trial + PrevType * meanWeight + 
#'     Frequency + NativeLanguage * Length + (1|Subject) + (1|Word), data = lexdec)
#' 
#' m1
#' # gives:
#' ##                   Effect df1       df2      Fstat p.value
#' ## 1            (Intercept)   1   96.6379 13573.1410  0.0000
#' ## 2                Correct   1 1627.7303     8.1452  0.0044
#' ## 3                  Trial   1 1592.4301     7.5738  0.0060
#' ## 4               PrevType   1 1605.3939     0.1700  0.6802
#' ## 5             meanWeight   1   75.3919    14.8545  0.0002
#' ## 6              Frequency   1   76.0821    56.5348  0.0000
#' ## 7         NativeLanguage   1   27.1213     0.6953  0.4117
#' ## 8                 Length   1   75.8259     8.6959  0.0042
#' ## 9    PrevType:meanWeight   1 1601.1850     6.1823  0.0130
#' ## 10 NativeLanguage:Length   1 1555.4858    14.2445  0.0002
#' 
#' # Fitting a GLMM using parametric bootstrap:
#' require("mlmRev") # for the data, see ?Contraception
#' 
#' gm1 <- mixed(use ~ age + I(age^2) + urban + livch + (1 | district), 
#'  family = binomial, data = Contraception, args.test = list(nsim = 10), method = "PB")
#' }

mixed <- function(formula, data, type = 3, method = c("KR", "PB"), args.test = list(), ...) {
	if ((type == 3 | type == "III") & options("contrasts")[[1]][1] != "contr.sum") warning(str_c("Calculating Type 3 sums with contrasts = ", options("contrasts")[[1]][1], ".\n  Use options(contrasts=c('contr.sum','contr.poly')) instead"))
	# browser()
	# prepare fitting
	mc <- match.call()
	formula.f <- as.formula(formula)
	dv <- as.character(formula.f)[[2]]
	all.terms <- attr(terms(formula.f), "term.labels")
	effect.order <- attr(terms(formula.f), "order")
	effect.order <- effect.order[!grepl("\\|", all.terms)]
	max.effect.order <- max(effect.order)
	random <- str_c(str_c("(", all.terms[grepl("\\|", all.terms)], ")"), collapse = " + ")
	rh2 <- lme4:::nobars(formula.f)
	rh2[[2]] <- NULL
	m.matrix <- model.matrix(rh2, data = data)
	fixed.effects <- attr(terms(rh2, data = data), "term.labels")
	mapping <- attr(m.matrix, "assign")
	# obtain the lmer fits
	#browser() 
	mf <- mc[!names(mc) %in% c("type", "method", "args.test")]
	mf[[1]] <- as.name("lmer")
    mf[["data"]] <- as.name("data")
    if (method[1] == "PB") if ((!"REML" %in% names(mf)) || mf[["REML"]]) {
        warning("REML argument to lmer() set to FALSE for method = 'PB'")
        mf[["REML"]] <- FALSE
    }
	if (type == 3 | type == "III") {
		if (attr(terms(rh2, data = data), "intercept") == 1) fixed.effects <- c("(Intercept)", fixed.effects)
		# prepare lmer call:
		cat(str_c("Fitting ", length(fixed.effects) + 1, " lmer() models:\n["))
        full.model <- eval(mf)
		#full.model <- eval(mf, data, parent.frame())
		cat(".")
		fits <- vector("list", length(fixed.effects))
		for (c in seq_along(fixed.effects)) {
			tmp.columns <- str_c(deparse(-which(mapping == (c-1))), collapse = "")
			mf[[2]] <- as.formula(str_c(dv, "~ 0 + m.matrix[,", tmp.columns, "] +", random))
			fits[[c]] <- eval(mf)
			cat(".")
		}
		cat("]\n")
		names(fits) <- fixed.effects
	} else if (type == 2 | type == "II") {
        warning("Type 2 Methods are incorrectly implemented.\nUse car::Anova instead!")
		cat(str_c("Fitting ", length(fixed.effects) + max.effect.order, " lmer() models:\n["))
		full.model <- vector("list", max.effect.order)
		fits <- vector("list", length(fixed.effects))
		full.model[[length(full.model)]] <- eval(mf)
		cat(".")
		for (c in seq_len(max.effect.order)) {
			if (c == max.effect.order) next 
			tmp.columns <- str_c(deparse(-which(mapping == which(effect.order > c))), collapse = "")
			mf[[2]] <- as.formula(str_c(dv, "~ 0 + m.matrix[,", tmp.columns, "] +", random))
			full.model[[c]] <-  eval(mf)
			cat(".")
		}
		for (c in seq_along(fixed.effects)) {
			order.c <- effect.order[c]
			tmp.columns <- str_c(deparse(-which(mapping == (c) | mapping == if (order.c == max.effect.order) -1 else which(effect.order > order.c))), collapse = "")
			mf[[2]] <- as.formula(str_c(dv, "~ 0 + m.matrix[,", tmp.columns, "] +", random))
			fits[[c]] <- eval(mf)
			cat(".")
		}
		cat("]\n")
		names(fits) <- fixed.effects
	} else stop('Only type 3 and type 2 tests implemented.')
	# obtain p-values:
	cat(str_c("Obtaining ", length(fixed.effects), " p-values:\n["))
    #browser()
	if (method[1] == "KR") {
		tests <- vector("list", length(fixed.effects))
		for (c in seq_along(fixed.effects)) {
			if (type == 3 | type == "III") tests[[c]] <- KRmodcomp(full.model, fits[[c]])
			else if (type == 2 | type == "II") tests[[c]] <- KRmodcomp(full.model[[order.c]], fits[[c]])
			cat(".")
		}
		cat("]\n")
		names(tests) <- fixed.effects
		df.out <- data.frame(Effect = fixed.effects, stringsAsFactors = FALSE)
		df.out <- cbind(df.out, t(vapply(tests, function(x) unlist(x[["test"]][1,]), unlist(tests[[1]][["test"]][1,]))))
		FtestU <- vapply(tests, function(x) unlist(x[["test"]][2,]), unlist(tests[[1]][["test"]][2,]))
		row.names(FtestU) <- str_c(row.names(FtestU), ".U")
		df.out <- cbind(df.out, t(FtestU))
		rownames(df.out) <- NULL
		#browser()
	} else if (method[1] == "PB") {
        tests <- vector("list", length(fixed.effects))
		for (c in seq_along(fixed.effects)) {
			if (type == 3 | type == "III") tests[[c]] <- do.call(PBmodcomp, args = c(largeModel = full.model, smallModel = fits[[c]], args.test))
			else if (type == 2 | type == "II") tests[[c]] <- do.call(PBmodcomp, args = c(largeModel = full.model[[order.c]], smallModel = fits[[c]], args.test))
			cat(".")
        }
        cat("]\n")
        names(tests) <- fixed.effects
        df.out<- data.frame(Effect = fixed.effects, stringsAsFactors = FALSE)
        df.out <- cbind(df.out, t(vapply(tests, function(x) unlist(x[["test"]][2,]), unlist(tests[[1]][["test"]][2,])))[,-2])
        LRT <- vapply(tests, function(x) unlist(x[["test"]][1,]), unlist(tests[[1]][["test"]][1,]))
        row.names(LRT) <- str_c(row.names(LRT), ".LRT")
        df.out <- cbind(df.out, t(LRT))
        rownames(df.out) <- NULL
    } else stop('Only methods "KR" or "PB" currently implemented.')
	#prepare output object
	list.out <- list(anova.table = df.out, full.model = full.model, restricted.models = fits, tests = tests, type = type, method = method[[1]])
	class(list.out) <- "mixed"
	list.out
}

print.mixed <- function(x, ...) {
    if (x[["method"]] == "KR") {
        tmp <- x[[1]][,1:6]
        tmp[,2:6] <- apply(tmp[,2:6], c(1,2), round, digits = 4)
        print(tmp)
    } else if (x[["method"]] == "PB") {
        tmp <- x[[1]][,1:3]
        tmp[,2:3] <- apply(tmp[,2:3], c(1,2), round, digits = 4)
        print(tmp)
    }
    invisible(tmp)
}

summary.mixed <- function(object, ...) object[[1]]

anova.mixed <- function(object, ...) object[[1]]

# is.mixed <- function(x) inherits(x, "mixed")
