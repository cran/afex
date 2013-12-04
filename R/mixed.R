#' Obtain p-values for a mixed-model from lmer().
#'
#' Fits and calculates p-values for all effects in a mixed model fitted with \code{\link[lme4]{lmer}}. The default behavior calculates type 3 like p-values using the Kenward-Rogers approximation for degrees-of-freedom implemented in \code{\link[pbkrtest]{KRmodcomp}} (for LMMs only), but also allows for parametric bootstrap (\code{method = "PB"}) (for LMMs and GLMMs). \code{print}, \code{summary}, and \code{anova} methods for the returned object of class \code{"mixed"} are available (all return the same data.frame).
#'
#' @usage mixed(formula, data, type = 3, method = c("KR", "PB", "LRT"), per.parameter = NULL, args.test = list(), check.contrasts = TRUE, progress = TRUE, cl = NULL, ...)
#'
#' @param formula a formula describing the full mixed-model to be fitted. As this formula is passed to \code{lmer}, it needs at least one random term.
#' @param data data.frame containing the data. Should have all the variables present in \code{fixed}, \code{random}, and \code{dv} as columns.
#' @param type type of test on which effects are based. Only type 3 tests (\code{3} or \code{"III"}) are correctly implemented (see Details).
#' @param method character vector indicating which methods for obtaining p-values should be used. \code{"KR"} (the default) corresponds to the Kenward-Rogers approximation for degrees of freedom (only working with linear mixed models). \code{"PB"} calculates p-values based on parametric bootstrap. \code{"LRT"} calculates p-values via the likelihood ratio tests implemented in the \code{anova} method for \code{lmerMod} objects (and is not recommended).
#' @param per.parameter \code{character} vector specifying for which variable tests should be run for each parameter (instead for the overall effect). Can be useful e.g., for testing ordered factors. Relatively untested so results should be compared with a second run without setting this argument. Uses \code{\link{grep}} for selecting parameters among the fixed effects so regular expressions (\code{\link{regex}}) are possible. See Examples.
#' @param args.test \code{list} of arguments passed to the function calculating the p-values. See details.
#' @param check.contrasts \code{logical}. Should contrasts be checked and (if necessary) changed to be \code{"contr.sum"}. See details.
#' @param progress  if \code{TRUE}, shows progress with a text progress bar
#' @param cl  A vector identifying a cluster; used for distributing the estimation of the different models using several cores. See examples. If \code{ckeck.contrasts}, mixed sets the current contrasts (\code{getOption("contrasts")}) at the nodes.
#' @param ... further arguments (such as \code{weights}) passed to \code{\link{lmer}}.
#'
#' @return An object of class \code{"mixed"} (i.e., a list) with the following elements:
#'
#' \enumerate{
#' \item \code{anova.table} a data.frame containing the statistics returned from \code{\link[pbkrtest]{KRmodcomp}}.
#' \item \code{full.model} the \code{"lmerMod"} object returned from fitting the full mixed model.
#' \item \code{restricted.models} a list of \code{"lmerMod"} objects from fitting the restricted models (i.e., each model lacks the corresponding effect)
#' \item \code{tests} a list of objects returned by the function for obtaining the p-values.
#' \item \code{type} The \code{type} argument used when calling this function.
#' \item \code{method} The \code{method} argument used when calling this function.
#' }
#'
#' The following methods exist for objects of class \code{"mixed"}: \code{print} (which uses rounding and invisibly returns the output), \code{summary}, and \code{anova} (the latter two return the same data.frame).
#'
#' @details Type 3 tests are obtained by comparing a model in which only the corresponding effect is missing with the full model (containing all effects). This corresponds to the (type 3) Wald tests given by \code{car::Anova} for \code{"lmerMod"} models.
#'
#' Type 2 tests are obtained by comparing a model in which the corresponding effect and all higher oder effect (e.g., all three-way interactions for a two-way interaction) are missing with a model in which all effects of the relevant order are present and all higher order effects absent. Consequently, the results for lower order effects are identical of wether or not higher order effects are part of the model or not, which is rather dubious (but \href{https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q3/018992.html}{I didn't find a better way} of implementing the Type 2 tests). This \strong{does not} correspond to the (type 2) Wald Test reported by \code{car::Anova}. If you want type 2 Wald tests, use \code{car::Anova} with \code{test = "F"} (from version 2.0-13) instead of this function.
#'
#' For an introduction to mixed-modeling for experimental designs using p-values see Judd, Westfall, and Kenny (2012). Further introductions to mixed-modeling for experimental designs are given by Baayen and colleagues (Baayen, 2008; Baayen, Davidson & Bates, 2008; Baayen & Milin, 2010). Recommendations on how to specify the random effects structure for experimental designs can be found in Barr and colleagues (2013).
#'
#' p-values are per default calculated via methods from \pkg{pbkrtest}. When \code{method = "KR"}, the Kenward-Rogers approximation for degrees-of-freedom is calculated using \code{\link[pbkrtest]{KRmodcomp}}, which is only applicable to linear-mixed models. The possible argument to \code{pbkrtest} via \code{args.test} is \code{details}.
#'
#' \code{method = "PB"} calculates p-values using parametric bootstrap using \code{\link[pbkrtest]{PBmodcomp}}. This can be used for linear and also generalized linear mixed models (GLMM) by specifiying a \code{\link[stats]{family}} argument to \code{mixed}. Note that you should specify further arguments to \code{PBmodcomp} via \code{args.test}, especially \code{nsim} (the number of simulations to form the reference distribution. For other arguments see \code{\link[pbkrtest]{PBmodcomp}}. Note that \code{REML} (argument to \code{lmer}) will be set to \code{FALSE} if method is \code{PB}.
#'
#' \code{method = "LRT"} calculates p-values via likelihood ratio tests implemented in the \code{anova} method for \code{"merMod"} objects. This is recommended by Barr et al. (2013) which did not test the other methods implemented here. Furthermore, this is not recommended by the \href{http://glmm.wikidot.com/faq}{lme4 faq} which instead recommends the other methods implemented here.
#' 
#' If \code{check.contrasts = TRUE}, contrasts will be set to \code{"contr.sum"} for all factors in the formula if default contrasts are not equal to \code{"contr.sum"} or \code{attrib(factor, "contrasts") != "contr.sum"}. Furthermore, the current contrasts (obtained via \code{getOption("contrasts")}) will be set at the cluster nodes if \code{cl} is not \code{NULL}.
#'
#' @note Please report all bugs to henrik.singmann (at) psychologie.uni-freiburg.de \cr
#' There might be problems with rather big models when constructing the model matrix to fit the \code{lmer} models (potentially problematic with Type 2 tests). If you find any such bug, please send an example including code and data!
#'
#' This functions needs a lot of RAM and rather long time especially with complex random structures (when \code{method = "KR"}). The RAM demand is a problem especially on 32 bit Windows which only supports up to 2 or 3GB RAM (see \href{http://cran.r-project.org/bin/windows/base/rw-FAQ.html}{R Windows FAQ}).
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
#' Barr, D. J., Levy, R., Scheepers, C., & Tily, H. J. (2013). Random effects structure for confirmatory hypothesis testing: Keep it maximal. \emph{Journal of Memory and Language}, 68(3), 255-278. doi:10.1016/j.jml.2012.11.001
#'
#' Judd, C. M., Westfall, J., & Kenny, D. A. (2012). Treating stimuli as a random factor in social psychology: A new and comprehensive solution to a pervasive but largely ignored problem. \emph{Journal of Personality and Social Psychology}, 103(1), 54-69. doi:10.1037/a0028347
#'
#' @export mixed
#' @S3method print mixed
#' @S3method summary mixed
#' @S3method anova mixed
#' @import pbkrtest
#' @importFrom lme4 lmer glmer nobars
#' @examples
#' \dontrun{
#' 
#' # use the obk.long data (not reasonable, no random slopes)
#' data(obk.long)
#' mixed(value ~ treatment * phase + (1|id), obk.long)
#'
#' # Examples for using the per.parammeter argument:
#' data(obk.long, package = "afex")
#' obk.long$hour <- ordered(obk.long$hour)
#' 
#' # tests only the main effect parameters of hour individually per parameter.
#' mixed(value ~ treatment*phase*hour +(1|id), per.parameter = "^hour$", data = obk.long)
#' 
#' # tests all parameters including hour individually
#' mixed(value ~ treatment*phase*hour +(1|id), per.parameter = "hour", data = obk.long)
#' 
#' # tests all parameters individually
#' mixed(value ~ treatment*phase*hour +(1|id), per.parameter = ".", data = obk.long)
#'
#' # example data from package languageR:
#' # Lexical decision latencies elicited from 21 subjects for 79 English concrete nouns, 
#' # with variables linked to subject or word. 
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
#'  
#' ## using multicore
#' require(parallel)
#' (nc <- detectCores()) # number of cores
#' cl <- makeCluster(rep("localhost", nc)) # make cluster
#' # to keep track of what the function is doind, redirect output to outfile:
#' # cl <- makeCluster(rep("localhost", nc), outfile = "cl.log.txt")
#' 
#' # obtain fits with multicore:
#' mixed(value ~ treatment*phase*hour +(1|id), data = obk.long, method = "LRT", cl = cl)
#' 
#' }
#' 

mixed <- function(formula, data, type = 3, method = c("KR", "PB", "LRT"), per.parameter = NULL, args.test = list(), check.contrasts = TRUE, progress = TRUE, cl = NULL, ...) {
  if (check.contrasts) {
    #browser()
    vars.to.check <- all.vars(formula)
    resetted <- NULL
    for (i in vars.to.check) {
      if (is.factor(data[,i])) {
        if (is.null(attr(data[,i], "contrasts")) & (options("contrasts")[[1]][1] != "contr.sum")) {
          contrasts(data[,i]) <- "contr.sum"
          resetted  <- c(resetted, i)
        }
        else if (!is.null(attr(data[,i], "contrasts")) && attr(data[,i], "contrasts") != "contr.sum") {
          contrasts(data[,i]) <- "contr.sum"
          resetted  <- c(resetted, i)
        }
      }
    }
    if (!is.null(resetted)) message(str_c("Contrasts set to contr.sum for the following variables: ", str_c(resetted, collapse=", ")))
  }
  #warning(str_c("Calculating Type 3 sums with contrasts = ", options("contrasts")[[1]][1], ".\n  Use options(contrasts=c('contr.sum','contr.poly')) instead"))
  # browser()
  # prepare fitting (i.e., obtain model info)
  mc <- match.call()
  #browser()
  formula.f <- as.formula(formula)
  if (class(formula) != "formula") message("Formula (the first argument) converted to formula.")
  dv <- as.character(formula.f)[[2]]
  all.terms <- attr(terms(formula.f), "term.labels")
  effect.order <- attr(terms(formula.f), "order")
  effect.order <- effect.order[!grepl("\\|", all.terms)]
  max.effect.order <- max(effect.order)
  random <- str_c(str_c("(", all.terms[grepl("\\|", all.terms)], ")"), collapse = " + ")
  rh2 <- nobars(formula.f)
  rh2[[2]] <- NULL
  m.matrix <- model.matrix(rh2, data = data)
  fixed.effects <- attr(terms(rh2, data = data), "term.labels")
  mapping <- attr(m.matrix, "assign")
  fixed.vars <- all.vars(rh2)
  # check for missing values in variables used:
  if (nrow(m.matrix) != nrow(data)) {
    data <- model.frame(as.formula(str_c(vars.to.check[1], "~", str_c(vars.to.check[-1], collapse = "+"))), data = data)
    m.matrix <- model.matrix(rh2, data = data)
    warning(str_c("Due to missing values, reduced number of observations to ", nrow(data)))
  }
  
  # check if numerical variables are centered
  c.ns <- fixed.vars[vapply(data[, fixed.vars, drop = FALSE], is.numeric, TRUE)]
  if (length(c.ns) > 0) {
    non.null <- c.ns[!abs(vapply(data[, c.ns, drop = FALSE], mean, 0)) < .Machine$double.eps ^ 0.5]
    if (length(non.null) > 0) warning(str_c("Numerical variables NOT centered on 0 (i.e., likely bogus results if in interactions): ", str_c(non.null, collapse = ", ")))
  }
  # obtain the lmer fits
  mf <- mc[!names(mc) %in% c("type", "method", "args.test", "progress", "check.contrasts", "per.parameter", "cl")]
  mf[["formula"]] <- formula.f
  if ("family" %in% names(mf)) mf[[1]] <- as.name("glmer")
  else mf[[1]] <- as.name("lmer")
  mf[["data"]] <- as.name("data")
  if ((method[1] %in% c("PB", "LRT")) & !("family" %in% names(mf))) if ((!"REML" %in% names(mf)) || mf[["REML"]]) {
    message("REML argument to lmer() set to FALSE for method = 'PB' or 'LRT'")
    mf[["REML"]] <- FALSE
  }
  #browser()
  ## prepare (g)lmer formulas:
  if (type == 3 | type == "III") {
    if (attr(terms(rh2, data = data), "intercept") == 1) fixed.effects <- c("(Intercept)", fixed.effects)
    #per.parameter <- c("hour", "treatment")
    # The next part alters the mapping of parameters to effects/variables if
    # per.parameter is not NULL (this does the complete magic).
    if (!is.null(per.parameter)) {
      fixed.to.change <- c()
      for (parameter in per.parameter) {
        fixed.to.change <- c(fixed.to.change, grep(parameter, fixed.effects))
      }
      fixed.to.change <- fixed.effects[sort(unique(fixed.to.change))]
      if ("(Intercept)" %in% fixed.to.change) fixed.to.change <- fixed.to.change[-1]
      fixed.all <- dimnames(m.matrix)[[2]]
      #tf2 <- fixed.to.change[2]
      for (tf2 in fixed.to.change) {
        tf <- which(fixed.effects == tf2)
        fixed.lower <- fixed.effects[seq_len(tf-1)]
        fixed.upper <- if (tf < length(fixed.effects)) fixed.effects[(tf+1):length(fixed.effects)] else NULL
        fixed.effects <- c(fixed.lower, fixed.all[which(mapping == (tf-1))], fixed.upper)
        map.to.replace <- which(mapping == (tf-1))
        map.lower <- mapping[seq_len(map.to.replace[1]-1)]
        map.upper <- if (max(map.to.replace) < length(mapping)) mapping[(map.to.replace[length(map.to.replace)]+1):length(mapping)] else NULL
        mapping <- c(map.lower, seq_along(map.to.replace) + map.lower[length(map.lower)], map.upper + length(map.to.replace)-1)
      }
    }
    # make formulas
    formulas <- vector("list", length(fixed.effects) + 1)
    formulas[[1]] <- mf[["formula"]]
    for (i in seq_along(fixed.effects)) {
      tmp.columns <- str_c(deparse(-which(mapping == (i-1))), collapse = "")
      formulas[[i+1]] <- as.formula(str_c(dv, "~ 0 + m.matrix[,", tmp.columns, "] +", random))
    }
    names(formulas) <- fixed.effects
  } else if (type == 2 | type == "II") {
    warning("Implementation of Type 2 method not unproblematic.\n  Check documentation or use car::Anova (Wald tests).")
    if (!is.null(per.parameter)) stop("per.parameter argument only implemented for Type 3 tests.")
    full.model.formulas <- vector("list", max.effect.order)
    submodel.formulas <- vector("list", length(fixed.effects))
    full.model.formulas[[length(full.model.formulas)]] <- mf[["formula"]]
    for (c in seq_len(max.effect.order)) {
      if (c == max.effect.order) next 
      tmp.columns <- str_c(deparse(-which(mapping %in% which(effect.order > c))), collapse = "")
      full.model.formulas[[c]] <-  as.formula(str_c(dv, "~ 0 + m.matrix[,", tmp.columns, "] +", random))
    }
    for (c in seq_along(fixed.effects)) {
      order.c <- effect.order[c]
      tmp.columns <- str_c(deparse(-which(mapping == (c) | mapping %in% if (order.c == max.effect.order) -1 else which(effect.order > order.c))), collapse = "")
      submodel.formulas[[c]] <- as.formula(str_c(dv, "~ 0 + m.matrix[,", tmp.columns, "] +", random))
    }
    formulas <- c(full.model.formulas, submodel.formulas)
  } else stop('Only type 3 and type 2 tests implemented.')
  ## fit models
  # single core
  if (is.null(cl)) {
    if (progress) cat(str_c("Fitting ", length(formulas), " (g)lmer() models:\n["))
    fits <- vector("list", length(formulas))
    for (i in seq_along(formulas)) {
      mf[["formula"]] <- formulas[[i]]
      fits[[i]] <- eval(mf)
      if (progress) cat(".")
    }
    if (progress) cat("]\n")
  } else {  # multicore
    eval.cl <- function(formula, m.call, progress) {
      m.call[[2]] <- formula
      res <- eval(m.call)
      if (progress) cat(".")
      return(res)
    }
    if (progress) cat(paste0("Fitting ", length(formulas), " (g)lmer() models.\n"))
    #junk <- clusterEvalQ(cl = cl, library("lme4", character.only = TRUE))
    junk <- clusterEvalQ(cl = cl, loadNamespace("lme4"))
    if (check.contrasts)  {
      curr.contrasts <- getOption("contrasts")
      clusterExport(cl = cl, "curr.contrasts", envir = sys.nframe())
      junk <- clusterEvalQ(cl = cl, options(contrasts=curr.contrasts))
    }
    if (progress) junk <- clusterEvalQ(cl = cl, cat("["))
    fits <- clusterApplyLB(cl = cl, x = formulas, eval.cl, m.call = mf, progress = progress)
    if (progress) junk <- clusterEvalQ(cl = cl, cat("]"))
  }
  ## prepare for p-values:
  if (type == 3 | type == "III") {
    full.model <- fits[[1]]
    fits <- fits[-1]
  } else if (type == 2 | type == "II") {
    full.model <- fits[1:max.effect.order]
    fits <- fits[(max.effect.order+1):length(fits)]
  }
  names(fits) <- fixed.effects
  # obtain p-values:
  #browser()
  if (method[1] == "KR") {
    if (progress) cat(str_c("Obtaining ", length(fixed.effects), " p-values:\n["))
    tests <- vector("list", length(fixed.effects))
    for (c in seq_along(fixed.effects)) {
      if (type == 3 | type == "III") tests[[c]] <- KRmodcomp(full.model, fits[[c]])
      else if (type == 2 | type == "II") {
        order.c <- effect.order[c]
        tests[[c]] <- KRmodcomp(full.model[[order.c]], fits[[c]])
      }
      if (progress) cat(".")
    }
    if (progress) cat("]\n")
    names(tests) <- fixed.effects
    df.out <- data.frame(Effect = fixed.effects, stringsAsFactors = FALSE)
    df.out <- cbind(df.out, t(vapply(tests, function(x) unlist(x[["test"]][1,]), unlist(tests[[1]][["test"]][1,]))))
    FtestU <- vapply(tests, function(x) unlist(x[["test"]][2,]), unlist(tests[[1]][["test"]][2,]))
    row.names(FtestU) <- str_c(row.names(FtestU), ".U")
    df.out <- cbind(df.out, t(FtestU))
    rownames(df.out) <- NULL
    #browser()
  } else if (method[1] == "PB") {
    if (progress) cat(str_c("Obtaining ", length(fixed.effects), " p-values:\n["))
    tests <- vector("list", length(fixed.effects))
    for (c in seq_along(fixed.effects)) {
      if (type == 3 | type == "III") tests[[c]] <- do.call(PBmodcomp, args = c(largeModel = full.model, smallModel = fits[[c]], args.test))
      else if (type == 2 | type == "II") {
        order.c <- effect.order[c]
        tests[[c]] <- do.call(PBmodcomp, args = c(largeModel = full.model[[order.c]], smallModel = fits[[c]], args.test))
      }
      if (progress) cat(".")
    }
    if (progress) cat("]\n")
    names(tests) <- fixed.effects
    df.out<- data.frame(Effect = fixed.effects, stringsAsFactors = FALSE)
    df.out <- cbind(df.out, t(vapply(tests, function(x) unlist(x[["test"]][2,]), unlist(tests[[1]][["test"]][2,])))[,-2])
    LRT <- vapply(tests, function(x) unlist(x[["test"]][1,]), unlist(tests[[1]][["test"]][1,]))
    row.names(LRT) <- str_c(row.names(LRT), ".LRT")
    df.out <- cbind(df.out, t(LRT))
    rownames(df.out) <- NULL
  } else if (method[1] == "LRT") {
    tests <- vector("list", length(fixed.effects))
    for (c in seq_along(fixed.effects)) {
      if (type == 3 | type == "III") tests[[c]] <- anova(full.model, fits[[c]])
      else if (type == 2 | type == "II") {
        order.c <- effect.order[c]
        tmpModel  <- full.model[[order.c]] 
        tests[[c]] <- anova(tmpModel, fits[[c]])
      }
    }
    names(tests) <- fixed.effects
    df.large  <- vapply(tests, function(x) x[["Df"]][2], 0)
    df.small  <- vapply(tests, function(x) x[["Df"]][1], 0)
    chisq  <- vapply(tests, function(x) x[["Chisq"]][2], 0)
    df  <- vapply(tests, function(x) x[["Chi Df"]][2], 0)
    p.value  <- vapply(tests, function(x) x[["Pr(>Chisq)"]][2], 0)
    df.out <- data.frame(Effect = fixed.effects, df.large, df.small, chisq, df, p.value, stringsAsFactors = FALSE)
    rownames(df.out) <- NULL
  } else stop('Only methods "KR", "PB" or "LRT" currently implemented.')
  #prepare output object
  list.out <- list(anova.table = df.out, full.model = full.model, restricted.models = fits, tests = tests, type = type, method = method[[1]])
  class(list.out) <- "mixed"
  list.out
}

round.ps <- function(x) {
  substr(as.character(ifelse(x < 0.0001, " <.0001", ifelse(x < 0.001, formatC(x, digits = 4, format = "f"), ifelse(x < 0.01, formatC(x, digits = 3, format = "f"), ifelse(round(x, 2) == 1, " >.99", formatC(x, digits = 2, format = "f")))))), 2, 7)
}


print.mixed <- function(x, ...) {
  if (x[["method"]] == "KR") {
    tmp <- x[[1]][,1:6]
    tmp[,"stat"] <- formatC(tmp[,"stat"], format = "f", digits = 2)
    tmp[,"ddf"] <- prettyNum(tmp[,"ddf"], digits = 2)
    tmp[,"F.scaling"] <- prettyNum(tmp[,"F.scaling"], digits = 2)
    
  } else if (x[["method"]] == "PB") {
    tmp <- x[[1]][,1:3]
    tmp[,2] <- formatC(tmp[,2], format = "f", digits = 2)
  } else if (x[["method"]] == "LRT") {
    tmp <- x[[1]]
    tmp[,"chisq"] <- formatC(tmp[,"chisq"], format = "f", digits = 2)
  }
  tmp[,"p.value"] <- round.ps(tmp[,"p.value"])
  warnings <- lapply(x[[3]], function(y) y@optinfo$warnings)
  warn <- vapply(warnings, function(y)  !length(y)==0, NA)
  if (any(warn)) warning("At least the following warnings were obtained when fitting via lme4:\n", paste(paste(names(which(warn)), vapply(warnings[warn], function(x) x[[1]][1], ""), sep = ": "), collapse = "\n"))
  print(tmp)
  invisible(tmp)
}


summary.mixed <- function(object, ...) object[[1]]

anova.mixed <- function(object, ...) object[[1]]

# is.mixed <- function(x) inherits(x, "mixed")
