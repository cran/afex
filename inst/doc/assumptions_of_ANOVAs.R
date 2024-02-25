## ----echo=FALSE-------------------------------------------------------------------------
req_suggested_packages <- c("see", "performance", "ggplot2")
pcheck <- lapply(req_suggested_packages, requireNamespace, 
                 quietly = TRUE)
if (any(!unlist(pcheck))) {
   message("Required package(s) for this vignette are not available/installed and code will not be executed.")
   knitr::opts_chunk$set(eval = FALSE)
}

## ----set-options, echo=FALSE, cache=FALSE-----------------------------------------------
options(width = 90)
knitr::opts_chunk$set(dpi=72)

## ----message=FALSE----------------------------------------------------------------------
library(afex)
library(performance) # for assumption checks

## ---------------------------------------------------------------------------------------
data(obk.long, package = "afex")

o1 <- aov_ez("id", "value", obk.long, 
             between = c("treatment", "gender"))

check_homogeneity(o1)

## ---------------------------------------------------------------------------------------
data("fhch2010", package = "afex")

a1 <- aov_ez("id", "log_rt", fhch2010,
             between = "task", 
             within = c("density", "frequency", "length", "stimulus"))

## ---------------------------------------------------------------------------------------
check_sphericity(a1)

## ----eval = FALSE-----------------------------------------------------------------------
#  afex_options(
#    correction_aov = "GG", # or "HF"
#    emmeans_model = "multivariate"
#  )

## ---------------------------------------------------------------------------------------
data("stroop", package = "afex")

stroop1 <- subset(stroop, study == 1)
stroop1 <- na.omit(stroop1)

s1 <- aov_ez("pno", "rt", stroop1,
             within = c("condition", "congruency"))

is_norm <- check_normality(s1)

plot(is_norm)

plot(is_norm, type = "qq")

## ---------------------------------------------------------------------------------------
plot(is_norm, type = "qq", detrend = TRUE)

## ---------------------------------------------------------------------------------------
s2 <- aov_ez("pno", "rt", stroop1,
             transformation = "log",
             within = c("condition", "congruency"))

is_norm <- check_normality(s2)

plot(is_norm, type = "qq", detrend = TRUE)

