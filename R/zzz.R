.onAttach <- function(libname, pkgname) {
	#assign(".oldContrasts", options("contrasts"), envir = .GlobalEnv)
  packageStartupMessage("************\nWelcome to afex created by Henrik Singmann. Important notes:\n")
	if (options("contrasts")[[1]][1] != "contr.sum") {
		packageStartupMessage("Setting contrasts to effects coding: options(contrasts=c('contr.sum', 'contr.poly'))\nThis affects all functions using contrasts (e.g., lmer, lm, aov, ...).\nTo reset default settings run: options(contrasts=c('contr.treatment', 'contr.poly')) (all afex functions should be unaffected by this)\n")
    # \nPrevious contrasts saved in '.oldContrasts'.
		options(contrasts=c('contr.sum', 'contr.poly'))
	} else packageStartupMessage("Contrasts already set to effects coding: options(contrasts=c('contr.sum', '...'))\n")
  packageStartupMessage("afex loads the required packages (e.g., lme4, coin, car, pbkrtest) in an order that does not lead to problems.\nLoading any of the packages (specifically lme4) beforehand usually leads to problems.\nLoading nlme in addition to afex (before or after loading it), also leads to problems.\n************")
}

