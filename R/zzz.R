
.onAttach <- function(libname, pkgname) {
	assign(".oldContrasts", options("contrasts"), envir = .GlobalEnv)
	if (.oldContrasts[[1]][1] != "contr.sum") {
		packageStartupMessage("Setting contrasts to effects coding: options(contrasts=c('contr.sum', 'contr.poly'))\nPrevious contrasts saved in '.oldContrasts'.")
		options(contrasts=c('contr.sum', 'contr.poly'))
	} else packageStartupMessage("Contrasts already set to effects coding: options(contrasts=c('contr.sum', '...'))")
}

#' @export
.Last.lib <- function(libpath) {
	if (exists(".oldContrasts")) {
		if (.oldContrasts[[1]][1] == "contr.sum") message("Unloading afe and removing '.oldContrasts'\nContrasts not restored, was unnecessary.")
		else {
			message("Unloading afe, restoring previous contrasts, and removing '.oldContrasts'.")
			options(.oldContrasts)
		}
		rm(.oldContrasts)
	} else message("Unloading afe but '.oldContrasts' not found.\nContrasts not restored, see 'options('contrasts')'")
}
