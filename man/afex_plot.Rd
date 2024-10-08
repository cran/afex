% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/afex_plot.R, R/afex_plot_plotting_functions.R
\name{afex_plot}
\alias{afex_plot}
\alias{afex_plot.afex_aov}
\alias{afex_plot.mixed}
\alias{afex_plot.merMod}
\alias{afex_plot.default}
\alias{interaction_plot}
\alias{oneway_plot}
\title{m-way Plot with Error Bars and Raw Data}
\usage{
afex_plot(object, ...)

\method{afex_plot}{afex_aov}(
  object,
  x,
  trace,
  panel,
  mapping,
  error = "model",
  error_ci = TRUE,
  error_level = 0.95,
  error_arg = list(width = 0),
  data_plot = TRUE,
  data_geom,
  data_alpha = 0.5,
  data_color = "darkgrey",
  data_arg = list(),
  point_arg = list(),
  line_arg = list(),
  emmeans_arg = list(),
  dodge = 0.5,
  return = "plot",
  factor_levels = list(),
  plot_first = NULL,
  legend_title,
  ...
)

\method{afex_plot}{mixed}(
  object,
  x,
  trace,
  panel,
  mapping,
  id,
  error = "model",
  error_ci = TRUE,
  error_level = 0.95,
  error_arg = list(width = 0),
  data_plot = TRUE,
  data_geom,
  data_alpha = 0.5,
  data_color = "darkgrey",
  data_arg = list(),
  point_arg = list(),
  line_arg = list(),
  emmeans_arg = list(),
  dodge = 0.5,
  return = "plot",
  factor_levels = list(),
  plot_first = NULL,
  legend_title,
  ...
)

\method{afex_plot}{merMod}(
  object,
  x,
  trace,
  panel,
  mapping,
  id,
  error = "model",
  error_ci = TRUE,
  error_level = 0.95,
  error_arg = list(width = 0),
  data_plot = TRUE,
  data_geom,
  data_alpha = 0.5,
  data_color = "darkgrey",
  data_arg = list(),
  point_arg = list(),
  line_arg = list(),
  emmeans_arg = list(),
  dodge = 0.5,
  return = "plot",
  factor_levels = list(),
  plot_first = NULL,
  legend_title,
  ...
)

\method{afex_plot}{default}(
  object,
  x,
  trace,
  panel,
  mapping,
  id,
  dv,
  data,
  within_vars,
  between_vars,
  error = "model",
  error_ci = TRUE,
  error_level = 0.95,
  error_arg = list(width = 0),
  data_plot = TRUE,
  data_geom,
  data_alpha = 0.5,
  data_color = "darkgrey",
  data_arg = list(),
  point_arg = list(),
  line_arg = list(),
  emmeans_arg = list(),
  dodge = 0.5,
  return = "plot",
  factor_levels = list(),
  plot_first = NULL,
  legend_title,
  ...
)

interaction_plot(
  means,
  data,
  mapping = c("shape", "lineytpe"),
  error_plot = TRUE,
  error_arg = list(width = 0),
  data_plot = TRUE,
  data_geom = ggplot2::geom_point,
  data_alpha = 0.5,
  data_color = "darkgrey",
  data_arg = list(),
  point_arg = list(),
  line_arg = list(),
  dodge = 0.5,
  plot_first = NULL,
  legend_title,
  col_x = "x",
  col_y = "y",
  col_trace = "trace",
  col_panel = "panel",
  col_lower = "lower",
  col_upper = "upper"
)

oneway_plot(
  means,
  data,
  mapping = "",
  error_plot = TRUE,
  error_arg = list(width = 0),
  data_plot = TRUE,
  data_geom = ggbeeswarm::geom_beeswarm,
  data_alpha = 0.5,
  data_color = "darkgrey",
  data_arg = list(),
  point_arg = list(),
  plot_first = NULL,
  legend_title,
  col_x = "x",
  col_y = "y",
  col_panel = "panel",
  col_lower = "lower",
  col_upper = "upper"
)
}
\arguments{
\item{object}{\code{afex_aov}, \code{mixed}, \code{merMod} or other model
object supported by \pkg{emmeans} (for further examples see:
\code{vignette("afex_plot_supported_models")}).}

\item{...}{currently ignored.}

\item{x}{A \code{character} vector or one-sided \code{formula} specifying the
factor names of the predictors displayed on the x-axis. \code{mapping}
specifies further mappings for these factors if \code{trace} is missing.}

\item{trace}{An optional \code{character} vector or one-sided \code{formula}
specifying the factor names of the predictors connected by the same line.
\code{mapping} specifies further mappings for these factors.}

\item{panel}{An optional \code{character} vector or one-sided \code{formula} 
specifying the factor names of the predictors shown in different panels.}

\item{mapping}{A \code{character} vector specifying which aesthetic mappings 
should be applied to either the \code{trace} factors (if \code{trace} is 
specified) or the \code{x} factors. Useful options are any combination of 
\code{"shape"}, \code{"color"}, \code{"linetype"}, or also \code{"fill"} 
(see examples). The default (i.e., missing) uses \code{c("shape", 
"linetype")} if \code{trace} is specified and \code{""} otherwise (i.e., no
additional aesthetic). If specific mappings should not be applied to
specific graphical elements, one can override those via the corresponding
further arguments. For example, for \code{data_arg} the default is
\code{list(color = "darkgrey")} which prevents that \code{"color"} is
mapped onto points in the background.}

\item{error}{A scalar \code{character} vector specifying on which standard 
error the error bars should be based. Default is \code{"model"}, which
plots model-based standard errors. Further options are: \code{"none"} (or 
\code{NULL}), \code{"mean"}, \code{"within"} (or \code{"CMO"}), and 
\code{"between"}. See details.}

\item{error_ci}{Logical. Should error bars plot confidence intervals
(=\code{TRUE}, the default) or standard errors (=\code{FALSE})?}

\item{error_level}{Numeric value between 0 and 1 determing the width of the
confidence interval. Default is .95 corresponding to a 95\% confidence
interval.}

\item{error_arg}{A \code{list} of further arguments passed to 
\code{\link[ggplot2]{geom_errorbar}}, which draws the errorsbars. Default 
is \code{list(width = 0)} which suppresses the vertical bars at the end of 
the error bar.}

\item{data_plot}{\code{logical}. Should raw data be plotted in the 
background? Default is \code{TRUE}.}

\item{data_geom}{Geom \code{function} or \code{list} of geom functions used
for plotting data in background. The default (missing) uses
\code{\link[ggplot2]{geom_point}} if \code{trace} is specified, otherwise
\code{\link[ggbeeswarm]{geom_beeswarm}} (a good alternative in case of many
data points is \code{ggbeeswarm::geom_quasirandom}) . See examples fo
further options.}

\item{data_alpha}{numeric \code{alpha} value between 0 and 1 passed to
\code{data_geom}. Default is \code{0.5} which correspond to semitransparent
data points in the background such that overlapping data points are plotted
darker. If \code{NULL} it is not passed to \code{data_geom}, and can be set
via \code{data_arg}.}

\item{data_color}{color that should be used for the data in the background.
Default is \code{"darkgrey"}. If \code{NULL} it is not passed to
\code{data_geom}, and can be set via \code{data_arg}. Ignored if
\code{"color"} or \code{"colour"} in \code{mapping}.}

\item{data_arg}{A \code{list} of further arguments passed to
\code{data_geom}. Can also be a \code{list} of \code{list}s, in case
\code{data_geom} is a \code{list} of multiple geoms, which allows having
separate argument lists per \code{data_geom}.}

\item{point_arg, line_arg}{A \code{list} of further arguments passed to 
\code{\link[ggplot2]{geom_point}} or \code{\link[ggplot2]{geom_line}} which
draw the points and lines in the foreground. Default is \code{list()}.
\code{line_arg} is only used if \code{trace} is specified.}

\item{emmeans_arg}{A \code{list} of further arguments passed to 
\code{\link[emmeans]{emmeans}}. Of particular importance for ANOVAs is 
\code{model}, see \code{\link{afex_aov-methods}}.}

\item{dodge}{Numerical amount of dodging of factor-levels on x-axis. Default 
is \code{0.5}.}

\item{return}{A scalar \code{character} specifying what should be returned. 
The default \code{"plot"} returns the \pkg{ggplot2} plot. The other option 
\code{"data"} returns a list with two \code{data.frame}s containing the 
data used for plotting: \code{means} contains the means and standard errors
for the foreground, \code{data} contains the raw data in the background.}

\item{factor_levels}{A \code{list} of new factor levels that should be used
in the plot. The name of each list entry needs to correspond to one of the
factors in the plot. Each list element can optionally be a named character
vector where the name corresponds to the old factor level and the value to
the new factor level. Named vectors allow two things: (1) updating only a
subset of factor levels (if only a subset of levels is specified) and (2)
reordering (and renaming) the factor levels, as order of names within a
list element are the order that will be used for plotting. If specified,
emits a \code{message} with \code{old -> new} factor levels.}

\item{plot_first}{A \code{ggplot2} geom (or a list of geoms) that will be
added to the returned plot as a first element (i.e., before any of the
other graphical elements). Useful for adding reference lines or similar
(e.g., using \code{\link[ggplot2]{geom_hline}}).}

\item{legend_title}{A scalar \code{character} vector with a new title for the
legend.}

\item{id}{An optional \code{character} vector specifying over which variables
the raw data should be aggregated. Only relevant for \code{mixed}, 
\code{merMod}, and \code{default} method. The default (missing) uses all 
random effects grouping factors (for \code{mixed} and \code{merMod} method)
or assumes all data points are independent. This can lead to many data
points. \code{error = "within"} or \code{error = "between"} require that
\code{id} is of length 1. See examples.}

\item{dv}{An optional scalar \code{character} vector giving the name of the
column containing the dependent variable for the \code{afex_plot.default}
method. If missing, the function attempts to take it from the \code{call}
slot of \code{object}. This is also used as y-axis label.}

\item{data}{For the \code{afex_plot.default} method, an optional
\code{data.frame} containing the raw data used for fitting the model and
which will be used as basis for the data points in the background. If
missing, it will be attempted to obtain it from the model via
\code{\link[emmeans]{recover_data}}. For the plotting functions, a
\code{data.frame} with the data that has to be passed and contains the
background data points.}

\item{within_vars, between_vars}{For the \code{afex_plot.default} method, an
optional \code{character} vector specifying which variables should be
treated as within-subjects (or repeated-measures) factors and which as
between-subjects (or independent-samples) factors. If one of the two
arguments is given, all other factors are assumed to fall into the other
category.}

\item{means}{\code{data.frame}s used for plotting of the plotting
functions.}

\item{error_plot}{\code{logical}. Should error bars be plotted? Only used in 
plotting functions. To suppress plotting of error bars use \code{error =
"none"} in \code{afex_plot}.}

\item{col_y, col_x, col_trace, col_panel}{A scalar \code{character} string 
specifying the name of the corresponding column containing the information
used for plotting. Each column needs to exist in both the \code{means} and
the \code{data} \code{data.frame}.}

\item{col_lower, col_upper}{A scalar \code{character} string specifying the 
name of the columns containing lower and upper bounds for the error bars. 
These columns need to exist in \code{means}.}
}
\value{
Returns a \pkg{ggplot2} plot (i.e., object of class \code{c("gg",
  "ggplot")}) unless \code{return = "data"}.
}
\description{
Plots results from factorial experiments. Estimated marginal 
  means and error bars are plotted in the foreground, raw data is plotted in 
  the background. Error bars can be based on different standard errors (e.g.,
  model-based, within-subjects, between-subjects). Functions described here
  return a \pkg{ggplot2} plot object, thus allowing further customization of
  the plot.
  
  \code{afex_plot} is the user friendly function that does data preparation
  and plotting. It also allows to only return the prepared data (\code{return
  = "data"}).
  
  \code{interaction_plot} does the plotting when a \code{trace} factor is
  present. \code{oneway_plot} does the plotting when a \code{trace} factor is
  absent.
}
\details{
\code{afex_plot} obtains the estimated marginal means via 
  \code{\link[emmeans]{emmeans}} and aggregates the raw data to the same 
  level. It then calculates the desired confidence interval or standard error
  (see below) and passes the prepared data to one of the two plotting
  functions: \code{interaction_plot} when \code{trace} is specified and 
  \code{oneway_plot} otherwise.

  \subsection{Error Bars}{Error bars provide a grahical representation of the
  variability of the estimated means and should be routinely added to results
  figures. However, there exist several possibilities which particular 
  measure of variability to use. Because of this, any figure depicting error 
  bars should be accompanied by a note detailing which measure the error bars
  shows. The present functions allow plotting of different types of
  confidence intervals (if \code{error_ci = TRUE}, the default) or standard 
  errors (if \code{error_ci = FALSE}).
  
  A further complication is that readers routinely misinterpret confidence
  intervals. The most common error is to assume that non-overlapping error
  bars indicate a significant difference (e.g., Belia et al., 2005). This is
  often too strong an assumption. (see e.g., Cumming & Finch, 2005; Knol et
  al., 2011; Schenker & Gentleman, 2005). For example, in a fully
  between-subjects design in which the error bars depict 95\% confidence
  intervals and groups are of approximately equal size and have equal
  variance, even error bars that overlap by as much as 50\% still correspond
  to \emph{p} < .05. Error bars that are just touching roughly correspond to
  \emph{p} = .01.
  
  In the case of designs involving repeated-measures factors the usual
  confidence intervals or standard errors (i.e., model-based confidence
  intervals or intervals based on the standard error of the mean) cannot be
  used to gauge significant differences as this requires knowledge about the
  correlation between measures. One popular alternative in the psychological
  literature are intervals based on within-subjects standard
  errors/confidence intervals (e.g., Cousineau & O'Brien, 2014). These
  attempt to control for the correlation across individuals and thereby allow
  judging differences between repeated-measures condition. As a downside,
  when using within-subjects intervals no comparisons across between-subjects
  conditions or with respect to a fixed-value are possible anymore.
  
  In the case of a mixed-design, no single type of error bar is possible that
  allows comparison across all conditions. Likewise, for mixed models
  involving multiple \emph{crossed} random effects, no single set of error
  bars (or even data aggregation) adequately represent the true varibility in
  the data and adequately allows for "inference by eye". Therefore, special
  care is necessary in such cases. One possiblity is to avoid error bars
  altogether and plot only the raw data in the background (with \code{error =
  "none"}). The raw data in the background still provides a visual impression
  of the variability in the data and the precision of the mean estimate, but
  does not as easily suggest an incorrect inferences. Another possibility is
  to use the model-based standard error and note in the figure caption that
  it does not permit comparisons across repeated-measures factors.
  
  The following "rules of eye" (Cumming and Finch, 2005) hold, when permitted
  by design (i.e., within-subjects bars for within-subjects comparisons;
  other variants for between-subjects comparisons), and groups are
  approximately equal in size and variance. Note that for more complex
  designs ususally analyzed with mixed models, such as designs involving
  complicated dependencies across data points, these rules of thumbs may be
  highly misleading.
  \itemize{
    \item  \emph{p} < .05 when the overlap of the 95\% confidence intervals
    (CIs) is no more than about half the average margin of error, that is,
    when proportion overlap is about .50 or less.
    \item \emph{p} < .01 when the two CIs do not overlap, that is, when
    proportion overlap is about 0 or there is a positive gap.
    \item  \emph{p} < .05 when the gap between standard error (SE) bars is at
    least about the size of the average SE, that is, when the proportion gap
    is about 1 or greater.
    \item \emph{p} < .01 when the proportion gap between SE bars is about 2
    or more.
  }   
  }
  \subsection{Implemented Standard Errors}{The following lists the 
  implemented approaches to calculate confidence intervals (CIs) and standard
  errors (SEs). CIs are based on the SEs using the \emph{t}-distribution with
  degrees of freedom based on the cell or group size. For ANOVA models,
  \code{afex_plot} attempts to warn in case the chosen approach is misleading
  given the design (e.g., model-based error bars for purely
  within-subjects plots). For \code{mixed} models, no such warnings are
  produced, but users should be aware that all options beside \code{"model"}
  are not actually appropriate and have only heuristic value. But then again,
  \code{"model"} based error bars do not permit comparisons for factors
  varying within one of the random-effects grouping factors (i.e., factors
  for which random-slopes should be estimated).
  \describe{
    \item{\code{"model"}}{Uses model-based CIs and SEs. For ANOVAs, the
    variant based on the \code{lm} or \code{mlm} model (i.e.,
    \code{emmeans_arg = list(model = "multivariate")}) seems generally
    preferrable.}
    \item{\code{"mean"}}{Calculates the standard error of the mean for
    each cell ignoring any repeated-measures factors.}
    \item{\code{"within"} or \code{"CMO"}}{Calculates within-subjects SEs
    using the Cosineau-Morey-O'Brien (Cousineau & O'Brien, 2014) method. This
    method is based on a double normalization of the data. SEs and CIs are
    then calculated independently for each cell (i.e., if the desired output
    contains between-subjects factors, SEs are calculated for each cell
    including the between-subjects factors).}
    \item{\code{"between"}}{First aggregates the data per participant and 
    then calculates the SEs for each between-subjects condition. Results in 
    one SE and \emph{t}-quantile for all conditions in purely within-subjects
    designs.}
    \item{\code{"none"} or \code{NULL}}{Suppresses calculation of SEs and
    plots no error bars.}
  }
  For \code{mixed} models, the within-subjects/repeated-measures factors are
  relative to the chosen \code{id} effects grouping factor. They are
  automatically detected based on the random-slopes of the random-effects
  grouping factor in \code{id}. All other factors are treated as
  independent-samples or between-subjects factors.
  }
}
\note{
Only the DV/response variable can be called \code{y}, but no
  factor/variable used for plotting.
}
\examples{
# note: use library("ggplot") to avoid "ggplot2::" in the following

##################################################################
##                2-factor Within-Subject Design                ##
##################################################################

data(md_12.1)
aw <- aov_ez("id", "rt", md_12.1, within = c("angle", "noise"))

##---------------------------------------------------------------
##                    Basic Interaction Plots                   -
##---------------------------------------------------------------

## all examples require emmeans and ggplot2:
if (requireNamespace("emmeans") && requireNamespace("ggplot2")) {

afex_plot(aw, x = "angle", trace = "noise") 
# or: afex_plot(aw, x = ~angle, trace = ~noise)

afex_plot(aw, x = "noise", trace = "angle")

### For within-subject designs, using within-subject CIs is better:
afex_plot(aw, x = "angle", trace = "noise", error = "within") 
(p1 <- afex_plot(aw, x = "noise", trace = "angle", error = "within"))

## use different themes for nicer graphs:
p1 + ggplot2::theme_bw()
}
\dontrun{
p1 + ggplot2::theme_light()
p1 + ggplot2::theme_minimal()
p1 + jtools::theme_apa()
p1 + ggpubr::theme_pubr()

### set theme globally for R session:
ggplot2::theme_set(ggplot2::theme_bw())

### There are several ways to deal with overlapping points in the background besides alpha
# Using the default data geom and ggplot2::position_jitterdodge
afex_plot(aw, x = "noise", trace = "angle", error = "within", dodge = 0.3,
          data_arg = list(
            position = 
              ggplot2::position_jitterdodge(
                jitter.width = 0, 
                jitter.height = 5, 
                dodge.width = 0.3  ## needs to be same as dodge
                )))

# Overlapping points are shown as larger points using geom_count
afex_plot(aw, x = "noise", trace = "angle", error = "within", dodge = 0.5,
                data_geom = ggplot2::geom_count)

# Using ggbeeswarm::geom_quasirandom (overlapping points shown in violin shape)
afex_plot(aw, x = "noise", trace = "angle", error = "within", dodge = 0.5,
          data_geom = ggbeeswarm::geom_quasirandom,
          data_arg = list(
            dodge.width = 0.5,  ## needs to be same as dodge
            cex = 0.8, 
            width = 0.05  ## small value ensure data points match means
            ))

# Using ggbeeswarm::geom_beeswarm (overlapping points are adjacent on y-axis)
afex_plot(aw, x = "noise", trace = "angle", error = "within", dodge = 0.5,
          data_geom = ggbeeswarm::geom_beeswarm,
          data_arg = list(
            dodge.width = 0.5,  ## needs to be same as dodge
            cex = 0.8))

# Do not display points, but use a violinplot: ggplot2::geom_violin
afex_plot(aw, x = "noise", trace = "angle", error = "within", 
          data_geom = ggplot2::geom_violin, 
          data_arg = list(width = 0.5))

# violinplots with color: ggplot2::geom_violin
afex_plot(aw, x = "noise", trace = "angle", error = "within", 
          mapping = c("linetype", "shape", "fill"),
          data_geom = ggplot2::geom_violin, 
          data_arg = list(width = 0.5))

# do not display points, but use a boxplot: ggplot2::geom_boxplot
afex_plot(aw, x = "noise", trace = "angle", error = "within", 
          data_geom = ggplot2::geom_boxplot, 
          data_arg = list(width = 0.3))
          
# combine points with boxplot: ggpol::geom_boxjitter
afex_plot(aw, x = "noise", trace = "angle", error = "within", 
          data_geom = ggpol::geom_boxjitter, 
          data_arg = list(width = 0.3))
## hides error bars!

# nicer variant of ggpol::geom_boxjitter
afex_plot(aw, x = "noise", trace = "angle", error = "within", 
          mapping = c("shape", "fill"),
          data_geom = ggpol::geom_boxjitter, 
          data_arg = list(
            width = 0.3, 
            jitter.params = list(width = 0, height = 10),
            outlier.intersect = TRUE),
          point_arg = list(size = 2.5), 
          error_arg = list(linewidth = 1.5, width = 0))

# nicer variant of ggpol::geom_boxjitter without lines
afex_plot(aw, x = "noise", trace = "angle", error = "within", dodge = 0.7,
          mapping = c("shape", "fill"),
          data_geom = ggpol::geom_boxjitter, 
          data_arg = list(
            width = 0.5, 
            jitter.params = list(width = 0, height = 10),
            outlier.intersect = TRUE),
          point_arg = list(size = 2.5), 
          line_arg = list(linetype = 0),
          error_arg = list(linewidth = 1.5, width = 0))

### we can also use multiple geoms for the background by passing a list of geoms
afex_plot(aw, x = "noise", trace = "angle", error = "within", 
          data_geom = list(
            ggplot2::geom_violin, 
            ggplot2::geom_point
          ))

## with separate extra arguments:

afex_plot(aw, x = "noise", trace = "angle", error = "within", 
          dodge = 0.5,
          data_geom = list(
            ggplot2::geom_violin, 
            ggplot2::geom_point
          ), 
          data_arg = list(
            list(width = 0.4),
            list(position = 
                   ggplot2::position_jitterdodge(
                     jitter.width = 0, 
                     jitter.height = 5, 
                     dodge.width = 0.5  ## needs to be same as dodge
                   )))
          )
}


##---------------------------------------------------------------
##                        One-Way Plots                         -
##---------------------------------------------------------------

\dontrun{
afex_plot(aw, x = "angle", error = "within") ## default

## with color we need larger points
afex_plot(aw, x = "angle", mapping = "color", error = "within", 
          point_arg = list(size = 2.5), 
          error_arg = list(linewidth = 1.5, width = 0.05)) 

afex_plot(aw, x = "angle", error = "within", data_geom = ggpol::geom_boxjitter)

## nicer
afex_plot(aw, x = "angle", error = "within", data_geom = ggpol::geom_boxjitter, 
          mapping = "fill", data_alpha = 0.7, 
          data_arg = list(
            width = 0.6, 
            jitter.params = list(width = 0.07, height = 10),
            outlier.intersect = TRUE
          ),
          point_arg = list(size = 2.5), 
          error_arg = list(linewidth = 1.5, width = 0.05))

## we can use multiple geoms with separate argument lists:
afex_plot(aw, x = "angle", error = "within", 
          data_geom = 
            list(ggplot2::geom_violin, ggplot2::geom_boxplot), 
          data_arg = 
            list(list(width = 0.7), list(width = 0.1)))

## we can add a line connecting the means using geom_point(aes(group = 1)):
afex_plot(aw, x = "angle", error = "within") +
  ggplot2::geom_line(ggplot2::aes(group = 1))

## we can also add lines connecting the individual data-point in the bg.
# to deal with overlapping points, we use geom_count and make means larger
afex_plot(aw, x = "angle", error = "within", 
          data_geom = list(ggplot2::geom_count, ggplot2::geom_line), 
          data_arg = list(list(), list(mapping = ggplot2::aes(group = id))), 
          point_arg = list(size = 2.5), 
          error_arg = list(width = 0, linewidth = 1.5)) +
  ggplot2::geom_line(ggplot2::aes(group = 1), linewidth = 1.5)

## One-way plots also supports panels:
afex_plot(aw, x = "angle", panel = "noise", error = "within")

## And panels with lines:
afex_plot(aw, x = "angle", panel = "noise", error = "within") +
  ggplot2::geom_line(ggplot2::aes(group = 1))


## For more complicated plots it is easier to attach ggplot2:
library("ggplot2")

## We can hide geoms by plotting them in transparent colour and add them
## afterward to use a mapping not directly supported. 
## For example, the next plot adds a line to a one-way plot with panels, but 
## with all geoms in the foreground having a colour conditional on the panel.

afex_plot(aw, x = "angle", panel = "noise", error = "within", 
          point_arg = list(color = "transparent"), 
          error_arg = list(color = "transparent")) +
  geom_point(aes(color = panel)) +
  geom_linerange(aes(color = panel, ymin = lower, ymax = upper)) + 
  geom_line(aes(group = 1, color = panel)) +
  guides(color = guide_legend(title = "NOISE"))
## Note that we need to use guides explicitly, otherwise the legend title would 
## be "panel". legend_title does not work in this case.

##---------------------------------------------------------------
##                      Other Basic Options                     -
##---------------------------------------------------------------

## relabel factor levels via factor_levels (with message)
afex_plot(aw, x = "noise", trace = "angle", 
          factor_levels = list(angle = c("0°", "4°", "8°"),
                               noise = c("Absent", "Present")))

## factor_levels allows named vectors which enable reordering the factor levels 
### and renaming subsets of levels:
afex_plot(aw, x = "noise", trace = "angle", 
          factor_levels = list(
            angle = c(X8 = "8°", X4 = "4°", X0 = "0°"),
            noise = c(present = "Present")
          )
)


## Change title of legend
afex_plot(aw, x = "noise", trace = "angle", 
          legend_title = "Noise Condition")

## Add reference line in the background
afex_plot(aw, x = "noise", trace = "angle", 
          plot_first = ggplot2::geom_hline(yintercept = 450, 
                                            colour = "darkgrey"))

## for plots with few factor levels, smaller dodge might be better:
afex_plot(aw, x = "angle", trace = "noise", dodge = 0.25) 

#################################################################
##                    4-factor Mixed Design                    ##
#################################################################

data(obk.long, package = "afex")
a1 <- aov_car(value ~ treatment * gender + Error(id/(phase*hour)), 
              data = obk.long, observed = "gender")

## too difficult to see anything
afex_plot(a1, ~phase*hour, ~treatment) +
  ggplot2::theme_light()

## better
afex_plot(a1, ~hour, ~treatment, ~phase) +
  ggplot2::theme_light()

## even better
afex_plot(a1, ~hour, ~treatment, ~phase, 
          dodge = 0.65, 
          data_arg = list(
            position = 
              ggplot2::position_jitterdodge(
                jitter.width = 0, 
                jitter.height = 0.2, 
                dodge.width = 0.65  ## needs to be same as dodge
                ),
            color = "darkgrey")) +
  ggplot2::theme_classic()

# with color instead of linetype to separate trace factor
afex_plot(a1, ~hour, ~treatment, ~phase, 
          mapping = c("shape", "color"),
          dodge = 0.65, 
          data_arg = list(
            position = 
              ggplot2::position_jitterdodge(
                jitter.width = 0, 
                jitter.height = 0.2, 
                dodge.width = 0.65  ## needs to be same as dodge
                ))) +
  ggplot2::theme_light()

# only color to separate trace factor
afex_plot(a1, ~hour, ~treatment, ~phase, 
          mapping = c("color"),
          dodge = 0.65, 
          data_color = NULL,  ## needs to be set to NULL to avoid error
          data_arg = list(
            position = 
              ggplot2::position_jitterdodge(
                jitter.width = 0, 
                jitter.height = 0.2, 
                dodge.width = 0.65  ## needs to be same as dodge
                ))) +
  ggplot2::theme_classic()


## plot involving all 4 factors:
afex_plot(a1, ~hour, ~treatment, ~gender+phase, 
          dodge = 0.65, 
          data_arg = list(
            position = 
              ggplot2::position_jitterdodge(
                jitter.width = 0, 
                jitter.height = 0.2, 
                dodge.width = 0.65  ## needs to be same as dodge
                ),
            color = "darkgrey")) +
  ggplot2::theme_bw()


##---------------------------------------------------------------
##              Different Standard Errors Available             -
##---------------------------------------------------------------

## purely within-design
cbind(
  afex_plot(a1, ~phase, ~hour, 
            error = "model", return = "data")$means[,c("phase", "hour", "y", "SE")],
  multivariate = afex_plot(a1, ~phase, ~hour,
                           error = "model", return = "data")$means$error,
  mean = afex_plot(a1, ~phase, ~hour, 
                    error = "mean", return = "data")$means$error,
  within = afex_plot(a1, ~phase, ~hour, 
                     error = "within", return = "data")$means$error,
  between = afex_plot(a1, ~phase, ~hour, 
                      error = "between", return = "data")$means$error)
## mixed design
cbind(
  afex_plot(a1, ~phase, ~treatment, 
            error = "model", return = "data")$means[,c("phase", "treatment", "y", "SE")],
  multivariate = afex_plot(a1, ~phase, ~treatment,
                           error = "model", return = "data")$means$error,
  mean = afex_plot(a1, ~phase, ~treatment, 
                    error = "mean", return = "data")$means$error,
  within = afex_plot(a1, ~phase, ~treatment, 
                     error = "within", return = "data")$means$error,
  between = afex_plot(a1, ~phase, ~treatment, 
                      error = "between", return = "data")$means$error)
}

##################################################################
##                         Mixed Models                         ##
##################################################################
if (requireNamespace("MEMSS") && 
    requireNamespace("emmeans") && 
    requireNamespace("ggplot2")) {

data("Machines", package = "MEMSS") 
m1 <- mixed(score ~ Machine + (Machine|Worker), data=Machines)
pairs(emmeans::emmeans(m1, "Machine"))
# contrast   estimate       SE df t.ratio p.value
# A - B     -7.966667 2.420850  5  -3.291  0.0481
# A - C    -13.916667 1.540100  5  -9.036  0.0007
# B - C     -5.950000 2.446475  5  -2.432  0.1253

## Default (i.e., model-based) error bars suggest no difference between Machines.
## This contrasts with pairwise comparisons above.
afex_plot(m1, "Machine")

## Impression from within-subject error bars is more in line with pattern of differences.
afex_plot(m1, "Machine", error = "within")
}

\dontrun{
data("fhch2010") # load 
fhch <- droplevels(fhch2010[ fhch2010$correct,]) # remove errors
### following model should take less than a minute to fit:
mrt <- mixed(log_rt ~ task*stimulus*frequency + (stimulus*frequency||id)+
               (task||item), fhch, method = "S", expand_re = TRUE)

## way too many points in background:
afex_plot(mrt, "stimulus", "frequency", "task") 

## better to restrict plot of data to one random-effects grouping variable
afex_plot(mrt, "stimulus", "frequency", "task", id = "id")
## when plotting data from a single random effect, different error bars are possible:
afex_plot(mrt, "stimulus", "frequency", "task", id = "id", error = "within")
afex_plot(mrt, "stimulus", "frequency", "task", id = "id", error = "mean")

## compare visual impression with:
pairs(emmeans::emmeans(mrt, c("stimulus", "frequency"), by = "task"))

## same logic also possible for other random-effects grouping factor
afex_plot(mrt, "stimulus", "frequency", "task", id = "item")
## within-item error bars are misleading here. task is sole within-items factor.
afex_plot(mrt, "stimulus", "frequency", "task", id = "item", error = "within")
## CIs based on standard error of mean look small, but not unreasonable given results.
afex_plot(mrt, "stimulus", "frequency", "task", id = "item", error = "mean")

### compare distribution of individual data for different random effects:
## requires package cowplot
p_id <- afex_plot(mrt, "stimulus", "frequency", "task", id = "id", 
                  error = "within", dodge = 0.7,
                  data_geom = ggplot2::geom_violin, 
                  mapping = c("shape", "fill"),
                  data_arg = list(width = 0.7)) +
  ggplot2::scale_shape_manual(values = c(4, 17)) +
  ggplot2::labs(title = "ID")

p_item <- afex_plot(mrt, "stimulus", "frequency", "task", id = "item", 
          error = "within", dodge = 0.7,
          data_geom = ggplot2::geom_violin, 
          mapping = c("shape", "fill"),
          data_arg = list(width = 0.7)) +
  ggplot2::scale_shape_manual(values = c(4, 17)) +
  ggplot2::labs(title = "Item")

### see: https://cran.r-project.org/package=cowplot/vignettes/shared_legends.html
p_comb <- cowplot::plot_grid(
  p_id + ggplot2::theme_light() + ggplot2::theme(legend.position="none"),
  p_item + ggplot2::theme_light() + ggplot2::theme(legend.position="none")
  )
legend <- cowplot::get_legend(p_id + ggplot2::theme(legend.position="bottom"))
cowplot::plot_grid(p_comb, legend, 
                   ncol = 1, 
                   rel_heights = c(1, 0.1))

##----------------------------------------------------------------
##                    Support for lme4::lmer                     -
##----------------------------------------------------------------

Oats <- nlme::Oats
## afex_plot does currently not support implicit nesting: (1|Block/Variety)
## Instead, we need to create the factor explicitly
Oats$VarBlock <- Oats$Variety:Oats$Block
Oats.lmer <- lmer(yield ~ Variety * factor(nitro) + (1|VarBlock) + (1|Block),
                        data = Oats)
afex_plot(Oats.lmer, "nitro", "Variety")
afex_plot(Oats.lmer, "nitro", panel = "Variety")

##################################################################
##     Default Method works for Models Supported by emmeans     ##
##################################################################

## lm
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
afex_plot(warp.lm, "tension")
afex_plot(warp.lm, "tension", "wool")

## poisson glm
ins <- data.frame(
    n = c(500, 1200, 100, 400, 500, 300),
    size = factor(rep(1:3,2), labels = c("S","M","L")),
    age = factor(rep(1:2, each = 3)),
    claims = c(42, 37, 1, 101, 73, 14))
ins.glm <- glm(claims ~ size + age + offset(log(n)), 
               data = ins, family = "poisson")
afex_plot(ins.glm, "size", "age")

## binomial glm adapted from ?predict.glm
ldose <- factor(rep(0:5, 2))
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- numdead/20  ## dv should be a vector, no matrix
budworm.lg <- glm(SF ~ sex*ldose, family = binomial, 
                  weights = rep(20, length(numdead)))
afex_plot(budworm.lg, "ldose")
afex_plot(budworm.lg, "ldose", "sex") ## data point is hidden behind mean!
afex_plot(budworm.lg, "ldose", "sex", 
          data_arg = list(size = 4, color = "red"))

## nlme mixed model
data(Oats, package = "nlme")
Oats$nitro <- factor(Oats$nitro)
oats.1 <- nlme::lme(yield ~ nitro * Variety, 
                    random = ~ 1 | Block / Variety,
                    data = Oats)
afex_plot(oats.1, "nitro", "Variety", data = Oats)
afex_plot(oats.1, "nitro", "Variety", data = Oats, id = "Block")
afex_plot(oats.1, "nitro", data = Oats)
afex_plot(oats.1, "nitro", data = Oats, id = c("Block", "Variety"))
afex_plot(oats.1, "nitro", data = Oats, id = "Block")

}
}
\references{
Belia, S., Fidler, F., Williams, J., & Cumming, G. (2005).
  Researchers Misunderstand Confidence Intervals and Standard Error Bars.
  \emph{Psychological Methods}, 10(4), 389-396.
  https://doi.org/10.1037/1082-989X.10.4.389
  
  Cousineau, D., & O'Brien, F. (2014). Error bars in within-subject designs:
  a comment on Baguley (2012). \emph{Behavior Research Methods}, 46(4),
  1149-1151. https://doi.org/10.3758/s13428-013-0441-z
  
  Cumming, G., & Finch, S. (2005). Inference by Eye: Confidence Intervals and
  How to Read Pictures of Data. \emph{American Psychologist}, 60(2), 170-180.
  https://doi.org/10.1037/0003-066X.60.2.170
  
  Knol, M. J., Pestman, W. R., & Grobbee, D. E. (2011). The (mis)use of
  overlap of confidence intervals to assess effect modification.
  \emph{European Journal of Epidemiology}, 26(4), 253-254.
  https://doi.org/10.1007/s10654-011-9563-8
  
  Schenker, N., & Gentleman, J. F. (2001). On Judging the Significance of
  Differences by Examining the Overlap Between Confidence Intervals.
  \emph{The American Statistician}, 55(3), 182-186.
  https://doi.org/10.1198/000313001317097960
}
