\name{HMR}

\alias{HMR}

\title{Trace gas flux estimation with static chamber data}

\description{HMR reads static chamber data from a semicolon/comma separated ASCII text file and analyzes the selected data series
by either user selected (with decision support) models or by user configured automatically selected models. Results are exported
to a semicolon/comma separated ASCII text file.}

\usage{HMR(filename, series = NA, dec = '.', sep = ';', SatPct = NA, SatTimeMin = NA,
pfvar = NA, pfalpha = 0.05, LR.always = FALSE, FollowHMR = FALSE,
IfNoValidHMR = 'No flux', IfNoFlux = 'No flux', IfNoSignal = 'No flux')}

\arguments{
  \item{filename}{The name of an HMR data file. It is assumed that the data file folder has previously been set by the \code{setwd}
  command, eg. \code{setwd('C:/My HMR applications')}, or through the \code{'File|Change dir...'} menu. The required fixed format
  of HMR data files is described in the details section below.}

  \item{series}{Data series in the data file to be analyzed by HMR. The default, \code{series=NA}, requests analysis of all available
  data series, whereas eg. \code{series=c('A','B')} requests analysis of series named 'A' and 'B' in the data file.}

  \item{dec}{The decimal separator used in the data file. Options are \code{dec='.'} and \code{dec=','} and the default is \code{dec='.'}.}

  \item{sep}{The column separator used in the data file. Options are \code{sep=';'} and \code{sep=','} and the default is \code{sep=';'}.}

  \item{SatPct}{The chamber saturation percentage assumed to be met no earlier than \code{SatTimeMin} after deployment; see details
  on flux limiting below. The percentage must be >0 and <100 or \code{NA} (no saturation assumption). The default is \code{SatPct=NA}.}

  \item{SatTimeMin}{The earliest time for chamber saturation \code{SatPct} assumed for for flux limiting; see details below. The
  earliest saturation time must be >0 or \code{NA} (no saturation assumption). The default is \code{SatTimeMin=NA}.}

  \item{pfvar}{The assumed variance of replicate measurements of the ambient trace gas concentration at the chamber site in case
  of no gas emission used for prefiltering of data series before analysis; see details below. The variance must be >0 or \code{NA}
  (no prefiltering). The default is \code{pfvar=NA}.}

  \item{pfalpha}{The statistical significance level used for the prefiltering test; see details below. Data series with prefiltering
  p-value less than \code{pfalpha} are classified as \code{'Noise'}, and \code{pfalpha} must be >0 and <1. The default is \code{pfalpha=0.05}.}

  \item{LR.always}{A logical indicating whether data series should be analysed by linear regression in addition to the selected
  analysis. The default is \code{LR.always=FALSE}.}

  \item{FollowHMR}{A logical indicating whether or not automatic model selection should be applied. This cancels all plotting
  and interactive user selection of analysis. The automatic selection can be configured by options \code{IfNoValidHMR},
  \code{IfNoSignal} and \code{IfNoFlux}, and interacts with flux limiting and prefiltering as described below. The default is
  \code{FollowHMR=FALSE}.}

  \item{IfNoValidHMR}{With automatic model selection, ie. \code{FollowHMR=TRUE}, the analysis specified by \code{IfNoValidHMR} is
  performed when no valid non-linear HMR model can be fitted; see further details below. Options are \code{'LR'} and \code{'No flux'}
  and the default is \code{IfNoValidHMR='No flux'}.}

  \item{IfNoFlux}{With automatic model selection, ie. \code{FollowHMR=TRUE}, the analysis specified by \code{IfNoFlux} is performed
  when the HMR function recommends \code{'No flux'}; see further details below. Options are \code{'LR'} and \code{'No flux'} and the
  default is \code{IfNoFlux='No flux'}.}

  \item{IfNoSignal}{With automatic model selection and prefiltering, ie. \code{FollowHMR=TRUE} and \code{pfvar} are specified,
  the analysis specified by \code{IfNoSignal} is performed when the prefiltering test classifies the data series as \code{'Noise'};
  see further details below. Options are \code{'LR'} and \code{'No flux'} and the default is \code{IfNoSignal='No flux'}.}
}

\details{
\subsection{HMR basic methodology}{
The HMR function implements the methods introduced in Pedersen et al. (2010) supplemented with methods for flux limiting, noise/signal
prefiltering and options for automatic model selection investigated by Pullens et al. (2023).

\eqn{\; \; \; \; \; } The HMR approach is based on the exponential model (\code{'HMR'}):
\deqn{C_t=\varphi+f_0 \frac{\exp(-\kappa t)}{-\kappa h}}
Here \eqn{C_t} denotes the chamber concentration at time \eqn{t}, \eqn{\varphi} denotes the chamber equilibrium concentration, \eqn{f_0}
denotes the initial flux to be estimated, \eqn{\kappa>0} is a non-linearity parameter, and \eqn{h=V/A} denotes the effective chamber height,
ie. the chamber volume divided by the cross-sectional area. The HMR model can be derived from the diffusion model suggested by Hutchinson
and Mosier (1981) (augmented to allow for horizontal gas transport by Pedersen et al. 2010) or from a standard two compartment model
(Seber and Wild 1989, ch. 8) with first order trace gas transport in and out of the chamber compartment (in: emission; out: leakage through
or horizontal transport below chamber walls). It can also simply be seen as a non-linear model with exponential curvature used for fitting
exponentially curved concentration data and with:
\deqn{f_0=h\frac{dC_t}{dt}|_{t=0}}
The HMR approach does, however, insist on parameter values that correspond to positive fitted values of the initial and the equilibrium
concentration, ie. \eqn{C_0} and \eqn{\varphi}. Further details can be found in Pedersen et al. (2010). For small values of \eqn{\kappa},
the HMR model is close to the linear model (\code{'LR'}):
\deqn{C_t=\varphi+f_0\frac{t}{h}}
For large values of \eqn{\kappa}, it is close to the constant model (\code{'No flux'}):
\deqn{C_t=C_0}
The HMR function estimates the model parameters by minimzing the mean squared error (MSE) criterion. If the optimal value of \eqn{\kappa} is
indeterminate by the MSE criterion, ie. the MSE continues to decrease for \eqn{\kappa} approaching either zero or infinity, the HMR function
instead recommends the corresponding limiting model for flux estimation, ie. \code{'LR'} and \code{'No flux'}, respectively. In these limiting
cases, the estimation uncertainty in the estimate of \eqn{\kappa} (zero or infinity) is not reflected in the standard error of the estimated
flux (nor the 95\%-confidence interval or the p-value) implying a discontinuity when going from very small/large estimated values of \eqn{\kappa}
in the non-linear model (\code{'HMR'}) to the limiting linear (\code{'LR'}) or constant model (\code{'No flux'}).
} % Afslutning paa 'HMR basic methodology'

\subsection{HMR flux limiting}{
The family of non-linear HMR models is clearly too large as all values of \eqn{\kappa>0} are permissible and large values of \eqn{\kappa} may
correspond to huge and unrealistic fluxes. With small data series, such excessive fluxes may appear optimal by the MSE criterion when random
patterns in data, erroneously, exhibit strong exponential curvature, so it may be important to limit the permissible values of \eqn{\kappa}.
A priori, \eqn{\kappa} can not be bounded by any fixed number, as it depends on the time unit, the chamber design, local characteristica and
more, but it may be limited by making assumptions about the chamber saturation. The non-linear model can be rewritten like:
\deqn{C_t=C_0+(\varphi-C_0)(1-\exp(-\kappa t))}
Hence, \eqn{1-\exp(-\kappa t)} denotes the saturation fraction at time \eqn{t>0} after deployment. According to the model, complete chamber
is not obtained until \eqn{t=\infty}, but for a given saturation fraction, \eqn{0<p<1}, and assuming that this saturation is not obtained
before time-point \eqn{T} after deployment, an upper limit for \eqn{\kappa} is implied by the equation:
\deqn{\kappa \leq \frac{1}{T} \log(\frac{1}{1-p})}
Even very crude assumptions about the chamber saturation may limit the permissible values of \eqn{\kappa} in an important way. \emph{Example}:
As complete chamber saturation does not happen within finite time in the theoretical model, one may use eg. 99\% saturation in the model as
a proxy for what is considered complete chamber saturation in real measurement situations. Assuming that this saturation will not occur until
after \eqn{T=1} or \eqn{T=2} limits the permissible values of \eqn{\kappa} to \eqn{\kappa\leq 4.6} and \eqn{\kappa\leq 2.3}, respectively.

\eqn{\; \; \; \; \; } Flux limiting with HMR is controlled by options \code{SatPct} and \code{SatTimeMin} that specifies the values of the
saturation fraction (in \%) and the assumed earliest time-point for the saturation percentage, respectively. If flux limiting is active, ie.
the globally optimal value of \eqn{\kappa} according to the MSE criterion is larger than the user specfified maximal value of \eqn{\kappa}
(whereby the MSE criterion recommends \code{'No flux'}), then \code{'LR'} should be considered for estimating the flux as the mathematical
justification for \code{'No flux'} is violated.
} % Afslutning paa 'HMR flux limiting'

\subsection{HMR data prefiltering}{
The MSE criterion is unable to detect if the variation in chamber concentrations is low relative to the variation in data that can be expected
merely from variation in ambient trace gas concentration measurements at the site in case of no trace gas emission. Hence, and particularly with
small samples, the MSE criterion may be fooled by random patterns in data and estimate a, potentially large, flux even if there is none. This is
particularly problematic with the non-linear HMR model. To avoid this, external information is needed about the expected variation in measured
chamber concentrations at the site in case of no flux. This variation depends on the analytical laboratory precision and on the natural variation
of true trace gas concentrations in replicate samples from the site in case of no flux. Although this abstract quantity is presumably unknown,
assessment of its value may prove useful for prefiltering of data series. For instance, the requested variability could be estimated from replicated
samples of ambient air before chamber deployment, although this may lead to an overestimation due to the potential presence of trace gas emission
at the site. Alternatively, one may consider the natural variation of true trace gas concentrations in replicate samples from the site to be
negligible, in which case the requested variability is simply the analytical laboratory precision. This may of course underestimate the variability.
In this way the abstract variability measure may be framed. If the frame is narrow, an intermediate value may be a good assessment of the required
abstract variability, whereas the two framing values can be used for, respectively, conservative (overestimated variability: classifies too few
data series as \code{'Signal'}) or liberal (underestimated variability: classifies too many data series as \code{'Signal'}) prefiltering.

\eqn{\; \; \; \; \; } Given that a reasonable measure of the variation of replicate concentrations measurements at the site without flux can be
provided, this can be used for prefiltering of data series, ie. the classification of series into \code{'Noise'} or \code{'Signal'}. The idea
is then to avoid erroneous flux estimation for data series identified as \code{'Noise'}. The HMR prefiltering test assumes that replicate
concentration measurements at the site without trace gas emission follow a normal distribution with variance \eqn{\sigma^2_0}. The prefiltering
classification is performed by a one-sided statistical test of the null hypothesis of no variation in data in excess of what can be expected from
concentration measurements from the site in case of no flux against the alternative hypothesis of larger variation in data than can be expected
at the site in case of no flux: The prefiltering p-value is given by:
\deqn{p=1-F_{\chi^2(n-1)}((n-1)\frac{s^2}{\sigma^2_0})}
Here \eqn{n} denotes the sample size, \eqn{s^2} denotes the sample variance of measured chamber concentrations, and \eqn{F_{\chi^2(n-1)}} denotes
the cumulative distribution function of the \eqn{\chi^2}-distribution with \eqn{n-1} degrees of freedom.

\eqn{\; \; \; \; \; } With prefiltering at statistical significance level \eqn{\alpha}, the HMR function classifies data series as \code{'Signal'}
if \eqn{p<\alpha} and otherwise as \code{'Noise'}. The level of prefiltering can be controlled by the selected values of \eqn{\sigma^2_0} and
\eqn{\alpha}. With statistical significance level \eqn{\alpha}, \eqn{100\alpha\%} of data series will, erroneously, be classified as \code{'Signal'}
due to statistical type I errors, so increasing \eqn{\alpha} will reduce the noise filtering rate for a given value of \eqn{\sigma^2_0}. Increasing
the assumed variation of concentration measurements at the site without flux, ie. \eqn{\sigma^2_0}, increases the noise filtering rate as more data
series will be classified as \code{'Noise'}.

\eqn{\; \; \; \; \; } Prefiltering with HMR is controlled by options \code{pfvar} and \code{pfalpha} that specify \eqn{\sigma^2_0} and
\eqn{\alpha}, respectively.
} % Afslutning paa 'HMR prefiltering'

\subsection{HMR user selection of analysis}{
To assist the user in selecting the appropriate method for flux estimation for a given data series, the HMR function displays -- organised
in a 2x2 matrix -- plots of the criterion function, the various models fits, decision support from the flux limiting and prefiltering methods,
and a panel with buttons for user selection of method for flux estimation. The upper left plot displays the criterion function over the range
of numerically feasible values of \eqn{\kappa} -- possibly further limited by the user through flux limiting, cf. sections above. Green parts
of the curve represent valid values of \eqn{\kappa}, red parts represent invalid values (cf. above), and the optimal value according to the MSE
criterion is indicated by a blue square. The upper right plot is a zoom into the upper left plot, and the lower left plot shows the fit of the
possible models with the model selected by the MSE criterion displayed in the headline. Moreover, decision support from the prefiltering and
flux limiting methods is written in red text in this plot. The lower right panel contains buttons for user selection of the method for flux
estimation (select by left mouse button click). Pressing the cancel button interrupts the HMR function.
} % HMR user selection of analysis

\subsection{HMR automatic selection of analysis}{
Although users are encouraged to do sequentially user selection of model for flux estimation, HMR does facilitate automatic selection of analysis
(\code{FollowHMR=TRUE}). The automatic selection is based on exactly the same decision support, however, with consistent unsubjective choices
configured by the user through options \code{IfNoValidHMR}, \code{IfNoFlux} and \code{IfNoSignal} (cf. above). Note that \code{IfNoSignal} has
precedence over \code{IfNoValidHMR} and \code{IfNoFlux} when relevant, and that \code{'LR'} is selected irrespective of the user selection for
\code{IfNoFlux} if flux limiting is active, ie. the globally optimal value of \eqn{\kappa} according to the MSE criterion is larger than the
user specfified maximal value of \eqn{\kappa} (whereby the MSE criterion recommends \code{'No flux'}). \emph{The automatic HMR decision tree with
option \code{FollowHMR=TRUE} is illustrated in the figure below (only visible with HTML help, \code{help(HMR,help_type='html')}, and in the PDF
manual).}
\if{html}{\figure{HMRtree.png}{options: width=1280 height=720 alt="Figure: HMRtree.png"}}
\if{latex}{\figure{HMRtree.pdf}{options: width=1280px height=720px}}
} % HMR automatic selection of analysis

\subsection{HMR data processing order}{
Firstly, HMR analyzes if the data series can be fitted by valid HMR models. If not, user selection between \code{'LR'} and \code{'No flux'} is
supported by results from the prefiltering test if selected, and automatic model selection is controlled by the user selected values of
\code{IfNoValidHMR} and \code{IfNoSignal}. If valid HMR models can be fitted, the range of values of \eqn{\kappa} is potentielly further limited
by user selected flux limiting assumptions. User selection of models is then supported by the prefiltering and flux limiting methods, if selected,
whereas automatic model selection is controlled by the user specified values of \code{IfNoFlux} and \code{IfNoSignal}.
} % Afslutning paa 'HMR data processing order'

\subsection{HMR data files}{
HMR data files are semicolon or comma separated files organised in five columns with, respectively, the \emph{data series names}, \emph{chamber
volumes}, \emph{chamber cross-sectional areas}, \emph{observation time-points}, and the \emph{observed chamber concentrations}. Semicolon/comma
separated files can for instance be created and edited by ASCII text editors or spreadsheet software. The following 36 lines shows an HMR data
file containing four data series of static chamber data:

\preformatted{
Series;V;A;Time;Concentration
k0d;140.625;0.5625;0;15.60
k0d;140.625;0.5625;10;15.62
k0d;140.625;0.5625;20;16.53
k0d;140.625;0.5625;30;16.90
k0d;140.625;0.5625;40;17.40
k0d;140.625;0.5625;50;17.69
k0d;140.625;0.5625;60;18.64
k0d;140.625;0.5625;70;18.36
k0d;140.625;0.5625;80;19.14
k0d;140.625;0.5625;110;18.83
k0d;140.625;0.5625;120;19.27
k10d;140.625;0.5625;0;0.3517
k10d;140.625;0.5625;10;0.3523
k10d;140.625;0.5625;20;0.3660
k10d;140.625;0.5625;30;0.3673
k10d;140.625;0.5625;40;0.3603
k10d;140.625;0.5625;50;0.3623
k10d;140.625;0.5625;60;0.3580
k10d;140.625;0.5625;70;0.3650
k10d;140.625;0.5625;80;0.3700
k10d;140.625;0.5625;90;0.3673
k10d;140.625;0.5625;110;0.3647
k10d;140.625;0.5625;120;0.3693
F2T2;2.0101;0.0201;0;10.87
F2T2;2.0101;0.0201;20;19.49
F2T2;2.0101;0.0201;54;24.99
F2T2;2.0101;0.0201;85;27.24
F2T2;2.0101;0.0201;119;33.13
F2T2;2.0101;0.0201;155;30.14
F2V2;2.0101;0.0201;0;9.940
F2V2;2.0101;0.0201;28;31.64
F2V2;2.0101;0.0201;60;48.88
F2V2;2.0101;0.0201;91;58.08
F2V2;2.0101;0.0201;123;76.16
F2V2;2.0101;0.0201;162;106.8
}

Apart from the required header (with optional ASCII character content), the five columns contain:
\describe{
  \item{Column 1}{Text labels that identify data series. Hence, labels must be identical within and different between data series. In the
  sample data above, the first column identifies four data series named \code{k0d}, \code{k10d}, \code{F2T2}, and \code{F2V2}.}
  \item{Column 2}{The chamber volumes, \eqn{V}. Chamber volumes must be identical within data series. In the sample data above, \eqn{V=140.625}
  [\eqn{L}] for data series \code{k0d} and \code{k10d}, and \eqn{V=2.0101} [\eqn{L}] for data series \code{F2T2} and \code{F2V2}.}
  \item{Column 3}{The chamber cross-sectional areas, \eqn{A}. Chamber cross-sectional areas must be identical within data series. In the sample
  data above, \eqn{A=0.5625} [\eqn{m^2}] for data series \code{k0d} and \code{k10d}, and \eqn{A=0.0201} [\eqn{m^2}] for data series \code{F2T2}
  and \code{F2V2}.}
  \item{Column 4}{The measurement time-points within data series in strictly increasing order and with the first time-point equal to zero. At
  least three observation time-points is required per data series. In the sample data above, time-points are in minutes and cover, approximately,
  two-four hour periods per data series.}
  \item{Column 5}{The measured chamber concentrations corresponding to the time-points in the fourth column. Chamber concentations must be
  strictly positive. In the sample data above, the fifth column contains measured nitrous oxide concentrations [\eqn{\mu g/L}].}
}
Missing values (\code{NA}'s) or empty lines or columns are not allowed in HMR data files.
} % Afslutning paa 'HMR data files'

\subsection{HMR data and flux physical units}{
For maximal flexibility, HMR has no requirements for the physical units of input data. The chosen units do, however, determine the unit of
the estimated flux, which has the physical unit of \eqn{(VC)/ (At)}, where \eqn{t} and \eqn{C} denote, respectively, time and concentration.
Some examples:
\tabular{lcl}{
\eqn{V} \eqn{[L]},   \eqn{A} \eqn{[m^2]},  \eqn{t} \eqn{[h]},   \eqn{C} \eqn{[\mu g/L]} \tab \eqn{\Rightarrow} \tab \eqn{f_0} \eqn{[\mu g/m^2/h]} \cr
\eqn{V} \eqn{[L]},   \eqn{A} \eqn{[m^2]},  \eqn{t} \eqn{[min]}, \eqn{C} \eqn{[\mu L/L]} \tab \eqn{\Rightarrow} \tab \eqn{f_0} \eqn{[\mu L/m^2/min]} \cr
\eqn{V} \eqn{[m^2]}, \eqn{A} \eqn{[km^2]}, \eqn{t} \eqn{[s]},   \eqn{C} \eqn{[kg/m^3] } \tab \eqn{\Rightarrow} \tab \eqn{f_0} \eqn{[kg/km^2/s]}
}
} % Afslutning paa 'HMR data and flux units'
} %% Afslutningen paa \details

\value{
HMR first examines the specified function arguments. If problems are detected, HMR halts and issues the comment \code{'Error in input parameters'}. Rules
for the HMR function arguments are provided above in the \bold{Arguments} section. If specified function arguments are valid, HMR continues to examine the
provided data file. If fatal problems are detected, HMR halts and issues the comment \code{'Data file could not be read'}. Typical problems with data files are:
\item{(i)}{Inconsistency between the decimal and column separators used in the data file and the standards for these separators on the computer. The
easiest way to resolve this issue is to specify the separators used in the data file as function arguments to HMR, ie. \code{dec} and \code{sep},
respectively.}
\item{(ii)}{Empty rows below or columns beside valid HMR data are not allowed in HMR data files.}
If specified function arguments are valid, and the data file contains at least one valid data series, HMR continues to analyze data. The output is then
a data frame with one row per analysed data series and variables:
  \item{Series}{Name of the data series.}
  \item{f0}{The estimated flux.}
  \item{f0.se}{The standard error of the estimated flux.}
  \item{f0.p}{The p-value for the null hypothesis of zero flux.}
  \item{f0.lo95}{The lower end-point of the 95\%-confidence interval for the flux.}
  \item{f0.up95}{The upper end-point of the 95\%-confidence interval for the flux.}
  \item{Method}{The method used for estimating the flux (\code{'HMR'}, \code{'LR'}, \code{'No flux'} or \code{'None'}).}
  \item{Warning}{A character string with a warning message in case of estimation problems.}
  \item{Prefilter}{The prefiltering classification (\code{'Signal'}, \code{'Noise'} or \code{'None'}).}
  \item{Prefilter.p}{The prefiltering p-value.}
  \item{SatCrit.Warning}{A character string with a warning message if flux limiting is active.}
  \item{LR.f0}{The flux estimated by linear regression. (Only present if \code{LR.always=TRUE}.)}
  \item{LR.f0.se}{The standard error of the flux estimated by linear regression. (Only present if \code{LR.always=TRUE}.)}
  \item{LR.f0.p}{The p-value for the null hypothesis of zero flux calculated by linear regression. (Only present if \code{LR.always=TRUE}.)}
  \item{LR.f0.lo95}{The lower end-point of the 95\%-confidence interval for the flux calculated by linear regression. (Only present if \code{LR.always=TRUE}.)}
  \item{LR.f0.up95}{The upper end-point of the 95\%-confidence interval for the flux calculated by linear regression. (Only present if \code{LR.always=TRUE}.)}
  \item{LR.Warning}{A character string with a warning message if linear regression estimated a negative predeployment concentration. (Only present if \code{LR.always=TRUE}.)}
The data frame is also exported to a semicolon/comma separated file with the name of the data file preceded by 'HMR - ' and located in the data file folder.
The exported file uses the same column and decimal separators as the data file.
}

\author{Asger R. Pedersen, Ph.D. in statistics, SEGES Innovation, Aarhus, Denmark}

\references{
Hutchinson, G.L. and Mosier, A.R. (1981). Improved soil cover method for field measurement of nitrous oxide fluxes.
\emph{Soil Science Society of America Journal}, \bold{45}, pp. 311-316

Seber, G.A.F. and Wild, C.J. (1989). \emph{Nonlinear regression}. Wiley, New York

Pedersen, A.R., Petersen, S.O. and Schelde, K. (2010): A comprehensive approach to soil-atmosphere trace-gas flux estimation with static chambers.
\emph{European Journal of Soil Science}, \bold{61}, pp. 888-902

Pullens, J.W.M., Abalos, D., Petersen, S.O. and Pedersen, A.R. (2023). Identifying criteria for greenhouse gas flux estimation with automatic and manual chambers: A case study for N2O.
\emph{European Journal of Soil Science}, \bold{74}, e13340. https://doi.org/10.1111/ejss.13340
}

\examples{
\dontrun{
# Suppose the sample data above are located on a Windows machine in the file
# 'C:\My HMR applications\N2O.csv'

# Start by setting the data file folder:
setwd('C:/My HMR applications')

# Notice that R uses '/' in folder declarations, whereas Windows uses '\'.

# Analyse all data series with default settings:
HMR(filename='N2O.csv')

# Produces the following (slightly edited) output when HMR recommendations are followed:
  Series        f0     f0.se      f0.p    f0.lo95   f0.up95
1    k0d 1.872e+01 4.188e+00 2.084e-03  9.062e+00 2.838e+01
2   k10d 1.797e-01 1.385e-01 2.269e-01 -1.337e-01 4.930e-01
3   F2T2 4.509e+01 1.470e+01 5.469e-02 -1.697e+00 9.188e+01
4   F2V2 5.584e+01 3.586e+00 9.926e-05  4.588e+01 6.579e+01
       Method Warning Prefilter Prefilter.p SatCrit.Warning
     1    HMR    None      None          NA              NA
     2    HMR    None      None          NA              NA
     3    HMR    None      None          NA              NA
     4     LR    None      None          NA              NA
# The non-linear 'HMR' analysis was recommended by the MSE criterion for three data series,
# whereas the 'LR' analysis was recommended for the fourth.

# The output was also exported to the semicolon-separated file:
# 'C:\My HMR applications\HMR - N2O.csv'

# Analyse all data series with flux limiting assuming 90\% chamber saturation not before 60 minutes:
HMR(filename='N2O.csv',SatPct=90,SatTimeMin=60)

# Produces the following (slightly edited) output when HMR recommendations are followed:
  Series        f0     f0.se      f0.p    f0.lo95   f0.up95
1    k0d 1.872e+01 4.188e+00 2.084e-03  9.061e+00 2.838e+01
2   k10d 2.689e-02 9.257e-03 1.571e-02  6.261e-03 4.751e-02
3   F2T2 4.509e+01 1.470e+01 5.469e-02 -1.697e+00 9.188e+01
4   F2V2 5.584e+01 3.586e+00 9.926e-05  4.588e+01 6.579e+01
       Method Warning Prefilter Prefilter.p
     1    HMR    None      None          NA
     2     LR    None      None          NA
     3    HMR    None      None          NA
     4     LR    None      None          NA
                                                SatCrit.Warning
          1                                                None
          2 Flux limited by saturation assumption - consider LR
          3                                                None
          4                                                None
# The chamber saturation assumption excluded the MSE optimal value of 'kappa', and HMR
# therefore recommends 'LR'.

# The output was in both HMR analyses above exported to the semicolon-separated file
# named 'C:\My HMR applications\HMR - N2O.csv'. Hence, several analyses of the same
# data file overwrites the output file, so to save a particular output file it has to
# be renamed before the next HMR analysis of the same data file.
}
}

\keyword{htest}
\keyword{regression}
\keyword{nonlinear}
