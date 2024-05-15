# Time Series Methods for the Analysis of Soundscapes and Other Cyclical Ecological Data

In this notebook, we provide code to produce the figures for the title journal article, published to Methods in Ecology and Evolution (2024).

## Preliminaries

Begin by loading packages and setting up directories in which to find the data and save the figures. \### Packages

``` r
## For spectrum analysis
library(astsa)

## For ordinary least squares 
library(pracma)

## To deal with the dates
library(lubridate)

## For minor ticks on base R plots
library(Hmisc)

## PCA requires multitaper
library(multitaper)

## For better plots
library(ggplot2)
library(dplyr)
library(patchwork) # To display 2 charts together
library(hrbrthemes)
extrafont::loadfonts()

## Data manipulation
library(reshape)

library(forcats)

## To export data.frames to LaTeX-formatted tables
library(xtable)
```

### Plot settings

``` r
## Use this to make ggplot2 look like base
theme_set(theme_bw())
theme_update(text = element_text(size=12),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
strip.background = element_blank()
)

## Set plotting window size default
options(jupyter.plot_scale=1)
```

## Import the data

The data represents soundscape sampling data collected from eight tropical forest sites across the Ogooué-Ivindo province of Gabon. Data was collected during the dry season at four sites (July 17th and July 23rd 2021) and during the wet season at four additional sites (February 19th and March 2nd 2021). Refer to the paper for further information attaining to the methods.

The soundscape has been characterized using Power Minus Noise (PMN; **column pmn\$TotPMN**). PMN was calculated seperately for 256 frequency bins between 0-11 kHz (43 Hz bandwidth each) and summed. PMN was calculated seperately for each minute of the day, yielding 1440 data points per day per site.

**pmn\$DatetimeFinal** is the time converted to local time.

### Specify your working directories as appropriate and import

``` r
# Load the data 
load(paste(pwd(),"/pmn.Rdata",sep=""))
colnames(pmn)

# Specify directory to save the figures (optional)
figdir = pwd()

# Change from UTC to local time, that is, West Africa time (Gabon)
pmn$DatetimeFinal = pmn$DatetimeFinal - hours(6)
```

The data is saved under two seasons "dry" and "rainy", we rename them here for clarity

``` r
# Rename factor levels for plotting
pmn$Season <-as.factor(pmn$Season)
levels(pmn$Season)
levels(pmn$Season) <-c("Dry", "Wet")
levels(pmn$Season)
```

```{=html}
<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
```
<ol class="list-inline">

<li>'dry'</li>

<li>'rainy'</li>

</ol>

```{=html}
<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
```
<ol class="list-inline">

<li>'Dry'</li>

<li>'Wet'</li>

</ol>

## Visualize the raw data

``` r
# Get the unique names of the sites 
usite = paste(unique(pmn$Site))
nsites = length(usite)
print("Names of the sites: ")
print(usite)
```

    [1] "Names of the sites: "
    [1] "IPA10ST"  "IPA11ST"  "IPA14ST"  "IPA15ST"  "TBNI21ST" "TBNI22ST" "TBNI24ST"
    [8] "TBNI25ST"

``` r
# Head the data
print("First few rows of the dataframe: ")
head(pmn)
```

    [1] "First few rows of the dataframe: "

+---+---------+---------------------+----------+---------+
|   | Site    | DatetimeFinal       | TotPMN   | Season  |
+===+=========+=====================+==========+=========+
|   | \<fct\> | \<dttm\>            | \<dbl\>  | \<fct\> |
+---+---------+---------------------+----------+---------+
| 1 | IPA10ST | 2021-02-20 00:00:00 | 746.4980 | Wet     |
+---+---------+---------------------+----------+---------+
| 2 | IPA10ST | 2021-02-20 00:01:00 | 554.3752 | Wet     |
+---+---------+---------------------+----------+---------+
| 3 | IPA10ST | 2021-02-20 00:02:00 | 717.8012 | Wet     |
+---+---------+---------------------+----------+---------+
| 4 | IPA10ST | 2021-02-20 00:03:00 | 717.0829 | Wet     |
+---+---------+---------------------+----------+---------+
| 5 | IPA10ST | 2021-02-20 00:04:00 | 641.4947 | Wet     |
+---+---------+---------------------+----------+---------+
| 6 | IPA10ST | 2021-02-20 00:05:00 | 655.8764 | Wet     |
+---+---------+---------------------+----------+---------+

: A data.frame: 6 × 4

``` r
# View summary statistics
print("Summary statistics: ")
summary(pmn)
```

    [1] "Summary statistics: "



           Site       DatetimeFinal                     TotPMN       Season     
     IPA14ST :11520   Min.   :2021-02-20 00:00:00   Min.   :  85.8   Dry:37440  
     IPA15ST :11520   1st Qu.:2021-02-25 11:59:45   1st Qu.: 680.8   Wet:34560  
     TBNI24ST:10080   Median :2021-07-17 04:59:30   Median : 861.0              
     TBNI25ST:10080   Mean   :2021-05-11 18:11:30   Mean   : 892.9              
     TBNI21ST: 8640   3rd Qu.:2021-07-20 07:59:15   3rd Qu.:1056.3              
     TBNI22ST: 8640   Max.   :2021-07-23 22:59:00   Max.   :7123.6              
     (Other) :11520                                                             

## Visualise when each survey occurred and the data timespans

We will make a figure that contains the beginning and end times of each site's data and shows the season. To do this, we first need to calculate the duration of each time series.

``` r
# Helper function to get the first and last elements of an array
firstlast <- function(x) { c(x[1], tail(x, n = 1)) }

# The ind matrix will hold the first and last indices of each site
ind = matrix(0,length(usite),2)
for (i in 1:length(usite)){
  ind[i,] = firstlast(which(pmn$Site == usite[i]))
}

# Create vectors
begin_time = pmn$DatetimeFinal[ind[,1]]
end_time   = pmn$DatetimeFinal[ind[,2]]
duration_mins = difftime(end_time, begin_time, units="mins")
```

Create a key which specifies whether each site represents dry or wet season data.

``` r
key = data.frame(Site=usite, Season = pmn$Season[ind[,1]])
key
```

+----------+---------+
| Site     | Season  |
+==========+=========+
| \<chr\>  | \<fct\> |
+----------+---------+
| IPA10ST  | Wet     |
+----------+---------+
| IPA11ST  | Wet     |
+----------+---------+
| IPA14ST  | Wet     |
+----------+---------+
| IPA15ST  | Wet     |
+----------+---------+
| TBNI21ST | Dry     |
+----------+---------+
| TBNI22ST | Dry     |
+----------+---------+
| TBNI24ST | Dry     |
+----------+---------+
| TBNI25ST | Dry     |
+----------+---------+

: A data.frame: 8 × 2

Create the data frame used for plotting

``` r
# Make a data frame, tabl, that contains the beginning and ending time data for each site, which is suitable for ggplot2
tabl = data.frame(Site=rep(usite,each=2),beginend_time=pmn$DatetimeFinal[t(ind[,c(1,2)])])
tabl <- arrange(mutate(tabl,Site=factor(Site,levels=usite)),Site)
tabl

# Date tick information
datetick = pretty_dates(c(min(begin_time),add_with_rollback(max(end_time), months(2))), 6)
datetick2 = c(-Inf, ymd_hms("2021-05-01 00:00:01"), Inf)
ytick = c(-Inf,Inf,-Inf,Inf,-Inf)

poly_x = c(1,1,2,2)
poly_y = c(1,2,2,1)

# Labels
ids <- factor(c("Wet", "Dry","Wet", "Dry"))

# The second data frame shows the polygons we want to shade in dry and wet colors
positions <- data.frame(
  id = rep(ids, each = 2),
  beginend_time = as_datetime(c(cbind(datetick2[1:2][poly_x],datetick2[2:3][poly_x]))),
  Site = c(cbind(ytick[1:2][poly_y],ytick[2:3][poly_y]))
)
```

+----------+---------------------+
| Site     | beginend_time       |
+==========+=====================+
| \<fct\>  | \<dttm\>            |
+----------+---------------------+
| IPA10ST  | 2021-02-20 00:00:00 |
+----------+---------------------+
| IPA10ST  | 2021-02-23 23:59:00 |
+----------+---------------------+
| IPA11ST  | 2021-02-23 00:00:00 |
+----------+---------------------+
| IPA11ST  | 2021-02-26 23:59:00 |
+----------+---------------------+
| IPA14ST  | 2021-02-23 00:00:00 |
+----------+---------------------+
| IPA14ST  | 2021-03-02 23:59:00 |
+----------+---------------------+
| IPA15ST  | 2021-02-22 00:00:00 |
+----------+---------------------+
| IPA15ST  | 2021-03-01 23:59:00 |
+----------+---------------------+
| TBNI21ST | 2021-07-16 23:00:00 |
+----------+---------------------+
| TBNI21ST | 2021-07-22 22:59:00 |
+----------+---------------------+
| TBNI22ST | 2021-07-16 23:00:00 |
+----------+---------------------+
| TBNI22ST | 2021-07-22 22:59:00 |
+----------+---------------------+
| TBNI24ST | 2021-07-16 23:00:00 |
+----------+---------------------+
| TBNI24ST | 2021-07-23 22:59:00 |
+----------+---------------------+
| TBNI25ST | 2021-07-16 23:00:00 |
+----------+---------------------+
| TBNI25ST | 2021-07-23 22:59:00 |
+----------+---------------------+

: A data.frame: 16 × 2

Generate figure (export optional)

``` r
q = ggplot(data = positions, aes(x=beginend_time,y=Site)) + 
     # shade the seasons
     geom_polygon(data = positions, aes(fill = id, group = id, alpha=0.01), show.legend = FALSE) +
     # line segments representing the data spans in time
     geom_line(data = tabl, aes(colour = Site), lwd=1.2, show.legend = FALSE) +
     # label the site names beside the line segments in the same colors
     geom_text(data = data.frame(beginend_time=pmn$DatetimeFinal[ind[,2]] + days(1), Site = usite), 
                aes(label=Site,hjust=0, colour = Site), show.legend = FALSE, check_overlap = TRUE) +
     # Remove y-axis labels, set x-axis labels (dates) at 45 degrees, and add a few other beautifications
     theme(axis.text.y = element_blank(), axis.text.x = element_text(angle=45, hjust = 1)) +
     scale_x_continuous(limits=c(ymd_hms("2021-02-01 00:00:01"), datetick2[2]+months(4)),breaks=datetick) +
     labs(
        x = "Date",  
        y = "",
        colour = "Site Type",
        fill = "Season"
      ) +
      scale_fill_manual(values=c("beige", "lightblue")) +
      guides(alpha = "none")

#ggsave(paste(figdir,"Datestimes_large_sitenameson.png",sep=""), plot=q, width=7, height=7)

q
```

![png](output_20_0.png)

## Plot raw PMN for each site together

Visualise how PMN varies across and between days for each site. I've created a ggplot2 function so that it is easier to do this for all of the variables in the pmn dataset, not just TotPMN column, which I have done here. Line 15 in the below cell is all that's required in subsequent cells.

``` r
ctfun <- function(data, var) {
  ggplot(data, aes(x=DatetimeFinal, y = {{ var }})) +
    geom_rect(aes(fill = Season),xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf) +
  geom_line(aes(colour = factor(Site))) +
  labs(
    x = "Date",  
    colour = "Site Name",
    fill="Season",
    y = "Power Minus Noise"
  ) +
    scale_fill_manual(values=c("beige", "lightblue")) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  facet_wrap(vars(Site), nrow = 2, scales = "free_x")
}

pmn |> ctfun(.data[["TotPMN"]]) 
```

![png](output_22_0.png)

## Boxplot figure showing data at different granularities

This plot demonstrates how binning time can impact our interpretation of the data. In the first plot, we lump together all of the time series for the dry and wet seasons. The second plot shows the dry and rainy seasons separated into night and day. The third plot breaks the data down further by hour, the fourth by fifteen minute. And finally, the fifth shows a scatterplot of the raw data by minute.

### Plot 1: Dry VS Rainy

``` r
# Specify the lower limit for y-axis, upper limit for y-axis
lly = 300
uly = 3000

# Generate plot
p1 = ggplot(subset(pmn, Site %in% c("IPA10ST","TBNI21ST")), aes(x=Season, y = TotPMN)) + 
    geom_boxplot(col=c("black","blue"), outlier.size=0.05) + ylim(lly,uly) + labs(x="",y="")
```

### Plot 2: Day VS Night

``` r
# First of all, we create a column of the dataframe to determine if the timestamp is during the day or night
# we assume day runs from 5am to 5pm, and night is 5pm to 5am
if(!("daynight" %in% colnames(pmn))){
pmn = pmn %>% 
  mutate(daynight = case_when(
    hour(DatetimeFinal) < 5 ~ 'Night',
    hour(DatetimeFinal) >= 17 ~ 'Night',
    TRUE ~ 'Day'
  ))
}

# Then the wet and dry season data each get their own axes and we plot these side by side
prainy = ggplot(data = subset(pmn, Site=="IPA10ST"), aes(x=daynight, y = TotPMN)) + 
            geom_boxplot(col="blue", outlier.size=0.05) + 
            ylim(lly,uly) + labs(x="",y="")

pdry = ggplot(data = subset(pmn, Site=="TBNI21ST"), aes(x=daynight, y = TotPMN)) + 
            geom_boxplot(col="black", outlier.size=0.05) + 
            ylim(lly,uly) + labs(x="",y="")

# Plot together
p2 = pdry + prainy 
```

### Plot 3: Hourly

``` r
# Bins and labels for boxplots to come - to avoid crowding the time axis, we place hour labels every 6 hours.
repr = 5 
labl6 = c("0", rep("", repr), "6", rep("", repr), "12", rep("", repr), "18", rep("", repr), "24")

# Similar to the last plot, we simply use the hour to bin the data and have two axes, one for dry and 
# one for wet seasons
prainy = ggplot(data = subset(pmn, Site=="IPA10ST"), aes(x=sprintf("%02s", hour(DatetimeFinal)), y = TotPMN)) + 
            geom_boxplot(col="blue", outlier.size=0.05) + 
            scale_x_discrete(labels=labl6) + ylim(lly,uly) + labs(x="",y="")

pdry = ggplot(data = subset(pmn, Site=="TBNI21ST"), aes(x=sprintf("%02s", hour(DatetimeFinal)), y = TotPMN)) + 
            geom_boxplot(col="black", outlier.size=0.05) + 
            scale_x_discrete(labels=labl6) + ylim(lly,uly) + labs(x="",y="Power Minus Noise")

# Plot together
p3 = pdry + prainy  
```

### Plot 4: Fifteen minutes

``` r
# Bin the data in the dataframe
if(!("quarterhr" %in% colnames(pmn))){
pmn = pmn %>% 
  mutate(quarterhr = case_when(
    (minute(DatetimeFinal) > 7)&(minute(DatetimeFinal) <= 22) ~ paste(sprintf("%02s", hour(DatetimeFinal)),':15', sep=""),
    (minute(DatetimeFinal) > 22)&(minute(DatetimeFinal) <= 37) ~ paste(sprintf("%02s", hour(DatetimeFinal)),':15', sep=""),
    (minute(DatetimeFinal) > 37)&(minute(DatetimeFinal) <= 52) ~ paste(sprintf("%02s", hour(DatetimeFinal)),':15', sep=""),
    TRUE ~ paste(sprintf("%02s", hour(DatetimeFinal)),':00', sep="")
  ))
}

# Set the x-axis labels
repr = 11 
labl24 = c("0", rep("", repr), "6", rep("", repr), "12", rep("", repr), "18", rep("", repr), "24")

# Compute wet and dry boxplots as before
prainy = ggplot(data = subset(pmn, Site=="IPA10ST"), aes(x=quarterhr, y = TotPMN)) + 
           geom_boxplot(col="blue", outlier.size=0.05) +
           scale_x_discrete(labels=labl24) + ylim(lly,uly) + labs(x="",y="")

pdry = ggplot(data = subset(pmn, Site=="TBNI21ST"), aes(x=quarterhr, y = TotPMN)) + 
            geom_boxplot(col="black", outlier.size=0.05) +
           scale_x_discrete(labels=labl24) + ylim(lly,uly) + labs(x="",y="")

# Plot together
p4 = pdry + prainy 
```

### Plot 5: Raw data by minute

``` r
# ONE MINUTE DATA - p5
if(!("onemin" %in% colnames(pmn))){
pmn = pmn %>% mutate(onemin = hour(DatetimeFinal) + (minute(DatetimeFinal)%%60)/60)
}


# In this plot we use a representative rainy site, IPA10ST, so as not to crowd the plot with points
prainy = ggplot(data = subset(pmn, Site=="IPA10ST"), aes(x=onemin, y = TotPMN)) + 
                  geom_point(col="blue", size=0.05) +
                  ylim(lly,uly) + xlim(0,24) + 
                  scale_x_continuous(breaks=seq(0,24,by=6)) + 
                  labs(x="Hour of day",y="")

# Here we use the representative dry site, TBNI21ST, for illustration
pdry = ggplot(data = subset(pmn, Site=="TBNI21ST"), aes(x=onemin, y = TotPMN)) + 
                  geom_point(col="black", size=0.05) +
                  ylim(lly,uly) + 
                  xlim(0,24) + 
                  scale_x_continuous(breaks=seq(0,24,by=6)) + 
                  labs(x="Hour of day",y="")

# Plot together
p5 = pdry + prainy 
```

    [1m[22mScale for [32mx[39m is already present.
    Adding another scale for [32mx[39m, which will replace the existing scale.
    [1m[22mScale for [32mx[39m is already present.
    Adding another scale for [32mx[39m, which will replace the existing scale.

### Assemble final plot

``` r
# Assemble the completed stacked plot
pscales = p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 1)

# It helps to use a skinny aspect ratio
options(repr.plot.width=10, repr.plot.height=15)
# Fig 3 in 1st version of ms, Fig 3 in 2nd version as well
#ggsave(paste(figdir,"DryRainy_Scales.png",sep=""), plot = pscales, width=5, height=9)

pscales
```

    Warning message:
    "[1m[22mRemoved 27 rows containing non-finite values (`stat_boxplot()`)."
    Warning message:
    "[1m[22mRemoved 14 rows containing non-finite values (`stat_boxplot()`)."
    Warning message:
    "[1m[22mRemoved 13 rows containing non-finite values (`stat_boxplot()`)."
    Warning message:
    "[1m[22mRemoved 14 rows containing non-finite values (`stat_boxplot()`)."
    Warning message:
    "[1m[22mRemoved 13 rows containing non-finite values (`stat_boxplot()`)."
    Warning message:
    "[1m[22mRemoved 14 rows containing non-finite values (`stat_boxplot()`)."
    Warning message:
    "[1m[22mRemoved 13 rows containing non-finite values (`stat_boxplot()`)."
    Warning message:
    "[1m[22mRemoved 14 rows containing missing values (`geom_point()`)."
    Warning message:
    "[1m[22mRemoved 13 rows containing missing values (`geom_point()`)."

![png](output_34_1.png)

## Boxplot of Total PMN for each site

This figure will show, for each time series, box plots representing each hour of the day. We will first create the label for the x-axis in the next cell and then create the ggplot2 function that computes these generically for all of the variables in the pmn dataset.

We will return to this boxplot later to add our waveform to it. So it is useful to give it a name, pbox.

``` r
boxfun <- function(data, var, nr=2) {
  # sprintf etc. creates a 2-digit hour
  ggplot(data, aes(x=sprintf("%02s", hour(DatetimeFinal)), y = {{var}}, color=fct_rev(factor(Season)))) +
  geom_boxplot() + scale_colour_manual(values=c("blue", "black")) +
  scale_x_discrete(labels=labl6) +
  scale_y_continuous(limits=c(0,3200)) +
  labs(
    x = "Hour of day", 
    y = "Power Minus Noise", 
    colour = "Site Name",
    lty = "Season"
  ) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  facet_wrap(vars(Site), nrow = nr, scales = "free_x")
}

#pmn %>% dplyr::filter(Season %in% c('rainy', 'dry'))
              
pbox = pmn |> boxfun(TotPMN)
pbox
```

    Warning message:
    "[1m[22mRemoved 76 rows containing non-finite values (`stat_boxplot()`)."

![png](output_36_1.png)

## Compute the periodogram

For the analysis of periodic components of time series, it is helpful to consider the Fourier series of the time series as well as its power spectrum estimate. Denote by $x_t$ the time series of interest, where $t = 1, \ldots, N$ and assume $N$ is odd, for simplicity of the following expressions. We will compute the periodogram estimate of the spectrum by way of the Discrete Fourier Transform of $x_t$,

$$ d(j/N) = \frac{1}{\sqrt{N}} \sum_{t = 1}^N x_t \exp (-2 \pi i t j/N). $$

The series $|d(j/N)|^2$ is referred to as the periodogram which estimates a quantity called the power spectrum of a time series.

``` r
# Helper function to get the scaled Fourier transformed data
ftgram <- function(x, M) { fft(c(x, matrix(0, M-length(x))))*sqrt(1/length(x))}

# Helper function to get the frequency grid
fgrid <- function(N, dt) { dt*seq(from = 0, to = 1, length.out = N+1)[seq(1,N)]}
```

``` r
# Set up to save the magnitudes of each pgram peak
N = array(0, nsites) # Number of samples per series
pad = 4 # padding factor (interpolation in frequency)
M = array(0, nsites) # Number of frequency bins (4*N here)
xbar = array(0, nsites) # sample mean of each series
maxcycd = 12 # maximum number of harmonics of 1 cyc/d to use
cycind = zeros(maxcycd, nsites) # Frequency bin for which 1 cyc/d, 2 cyc/d, ..., 12 cyc/d happens 

# This array will hold the periodograms
simplespec1 = array(0, pad*length(pmn$TotPMN))

# This array will hold the frequency grid
freq1 = array(0, pad*length(pmn$TotPMN))

# This array will contain the name of the site for which we've calculated the periodogram
names = array("", pad*length(pmn$TotPMN))

# Sample means of the series
xbar = array(0, nsites)

# Variance
ssq = array(0, nsites)

# For the reconstruction, Rsq is the R-squared value for the regression, 
# coeffs_pgram is the coefficent in front of the cosine for 1, 2,..., 12 
# cyc/d, and phase_pgram is the corresponding phase
coeffs_pgram = zeros(maxcycd, nsites)
phase_pgram = zeros(maxcycd, nsites)
soln = rep(0,length(pmn$TotPMN))

# Compute the periodogram and do the regression
for (i in seq(from = 1, to = nsites)){
    N[i] = length(seq(ind[i, 1], ind[i, 2]))
    M[i] = pad*N[i]
    # sind is the first index of each of the series after padding 
    sind = cumsum(c(0, M))
    xbar[i] = mean(pmn$TotPMN[seq(ind[i,1], ind[i,2])])
    ssq[i] = var(pmn$TotPMN[seq(ind[i,1], ind[i,2])])
    # Compute the fourier transform of x-xbar
    before_squaring = ftgram(pmn$TotPMN[seq(ind[i,1], ind[i,2])] - xbar[i], M[i])
    # Square it to get the periodogram
    simplespec1[seq(sind[i]+1, sind[i+1])] = abs(before_squaring)**2
    # Frequency grid and names for the i-th series
    freq1[seq(sind[i]+1, sind[i+1])] = fgrid(M[i], 60*24)
    names[seq(sind[i]+1, sind[i+1])] = usite[i]
    # Find the index in frequency at which 1 cyc/d, 2 cyc/d, ... 12 cyc/d lives
    for(ii in seq(1, maxcycd)){
        cycind[ii, i]= sum(freq1[seq(sind[i]+1, sind[i+1])]<=ii)
    }
    # Find the magnitude and phase
    coeffs_pgram[,i] = 2*sqrt(simplespec1[seq(sind[i]+1, sind[i+1])][cycind[,i]])/sqrt(N[i])
    phase_pgram[,i] = -atan2(-Im(before_squaring[cycind[,i]]), Re(before_squaring[cycind[,i]]))
    # Compute the waveform using the coefficients, phases, and the mean value of the time series
    soln[seq(ind[i, 1], ind[i, 2])] = cos(2*pi*kronecker(matrix(1,N[i],1), t(seq(1,12)))*kronecker(matrix(1,1,12), seq(1,N[i]))/(24*60) + kronecker(matrix(1,N[i],1), t(phase_pgram[,i])))  %*% coeffs_pgram[,i] + xbar[i]
}
```

## Compute the R-squared value for the model with 12 cosines

That is, given the model equation for the waveform below:

$$  w_t = \bar{x} + \sum_{j = 1}^{12} c_{j} \cos (2\pi t j + \phi_j) $$

We compute the difference between $x_t$ and $w_t$ and compare the fit.

``` r
# Get the R^2s - here we just use the linear model (lm) function to determine the significance of the
# variance remaining after computing the model
Rsq = array(0, nsites)
for (i in seq(from = 1, to = nsites)){
    
    # Get the indices of 12.5 cyc/d, 1 cyc/d, ... 12 cyc/d, etc
    indfreqbin1 = sum(freq1[seq(sind[i], sind[i+1])] <= 12.5)
    
    print(paste("A regression of site", usite[i], "on the reconstruction of its first", maxcycd, "harmonics of 1 cyc/d yields:"))
    print(summary(fit <- lm(pmn$TotPMN[seq(ind[i,1], ind[i,2])]~soln[seq(ind[i, 1], ind[i, 2])], na.action=NULL)))
    Rsq[i] = summary(fit)$r.squared

}
```

    [1] "A regression of site IPA10ST on the reconstruction of its first 12 harmonics of 1 cyc/d yields:"

    Call:
    lm(formula = pmn$TotPMN[seq(ind[i, 1], ind[i, 2])] ~ soln[seq(ind[i, 
        1], ind[i, 2])], na.action = NULL)

    Residuals:
       Min     1Q Median     3Q    Max 
    -639.8 -158.2  -40.6   89.2 3564.4 

    Coefficients:
                                    Estimate Std. Error t value Pr(>|t|)    
    (Intercept)                      0.16757   26.93085   0.006    0.995    
    soln[seq(ind[i, 1], ind[i, 2])]  0.99981    0.03081  32.456   <2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 302.2 on 5758 degrees of freedom
    Multiple R-squared:  0.1547,    Adjusted R-squared:  0.1545 
    F-statistic:  1053 on 1 and 5758 DF,  p-value: < 2.2e-16

    [1] "A regression of site IPA11ST on the reconstruction of its first 12 harmonics of 1 cyc/d yields:"

    Call:
    lm(formula = pmn$TotPMN[seq(ind[i, 1], ind[i, 2])] ~ soln[seq(ind[i, 
        1], ind[i, 2])], na.action = NULL)

    Residuals:
       Min     1Q Median     3Q    Max 
    -811.5 -149.7  -37.0   89.9 4038.5 

    Coefficients:
                                    Estimate Std. Error t value Pr(>|t|)    
    (Intercept)                       0.1355    25.9133   0.005    0.996    
    soln[seq(ind[i, 1], ind[i, 2])]   0.9998     0.0315  31.742   <2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 291.7 on 5758 degrees of freedom
    Multiple R-squared:  0.1489,    Adjusted R-squared:  0.1488 
    F-statistic:  1008 on 1 and 5758 DF,  p-value: < 2.2e-16

    [1] "A regression of site IPA14ST on the reconstruction of its first 12 harmonics of 1 cyc/d yields:"

    Call:
    lm(formula = pmn$TotPMN[seq(ind[i, 1], ind[i, 2])] ~ soln[seq(ind[i, 
        1], ind[i, 2])], na.action = NULL)

    Residuals:
       Min     1Q Median     3Q    Max 
    -785.7 -167.4  -32.5  113.8 5794.9 

    Coefficients:
                                    Estimate Std. Error t value Pr(>|t|)    
    (Intercept)                      0.35928   42.33900   0.008    0.993    
    soln[seq(ind[i, 1], ind[i, 2])]  0.99962    0.04486  22.283   <2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 328.6 on 11518 degrees of freedom
    Multiple R-squared:  0.04133,   Adjusted R-squared:  0.04125 
    F-statistic: 496.5 on 1 and 11518 DF,  p-value: < 2.2e-16

    [1] "A regression of site IPA15ST on the reconstruction of its first 12 harmonics of 1 cyc/d yields:"

    Call:
    lm(formula = pmn$TotPMN[seq(ind[i, 1], ind[i, 2])] ~ soln[seq(ind[i, 
        1], ind[i, 2])], na.action = NULL)

    Residuals:
        Min      1Q  Median      3Q     Max 
    -901.54 -148.63  -20.44  104.63 2636.22 

    Coefficients:
                                    Estimate Std. Error t value Pr(>|t|)    
    (Intercept)                      0.20854   19.54759   0.011    0.991    
    soln[seq(ind[i, 1], ind[i, 2])]  0.99977    0.02164  46.199   <2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 258.6 on 11518 degrees of freedom
    Multiple R-squared:  0.1563,    Adjusted R-squared:  0.1563 
    F-statistic:  2134 on 1 and 11518 DF,  p-value: < 2.2e-16

    [1] "A regression of site TBNI21ST on the reconstruction of its first 12 harmonics of 1 cyc/d yields:"

    Call:
    lm(formula = pmn$TotPMN[seq(ind[i, 1], ind[i, 2])] ~ soln[seq(ind[i, 
        1], ind[i, 2])], na.action = NULL)

    Residuals:
       Min     1Q Median     3Q    Max 
    -538.5 -140.4  -37.1   97.3 5079.5 

    Coefficients:
                                    Estimate Std. Error t value Pr(>|t|)    
    (Intercept)                      0.05986    9.29767   0.006    0.995    
    soln[seq(ind[i, 1], ind[i, 2])]  0.99993    0.01049  95.298   <2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 228.8 on 8638 degrees of freedom
    Multiple R-squared:  0.5125,    Adjusted R-squared:  0.5125 
    F-statistic:  9082 on 1 and 8638 DF,  p-value: < 2.2e-16

    [1] "A regression of site TBNI22ST on the reconstruction of its first 12 harmonics of 1 cyc/d yields:"

    Call:
    lm(formula = pmn$TotPMN[seq(ind[i, 1], ind[i, 2])] ~ soln[seq(ind[i, 
        1], ind[i, 2])], na.action = NULL)

    Residuals:
       Min     1Q Median     3Q    Max 
    -555.9 -136.5  -33.8   88.9 4698.9 

    Coefficients:
                                    Estimate Std. Error t value Pr(>|t|)    
    (Intercept)                      0.05968   10.38731   0.006    0.995    
    soln[seq(ind[i, 1], ind[i, 2])]  0.99993    0.01234  81.001   <2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 229.7 on 8638 degrees of freedom
    Multiple R-squared:  0.4317,    Adjusted R-squared:  0.4316 
    F-statistic:  6561 on 1 and 8638 DF,  p-value: < 2.2e-16

    [1] "A regression of site TBNI24ST on the reconstruction of its first 12 harmonics of 1 cyc/d yields:"

    Call:
    lm(formula = pmn$TotPMN[seq(ind[i, 1], ind[i, 2])] ~ soln[seq(ind[i, 
        1], ind[i, 2])], na.action = NULL)

    Residuals:
       Min     1Q Median     3Q    Max 
    -998.2 -167.7  -37.8  123.6 5513.7 

    Coefficients:
                                    Estimate Std. Error t value Pr(>|t|)    
    (Intercept)                      0.07721   14.00092   0.006    0.996    
    soln[seq(ind[i, 1], ind[i, 2])]  0.99992    0.01381  72.402   <2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 274.4 on 10078 degrees of freedom
    Multiple R-squared:  0.3422,    Adjusted R-squared:  0.3421 
    F-statistic:  5242 on 1 and 10078 DF,  p-value: < 2.2e-16

    [1] "A regression of site TBNI25ST on the reconstruction of its first 12 harmonics of 1 cyc/d yields:"

    Call:
    lm(formula = pmn$TotPMN[seq(ind[i, 1], ind[i, 2])] ~ soln[seq(ind[i, 
        1], ind[i, 2])], na.action = NULL)

    Residuals:
       Min     1Q Median     3Q    Max 
    -786.1 -169.1  -46.1  103.2 6269.6 

    Coefficients:
                                    Estimate Std. Error t value Pr(>|t|)    
    (Intercept)                      0.07935   15.59682   0.005    0.996    
    soln[seq(ind[i, 1], ind[i, 2])]  0.99991    0.01712  58.420   <2e-16 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    Residual standard error: 323.8 on 10078 degrees of freedom
    Multiple R-squared:  0.253, Adjusted R-squared:  0.2529 
    F-statistic:  3413 on 1 and 10078 DF,  p-value: < 2.2e-16

Since we are about to plot the periodograms, it is nice to have a dataframe with some annotations for reference. Here we just collect the values of the periodogram at the 1 cyc/d, 2 cyc/d, etc. frequencies for ggplot2

``` r
freq_dots = c()
pgram_dots = c()
names_dots = c()
coeffs_dots = c()
phase_dots = c()
label_dots = c()
for (i in seq(1,nsites)){
    names_dots = c(names_dots, rep(usite[i], each=12))
    freq_dots = c(freq_dots, freq1[seq(sind[i]+1, sind[i+1])][cycind[,i]])
    pgram_dots = c(pgram_dots, simplespec1[seq(sind[i]+1, sind[i+1])][cycind[,i]])
    coeffs_dots = c(coeffs_dots, round(coeffs_pgram[,i], digits=0))
    phase_dots = c(phase_dots, round(phase_pgram[,i]*180/pi, digits=0))
    label_dots = c(label_dots, paste("(",round(coeffs_pgram[,i], digits=0),", ",
                           round(phase_pgram[,i]*180/pi, digits=0),"°)", sep=""))
}

ann_df = data.frame(frequency=freq_dots, periodogram=pgram_dots, Site = names_dots, coeffs=coeffs_dots, 
                    phase = phase_dots, labels = label_dots)
if(!("Rsq" %in% colnames(key))){
key = data.frame(key,Rsq)
    }
ann_df <- merge(ann_df, key, by = c("Site"))
head(ann_df)
```

+---+---------+-----------+-------------+---------+---------+-------------+---------+-----------+
|   | Site    | frequency | periodogram | coeffs  | phase   | labels      | Season  | Rsq       |
+===+=========+===========+=============+=========+=========+=============+=========+===========+
|   | \<chr\> | \<dbl\>   | \<dbl\>     | \<dbl\> | \<dbl\> | \<chr\>     | \<fct\> | \<dbl\>   |
+---+---------+-----------+-------------+---------+---------+-------------+---------+-----------+
| 1 | IPA10ST | 1         | 10900255.6  | 87      | -169    | (87, -169°) | Wet     | 0.1546511 |
+---+---------+-----------+-------------+---------+---------+-------------+---------+-----------+
| 2 | IPA10ST | 2         | 6737543.6   | 68      | 68      | (68, 68°)   | Wet     | 0.1546511 |
+---+---------+-----------+-------------+---------+---------+-------------+---------+-----------+
| 3 | IPA10ST | 3         | 13556384.7  | 97      | 135     | (97, 135°)  | Wet     | 0.1546511 |
+---+---------+-----------+-------------+---------+---------+-------------+---------+-----------+
| 4 | IPA10ST | 4         | 7703796.0   | 73      | -36     | (73, -36°)  | Wet     | 0.1546511 |
+---+---------+-----------+-------------+---------+---------+-------------+---------+-----------+
| 5 | IPA10ST | 5         | 857290.9    | 24      | 70      | (24, 70°)   | Wet     | 0.1546511 |
+---+---------+-----------+-------------+---------+---------+-------------+---------+-----------+
| 6 | IPA10ST | 6         | 563694.6    | 20      | 99      | (20, 99°)   | Wet     | 0.1546511 |
+---+---------+-----------+-------------+---------+---------+-------------+---------+-----------+

: A data.frame: 6 × 8

This dataframe will contain the periodograms for ggplot2 to read

``` r
df = data.frame(frequency = freq1, periodogram=simplespec1, Site=names, prefix = gsub("^\\d+|\\d+$","",substr(names,1,3)))
df <- arrange(mutate(df,Site=factor(Site,levels=usite)),Site)
df <- merge(df, key, by = c("Site"))
head(df)
```

+---+---------+-----------+--------------+---------+---------+-----------+
|   | Site    | frequency | periodogram  | prefix  | Season  | Rsq       |
+===+=========+===========+==============+=========+=========+===========+
|   | \<fct\> | \<dbl\>   | \<dbl\>      | \<chr\> | \<fct\> | \<dbl\>   |
+---+---------+-----------+--------------+---------+---------+-----------+
| 1 | IPA10ST | 0.0000    | 4.830965e-25 | IPA     | Wet     | 0.1546511 |
+---+---------+-----------+--------------+---------+---------+-----------+
| 2 | IPA10ST | 0.0625    | 2.612780e+06 | IPA     | Wet     | 0.1546511 |
+---+---------+-----------+--------------+---------+---------+-----------+
| 3 | IPA10ST | 0.1250    | 7.403034e+06 | IPA     | Wet     | 0.1546511 |
+---+---------+-----------+--------------+---------+---------+-----------+
| 4 | IPA10ST | 0.1875    | 9.053150e+06 | IPA     | Wet     | 0.1546511 |
+---+---------+-----------+--------------+---------+---------+-----------+
| 5 | IPA10ST | 0.2500    | 6.218598e+06 | IPA     | Wet     | 0.1546511 |
+---+---------+-----------+--------------+---------+---------+-----------+
| 6 | IPA10ST | 0.3125    | 2.166938e+06 | IPA     | Wet     | 0.1546511 |
+---+---------+-----------+--------------+---------+---------+-----------+

: A data.frame: 6 × 6

``` r
# Print summary of data in a table - available in supplementary materials
begin_time = pmn$DatetimeFinal[ind[,1]]
end_time = pmn$DatetimeFinal[ind[,2]]
table = cbind(usite, paste(begin_time), paste(end_time)) # Minute(end_time-begin_time))
table = cbind(table, round(xbar,digits=2), round(sqrt(ssq),digits=2), round(Rsq,digits=3)) 
print(table)

#fn = paste(figdir,"Rsq_table.tex",sep="")
#tab<-xtable(table, caption= "\\label{tab:eightsitesRsq} R-squared, mean, variance for the total power minus noise for the eight sites. Note that the periodic waveform captures about $15$ percent of the variance for the rainy time series (excepting IPA14ST) and $25-50$ percent of the dry time series.", 
#            align=c("|c","|c","|c","|c","|c","|c","|c|"))
#if (file.exists(fn)) {
  # Delete file if it exists
#  file.remove(fn)
#}
#names(tab) <- c('Site ID','Start Date', 'End Date', 'Sample Mean', 'Sample SD', 'R-sq.' )
# Uses latex format so as to drop into our document
#print(tab,file=fn,append=T,table.placement = "h", include.rownames=FALSE,
# caption.placement="bottom", hline.after=seq(from=-1,to=nrow(tab),by=1))
```

         usite                                                                   
    [1,] "IPA10ST"  "2021-02-20 00:00:00" "2021-02-23 23:59:00" "864.63" "328.63"
    [2,] "IPA11ST"  "2021-02-23 00:00:00" "2021-02-26 23:59:00" "813.57" "316.17"
    [3,] "IPA14ST"  "2021-02-23 00:00:00" "2021-03-02 23:59:00" "941.34" "335.57"
    [4,] "IPA15ST"  "2021-02-22 00:00:00" "2021-03-01 23:59:00" "896.4"  "281.55"
    [5,] "TBNI21ST" "2021-07-16 23:00:00" "2021-07-22 22:59:00" "854.51" "327.63"
    [6,] "TBNI22ST" "2021-07-16 23:00:00" "2021-07-22 22:59:00" "817.28" "304.71"
    [7,] "TBNI24ST" "2021-07-16 23:00:00" "2021-07-23 22:59:00" "994.27" "338.29"
    [8,] "TBNI25ST" "2021-07-16 23:00:00" "2021-07-23 22:59:00" "891.55" "374.6" 
                
    [1,] "0.155"
    [2,] "0.149"
    [3,] "0.041"
    [4,] "0.156"
    [5,] "0.513"
    [6,] "0.432"
    [7,] "0.342"
    [8,] "0.253"

## Plot all of the periodograms together

Plot the periodogram generated for each site along with the R\^2 value for each model fit

``` r
# Plot specifications
options(repr.plot.width=15, repr.plot.height=10)

# Plot all of the periodograms together
p = ggplot(df, aes(frequency, periodogram)) + 
    geom_line(aes(colour = Site)) +
    scale_y_log10(name = "Periodogram",limits = c(1e4,1e8)) +
    # The value of R^2 should be in the top right of each periodogram
    geom_text(aes(x=10,y=5e7,label=round(Rsq,digits=2))) +
    scale_x_continuous(name = "Frequency (cyc/d)", limits = c(0,13), breaks=seq(0,13, by=2)) +
    geom_point(data=ann_df, aes(x=frequency,y=periodogram), col=4) +
    geom_text(data=ann_df, aes(x = frequency, y = periodogram, label=coeffs, angle=90), col=4, hjust=-0.5) +
    geom_text(data=ann_df, aes(x = 8, y = 5e8, label=Site)) +
    facet_wrap(~factor(Site, levels=usite),nrow=2)

# ggsave(paste(figdir, "pgram_All.png", sep=""), width=25,height=15,plot=p)
p
```

    Warning message:
    "[1m[22mRemoved 285406 rows containing missing values (`geom_line()`)."
    Warning message:
    "[1m[22mRemoved 5 rows containing missing values (`geom_point()`)."
    Warning message:
    "[1m[22mRemoved 5 rows containing missing values (`geom_text()`)."
    Warning message:
    "[1m[22mRemoved 96 rows containing missing values (`geom_text()`)."

![png](output_48_1.png)

## Visualize the waveforms

Now visualise the waveforms themselves. All of the sites have start times starting at 1 minute past midnight, so no shifting is required.

``` r
# Confirm start times
pmn$DatetimeFinal[ind[,1]]
```

    [1] "2021-02-20 00:00:00 GMT" "2021-02-23 00:00:00 GMT"
    [3] "2021-02-23 00:00:00 GMT" "2021-02-22 00:00:00 GMT"
    [5] "2021-07-16 23:00:00 GMT" "2021-07-16 23:00:00 GMT"
    [7] "2021-07-16 23:00:00 GMT" "2021-07-16 23:00:00 GMT"

``` r
# Arrange the waveforms in a dataframe for ggplot to read - we'll take only the first 24 hours for plotting
soln24h = array(0,60*24*nsites)
hour = array(0,60*24*nsites)
for (i in 1:nsites){
  soln24h[seq((1+60*24*(i-1)),(60*24*i))] = soln[seq(ind[i, 1], ind[i, 2])][seq(1,60*24)] #- xbar[i]
  hour[seq((1+60*24*(i-1)),(60*24*i))] = seq(0,(60*24)-1)/60
}
soln24h = data.frame(Waveform=soln24h, hour = hour, Site = rep(usite, each=1440))
soln24h <- merge(soln24h, key, by = c("Site"))

head(soln24h)
```

+---+---------+----------+------------+---------+-----------+
|   | Site    | Waveform | hour       | Season  | Rsq       |
+===+=========+==========+============+=========+===========+
|   | \<chr\> | \<dbl\>  | \<dbl\>    | \<fct\> | \<dbl\>   |
+---+---------+----------+------------+---------+-----------+
| 1 | IPA10ST | 789.2587 | 0.00000000 | Wet     | 0.1546511 |
+---+---------+----------+------------+---------+-----------+
| 2 | IPA10ST | 787.8025 | 0.01666667 | Wet     | 0.1546511 |
+---+---------+----------+------------+---------+-----------+
| 3 | IPA10ST | 786.3649 | 0.03333333 | Wet     | 0.1546511 |
+---+---------+----------+------------+---------+-----------+
| 4 | IPA10ST | 784.9459 | 0.05000000 | Wet     | 0.1546511 |
+---+---------+----------+------------+---------+-----------+
| 5 | IPA10ST | 783.5454 | 0.06666667 | Wet     | 0.1546511 |
+---+---------+----------+------------+---------+-----------+
| 6 | IPA10ST | 782.1631 | 0.08333333 | Wet     | 0.1546511 |
+---+---------+----------+------------+---------+-----------+

: A data.frame: 6 × 5

## Compare waveform with boxplot

Let's superimpose the waveform obtained using Fourier series methods on the boxplots made with hourly binning. This comparison should be helpful in assessing outliers and such.

``` r
# Plot specifications
options(repr.plot.width=15, repr.plot.height=12)

# Plot
pbox2 = pbox + 
    geom_line(data=subset(soln24h,hour%%1==0), aes(x = hour,y=Waveform), lwd=1.5, colour="red",alpha=0.5)

#ggsave(paste(figdir,"TOTPMN_HourBox_waveform.png",sep=""), width=10, height=6)
pbox2
```

    Warning message:
    "[1m[22mRemoved 76 rows containing non-finite values (`stat_boxplot()`)."
    Warning message:
    "[1m[22mRemoved 76 rows containing non-finite values (`stat_boxplot()`)."

![png](output_54_1.png)

## Plot waveforms between seasons

With dawn and dusk shaded, let's plot the waveforms altogether, grouping dry and wet season time series to identify any similarities within seasons.

``` r
unique(soln24h$Site)
```

```{=html}
<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
```
<ol class="list-inline">

<li>'IPA10ST'</li>

<li>'IPA11ST'</li>

<li>'IPA14ST'</li>

<li>'IPA15ST'</li>

<li>'TBNI21ST'</li>

<li>'TBNI22ST'</li>

<li>'TBNI24ST'</li>

<li>'TBNI25ST'</li>

</ol>

``` r
# Plot specifications
options(repr.plot.width=12, repr.plot.height=16)

# all blue for rainy, all black for dry
ggplot(soln24h,aes(x = hour, y = Waveform)) +
    # Shade dawn and dusk
    geom_rect(aes(fill = Site), xmin = 4, xmax = 6, ymin = -Inf,ymax = Inf,alpha = 0.01) +
    geom_rect(aes(fill = Site), xmin = 16, xmax = 18, ymin = -Inf,ymax = Inf,alpha = 0.01) +
    scale_fill_manual(values=rep(c("grey"),nsites)) +
    geom_line(aes(x = hour, col = Site)) +  
    scale_colour_manual(values=c('dodgerblue','dodgerblue1','dodgerblue2','dodgerblue3','grey4','grey5','grey6','grey7')) +
    scale_x_continuous(name = "Hour of day", limits = c(0,24), breaks = seq(0, 24, by=6)) +
    scale_y_continuous(name = "Power Minus Noise") +
    facet_grid(vars(Season)) + 
    theme(legend.position = "none")

#ggsave(paste(figdir,"Waveforms_all_RD.png",sep="/"), width=6, height=6)
#ggsave(paste(figdir,"Waveforms_all_RD.pdf",sep="/"), width=6, height=6)
```

![png](output_57_0.png)

## Compare waveform with the raw data

``` r
# Put the waveform together with the original dataframe and make a plot with waveforms and data 
# as well as periodograms and annotations
pmn = cbind(pmn, waveform = soln)
```

``` r
substr(usite[1], start = 1, stop = 1) == "I"
```

TRUE

``` r
# Plot specifications
options(repr.plot.width=14, repr.plot.height=10)

# Plot waveform on left, periodogram on right for each site (Figure 3 in the paper)
for (st in usite) {

    # Sites beginning with I should be colored blue, else black (wet/dry)
    if (substr(st, start = 1, stop = 1) == "I"){
        cc = "blue"}
    else {cc = "black" }
    
    # Left plot is the waveform (red) over the data (black)
p_left = ggplot(subset(pmn, Site==st), aes(x=DatetimeFinal, y = TotPMN)) +
  geom_line(col=cc) +
  labs(x = "Date",  
    colour = "Site Name",
    fill="Season",
    y = "Power Minus Noise"
  ) + ggtitle(st) +
  scale_fill_manual(values=c("beige", "lightblue")) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  geom_line(aes(x=DatetimeFinal,y=waveform),col="red",lwd=1.2)

p_right = ggplot(subset(df, Site==st), aes(frequency, periodogram)) + 
    geom_line(col=cc) +
    scale_y_log10(name = "Periodogram",limits = c(1e4,1e9)) +
    scale_x_continuous(name = "Frequency (cyc/d)", limits = c(0,13), breaks=seq(0,13, by=2),
                       sec.axis = sec_axis(~ 24/., name="Period (hrs)", breaks=c(48, 24,12,8,6,4,2))) +
    geom_point(data=subset(ann_df, Site==st), aes(x=frequency,y=periodogram), col=4) +
    geom_text(data=subset(ann_df, Site==st), aes(x = frequency, y = periodogram, label=coeffs, angle=90), col=4, hjust=-0.5) +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    theme(legend.position = "none")

p_left + p_right

#ggsave(paste(figdir,"reconstructionpgram",st,".png",sep=""), p_left + p_right, width=8, height=4)
}

p_left + p_right
```

    Warning message:
    "[1m[22mRemoved 39956 rows containing missing values (`geom_line()`)."

![png](output_61_1.png)

## Magnitude squared coherence (MSC)

### Comparing multiple sites at a time

Suppose we want to compare sites i and j. This is similar to example 4.21 in Shumway and Stoffer. We use the astsa package to simply compute the coherences. We have used the default smoothing from the R *astsa* package (which matches *spec.pgram* from *baseR*, a modified Daniell smoother with parameter 6, to compute the MSC. Finally, the MSC was computed by using the ratio

$$ c^2(f) = \frac{|S_{1,2}(f)|^2}{\sqrt{S_{1,1}(f) S_{2,2}(f)}} $$

the \emph{phase} was computed as

$$ \phi(f) = \angle S_{1,2}(f) $$

where the angle notation in the above is accomplished in practice using the arctangent function that takes two arguments, the first being the negative imaginary part of $S_{1,2}(f)$ and the second being the real part of $S_{1,2}(f)$.

``` r
# List all of the pairs of sites to be compared 
pairs = t(combn(seq(1,nsites), 2))
Np = (nsites-1)*nsites/2
```

For MSC to be computed, the time series needs to be the same length. Therefore, we clip the data to the first n samples of the data.

``` r
# Calculate duration for each site 
duration = ind[,2]-ind[,1]
duration

# Subset data to the first 5760 samples* from each time series
mod_ind = cbind(ind[,1],ind[,1]+5760)
```

```{=html}
<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
```
<ol class="list-inline">

<li>5759</li>

<li>5759</li>

<li>11519</li>

<li>11519</li>

<li>8639</li>

<li>8639</li>

<li>10079</li>

<li>10079</li>

</ol>

\*Note: This could potentially be relaxed for some pairs, as 5760 is actually just the length of the shortest series, but for looping, we simply use 5760

``` r
# Print table of the 1 cyc/d 2 cyc/d etc harmonics, and values of the squared coherence and phase, for reference
# harmonics <- data.frame (i, j, fcyc, sr$coh[cycind], sr$phase[cycind]*180/pi)
```

``` r
# Compute the MSC
site1 = ts(pmn$TotPMN[seq(mod_ind[1,1], mod_ind[1,2])], frequency = (60*24))
site2 = ts(pmn$TotPMN[seq(mod_ind[2,1], mod_ind[2,2])], frequency = (60*24))
site3 = ts(pmn$TotPMN[seq(mod_ind[3,1], mod_ind[3,2])], frequency = (60*24))
site4 = ts(pmn$TotPMN[seq(mod_ind[4,1], mod_ind[4,2])], frequency = (60*24))

site5 = ts(pmn$TotPMN[seq(mod_ind[5,1], mod_ind[5,2])], frequency = (60*24))
site6 = ts(pmn$TotPMN[seq(mod_ind[6,1], mod_ind[6,2])], frequency = (60*24))
site7 = ts(pmn$TotPMN[seq(mod_ind[7,1], mod_ind[7,2])], frequency = (60*24))
site8 = ts(pmn$TotPMN[seq(mod_ind[8,1], mod_ind[8,2])], frequency = (60*24))

sr = mvspec(cbind(site1,site2,site3,site4,site5,site6,site7,site8), kernel("daniell",6), plot=FALSE, pad = 4)
```

``` r
# The following code helps us to assign significance to the coherence as described in Shumway and Stoffer text.
logRange <- function(logmin, logmax, N) { 1 - 10**(-seq(from=logmin, to=logmax, len=N)) }
sigMax = 4
sigs = logRange(1,sigMax,sigMax)
L=13 # daniell 6 means length is 2*6+1
print(sigs)
labs = qf(sigs, 2, 2*L-2)/(L-1+qf(sigs, 2, 2*L-2))
print(labs)
```

    [1] 0.9000 0.9900 0.9990 0.9999
    [1] 0.1745958 0.3187079 0.4376587 0.5358411

``` r
# Find the {one, two, three} cyc/day frequency index - 
# We need to compute these again because they are not the same indices as for the spectra (as these were subsetted)
fcyc = seq(1, 12, by=1)
cycind = fcyc
for(ii in fcyc){
    cycind[ii]= sum(sr$freq<=ii)
}
cycind
```

```{=html}
<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
```
<ol class="list-inline">

<li>20</li>

<li>40</li>

<li>60</li>

<li>81</li>

<li>101</li>

<li>121</li>

<li>141</li>

<li>162</li>

<li>182</li>

<li>202</li>

<li>222</li>

<li>243</li>

</ol>

``` r
# Calculate the mean of the squared coherence at the harmonics of 1 cyc/d
ddf = data.frame(first=usite[pairs[,1]],  
      second=usite[pairs[,2]], 
      avecoh = colMeans(sr$coh[cycind,]))

# Sort in descending order
df2 <- ddf[order(ddf$avecoh,decreasing=TRUE),]
head(df2)
```

+----+----------+----------+-----------+
|    | first    | second   | avecoh    |
+====+==========+==========+===========+
|    | \<chr\>  | \<chr\>  | \<dbl\>   |
+----+----------+----------+-----------+
| 12 | IPA11ST  | TBNI24ST | 0.5361143 |
+----+----------+----------+-----------+
| 27 | TBNI22ST | TBNI25ST | 0.5187248 |
+----+----------+----------+-----------+
| 11 | IPA11ST  | TBNI22ST | 0.5111882 |
+----+----------+----------+-----------+
| 21 | IPA15ST  | TBNI24ST | 0.5045685 |
+----+----------+----------+-----------+
| 7  | IPA10ST  | TBNI25ST | 0.5034458 |
+----+----------+----------+-----------+
| 26 | TBNI22ST | TBNI24ST | 0.5000783 |
+----+----------+----------+-----------+

: A data.frame: 6 × 3

## Magnitude squared coherence plot

The msc is to be interpreted as a frequency dependent correlation coefficient, that is, say, at 1 cycle per day we see that (as labeled) the top plot shows that the squared coherence is 0.655 which means that at one cycle per day the two series, IPA10ST and IPA11ST are 0.655 correlated. Squared coherence is bounded between 0 and 1 and thus 0.655 is unlikely to happen if the true correlation at 1 cycle per day is actually zero. Ditto for three cycles per day (0.929), four cycles per day (0.9), etc.

The second number in the annotated tuple on the plot is the phase difference between the two series as a function of frequency.

``` r
# Condensed msc for all pairs (doesn't show phase)
for (k in seq(1,15)){
    i=pairs[k,1]
    j=pairs[k,2]
    
# Plot specifications
options(repr.plot.width=14, repr.plot.height=10)

# Figures 8
#png(paste(figdir, "msc_",usite[i],usite[j],"_condensed.png", sep=""))
    par(mar = c(3, 3, 1, 1), mgp = c(1.6, 0.6, 0), xpd=FALSE) #, cin=c(3.5,7))
    plot(sr$freq, sr$coh[,k], type='l', 
         main = paste("Totpmn Coherence between", usite[i], "and", usite[j]),
         xlim = c(0,12), col = 1, ylab="Sq. Coherence", xlab="", ylim=c(0,1.4))
    minor.tick(nx = 2, ny = 2)
    points(sr$freq[cycind], sr$coh[cycind,k], col=4)
    text(x=sr$freq[cycind], y=sr$coh[cycind,k], 
         labels=paste("(",round(sr$coh[cycind,k], digits=3),", ",round(sr$phase[cycind,k]*180/pi, digits=1),"°)", sep=""), 
         adj=-0.3, srt=90, col=4)
    for (ii in fcyc){
      lines(x=c(ii,ii),y=c(0,sr$coh[cycind[ii],k]), lty = 3)
      }
    for (ll in labs){
      lines(x=c(-0.25,12.25), y=c(ll,ll), lty=2, col='red')
    }
#dev.off()
}
```

![png](output_73_0.png)

![png](output_73_1.png)

![png](output_73_2.png)

![png](output_73_3.png)

![png](output_73_4.png)

![png](output_73_5.png)

![png](output_73_6.png)

![png](output_73_7.png)

![png](output_73_8.png)

![png](output_73_9.png)

![png](output_73_10.png)

![png](output_73_11.png)

![png](output_73_12.png)

![png](output_73_13.png)

![png](output_73_14.png)

``` r
# Condensed msc for all pairs (doesn't show phase)
for (k in seq(1,28)){
    i=pairs[k,1]
    j=pairs[k,2]
    
# Plot specifications
options(repr.plot.width=40, repr.plot.height=13)


#png(paste(figdir, "/msc_",usite[i],usite[j],"_condensed.png", sep=""))
    par(mar = c(4, 4, 2, 2), mgp = c(2.5, 0.6, 0), xpd=FALSE, las=1, cex=1) #, cin=c(3.5,7))
    plot(sr$freq, sr$coh[,k], type='l', 
         main = paste("Totpmn Coherence between", usite[i], "and", usite[j]),
         xlim = c(0,12), col = 1, ylab="Sq. Coherence", xlab="Frequency (cyc/d)", ylim=c(0,1.4))
   
    # Set custom x-axis labels with horizontal orientation
axis(1, at = 0:12, labels = 0:12, las = 1) # las = 1 for horizontal labels

# Increase y-axis label size
#axis(2, cex.axis = 5)
    minor.tick(nx = 2, ny = 2)
    points(sr$freq[cycind], sr$coh[cycind,k], col=4)
    text(x=sr$freq[cycind], y=sr$coh[cycind,k], 
         labels=paste("(",round(sr$coh[cycind,k], digits=3),", ",round(sr$phase[cycind,k]*180/pi, digits=1),"°)", sep=""), 
         adj=-0.3, srt=90, col=4)
    for (ii in fcyc){
      lines(x=c(ii,ii),y=c(0,sr$coh[cycind[ii],k]), lty = 3)
      }
    for (ll in labs){
      lines(x=c(-0.25,12.25), y=c(ll,ll), lty=2, col='red')
    }

  #  dev.off()
}
```

![png](output_74_0.png)

![png](output_74_1.png)

![png](output_74_2.png)

![png](output_74_3.png)

![png](output_74_4.png)

![png](output_74_5.png)

![png](output_74_6.png)

![png](output_74_7.png)

![png](output_74_8.png)

![png](output_74_9.png)

![png](output_74_10.png)

![png](output_74_11.png)

![png](output_74_12.png)

![png](output_74_13.png)

![png](output_74_14.png)

![png](output_74_15.png)

![png](output_74_16.png)

![png](output_74_17.png)

![png](output_74_18.png)

![png](output_74_19.png)

![png](output_74_20.png)

![png](output_74_21.png)

![png](output_74_22.png)

![png](output_74_23.png)

![png](output_74_24.png)

![png](output_74_25.png)

![png](output_74_26.png)

![png](output_74_27.png)

# Multitaper Principal Component Analysis (PCA) or Singular Value decomposition (SVD)

If we want to compare more than two sites simultaneously, we can use PCA. For example, if we utilize all the data (for which we have a full week) together to project the complex version of the data onto a principal component, then what is the magnitude of that vector for each frequency?

In the PCA method, we denote by $x_{j,t}$ where $j = 0, \ldots, J-1$, $n = 0, \ldots, N-1$ the $j$th time series, where we have first subtracted the sample mean and divided by the sample variance of each time series. Compute first the Fourier transformed, windowed data as $$ X^{(j,k)}(f) = \sum_{n=0}^{N-1} x_{j,n} v_n^{(k)} e^{-2 \pi i f n} $$ and collect the $X^{(j,k)}(f)$'s into a $J \times K$ matrix ${X}(f)$. The spectral matrix ${X}(f)^T {X}(f)$ contains all of the spectra and cross spectra at the frequency $f$. We decompose the complex-valued ${X}(f)$ as $$ {X}(f) = {U}(f) {D}(f) {V}(f)$$ The complex-valued matrices ${U}(f)$ and ${V}(f)$ are $J \times J$ and $K\times K$, respectively and contains the left and right principal components of the matrix ${X}(f)$. The matrix ${D}(f)$ is a real-valued $J\times K$ diagonal matrix containing the principal components of the spectral matrix. The principal components are to be interpreted as one interprets multitaper spectrum estimates.

``` r
# First reshape a vector of data with indices in a n x 2 array into an N x n array
reshapedat <- function(vecdat,ind){
    # only works if all datasets have the same length
    Nf = size(ind,1)
    NN = array(0, Nf)
    for (i in seq(1, Nf)){
        NN[i] = length(vecdat[seq(ind[i, 1], ind[i, 2])])
        if ((i > 1)&&(NN[i] != NN[1])){
            stop("Cannot concatenate datasets with different lengths.")
        }
    }
    dat = array(0,c(NN[1], size(ind,1)))
    for (i in seq(1, Nf)){
        dat[,i] = vecdat[seq(ind[i, 1], ind[i, 2])]
    }
    dat
}

# These are the time series that are not missing any data in the middle
fulldat = c(1,2,3,4)
dat = reshapedat(pmn$TotPMN,mod_ind)
```

Now, compute the PCA

``` r
# Set some parameters for the PCA
nw = 1
k = 2

# Function that computes the MTM-PCA
mtmsvd <- function(dat,nw,k){
    nft = 500 # shortened to preserve memory
    ncol = size(dat,2)
    mts = array(0,c(nft,k,ncol))
    mtspec = array(0,c(nft,ncol))
    for (i in seq(1,ncol)){
        temp = spec.mtm(as.ts(dat[,i], deltat=60*24), 
                    returnInternals = TRUE, nw=nw, k = k, plot=FALSE)
        mts[,,i] = temp$mtm$eigenCoefs[1:nft,]
        mtspec[,i] = temp$spec[1:nft]
    }
    d = array(0, dim=c(nft, min(ncol,k)))

    # Compute the SVD at each frequency
    for (ff in seq(1, nft)){
      d[ff,] = svd(mts[ff,,])$d  
    }
    data.frame(freq = temp$freq[seq(1,nft)], spec = mtspec, sv = d)
}
```

Let's run the PCA analysis on the two groups: wet and dry seasons, and look at the first eigenvalue.

``` r
dat_cw = reshapedat(pmn$TotPMN, mod_ind[c(1,2,3,4),])
svd_cw = mtmsvd(dat_cw,nw,k)
                                   
dat_td = reshapedat(pmn$TotPMN, mod_ind[c(5,6,7,8),])
svd_td = mtmsvd(dat_td,nw,k)

singvals = rbind(data.frame(names = rep("Wet", each=500), svd_cw), 
                 data.frame(names = rep("Dry", each=500), svd_td))
head(singvals,5)
```

+---+---------+--------------+---------+-----------+----------+----------+----------+--------------+
|   | names   | freq         | spec.1  | spec.2    | spec.3   | spec.4   | sv.1     | sv.2         |
+===+=========+==============+=========+===========+==========+==========+==========+==============+
|   | \<chr\> | \<dbl\>      | \<dbl\> | \<dbl\>   | \<dbl\>  | \<dbl\>  | \<dbl\>  | \<dbl\>      |
+---+---------+--------------+---------+-----------+----------+----------+----------+--------------+
| 1 | Wet     | 0.000000e+00 | 7251770 | 36608.21  | 122802.4 | 322020.1 | 3964.173 | 1.127340e-10 |
+---+---------+--------------+---------+-----------+----------+----------+----------+--------------+
| 2 | Wet     | 6.103516e-05 | 5186772 | 55283.04  | 155121.1 | 289563.3 | 3390.295 | 7.436099e+01 |
+---+---------+--------------+---------+-----------+----------+----------+----------+--------------+
| 3 | Wet     | 1.220703e-04 | 3026404 | 77252.85  | 248570.7 | 306940.5 | 2644.328 | 4.834914e+02 |
+---+---------+--------------+---------+-----------+----------+----------+----------+--------------+
| 4 | Wet     | 1.831055e-04 | 3752447 | 239825.49 | 365694.7 | 413710.2 | 2852.084 | 1.165587e+03 |
+---+---------+--------------+---------+-----------+----------+----------+----------+--------------+
| 5 | Wet     | 2.441406e-04 | 3886662 | 745953.41 | 370645.1 | 562983.6 | 2909.275 | 1.632477e+03 |
+---+---------+--------------+---------+-----------+----------+----------+----------+--------------+

: A data.frame: 5 × 8

## Plot the PCA results

``` r
# Plot PCA for the wet season
p1 = ggplot(subset(singvals, names=="Wet"), aes(x = freq*60*24)) +
    geom_line(aes(y=sv.1), colour = "blue") + 
    geom_line(aes(y=sv.2), colour = "blue", lty=2) +
    scale_y_log10(name = "Principal component, Wet",limits = c(8e2,2e4)) +
    scale_x_continuous(name = "Frequency (cyc/d)", limits = c(0,13), breaks=seq(0,13, by=1)) +
    theme(legend.position="none") 

# Plot PCA for the dry season
p2 = ggplot(subset(singvals, names=="Dry"), aes(x = freq*60*24)) +
    geom_line(aes(y=sv.1)) + 
    geom_line(aes(y=sv.2), lty=2) + 
    scale_y_log10(name = "Principal component, Dry",limits = c(8e2,2e4)) +
    scale_x_continuous(name = "Frequency (cyc/d)", limits = c(0,13), breaks=seq(0,13, by=1)) + 
    theme(legend.position="none")

# Combine
pp = p1 + p2 + plot_layout(ncol = 1)

# Plot PCA for both seasons together for comparison
psingular = ggplot(singvals, aes(freq*60*24, sv.1)) +
    geom_line(aes(colour=factor(names))) + scale_colour_manual(values=c("black","blue")) +
    scale_y_log10(name = "First principal component",limits = c(8e2,2e4)) +
    scale_x_continuous(name = "Frequency (cyc/d)", limits = c(0,13), breaks=seq(0,13, by=1))   + 
    theme(legend.position="none")

# Set plot settings
options(repr.plot.width=10, repr.plot.height=5)

psval = pp | psingular 

# Save
 #ggsave(paste(figdir, "FirstSingvals_RD.png", sep =""), width=10, height=5)

psval
```

    Warning message:
    "[1m[22mRemoved 352 rows containing missing values (`geom_line()`)."
    Warning message:
    "[1m[22mRemoved 355 rows containing missing values (`geom_line()`)."
    Warning message:
    "[1m[22mRemoved 352 rows containing missing values (`geom_line()`)."
    Warning message:
    "[1m[22mRemoved 358 rows containing missing values (`geom_line()`)."
    Warning message:
    "[1m[22mRemoved 704 rows containing missing values (`geom_line()`)."

![png](output_82_1.png)

For the rainy season, there are peaks present for the first four cyc/d but no discernible pattern at higher harmonics. Singular values are marginally higher for 3 cyc/d and 4 cyc/d compared to 1 cyc/d and 2 cyc/d. In contrast, the SVD for the dry season shows peaks at 1 cyc/d, 3 cyc/d, 5 cyc/d, 7 cyc/d, 9 cyc/d, and 11 cyc/d, though the peak at 1 cyc/d is much greater than at other harmonics. When we overlay the first singular value for both seasons, we see the peak at 1 cyc/d (24 hours) is a factor of three greater in the dry season compared to the rainy season, whereas the singular value at 3 cyc/d (8 hours) are very similar. The peak at 2 cyc/d (12 hours) is missing from the dry season as the pseudo peak is less than $2W$ wide. So too the peak at 4 cyc/d (4 hours) present in the rainy season is absent from the dry season. As a whole, we observe that the dry season has a larger number of higher harmonics required to capture the variance reflective of the square nature of the waveforms.
