# BES Palaeo SIG: Masterclass in R
# "tapas, and R version of CharAnalysis"
# author: Walter Finsinger
# date: 2022
#
### Introduction --------------------------------------------------------------

# Hi,
# many thanks for your interest in this R package.
# The following R script is meant to illustrate a simple workflow to perform
# trend and peak-detection analysis of macroscopic charcoal records.

# The code builds on CharAnalysis (https://github.com/phiguera/CharAnalysis),
# a software for analyzing sediment-charcoal records written in and compiled
# with Matlab 7.0 by Phil Higuera (Higuera et al., 2009), with significant input
# by (amongst others) Patrick Bartlein (U of OR), Daniel Gavin (U of OR),
# Jennifer Marlon, and Ryan Kelly.
# Most functions that are included in the 'tapas' package were translated
# verbatim from CharAnalysis.


# In this tutorial, you'll be using the macroscopic charcoal record from
# Code Lake (Higuera et al., 2009). Please be aware that the tutorial shows
# how the tapas package performs the peak-detection analysis with this data set
# using different options than those chosen by Higuera et al. (2009).

# Specifically, the tutorial uses the so-called 'global' threshold rather than
# the so-called 'local threshold'. However, the tapas package also allows
# performing the anlysis with the so-called 'local' threshold, as does
# CharAnalysis.



# 1. Download, install and load the R 'tapas' package -------------------------

## Run these lines of code to install the latest build on your local computer.
# install.packages("devtools")
devtools::install_github("wfinsinger/tapas")



## Load packages into the local R Environment
library(tapas)   # Load the 'tapas' package
library(ggplot2) # additionally load the ggplot package



# 2. Load the toy data that comes with the 'tapas' package --------------------
# This is a dataframe with the macroscopic charcoal record from Code Lake
#       (Higuera et al., 2009), including:
#       - sample depths (cmTop and cmBot),
#       - sample ages (AgeTop and AgeBot), and
#       - the number of charcoal pieces in each sample (char).
# There is an additional variable (dummy_chararea) that was added to illustrate
#       the functioning of some of the tapas-package functions, which treat
#       several variables in one single run.
co <- tapas::co_char_data

# NB: if you want to run the analysis with your own data, you should load your
# data into the R Environment, for instance using 'read.csv()'.
# The format of the input data is described in the readme file available on
# https://github.com/wfinsinger/tapas, or can be viewed by running:
head(co, 13)


#### 2.1. Plot raw data
tapas::plot_raw(co)


# 3. Check the data  ----------------------------------------------------------
# Are there any gaps?
# contiguous sampling?
# Any duplicates in the sample-depth and sample-ages columns?

co <- tapas::check_pretreat(co)


co_short <- co[1:30, ]
tapas::plot_raw(co_short, bars = TRUE)

co_gaps <- co_short[-c(11, 12), ]
tapas::plot_raw(co_gaps, bars = TRUE)

co_gaps_checked <- tapas::check_pretreat(co_gaps)
head(co_gaps_checked, 13)
head(co_gaps, 13)
plot_raw(co_gaps_checked, bars = TRUE)
rm(co_gaps, co_gaps_checked)


co_dupl <- co_short
co_dupl$cmTop[1] <- co_dupl$cmTop[2]
head(co_dupl)
tapas::check_pretreat(co_dupl)
rm(co_dupl)
rm(co_short)


# 4. Bin the data (re-sample the data to equal sampling intervals) ------------

#### 4.1. Bin the data using the default options ------------------------------

co_i <- tapas::pretreatment_data(co, yrInterp = NULL, series.name = "CO")

# Question: Are default options ok?
# Specifically: the argument 'yrInterp', which defines the resolution of the
#    binned series (by default, yrInterp = the median sampling resolution).
# To check the value that was chosen for the 'yrInterp' argument:
co_i$int$yr.interp



#### 4.2. Explore the sampling resolution of the record -----------------------

# To plot the depth-age relationship,
# first set x-axis and y-axis limits, and the direction of the axes:
x_lim <- c(max(co$AgeTop), min(co$AgeTop))
y_lim <- c(max(co$cmTop), min(co$cmTop))

# then plot the depth-age relationship:
plot(co$AgeTop, co$cmTop, type = "l",
     xlim = x_lim, ylim = y_lim,
     xlab = "Age (cal BP; AgeTop)",
     ylab = "Depth (cm; cmTop)")


# Plot the raw data as simple bar plot: counts against age
plot_raw(co, bars = TRUE)


# Calculate and plot the sampling resolution (sample integration time)
# for the record (yr/sample),
# and the median sampling resolution (red dashed line)
yr_smpl <- co$AgeBot - co$AgeTop

plot(co$AgeTop, yr_smpl, type = "l",
     xlab = "Age (cal BP)",
     ylab = "Sample integration time\n(yr/sample)")
abline(h = median(yr_smpl), col = "red", lty = 2)

#### Check the distribution of sampling resolution values
boxplot(yr_smpl, ylab = "Sample integration time\n(yr/sample)")
summary(yr_smpl)

rm(yr_smpl)


#### 4.3. Bin the data --------------------------------------------------------
# Given the results above, one may set, for instance, yrInterp = 16 years.
co_i <- tapas::pretreatment_data(series = co, out = "accI",
                                 first = -51, last = 7500,
                                 yrInterp = 16)


# 5. Detrend the data ---------------------------------------------------------
## To decompose the 'background' from the peaks

co_detr <- tapas::SeriesDetrend(co_i, detr.type = "mov.median",
                                smoothing.yr = 500)
# The SeriesDetrend() function sends plots (one for each variable included in
# the input data.frame) to the device.
# To save plots directly to the hard disk, one can use the pdf() and dev.off()
# base R functions, and specify a sub-directory (here the folder "figures")
pdf("./figures/co_i.pdf")
co_detr <- tapas::SeriesDetrend(co_i, detr.type = "mov.median",
                                smoothing.yr = 500)
dev.off()


# 6. Decompose noise & potential peaks ----------------------------------------
## To remove the 'noise' from the peaks
co_glob <- tapas::global_thresh(co_detr, proxy = "charAR",
                                thresh.value = 0.95)

# The global_thresh() function produces two figures.
# The first one shows the results of the classification, which is based on a
# Gaussian Mixture Model analysis.
# The second figure shows the detrended series with the thresholds
# and the peaks (grey circles: peaks that did not pass the minimum-count test,
# red crosses: those that passed the test).


# 7. Evaluate results ---------------------------------------------------------
# - check SNI plot produced with previous function (the 2nd figure)
#
# NB: SNI is calculated on moving windows (as of the 'smoothing.yr' argument).
# It quantifies the separation of the peaks from the ‘noise’.
# It tells if a record that was treated with a given set of user-determined
# parameters is suitable for peak-detection analysis.
# SNI > 3 : samples of Signal are on average > 3 SD above the mean of Noise



# 8. Evaluate sensitivity to smoothing-window widths --------------------------
# using the peak_detection() wrapper function (which runs the functions used in
# steps 3.-6. in one go):
co_glob2 <- tapas::peak_detection(series = co, proxy = "char",
                                  first = -51, last = 7500,
                                  yrInterp = 16,
                                  thresh_type = "global",
                                  detr_type = "mov.median",
                                  smoothing_yr = 500,
                                  thresh_value = 0.95, min_CountP = 0.05,
                                  sens = T,
                                  smoothing_yr_seq = c(500, 600, 700,
                                                       800, 900, 1000))
# and check the 2nd to last plot...

# 9. ...and if you are happy with it, plot results & export data -------
par(mfrow = c(2,1))
tapas::Plot.Anomalies(co_glob, plot.neg = FALSE)
tapas::Plot_ReturnIntervals(co_glob, plot.x = TRUE)

# extract data from the package's output as a data.frame.
co_glob_exp <- tapas::tapas_export(co_glob)

# NB: This can then be saved as *.csv file to the local hard disk
# with write.csv().
# For instance,
write.csv(co_glob_exp, "./data_out/co_glob_exp.csv", row.names = FALSE)





# 10. Influence of the 'yrInterp' argument on the SNI? ------------------------

# Here the code will loop analyses with different yrInterp values, then gather
# data obtained in each of the loops, and finally plot diagnostic figures:

# First define which binning-year values you want to use
values_want <- c(5, 15, 30, 50, 100)

# prepare empty lists where data will be stored
a_list <- list()
co_loc_list <- list()

# Run the peak_detection() analysis with different yrInterp values and store
#       the results of each loop in the lists:
for (i in 1:length(values_want)) {
        values_i <- values_want[i]
        co_loc_i <- tapas::peak_detection(series = co, proxy = "char",
                                          first = -51, last = 7500,
                                          yrInterp = values_i,
                                          thresh_type = "global",
                                          detr_type = "mov.median",
                                          smoothing_yr = 500,
                                          thresh_value = 0.95,
                                          sens = FALSE, plotit = FALSE)
        co_loc_list[[i]] <- co_loc_i  # stores output of analysis

        # Write the summary data.frame generated by tapas_export() in the list
        co_loc_i_exp <- tapas::tapas_export(co_loc_i)
        co_loc_i_exp$values_i <- values_i
        a_list[[i]] <- co_loc_i_exp
}

# Clean environment
rm(co_loc_i, co_loc_i_exp)

# Gather the data from the list into one data.frame
arg_sens <- dplyr::bind_rows(a_list)

# Plot results
ggplot(data = arg_sens, aes(x = age_top_i, y = sni_smooth)) +
        geom_line(aes(group = values_i, colour = factor(values_i))) +
        scale_x_reverse() + #scale_y_continuous(limits = c(0, 20)) +
        geom_hline(yintercept = 3) +
        ggtitle(label = "Code Lake - SNI records for different yrInterp",
                subtitle = "local GMMs, mov.median 500yrs, thresh.value 0.95")

# Question: what emerges from this comparison ?



# We can also visually compare the results
par(mfrow = c(length(values_want), 1))
for (i in 1:length(values_want)) {
        tapas::Plot.Anomalies(co_loc_list[[i]], plot.neg = FALSE)
        mtext(paste("thresh.value = ", values_want[i]), side = 3)
}


# 11. - n: Assignment ---------------------------------------------------------
# You may copy/paste the code under
# paragraph "10. Influence of the 'yrInterp'...",
# and modify it to explore the influence of other user-determined choices,
# such as:


#### 11.1 The threshold value (thresh.value) ----------------------------------

# Here the code will loop analyses with different yrInterp values, then gather
# data obtained in each of the loops, and finally plot diagnostic figures:

# First define which binning-year values you want to use
values_want <- c(0.95, 0.99, 0.999)

# prepare empty lists where data will be stored
a_list <- list()
co_loc_list <- list()

# Run the peak_detection() analysis with different yrInterp values and store
#       the results of each loop in the lists:
for (i in 1:length(values_want)) {
        values_i <- values_want[i]
        co_loc_i <- tapas::peak_detection(series = co, proxy = "char",
                                          first = -51, last = 7500,
                                          yrInterp = 16,
                                          thresh_type = "global",
                                          detr_type = "mov.median",
                                          smoothing_yr = 500,
                                          thresh_value = values_i,
                                          sens = FALSE, plotit = FALSE)
        co_loc_list[[i]] <- co_loc_i  # stores output of analysis

        # Write the summary data.frame generated by tapas_export() in the list
        co_loc_i_exp <- tapas::tapas_export(co_loc_i)
        co_loc_i_exp$values_i <- values_i
        a_list[[i]] <- co_loc_i_exp
}

# Clean environment
rm(co_loc_i, co_loc_i_exp)

# Gather the data from the list into one data.frame
arg_sens <- dplyr::bind_rows(a_list)

# Plot results
ggplot(data = arg_sens, aes(x = age_top_i, y = sni_smooth)) +
        geom_line(aes(group = values_i, colour = factor(values_i))) +
        scale_x_reverse() + #scale_y_continuous(limits = c(0, 20)) +
        geom_hline(yintercept = 3) +
        ggtitle(label = "Code Lake - SNI records for different yrInterp",
                subtitle = "local GMMs, mov.median 500yrs")

# Question: what emerges from this comparison ?



# We can also visually compare the results
par(mfrow = c(length(values_want), 1))
for (i in 1:length(values_want)) {
        tapas::Plot.Anomalies(co_loc_list[[i]], plot.neg = FALSE)
        mtext(paste("thresh.value = ", values_want[i]), side = 3)
}




#### 11.2 The minimum-count test (min_CountP) ---------------------------------





#### 11.3 The detrending method (detr_type) -----------------------------------










# 12. Model the 'background' trend with a Generalised Additive Model (GAM) ----
#
# NB:
# The functions that are illustrated here are still at a 'beta' stage.
# For more details on GAMs you can refer to the following paper:
# Simpson GL (2018) Modelling Palaeoecological Time Series Using Generalised
#  Additive Models. Frontiers in Ecology and Evolution 6: 149:
#  doi:10.3389/fevo.2018.00149.


#### 12.1. With binned data ---------------------------------------------------

#### Reshape the binned data that is stored in a list such that the
## 'mgcv' package can read it:
co_i_mgcv <- tapas::tapas2mgcv(series = co_i)

#### Determine trend with a GAM and let the gam() function choose the degree
## of the whiggliness of the smoothed curve:
co_i_gam1 <- mgcv::gam(charAR ~ s(age_top),
                       data = co_i_mgcv,
                       family = gaussian(link = "identity"),
                       method = "REML")

#### Check GAM model
par(mfrow = c(4,1), mar = c(2,5,2,2))
mgcv::gam.check(co_i_gam1)
summary(co_i_gam1)

par(mfrow = c(2,1), mar = c(2,5,2,2))
mgcv::plot.gam(co_i_gam1, shade = T)
plot(co_i_mgcv$age_top, co_i_mgcv$charAR, type = "l")
lines(co_i_mgcv$age_top, co_i_gam1$fitted.values, col = "red", lwd = 2)


#### Reshape the GAM output for tapas:: and detrend data
co_i_gam1_detr <- tapas::mgcv2tapas(series = co_i_gam1,
                                    data_type = "accI")

#### and finally move on to the 'threshold analysis'
co_gam_thresh <- global_thresh(co_i_gam1_detr, proxy = "charAR",
                               smoothing.yr = 500)
tapas::Plot.Anomalies(co_gam_thresh, plot.neg = F, plot.x = T)




#### 12.2. With non-binned data -----------------------------------------------

## Calculate the sample interval (= number of years included in each sample)
depostion_time <- co[ ,4] - co[ ,3]

## Calculate charcoal-accumulation rate manually as we don't use the
## pretreatment_data() function: (AR = Concentration / deposition_time;
##  AR = (counts / sample_volume) / deposition_time).
co$charAR <- (co[ ,6] / co[ ,5]) / depostion_time

# Determine trend with a GAM
d_gam1 <- mgcv::gam(charAR ~ s(AgeTop, k = 20), data = co,
                    family = gaussian(link = "identity"),
                    method = "REML",
                    weights = depostion_time / mean(depostion_time))

# Check GAM model
mgcv::gam.check(d_gam1)
summary(d_gam1)
par(mfrow = c(2,1), mar = c(2,5,2,2))
mgcv::plot.gam(d_gam1, n = 200, shade = T)
plot(co$AgeTop, co$charAR, type = "l")
lines(co$AgeTop, d_gam1$fitted.values, col = "red", lwd = 2)


# Reshape GAM output for tapas:: and detrend data
d_gam1_detr <- mgcv2tapas(series = d_gam1, data_type = "accI",
                          series.name = "CO")

rm(d_gam1, d_gam1_detr)
