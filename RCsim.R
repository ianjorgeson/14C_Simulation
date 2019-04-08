################################################################################
###----------- SIMULATED RADIOCARBON SAMPLING DISTRIBUTIONS FOR ------------###
###------------- MULTIPLE RADIOCARBON DATES OF A SINGLE EVENT -------------###
#############################################################################

# Last edited November 19, 2018, R v.3.5.1
# Ryan Breslawski, rbreslawski@smu.edu

# **INTRODUCTION**
# ****************
# This script simulates radiocarbon dates from synchronous events over a user-
# specified range of years. It evaluates the probability of observing a set of 
# radiocarbon dates with at least as much dispersion as that observed in a real,
# user supplied sample of radiocarbon dates. Due to the potentially intensive
# nature of these simulations, the script performs this task on a cluster. After
# about 20-40 seconds, it will generate daughter scripts to be executed on cluster
# nodes. This script will pause and wait for the daughter scripts to complete
# before it finishes. After the daughter scripts are generated, a shell script
# is also generated that the user can execute at the command line to run every
# daughter script. Currently, this shell script is specific to the Slurm workload
# manager as set up on the ManeFrame cluster at Southern Methodist University. If
# running this code on a different cluster, you may wish to modify the internal
# code here to specify a shell command that will execute the daughter scripts, or,
# you may simply type the necessary shell commands rather than run the output
# script. If the daughter scripts are never executed, this script will wait
# indefinitely, as it requires RData output from each daughter script to complete.
# 
# Before running the script, read the 'USER ARGUMENTS' section and modify as 
# necessary to set the simulation parameters, plotting options, and cluster settings.
# 
# It is recommended that this script is sourced non-interactively. After beginning the
# script, monitor your working directory for the daughter scripts, which should appear
# after 20-40 seconds. Run these daughter scripts, one per node, to generate the data
# that. The main script will aggregate the data into list objects and plots.
# 
# Simulation execution time varies dramatically based on the number of iterations.
# Ideally, the cluster is used such that every year is assigned to a core. In this
# scenario, it takes roughly 1.5 days for the simulation to complete 10,000 iterations
# on the SMU ManeFrame cluster. 1000 iterations takes about 4 hours and 100 iterations
# completes in under half an hour. Avoid using fewer cores than the number of years.
# If any core has to complete multiple years, this serializes part of the simulation.
# In this scenario, multiply the estimated completion times by the maximum number of
# years to by handled by any single core.
# 
# 
# **REQUIREMENTS**
# ****************
# -- Libraries (see below 'LOAD LIBRARIES' section)
# -- Rstan for first run of simulation; used for inter- and intra-laboratory
#    error model.
# -- csv files
#    - ydb.csv: csv file of sample radiocarbon dates, see 'USER ARGUMENTS' for 
#      formatting details.
#    - IRI.csv: Only required for the first run, which will create an RData file
#      in the working directory with simulation parameter values ('IRI#.RData', 
#      where # is the values of setsn). This RData file will be called in the 
#      simulation. If it already exists in the working directory, the script 
#      will not attempt to recreate it. This csv consists of results from the 
#      International Radiocarbon Intercomparison, and it is used to fit a model 
#      that describes inter- and intra-laboratory error in the measurement of 
#      C-14 values.
#      
##--- LOAD LIBRARIES ---##
# In CRAN:
library(ggplot2)
library(parallel)
library(reshape2)
library(rcarbon)
library(matrixStats)
# Not in CRAN:
library(patchwork) # Install with devtools: 
# "devtools::install_github("thomasp85/patchwork")"
library(rethinking) # Install with devtools: 
# "devtools::install_github("rmcelreath/rethinking")" 
#     
# 
# **OUTPUT FILES**
# ****************
# -- NodeSim#.R: 
#  Daughter script that simulates over a chunk of the year range on a node. # is 
#  an integer. The number of daughter scripts depends on the number of nodes.
# -- NodeDat#.RData: 
#  Output from NodeSim#.R. Aggregated by the main script. The number of daughter
#  scripts determines the number of data files.
# -- nodesim.sh:
#  Shell script that submits the daughter scripts to the SMU ManeFrame cluster.
# -- IRI#.RData (where # = setsn):
#  Data file containing samples from the laboratory bias and repetability model.
#  Required for the simulation to run. If this file already exists in the working
#  directory, the script will not recreate it.
# -- SimDatIntermediate.RData:
#  Simulation parameters generated from the main script and used in each daughter
#  script. Do not delete this file prior to running the daughter scripts.
# -- SimDat#.RData:
#  Aggregated results once all nodes have completed. # is the number of user-
#  specified iterations.
#  
#  
# ** SIMULATION RESULTS OBJECTS STORED IN THE R ENVIRONMENT **
# ***********************************************************
# These are the objects saved in SimDat#.RData. 
# 
# Data frames and lists:
#   datelist: 
#     Data frame of sample dates imported from a csv file.
#   event_index_list: 
#     Vectors of event index values for the daughter scripts run by each 
#     node during the simulation.
#   FailedIterations:
#     This data frame is only generated if the calibrate function attempts
#     to calibrate a 14C value that is too old for IntCal13 during the
#     simulation. This is an improbable event for most simulations, and in
#     most cases this data frame will not be generated for simulations over
#     date ranges within the most recent 40,000 calendar years. If this data
#     frame is generated, it indicates that a date too old for calibration
#     was attempted, leading to a failed iteration. The data frame stores
#     the information for this failed iteration. The failed iteration will
#     be replaced with a random successful iteration drawn from the same
#     simulation (i.e., using the same old wood and laboratory model
#     parameters).
#   IRI_diagnostics:
#     RStan diagnostics for the laboratory bias and repeatability model.
#   IRI_model:
#     map2stan model of inter- and intra-lab error (see Rethinking
#     documentation for details of object structure)
#   IRI_dates:
#     Data frame of International Radiocarbon Intercomparison data,
#     imported from csv file.
#   IRIplots:
#     ggplot2 object of IRI data and posterior predictions from IRI_model.
#     Each panel is a sample, with labs listed along the y-axis. Black
#     segments show measured C-14 means with one standard deviation. Orange
#     bars show 95% prediction intervals for AMS labs, and purple bars show
#     95% prediction intervals for non AMS labs. Red lines and bands show the
#     posterior mean and 95% highest density interval for the sample mean.
#   List_observed_master:
#     List of [1] observed C-14 standard deviations and [2] calibrated age
#     distances for each event. List contains as many positions as there
#     are events.
#   List_plots_ages:
#     List of ggplot2 objects illustrating calibrated dates, with one object
#     per event. Grey bands mark hypothesized date ranges for the event.
#   List_plots_simulation:
#     List of ggplot2 objects with panels showing the simulated results, with
#     one object per event. 
#   List_stats_master:
#     Simulated results, with events indexed by list position. Within each 
#     indexed event:
#       -Stats (raw simulated values)
#         --Uncalibrated_mu_sd
#           ---No_lab_error
#             [[x positions]] Matrices of results for x outlier models, rows
#             are simulated years and columns are simulation iterations.
#           ---With_lab_error
#             [[x positions]] Matrices of results for x outlier models, rows
#             are simulated years and columns are simulation iterations.
#         --Calibrated_mu_dist
#           ---No_lab_error
#             [[x positions]] Matrices of results for x outlier models, rows
#             are simulated years and columns are simulation iterations.
#           ---With_lab_error
#             [[x positions]] Matrices of results for x outlier models, rows
#             are simulated years and columns are simulation iterations.
#       -P_vals (proportions of sim iterations exceeding observed values)
#         --Uncalibrated_mu_sd
#           ---No_lab_error: Matrix of results, rows are simulated years and 
#             columns are different outlier models.
#           ---With_lab_error: Matrix of results, rows are simulated years and 
#             columns are different outlier models.
#         --Calibrated_mu_dist
#           ---No_lab_error: Matrix of results, rows are simulated years and 
#             columns are different outlier models.
#           ---With_lab_error: Matrix of results, rows are simulated years and 
#             columns are different outlier models.
#  List_years_master: Vectors of years for each event to be simulated. List
#    positions correspond to each event.
#  ml: List of event specific parameters for the simulation, list positions
#    correspond to each event.
#  Plot_YDB_LS_z: GGplot2 object that displays the the distance of the observed
#    Manhattan distances (logit transformed) and C-14 stdevs (log transformed)
#    from their respective simulated logit and log mean values, as represented
#    by z-scores. For the Younger Dryas Boundary (YDB) and Laacher See Tephra 
#    (LST). 
#  Plot_YDB_LS_sim_D: GGplot2 object that displays simulated Manhattan distance
#    results for the Younger Dryas Boundary and Laacher See Tephra.
#  Plot_YDB_LS_sim_sd: GGplot2 object that displays simulated standard deviation
#    results for the Younger Dryas Boundary and Laacher See Tephra.
#  SimList: List of simulation results matrices aggregated from the daughter 
#    R scripts. These data are reformatted into the List_stats_master list.
#  Years_list: List of year vectors, where each list position is read by a
#    a different daughter R script for the simulation.
#    
# Values and vectors:
#   CoreNum: Number of cores per cluster node.
#   DataFiles: Filenames for nodes pecific RData output.
#   event_index: Vector of event index values for the concatenated vector 
#     years across events.
#   events: Names of each event to be simulated.
#   IRIsd: Standard deviation of centered dates reported in the International
#     Radiocarbon Intercomparison study.
#   mpt: Minimum probability thresholds for radiocarbon calibration. Used in 
#     the rcarbon calibration function.
#   nNodes: Number of cluster nodes.
#   OM_lambda: Lambda parameter values for outlier models.
#   oml: Number of outlier models (excluding the default no outlier model)
#   Plot_posterior: Variable controlling whether the posterior distribution
#     is plotted in figures displayed calibrated observed dates.
#   setsn: Number of iterations for each simulated years.
#   ybuff: Buffer, in years, on each side of the hypothesized date ranges
#     for the simulation.
  


# Fix randomization for repeatability
set.seed(100)

######################################################################
########---------------- USER ARGUMENTS ----------------#############
####################################################################

# Read csv table of radiocarbon data. Must have following headers:
# Site: Name of site. String.
# lID: Laboratory sample ID (e.g., "Beta-1234"). Must contain a "-" that separates the
#      laboratory code from the sample number. Laboratory codes must be consistent
#      between entries for valid results in the simulations that include modelled intra-
#      and inter-laboratory error. String containing "-".
# Mus: Reported C-14 means. Numeric, only positive values within domain of IntCal13.
# SEs: Reported C-14 errors, in 1 standard deviation. Numeric, only positive values.
# Event: Name of the hypothesized synchronous event. All samples associated with the same
#     event must have the same name. String.
# Method: AMS (0) or non AMS (1). This is used in the intra- and inter-laboratory error 
#     simulations. Numeric, only 0 and 1.
# Outlier: No (0) or yes (1). If the sample may be associated with an event that predates
#    the hypothesized synchronous event (e.g., old wood), enter 1. This will be used in
#    the outlier model simulations. Numeric, only 0 and 1.
datelist <- read.csv("ydb.csv", stringsAsFactors=FALSE, header=TRUE)

# Number of sets of radiocarbon dates to sample for each simulated year. The 
# default value of 10 will only provide very noisy estimates of the simulated
# distributions. A value of 1000+ is preferable.
# Default: 1e4 sets
setsn <- 1e4

# Ranges of years for hypothesized events. Each event should have a start
# and end date. For multiple events, start and end should be vectors of
# dates. Default values are set for two events: YDB and Laacher See.
Range_years <- data.frame(start=c(12835, 12966), end=c(12735, 12866))

# Buffer each hypothesized range with a fixed number of years. This will
# add ybuff (a positive integer) years to either side of each hypothesized
# range. Set to zero if you only want to simulate over the hypothesized
# ranges.
# Default: 25
ybuff <- 25

# Plot posterior date for the observed sample of dates? Yes = "y". If
# "y" is marked, the posterior distribution will be displayed on the
# plot of calibrated dates. This posterior will only display if a 
# valid posterior can be obtained from the date sample under the user
# defined mpt threshold (see next argument).
# Default value: "n"
Plot_posterior <- "n"

# Minimum probabilities for radiocarbon calibration. This is the minimum
# threshold for calendar year probabilities generated through the calibrate
# function. If the current value does not generate a valid posterior, try 
# reducing it in increments x of 1e-x. It is generally recommended that 
# minimum probabilities remain below 1e-4. The first position of the mpt 
# variable is for plotting the observed dates and their posterior. The 
# second vector position is for the simulated distributions.
# Default values: 1e-5, 1e-5
mpt <- c(1e-5, 1e-5)

# Set lambda parameters for exponential distributions used in outlier models.
# The simulation will automatically generate a non-outlier version. If you
# do not wish to view outlier model output, replace the vector with 'NA'. For
# practical purposes, it is recommended that no more than three values are used.
# Smaller values will raise the probability of older dates. To get a graphical
# idea of the amount of potential offset for each lambda value, you can run
# "curve(dexp(x,lambda), from=0, to=300)" on your R console, replacing 'lambda'
# with the parameter value of interest. Values on the x-axis represent the number
# of calendar years by which a sample predates the event of interest, while the
# y-axis shows the probability density across these years.
# Default values: 0.04, 0.01
# For these defaults, outlier offsets have the following expected values:
# Lambda = 0.04 -- mean=25, 95% most probable=0-75
# Lambda = 0.01 -- mean=100, 95% most probable=0-300
OM_lambda <- c(0.04, 0.01)


# ** CLUSTER ARGUMENTS**

# Define number of cores per node. Optimally, one core should be used per 
# simulated year
# Default value: 34
CoreNum <- 34

# Define numer of nodes, must be a positive integer greater than 1. This will
# set the number of daughter scripts, where the range of simulated years is
# chunked up over nodes.
# Default value: 9
nNodes <- 9




#######################################################################
###- FIT LABORATORY MEASUREMENT BIAS AND REPEATABILITY MODEL -########
#####################################################################

# Fit model if it does not already exist in the working directory.
# Otherwise, read model parameters from working directory.
if(!file.exists(paste0("IRI",setsn,".RData"))){
    
  # Read data
  IRIdates <- read.csv("IRI.csv", header=TRUE, stringsAsFactors=FALSE)
  
  # Create dummy variable for counting and AMS methods, where
  # AMS is the default value (i.e., AMS = 0).
  IRIdates$method_d <- ifelse(IRIdates$method=="AMS", 0, 1)
  
  # Center sample means within samples
  IRIdates$mean_d <- sapply(1:nrow(IRIdates), function(x){
    return(IRIdates$s_mean[x]-
             median(IRIdates$s_mean[which(IRIdates$sample==IRIdates$sample[x])]))
  })
  
  # Remove labs that generated outlier dates, here defined as dates that are
  # more than 5X the IQR distance from either quartile
  IRI_OutlierLabs <- sapply(unique(IRIdates$lab), function(x){
    bounds <- c((quantile(IRIdates$mean_d)[2]-6*IQR(IRIdates$mean_d)),
                (quantile(IRIdates$mean_d)[4]+6*IQR(IRIdates$mean_d)))
    lab_dates <- IRIdates$mean_d[which(IRIdates$lab==x)]
    ifelse(min(lab_dates) < bounds[1] | max(lab_dates) > bounds[2],
           return(x), return(NA))
  })
  IRI_OutlierLabs <- IRI_OutlierLabs[!is.na(IRI_OutlierLabs)]
  IRIdates <- IRIdates[!(IRIdates$lab %in% IRI_OutlierLabs),]
  
  # Create indices for labs
  IRIdates$lab_i <- sapply(IRIdates$lab, function(x){
    which(sort(unique(IRIdates$lab))==x)})
  # Create indices for samples
  IRIdates$sample_j <- sapply(IRIdates$sample, function(x){
    which(sort(unique(IRIdates$sample))==x)})
  
  # Create z-scores for centered sample means
  IRIdates$mean_z <- (IRIdates$mean_d-mean(IRIdates$mean_d))/sd(IRIdates$mean_d)
  # Log-log scale measurement errors
  IRIdates$logsd <- log(IRIdates$s_sd)
  
  # Set iteration values based on number of available threads. Total non-
  # warmup iterations is set to 10000 summed across HMC chains.
  cl <- detectCores()
  if(cl <= 4){
    chain_count <- 4
    c_iter <- 7500
  }else if(cl < 8){
    chain_count <- 5
    c_iter <- 7000
  }else{
    chain_count <- 8
    c_iter <- 6250
  }
  
  # Fit a multilevel model to the date sample.
  IRI_model <- map2stan(alist(
    
    # centered dates are normally distributed with mean 'mu' 
    # and std dev 'sigma'.
    mean_z ~ dnorm(mu, sigma),
    
    # Where mean 'mu' varies by sample and by systematic lab error.
    # The baseline for lab error is AMS dates. For non-AMS dates,
    # a multiplier term is included to adjust the error ('method_offset').
    # This multiplier is strictly positive, existing in (0,Inf).
    # In other words, it is possible for 'non AMS' to be an parameter 
    # value that either reduces or increases the spread of interlab 
    # dates. The expression of this parameter is controlled by method_d,
    # which is a dummy variable (AMS = 0, non AMS = 1).
    mu <- sample_mean[sample_j] + 
      lab_offset[lab_i]*(1+(method_offset-1)*method_d),
    
    # Prior for non AMS multiplier parameter. It is centered on 1,
    # which cancels out the multiplier parameter (1-1=0), corresponding
    # to non AMS labs with means that are distributed identically
    # to AMS labs. Hence, this functions as a regularizing prior
    # for the 'method_offset' parameter.
    method_offset ~ dgamma2(1, 0.5),
    
    # Prior for sample mean. Since within sample values are
    # centered, sample means are likely close to zero.
    sample_mean[sample_j] ~ dnorm(0, 1),
    
    # std dev 'sigma' varies by AMS/non AMS and by measurement
    # error. Since the model fits mean values between labs,
    # this linear model corresponds to within lab variability
    # conditioned by reported measurement error and AMS vs non 
    # AMS. The intercept 'a' can be interpreted as the global
    # within lab variability when reported measurement
    # error is zero and with the method as AMS. a_l is an
    # intercept offset that varies by lab i. Coefficients 
    # should be interpreted on the log scale, as a log link 
    # was selected to keep sigma positive.
    log(sigma) <- a + a_l[lab_i] + bSD*logsd + bM*method_d,
    
    # Priors for method and reported measurement error parameters
    # within sigma model.
    c(bM, bSD) ~ dnorm(0, 1),
    # Prior for global within lab variability
    a ~ dnorm(0, 1),
    
    # Priors for lab parameters
    lab_offset[lab_i] ~ dnorm(0, lsigma1),
    a_l[lab_i] ~ dnorm(0, lsigma2),
    c(lsigma1, lsigma2) ~ dexp(2)),
    
    # Stan arguments
    data=IRIdates, chains=chain_count, iter=c_iter, 
    warmup=5e3, cores=chain_count, 
    control=list(adapt_delta=0.99, max_treedepth=12),
    constraints=list(sigma="lower=0", lsigma1="lower=0",
                     lsigma2="lower=0", method_offset="lower=0"))
  
  # Save HMC diagnostics
  IRI_diagnostics <- precis(IRI_model, depth=3, prob=0.95)
  
  # Export sd of IRI mean date distances
  IRIsd <- sd(IRIdates$mean_d)
  
  # Jitter lab indices for plotting
  IRIdates$yjit <- IRIdates$lab - rep(0.5, nrow(IRIdates)) + 
    rbeta(n=nrow(IRIdates), 5, 5)
  
  ### EXTRACT POSTERIOR SAMPLES FROM LABORATORY MODEL PARAMETERS ###
  # For small iteration numbers, draw a minimum of 1000 samples,
  # otherwise, draw as many samples as there are iterations.
  ifelse(setsn < 1e3, nsamples <- 1e3, nsamples <- setsn)
  IRI_model_samples <- extract.samples(IRI_model, n=nsamples)
  
  # IRI dates plot
  for(j in 1:max(IRIdates$sample_j)){
    
    IRIsub <- IRIdates[which(IRIdates$sample_j==j),]
    
    # Data frame of lab specific parameters for sample j,
    # observed measurement error is the mean of observed
    # errors when labs have multiple dates.
    sub_data <- data.frame(lab_i=unique(IRIsub$lab_i))
    sub_data$lab <- sapply(sub_data$lab_i, function(x){
      IRIsub$lab[which(IRIsub$lab_i==x)][1]})
    sub_data$method_d <- sapply(sub_data$lab_i, function(x){
      IRIsub$method_d[which(IRIsub$lab_i==x)][1]})
    sub_data$logsd <- sapply(sub_data$lab_i, function(x){
      log(mean(IRIsub$s_sd[which(IRIsub$lab_i==x)]))})
    sub_data$sample_j <- rep(j, nrow(sub_data))
    
    # Generate posterior predictions for each lab in sub_data
    lab_pp <- t(sapply(1:nrow(sub_data), function(x){
      
      # Lab ID position
      lid <- sub_data$lab_i[x]
      
      # Mean lab value (mean of sample j + lab i offset)
      mlv <- IRI_model_samples$lab_offset[,lid] + 
        IRI_model_samples$sample_mean[,j]
      # Adjust lab values based on AMS or non AMS method
      mlv <- mlv*(1+(IRI_model_samples$method_offset-1)*
                    sub_data$method_d[x])
      
      # Standard deviation; global + lab i offset + measurement
      # error term + AMS vs non AMS term. Put on exponential scale
      sdlv <- exp(IRI_model_samples$a + 
                    IRI_model_samples$a_l[,lid] + 
                    IRI_model_samples$bSD*sub_data$logsd[x] + 
                    IRI_model_samples$bM*sub_data$method_d[x])
      
      # Sample predictions from mlv and sdlv
      pp <- rnorm(length(mlv), mean=mlv, sd=sdlv)
      
      # Put back on scale of radiocarbon years from z-scores and
      # Add the median value to uncenter the values
      pp <- (IRIsd*pp) + median(IRIsub$s_mean)
      
      # Return posterior predictions
      return(pp)
      
    }))
    
    # Get 95% intervals of posterior predictions for each lab
    sub_data$pp_l <- sapply(1:nrow(sub_data), function(x){
      quantile(lab_pp[x,], probs=0.025)})
    sub_data$pp_u <- sapply(1:nrow(sub_data), function(x){
      quantile(lab_pp[x,], probs=0.975)})
    
    # Obtain posterior samples from sample j
    sample_j_samples <- (IRI_model_samples$sample_mean[,j]*IRIsd) +
      median(IRIsub$s_mean)
    sample_j_mean <- mean(sample_j_samples)
    sample_j_interval <- quantile(sample_j_samples, probs=c(0.025, 0.975))
    
    # Set x-axis position for sample ID text on each panel
    sxmin <- min(c(min(IRIsub$s_mean)-
                     IRIsub$s_sd[which.min(IRIsub$s_mean)],
                   min(sub_data$pp_l)))
    sxmax <- max(c(max(IRIsub$s_mean)+
                     IRIsub$s_sd[which.max(IRIsub$s_mean)],
                   max(sub_data$pp_u)))
    sx <- quantile(seq(sxmin, sxmax, length.out=10), p=0.92)
    
    # Plot data and posterior predictions for sample j
    IRIplot <- ggplot(sub_data, aes(group=factor(lab_i)))+
      annotate("rect", xmin=sample_j_interval[1], 
               xmax=sample_j_interval[2], ymin=-Inf, ymax=Inf,
               alpha=0.3, fill="red")+
      annotate("segment", y=Inf, yend=-Inf, x=sample_j_mean,
               xend=sample_j_mean, color="red", size=0.8)+
      annotate("segment", y=seq(1, max(IRIdates$lab)), 
               yend=seq(1, max(IRIdates$lab)),
               x=Inf, xend=-Inf, color="gray88")+
      geom_rect(aes(xmin=pp_l, xmax=pp_u, ymin=lab-0.5,
                    ymax=lab+0.5, fill=factor(method_d)),
                alpha=0.5)+
      geom_segment(data=IRIsub, aes(y=yjit, yend=yjit, 
                                    x=s_mean-s_sd, xend=s_mean+s_sd), size=0.4)+
      geom_point(data=IRIsub, aes(x=s_mean, y=yjit), shape=108)+
      scale_fill_manual(values=c("orange", "blue"))+
      annotate("text", label=paste("Sample", IRIsub$sample[1]), 
               x=sx, y=79)+
      labs(y="Laboratory ID", x=expression(""^14*"C years BP"))+
      scale_y_continuous(breaks=seq(5, max(IRIdates$lab), 5), 
                         expand=c(0.01,0.01), limits=c(0.5, 82.5))
    
    if(j==1){
      # Theme for first panel of row 1
      IRIplot <- IRIplot + 
        theme(panel.grid=element_blank(),
              legend.position="none", 
              axis.title.x=element_blank(),
              axis.ticks.y=element_line(color="gray88"),
              panel.background=element_blank())
      
    }else if(j==2 | j==3){
      # Theme for right two panels in first row
      IRIplot <- IRIplot + 
        theme(panel.grid=element_blank(),
              legend.position="none",
              axis.ticks.y=element_line(color="gray88"),
              axis.title.y=element_blank(),
              axis.title.x=element_blank(),
              panel.background=element_blank())   
      
    }else if(j==4){
      # Theme for first panel of row 2
      IRIplot <- IRIplot + 
        theme(panel.grid=element_blank(),
              legend.position="none",
              axis.ticks.y=element_line(color="gray88"),
              panel.background=element_blank())
      
    }else if(j==5 | j==6){
      # Theme for right two panels in row 2
      IRIplot <- IRIplot + 
        theme(panel.grid=element_blank(),
              legend.position="none",
              axis.ticks.y=element_line(color="gray88"),
              axis.title.y=element_blank(),
              panel.background=element_blank())   
      
    }
    # Patch together plot panels in loop
    ifelse(j==1, IRIplots <- IRIplot, IRIplots <- IRIplots + IRIplot)
  }
  
  # Remove extraneous variables
  rm(IRIplot, IRIsub, lab_pp, sub_data, c_iter, chain_count, cl,
     IRI_OutlierLabs, j, nsamples, sample_j_interval, sample_j_mean,
     sample_j_samples, sx, sxmax, sxmin)
  
  # Save data file
  save(IRIdates, IRI_diagnostics, IRI_model_samples, IRIplots, 
       IRIsd, IRI_model, file=paste0("IRI",setsn,".RData"))
}else{
  load(paste0("IRI",setsn,".RData"))
}


#######################################################################
#######------ SAMPLE VALUES FOR INTRA-ANNUAL C-14 OFFSETS ------######
#####################################################################

# Reported means and std error of intra-annual date distances, from
# McDonald et al. 2018:7. They report 4 distances, one each for both
# increasing and decreasing curve slopes (ie, different atmospheric
# C-14 production scenarios), and across two different calculation 
# methods. I also add distance for a 0-slope situation, which is 
# simply the midpoint between the increasing and decreasing slope
# distance values
dist_mean <- c(18, 23, 20.5, 26, 22, 24)
dist_se <- c(13, 16, 14.5, 16, 13, 14.5)

# Sample a distance randomly from one of the six distance measurements.
# Use this sampled value to scale a beta distribution, sample from this beta 
# distribution, and recenter the sampled value from 0 - 1 to -0.5 - 0.5.
RCoffsets <- sapply(1:1e4, function(v){
  
  # Sample from measured distances
  p_sample <- sample(1:4, size=1) # random distance position
  s <- abs(rnorm(dist_mean[p_sample], dist_se[p_sample], n=1))
    
  return((rbeta(1, 2, 2)-0.5)*s)
  
})

# Remove extraneous variables from environment
rm(dist_mean, dist_se)


#######################################################################
#####------ PREPARE EVENT SPECIFIC PARAMETERS FOR SIMULATION ----#####
#####################################################################

# Apply log-log to datelist errors to put on scale of IRIdates errors
datelist$logSE <- log(datelist$SEs)

# Create list of events over which to simulate synchroneity
events <- unique(datelist$Event)

# Length of non NA outlier model lambda parameter vector.
oml <- length(which(!is.na(OM_lambda)))

# Create main list to store lists of simulated results for each event.
List_stats_master <- rep(list(), length(events))
# Create main list to store vectors of years for each event.
List_years_master <- list()
# Create main list to store observed std dev and distance statistics
# for each event.
List_observed_master <- list()
# Create main list to store calibrated age plots for each event.
List_plots_ages <- list()
# Create main list to store simulation plots for each event.
List_plots_simulations <- list()
# Create master list of event specific parameters for simulation
ml <- rep(list(), length(events))

# For each event in events vector
for(q in 1:length(events)){
  
  # Subset datelist for event q
  d <- datelist[which(datelist$Event==events[q]), ]
  
  # Find number of unique labs in d, then assign integer designations to
  # each sample marking which lab it was produced by.
  labnames <- unlist(lapply(strsplit(d$lID, "-"), '[[', 1))
  d$lab <- sapply(labnames, function(x) which(unique(labnames)==x))
  
  ########################################################################
  ###-------- CALCULATE STATISTICS FOR OBSERVED DATES --------##########
  ####################################################################
  
  # Combine list of observed calibrated dates into a single data frame
  obs <- calibrate(d$Mus, d$SEs, eps=mpt[1])
  obsDF <- do.call("rbind", obs$grids)
  
  # Create sequence of calendar years based on calendar dates in obsDF
  obs_years <- seq(from=max(obsDF$calBP), to=min(obsDF$calBP), by=-1)
  
  # Calculate unstandardized posterior probabilities for each calendar year
  # from the calibrated radiocarbon samples
  
    # Create matrix to store probabilities for every sample over the entire
    # range of years. 
    obs_probs1m <- sapply(1:length(obs_years), function(x){
      sapply(1:length(obs$grids), function(z){
        q <- obs$grids[[z]]
        ifelse(obs_years[x] %in% q$calBP, 
               return(q$PrDens[which(q$calBP==obs_years[x])]), return(0))
      })
    })
    
    # Calculate posterior from obs_probs1m, use log scale calculation.
    obs_posterior <- exp(apply(log(obs_probs1m), 2, sum) -
      logSumExp(apply(log(obs_probs1m), 2, sum)))
  
  # Standard deviation of raw date means
  obs_sd <- sd(d$Mus)
  
  # Mean distance of calibrated dates
  obs_D <- mean(dist(obs_probs1m, method = "manhattan"))/2
  
  # Store statistics in master list
  List_observed_master[[q]] <- c(obs_sd, obs_D)
  
  # Print observed statistics for user
  cat(paste("\n",  events[q], 
            ": Standard deviation of measured C-14 means = ", 
            round(obs_sd, 3), "\n", sep=""))
  cat(paste(events[q], 
            ": Mean distance between calibrated dates = ", 
            round(obs_D, 3), "\n\n", sep=""))
  
  # Rank each date by the age of its mode
  d$rank <- sapply(1:nrow(d), function(x) {
    obs$grids[[x]]$calBP[which.max(obs$grids[[x]]$PrDens)]})
  d$rank <- rank(d$rank, ties.method="random")
  # Apply rank to data frame of calibrated probabilities
  obsDF$rank <- do.call("c", lapply(1:nrow(d), function(x){
    rep(d$rank[x], nrow(obs$grids[[x]]))
  }))
  # Assign AMS vs non AMS indicator variable to obsDF
  obsDF$AMS <- do.call("c", lapply(1:nrow(d), function(x){
    ifelse(d$Method[x]==0, return(rep("AMS", nrow(obs$grids[[x]]))),
                        return(rep("nonAMS", nrow(obs$grids[[x]]))))
  }))
  # Assign old wood indicator variable to obsDF
  obsDF$OW <- do.call("c", lapply(1:nrow(d), function(x){
    ifelse(d$Outlier[x]==0, return(rep("N", nrow(obs$grids[[x]]))),
           return(rep("Y", nrow(obs$grids[[x]]))))
  }))
  
  # Find calibrated dates spanning over 1000 years and truncate them
  # to avoid plot distortion. For these distributions, this loop
  # returns the 1000 years with the highest probability densities.
  for(j in unique(obsDF$rank)){
    if(length(which(obsDF$rank==j))>1000){
      th <- sort(obsDF$PrDens[which(obsDF$rank==j)], 
                 decreasing=TRUE)[1001]
      obsDF <- obsDF[-which(obsDF$rank==j & obsDF$PrDens < th),]
    }
  }
  
  # Remove small probs from obsDF to prevent plot distortion
  obsDF <- obsDF[which(obsDF$PrDens>1e-5),]
  
  # Scale probs in obsDF and map them to a y-axis
  obsDF$probST <- obsDF$PrDens/max(obsDF$PrDens) + obsDF$rank
  
  # Expression with full date
  d$FullDate <- sapply(1:nrow(d), function(x){
    paste(round(d$Mus[x],0), "\u00B1", d$SEs[x])})
  
  # Create data frame for posterior
  d1 <- data.frame(calBP=obs_years, 
                   probST=obs_posterior/max(obs_posterior), 
                   rank=rep(0, length(obs_years)))
  
  # Retrieve year sequence for event q (generated based on the
  # posterior), if the posterior is valid.
  if(!anyNA(d1$probST) & Plot_posterior=="y"){
    # Retrieve approximate 95% HPDI of posterior.
    post_HPDI <- sapply(c(0.025, 0.975), function(x){
      cumdens <- cumsum(obs_posterior)
      return(obs_years[which.min(abs(cumdens-x))])
    })
  }else{
    # If no posterior, adjust values in obsDF for plotting
    obsDF$rank <- obsDF$rank - 1
    obsDF$probST <- obsDF$probST - 1
  }
  
  # Remove 0 probabilities from d1 to prevent excessively
  # dispersed dates from distorting the plot horizontally.
  d1 <- d1[which(d1$probST>0),]
  
  # x-axis min. If both YDB and LS are simulated, set to same
  # scale for comparability. Otherwise, scale according to the
  # observed data.
  if(("YDB" %in% events) & ("Laacher See" %in% events)){
    if((events[q]=="YDB") | (events[q]=="Laacher See")){
      sitex <- max(obsDF$calBP) - 3450
      idx <- max(obsDF$calBP) - 2800
      datex <- max(obsDF$calBP) - 2200
    }
  }else{
  axm <- 0.2*(max(obsDF$calBP)-min(obsDF$calBP))
  sitex <- rep(min(obsDF$calBP)-(axm*1.3)+5, nrow(d))
  idx <- rep(min(obsDF$calBP), nrow(d))
  datex <- rep(min(obsDF$calBP)+axm, nrow(d))
  }
  
  
  # Plot calibrated ages
  Calplot <- ggplot(NULL, aes(ymin=rank, ymax=probST, x=calBP))+
    annotate("rect", ymin=0, ymax=max(obsDF$rank)+1, xmin=Range_years$start[q],
             xmax=Range_years$end[q], fill="grey", alpha=0.3)
  
  # If there is a valid posterior, add posterior plot elements and sample date
  # text elements and geometry.
  if(nrow(d1)>0 & Plot_posterior=="y"){
    Calplot <- Calplot +
    annotate("rect", ymin=0, ymax=max(obsDF$rank)+1, xmin=post_HPDI[1],
             xmax=post_HPDI[2], fill="grey86")+
    geom_ribbon(data=d1, fill="grey45")+
    geom_ribbon(data=obsDF, aes(group=factor(rank), fill=AMS, color=AMS,
                                alpha=OW), size=0.3)+
    scale_alpha_manual(values=c("Y"=0, "N"=0.6))+
    annotate("text", label="Posterior", x=min(obsDF$calBP)+axm, y=0.4, hjust=0)+
    annotate("text", label=d$Site, x=sitex, y=d$rank+0.4, hjust=0)+
    annotate("text", label=d$lID, x=idx, y=d$rank+0.4, hjust=0)+
    annotate("text", label=d$FullDate, x=datex, y=d$rank+0.4, hjust=0)
    
    # Remove post_HPDI variable after plotting
    rm(post_HPDI)
  # If there is no valid posterior, plot text and geometry for sample dates.
  }else{
    Calplot <- Calplot +
      geom_ribbon(data=obsDF, aes(group=factor(rank), fill=AMS, color=AMS,
                                  alpha=OW), size=0.3)+
      scale_alpha_manual(values=c("Y"=0, "N"=0.6))+
      annotate("text", label=d$Site, x=sitex, y=d$rank-0.6, hjust=0)+
      annotate("text", label=d$lID, x=idx, y=d$rank-0.6, hjust=0)+
      annotate("text", label=d$FullDate, x=datex, y=d$rank-0.6, hjust=0)
  }
  
  # Set manual color and fill scale based on whether the event is YDB
  if(events[q]=="YDB"){
    Calplot <- Calplot +
      scale_fill_manual(values=c("AMS"="red4", "nonAMS"="red"))+
      scale_color_manual(values=c("AMS"="red4", "nonAMS"="red"))
  }else{
    Calplot <- Calplot +
      scale_fill_manual(values=c("AMS"="blue4", "nonAMS"="blue"))+
      scale_color_manual(values=c("AMS"="blue4", "nonAMS"="blue"))
  }
  # Add remaining plot elements that apply to situations where there either
  # is or is not a valid posterior.
  Calplot <- Calplot +
    labs(x="cal BP")+
    scale_y_continuous(expand=c(0.01,0.02), limits=c(0, max(obsDF$rank)+1.2))+
    scale_x_continuous(expand=c(0.01,0.01), breaks=seq(from=0, to=20000, by=500),
                       limits=c(sitex-5, max(obsDF$calBP)))+
    theme(axis.line.x=element_line(color="grey"), legend.position="none",
          axis.ticks.x=element_line(color="grey"), 
          panel.background=element_blank(), axis.ticks.y=element_blank(),
          axis.text.y=element_blank(), axis.title.y=element_blank())
  
  # Assign calplot to list of plots of calibrated dates for each event
  List_plots_ages[[q]] <- Calplot
  
  # Create vector of years for the simulation, based on the hypothesized
  # range and the user specified buffer.
  Years <- seq(Range_years$start[q]+ybuff, Range_years$end[q]-ybuff, by=-1)
  
  # Store in list of year sequences for each event q
  List_years_master[[q]] <- Years
  
  #### Event specific parameters ####
  ###################################
  
  # Initialize matrix to store samples non lab-specific parameter values 
  # that are on the scale of radiocarbon years. Each row corresponds to a 
  # parameter.
  params <- matrix(NA, 4, setsn)
  # params[1,] Baseline within lab std dev
  params[1,] <- IRI_model_samples$a[1:setsn]
  # params[2,] AMS vs non AMS within lab std dev term
  params[2,] <- IRI_model_samples$bM[1:setsn]
  # params[3,] Reported measurement error within lab std dev term
  params[3,] <- IRI_model_samples$bSD[1:setsn]
  # params[4,] Method offset for AMS vs non AMS labs
  params[4,] <- IRI_model_samples$method_offset[1:setsn]
  
  # Simulate setsn samples of varying lab effects for the number
  # of labs that contributed to the dates for this event
  paramsl <- lapply(1:setsn, function(x){
    # Initialize matrix to hold parameters
    mat <- matrix(NA, max(d$lab), 2)
    msample <- sample(size=max(d$lab), 1:setsn)
    # Sample parameter values. Put lab offsets on radiocarbon
    # scale, but leave within-lab sds on log scale.
    mat[,1] <- rnorm(max(d$lab), mean=rep(0, max(d$lab)),
                     sd=IRI_model_samples$lsigma1[msample])*IRIsd
    mat[,2] <- rnorm(max(d$lab), mean=rep(0, max(d$lab)),
                     sd=IRI_model_samples$lsigma2[msample])
    return(mat)
    })
  
  # Store model parameters specific to event q in list
  ml[[q]] <- list(d=d, params=params, paramsl=paramsl)
  
  
  ####################################################################
  ## FOR EACH EVENT, CREATE A LIST OBJECT TO STORE SIMULATED STATS ##
  ##################################################################
  
  # Basic matrix to store stats over outlier models.
  StMa <- rep(list(matrix(NaN, length(Years), setsn)), 1 + oml)
  
  # List of matrices for the various simulation versions.
  List_stats <- list(Stats=list(Uncalibrated_mu_sd=list(No_lab_error=StMa,
                                                        With_lab_error=StMa),
                                Calibrated_mu_dist=list(No_lab_error=StMa,
                                                        With_lab_error=StMa)),
                     P_vals=list(Uncalibrated_mu_sd=list(No_lab_error=matrix(NaN, 
                                                         length(Years), 1 + oml),
                                                         With_lab_error=matrix(NaN, 
                                                         length(Years), 1 + oml)),
                                 Calibrated_mu_dist=list(No_lab_error=matrix(NaN, 
                                                         length(Years), 1 + oml),
                                                         With_lab_error=matrix(NaN, 
                                                         length(Years), 1 + oml))))
  # Store in master list that hold s simulated results
  List_stats_master[[q]] <- List_stats
  
  # Remove extraneous variables from environment
  rm(Calplot, d, d1, obs, obs_probs1m, obsDF, j, q, th, labnames,
     obs_posterior, obs_years, List_stats, StMa, sitex, idx, datex)
  
}



# Concatenate vector of years across events
Years <- do.call("c", List_years_master)

# Create index of events for each year in Years
event_index <- do.call("c", lapply(1:length(events), function(g){
  rep(g, length(List_years_master[[g]]))
}))

# Create list of years and event indices with vectors for each event
Years_list <- list()
event_index_list <- list()
init <- 1
for(f in 1:nNodes){
  ifelse((init + CoreNum) < length(Years),
         last <- init + CoreNum, last <- length(Years)+1)
  Years_list[[f]] <- Years[init:(last-1)]
  event_index_list[[f]] <- event_index[init:(last-1)]
  init <- last
}
rm(init, last) # Remove extraneous variables


#################################################################
################## GENERATE SCRIPTS BASED ON THE ###############
################## NUMBER OF SPECIFIED NODES ##################
##############################################################

# Write data file to be used by every node
save.image("SimDatIntermediate.RData")

# Create r scripts for each node
for(sc in 1:nNodes){
  
  SimScript <- paste("##--- LOAD LIBRARIES ---##
  # In CRAN:
  library(parallel) # For parallel processing
  library(rcarbon) # For radiocarbon calibration
  # Fix randomization for repeatability
  set.seed(100)

  # Load data to be used by every file
  load('SimDatIntermediate.RData')
  
  #######################################################################
  ############## EXPORT VARIABLES AND RCARBON TO CLUSTER ###############
  #####################################################################
  
  # Subset Years and year specific event indices for node
  node <- ", sc,
  "
  # Remove imported existing variables
  rm(Years)  

  Years <- Years_list[[node]]
  event_i <- event_index_list[[node]]

  # Start cluster
  ifelse(CoreNum > length(Years), cl <- makeCluster(length(Years)), 
  cl <- makeCluster(CoreNum))
  
  # Export variables to cluster
  clusterExport(cl, varlist=c('IRIsd', 'RCoffsets', 'setsn', 'oml',
                              'OM_lambda', 'ml', 'mpt', 'Years',
                              'event_i'))
  
  # Export rcarbon package to cluster
  clusterCall(cl, function() library(rcarbon))
  
    
  #######################################################################
  #######------------ SYNCHRONEITY SIMULATION -----------###############
  #####################################################################
  
  # For each year in the concatenated vector of years across events...
  SimList <- parLapply(cl, 1:length(Years), function(y){
  
    # Retrieve event specific parameters
    d <- ml[[event_i[y]]]$d
    params <- ml[[event_i[y]]]$params
    paramsl <- ml[[event_i[y]]]$paramsl
    
    # Use sapply function to simulate statistics for year y
    Stats_y <- sapply(1:setsn, function(j){
      
      # Simulate across outlier models for set j
      j_matrix <- sapply(1:(oml + 1), function(v){
        
        # Create vector of true calendar years
        ndates <- rep(Years[y], nrow(d))
        
        # Use outlier model if v > 1
        if(v>1){
          ndates <- ndates + (d$Outlier*rexp(nrow(d), OM_lambda[v-1]))

          # If any outlier dates exceed 49,000 cal BP, nearing the end
          # of the calibration curve, replace those dates with 49,000.
          ndates[which(ndates > 49000)] <- 49000

        } 
        
        # Uncalibrate each unique date in ndates.
        sampled_dates1 <- sapply(unique(ndates), function(h){
          return(uncalibrate(h)$rCRA)})
        # Assign each unique uncalibrated age to each sample in ndates.
        sampled_dates1 <- sapply(ndates, function(r){
          sampled_dates1[which(unique(ndates)==r)]
        })
        # Add intra-annual C14 offset.
        sampled_dates1 <- sampled_dates1 + sample(RCoffsets, 
                                                  length(ndates))
  
        # Run sampled dates through laboratory bias and repeatability
        # model for an alternative set of observed lab dates.
        sampled_dates2 <- sapply(1:nrow(d), function(x){
          
          # Add lab offset scaled by AMS vs non AMS
          m <- sampled_dates1[x] + 
            paramsl[[j]][d$lab[x],1]*(1+(params[4,j]-1)*d$Method[x])
          
          # Calculate within lab sd based on baseline variability,
          # reported error, and AMS vs non AMS. Transform from z-score
          # to radiocarbon age scale.
          s <- exp(params[1,j] + paramsl[[j]][d$lab[x],2] + 
            params[2,j]*d$logSE + params[3,j]*d$Method[x])*IRIsd
          
          # Sample a date from lab with mean m and within lab
          # std dev s. Return this value.
          return(rnorm(1, m, s))
          
        })
        
        # Calibrate sampled uncalibrated dates
        calibrated_dates1 <- calibrate(sampled_dates1, d$SEs, 
                                      eps=mpt[2])$grids
        calibrated_dates2 <- calibrate(sampled_dates2, d$SEs, 
                                      eps=mpt[2])$grids
        
        # Check for NA values in the unlikely event of a 
        # calibrated value outside the range of the calibration 
        # curve. Count instances of calibrated dates that are 
        # NA using negative integers.
        EmptyCount <- c(0,0)
        for(t in 1:length(calibrated_dates1)){
          if(anyNA(calibrated_dates1[[t]]$calBP)){
            EmptyCount[1] <- EmptyCount[1] - 1
          }
          if(anyNA(calibrated_dates2[[t]]$calBP)){
            EmptyCount[2] <- EmptyCount[2] - 1
          }
        }
        
        # If no calibrated dates are NA...
        if(sum(EmptyCount)==0){
        
          # Bind calibrated dates into a data frame
          cd1 <- do.call('rbind',calibrated_dates1)
          cd2 <- do.call('rbind',calibrated_dates2)
          
          # Create sequence of calendar years based on calendar dates in cd
          year_seq1 <- seq(from=min(cd1$calBP), to=max(cd1$calBP), by=1)
          year_seq2 <- seq(from=min(cd2$calBP), to=max(cd2$calBP), by=1)
          
          # Create matrices to store probabilities for every sample over
          # the entire year ranges
          obs_probs1m1 <- sapply(1:length(year_seq1), function(x){
            sapply(1:length(calibrated_dates1), function(z){
              q <- calibrated_dates1[[z]]
              ifelse(year_seq1[x] %in% q$calBP, 
                return(q$PrDens[which(q$calBP==year_seq1[x])]), 
                return(0))
            })
          })

          obs_probs1m2 <- sapply(1:length(year_seq2), function(x){
            sapply(1:length(calibrated_dates2), function(z){
              q <- calibrated_dates2[[z]]
              ifelse(year_seq2[x] %in% q$calBP, 
                return(q$PrDens[which(q$calBP==year_seq2[x])]), 
                return(0))
            })
          })
          
          # sd for uncalibrated dates in set j and outlier model v
          sample_sd1 <- sd(sampled_dates1)
          sample_sd2 <- sd(sampled_dates2)
          
          # Mean distance of calibrated dates in set j and outlier model v
          sample_D1 <- mean(dist(obs_probs1m1, method='manhattan'))/2
          sample_D2 <- mean(dist(obs_probs1m2, method='manhattan'))/2
          
        }else{
          
          # Else, if one or more calibrated dates is NA, return
          # negative integer in place of summary statistics.
          sample_sd1 <- sample_D1 <- EmptyCount[1]
          sample_sd2 <- sample_D2 <- EmptyCount[2]
        
        }
        
        # Return sample statistics for set j and outlier model v
        return(c(sample_sd1, sample_sd2, sample_D1, sample_D2))
      })
      
      # Convert matrix of sample statistics for all outlier models in 
      # set j into vector and return
      return(c(j_matrix))
      
    })
    
    # Return matrix of statistics for Year y
    return(Stats_y)
  })
  
  
  ####--- End cluster ---###
  ##########################
  stopCluster(cl)
  rm(cl)

  # Write data
  save(SimList, file=paste('NodeDat', node, '.RData', sep=''))",
  sep="")
  
  # Name R file according to node
  filename <- file(paste("NodeSim", sc, ".R", sep=""))
  
  # Write R file
  writeLines(SimScript, filename)
  close(filename)
  
}



###################################################################################
#########------ CREATE BASH FILE FOR EXECUTING R SCRIPTS ON EACH NODE -------#####
#################################################################################

# Create Bash array shell script to run R scripts across nodes
SBatchCmd <- paste('#!/bin/bash
#SBATCH -p standard-mem-m              # Partition (queue)
#SBATCH --exclusive                    # Exclusivity
#SBATCH --mem=250G                     # Total memory required per node
#SBATCH -o nodesim_%A-%a.out     # Job output; %A is job ID and %a is array index
#SBATCH --array=1-', nNodes, '     # Range of indices to be executed

module purge
module load r

R --vanilla < NodeSim${SLURM_ARRAY_TASK_ID}.R', sep="")

# Save shell script
filename <- file("nodesim.sh")
writeLines(SBatchCmd, filename)
close(filename)



#########################################################################################
############--- PAUSE SCRIPT TO WAIT FOR ALL NODES TO FINISH SIMULATIONS ---############
#######################################################################################

# Create vector of data file names for node specific output
DataFiles <- sapply(1:nNodes, function(z){
  paste("NodeDat", z, ".RData", sep="")
})

# Scan file system every 30 seconds to check for simulation output
while(!all(file.exists(DataFiles))){Sys.sleep(30)}


# Read results from each node output and aggregate into one list
SimList <- lapply(DataFiles, function(a){
  
  # Load data from Node a
  load(a)
  # Return simulated results
  return(SimList)
  
})

# Concatenate node specific lists into one list (i.e., remove list level
# that was created by dividing the year vector into node specific chunks)
SimList <- do.call("c", SimList)



#############################################################################
#################-- CONVERT LIST OF SIMULATED MATRICES -----################
#################-- INTO A NESTED LIST WITH NAMED OBJECTS --###############
##########################################################################

# If there were any failed simulation iterations (i.e., calibrated dates
# containing NA values), store these results in a data frame and replace
# the failed iterations with a random successful iteration.
if(any(do.call("rbind", SimList)<0)){
  
  # scan years X for simulations with negative stats...
  FailedIterations <- lapply(1:length(SimList), function(x){
    # If any stats in year x are negative...
    if(any(SimList[[x]]<0)){
      # apply over iterations Y with negative stats for year x...
      Failedx <- lapply(which(apply(SimList[[x]], 2, min)<0), function(y){
        # Apply over old wood models within iteration y in year x
        Failedy <- lapply(seq(4, 4*(oml + 1), 4), function(z){
          # Return empty data frame if old wood model does not contain
          # negative stats.
          if(SimList[[x]][z,y]>0){
            return(data.frame(Event=NA,
                              Year_BP=NA,
                              Failed_Iteration=NA, 
                              OM_lambda=NA,
                              Count_NoLabError=NA,
                              Count_WiLabError=NA,
                              stringsAsFactors=FALSE))
          }else{
            # Return data frame containing failed iteration info 
            # if old wood model contains negative stats.
            return(data.frame(Event=events[event_index[x]],
                              Year_BP=Years[x],
                              Failed_Iteration=y, 
                              OM_lambda=ifelse(z<5, NA, OM_lambda[(z/4)-1]),
                              Count_NoLabError=abs(SimList[[x]][z-1,y]),
                              Count_WiLabError=abs(SimList[[x]][z,y]),
                              stringsAsFactors=FALSE))
          }
          
        })
        return(do.call("rbind", Failedy))
      })
      
      return(do.call("rbind", Failedx))
      
    }else{
      return(data.frame(Event=NA,
                        Year_BP=NA,
                        Failed_Iteration=NA, 
                        OM_lambda=NA,
                        Count_NoLabError=NA,
                        Count_WiLabError=NA,
                        stringsAsFactors=FALSE))
    }

  })
  
  # Bind list into data frame and remove NAs
  FailedIterations <- do.call("rbind", FailedIterations)
  FailedIterations <- FailedIterations[!is.na(FailedIterations$Event),]
  
  # Replace failed iterations with random successful iterations
  for(i in 1:length(SimList)){    # For each year i in the simulation...
    if(any(SimList[[i]]<0)){      # If year i contains an NA stat..
      for(u in seq(4, nrow(SimList[[i]]), 4)){ # Scan 4-row sets for year i...
        if(any(SimList[[i]][(u-3):u,]<0)){ # If 4-row set contains NA stat...
          tm <- SimList[[i]][(u-3):u,] # Isolate 4-row set...
          tm <- tm[,colSums(tm)>0] # ...and remove NA stat occurences.
          for(r in 1:setsn){ # For each column r of the 4-row set with NAs...
            if(any(SimList[[i]][(u-3):u,r]<0)){ # If column r has NAs..
              SimList[[i]][(u-3):u,r] <- # Replace column r with random column
                tm[,sample(1, 1:ncol(tm))] # from tm, which lacks NAs.
            }
          }
        }
      }
    }
  }
  rm(i, u, r, tm) # Remove extraneous variables
}


# Loop through each year of simulated statistics
for(y in 1:length(Years)){
  
  # Simulated results for year y
  Stats_y <- SimList[[y]]
  
  # Observed stats for event associated with year y
  os <- List_observed_master[[event_index[y]]]
  
  # Create specific index for each event
  i <- y - length(event_index[which(event_index < event_index[y])])
  
  # Assign simulated D and sd statistics to appropriate matrices for year y
  for(u in 1:(oml + 1)){
    
    # Store summary statistics for year j
    List_stats_master[[event_index[y]]]$Stats$Uncalibrated_mu_sd$No_lab_error[[u]][i,] <- 
      Stats_y[4*u-3,]
    List_stats_master[[event_index[y]]]$Stats$Uncalibrated_mu_sd$With_lab_error[[u]][i,] <- 
      Stats_y[4*u-2,]
    List_stats_master[[event_index[y]]]$Stats$Calibrated_mu_dist$No_lab_error[[u]][i,] <- 
      Stats_y[4*u-1,]
    List_stats_master[[event_index[y]]]$Stats$Calibrated_mu_dist$With_lab_error[[u]][i,] <- 
      Stats_y[4*u,]
    
    # Store asymptotic p values for year j
    List_stats_master[[event_index[y]]]$P_vals$Uncalibrated_mu_sd$No_lab_error[i,u] <- 
      length(which(Stats_y[4*u-3,] > os[1]))/setsn
    List_stats_master[[event_index[y]]]$P_vals$Uncalibrated_mu_sd$With_lab_error[i,u] <- 
      length(which(Stats_y[4*u-2,] > os[1]))/setsn
    List_stats_master[[event_index[y]]]$P_vals$Calibrated_mu_dist$No_lab_error[i,u] <- 
      length(which(Stats_y[4*u-1,] > os[2]))/setsn
    List_stats_master[[event_index[y]]]$P_vals$Calibrated_mu_dist$With_lab_error[i,u] <- 
      length(which(Stats_y[4*u,] > os[2]))/setsn
  }
  
  # Remove extraneous variables
  rm(Stats_y, os, i)
  
}


#############################################################################
#########################------------------#################################
################------------ PLOTTING ------------#########################
#########################------------------###############################
#########################################################################

####-- Create unique plots for each event --#####
################################################

for(q in 1:length(events)){
  
  # Extract results for event q
  event_data <- List_stats_master[[q]]
  # Extract sequence of years for event q
  Years <- List_years_master[[q]]
  # Extract observed sd and D statistics for event q
  obs_sd <- List_observed_master[[q]][1]
  obs_D <- List_observed_master[[q]][2]
  
  
  # Set uniform y-axis limits for each statistic for event q
  ymaxsd <- max(c(obs_sd, 
                  max(do.call("rbind", event_data$Stats$Uncalibrated_mu_sd$No_lab_error)), 
                  max(do.call("rbind", event_data$Stats$Uncalibrated_mu_sd$With_lab_error))))*1.01
  ymaxD <- max(c(obs_D, 
                 max(do.call("rbind", event_data$Stats$Calibrated_mu_dist$No_lab_error)), 
                 max(do.call("rbind", event_data$Stats$Calibrated_mu_dist$With_lab_error))))*1.01
  
  # Set x-axis breaks based on number of calendar years in Years.
  if(length(Years) < 80){
    xbreaks <- 20} else if(length(Years) < 160){
      xbreaks <- 40} else if(length(Years) < 400){
        xbreaks <- 100} else if(length(Years) < 800){
          xbreaks <- 200} else {xbreaks <- 200}
  xbreaks <- seq(from=0, to=20000, by=xbreaks)
  
  # For each outlier model output, generate four plots for event q.
  for(w in 1:(1 + oml)){
    
    # Extract matrices of results from event_data
    simMsdNLE <- event_data$Stats$Uncalibrated_mu_sd$No_lab_error[[w]]
    simMsdWLE <- event_data$Stats$Uncalibrated_mu_sd$With_lab_error[[w]]
    simMDNLE <- event_data$Stats$Calibrated_mu_dist$No_lab_error[[w]]
    simMDWLE <- event_data$Stats$Calibrated_mu_dist$With_lab_error[[w]]
    
    # Initialize summary data frames
    summaryDFsdNLE <- data.frame(BP=Years)
    summaryDFsdWLE <- data.frame(BP=Years)
    summaryDFDNLE <- data.frame(BP=Years)
    summaryDFDWLE <- data.frame(BP=Years)
    
    # Get 95% and 50% quantiles for each statistic
    summaryDFsdNLE <- cbind(summaryDFsdNLE, 
                            t(sapply(1:nrow(simMsdNLE), function(x){
                              quantile(simMsdNLE[x,], probs=c(0.975, 0.025, 0.75, 0.25, 0.5))})))
    
    summaryDFsdWLE <- cbind(summaryDFsdWLE, 
                            t(sapply(1:nrow(simMsdWLE), function(x){
                              quantile(simMsdWLE[x,], probs=c(0.975, 0.025, 0.75, 0.25, 0.5))})))
    
    summaryDFDNLE <- cbind(summaryDFDNLE, 
                           t(sapply(1:nrow(simMDNLE), function(x){
                             quantile(simMDNLE[x,], probs=c(0.975, 0.025, 0.75, 0.25, 0.5))})))
    
    summaryDFDWLE <- cbind(summaryDFDWLE, 
                           t(sapply(1:nrow(simMDWLE), function(x){
                             quantile(simMDWLE[x,], probs=c(0.975, 0.025, 0.75, 0.25, 0.5))})))
    
    
    # Rename quantile columns
    colnames(summaryDFsdNLE)[2:ncol(summaryDFsdNLE)] <- 
      colnames(summaryDFsdWLE)[2:ncol(summaryDFsdWLE)] <-
      colnames(summaryDFDNLE)[2:ncol(summaryDFDNLE)] <-
      colnames(summaryDFDWLE)[2:ncol(summaryDFDWLE)] <-
      sapply(colnames(summaryDFsdNLE)[2:ncol(summaryDFsdNLE)],
             function(x){paste("b", substr(x, 1, nchar(x)-1), sep="")})
    
    # Create long form quantile lists
    lfDFsdNLE <- lapply(1:2, function(x){
      return(data.frame(BP=Years, upper=summaryDFsdNLE[,x*2], 
                        lower=summaryDFsdNLE[,x*2+1],
                        QuantileID=rep(c("a","b")[x], 
                                       nrow(summaryDFsdNLE))))
    })
    
    lfDFsdWLE <- lapply(1:2, function(x){
      return(data.frame(BP=Years, upper=summaryDFsdWLE[,x*2], 
                        lower=summaryDFsdWLE[,x*2+1],
                        QuantileID=rep(c("a","b")[x], 
                                       nrow(summaryDFsdWLE))))
    })
    
    lfDFDNLE <- lapply(1:2, function(x){
      return(data.frame(BP=Years, upper=summaryDFDNLE[,x*2], 
                        lower=summaryDFDNLE[,x*2+1],
                        QuantileID=rep(c("a","b")[x], 
                                       nrow(summaryDFDNLE))))
    })
    
    lfDFDWLE <- lapply(1:2, function(x){
      return(data.frame(BP=Years, upper=summaryDFDWLE[,x*2], 
                        lower=summaryDFDWLE[,x*2+1],
                        QuantileID=rep(c("a","b")[x], 
                                       nrow(summaryDFDWLE))))
    })
    
    # Collapse long form quantile lists into data frames
    lfDFsdNLE <- do.call("rbind", lfDFsdNLE)
    lfDFsdWLE <- do.call("rbind", lfDFsdWLE)
    lfDFDNLE <- do.call("rbind", lfDFDNLE)
    lfDFDWLE <- do.call("rbind", lfDFDWLE)
    
    # Plot for standard deviation without lab error
    Plot_sdNLE <- ggplot(data=lfDFsdNLE, aes(x=BP, ymin=lower, ymax=upper))+
      annotate("segment", x=Range_years$start[q], xend=Range_years$end[q], y=obs_sd, 
               yend=obs_sd, color="red", size=1)+
      geom_ribbon(aes(group=QuantileID), fill="red", alpha=0.1)+
      annotate("line", x=summaryDFsdNLE$BP, y=summaryDFsdNLE$b50, color="red")+
      scale_x_reverse(expand=c(0.015,0.015))+
      scale_y_continuous(limits=c(0, to=ymaxsd))
    if(w==1){
      Plot_sdNLE <- Plot_sdNLE +
        labs(y=expression(sigma[mu][""^14*"C"]))+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_line(color="grey"), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.line.x=element_blank(),
              axis.ticks.x=element_blank(), axis.title.x=element_blank(),
              axis.text.x=element_blank())
    }else{
      Plot_sdNLE <- Plot_sdNLE +
        labs(title=bquote(lambda==.(OM_lambda[w-1])))+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_blank(), panel.grid=element_blank(),
              axis.ticks.y=element_blank(), axis.title.y=element_blank(),
              axis.text.y=element_blank(), axis.line.x=element_blank(),
              axis.ticks.x=element_blank(), axis.title.x=element_blank(),
              axis.text.x=element_blank(), plot.title = element_text(hjust = 0.5))
    }
    
    # Plot for standard deviation with lab error
    Plot_sdWLE <- ggplot(data=lfDFsdWLE, aes(x=BP, ymin=lower, ymax=upper))+
      annotate("segment", x=Range_years$start[q], xend=Range_years$end[q], y=obs_sd, 
               yend=obs_sd, color="blue", size=1)+
      geom_ribbon(aes(group=QuantileID), fill="blue", alpha=0.1)+
      annotate("line", x=summaryDFsdWLE$BP, y=summaryDFsdWLE$b50, color="blue")+
      scale_x_reverse(expand=c(0.015,0.015))+
      scale_y_continuous(limits=c(0, to=ymaxsd))
    if(w==1){
      Plot_sdWLE <- Plot_sdWLE +
        labs(y=expression(sigma[mu][""^14*"C"]))+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_line(color="grey"), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.line.x=element_blank(),
              axis.ticks.x=element_blank(), axis.title.x=element_blank(),
              axis.text.x=element_blank())
    }else{
      Plot_sdWLE <- Plot_sdWLE +
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_blank(), panel.grid=element_blank(),
              axis.ticks.y=element_blank(), axis.title.y=element_blank(),
              axis.text.y=element_blank(), axis.line.x=element_blank(),
              axis.ticks.x=element_blank(), axis.title.x=element_blank(),
              axis.text.x=element_blank())
    }
    
    # Plot for distance without lab error
    Plot_DNLE <- ggplot(data=lfDFDNLE, aes(x=BP, ymin=lower, ymax=upper))+
      annotate("segment", x=Range_years$start[q], xend=Range_years$end[q], y=obs_D, 
               yend=obs_D, color="red", size=1)+
      geom_ribbon(aes(group=QuantileID), fill="red", alpha=0.1)+
      annotate("line", x=summaryDFDNLE$BP, y=summaryDFDNLE$b50, color="red")+
      scale_x_reverse(expand=c(0.015,0.015))+
      scale_y_continuous(limits=c(0, to=ymaxD), breaks=seq(0, 1, 0.25))
    if(w==1){
      Plot_DNLE <- Plot_DNLE +
        labs(y="MPMD")+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_line(color="grey"), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.line.x=element_blank(),
              axis.ticks.x=element_blank(), axis.title.x=element_blank(),
              axis.text.x=element_blank())
    }else{
      Plot_DNLE <- Plot_DNLE +
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_blank(), panel.grid=element_blank(),
              axis.ticks.y=element_blank(), axis.title.y=element_blank(),
              axis.text.y=element_blank(), axis.line.x=element_blank(),
              axis.ticks.x=element_blank(), axis.title.x=element_blank(),
              axis.text.x=element_blank())
    }
    
    # Plot for distance with lab error
    Plot_DWLE <- ggplot(data=lfDFDWLE, aes(x=BP, ymin=lower, ymax=upper))+
      annotate("segment", x=Range_years$start[q], xend=Range_years$end[q], y=obs_D, 
               yend=obs_D, color="blue", size=1)+
      geom_ribbon(aes(group=QuantileID), fill="blue", alpha=0.1)+
      annotate("line", x=summaryDFDWLE$BP, y=summaryDFDWLE$b50, color="blue")+
      scale_x_reverse(expand=c(0.015,0.015), breaks=xbreaks)+
      scale_y_continuous(limits=c(0, to=ymaxD), breaks=seq(0, 1, 0.25))
    if(w==1){
      Plot_DWLE <- Plot_DWLE +
        labs(x="cal BP", y="MPMD")+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_line(color="grey"), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.line.x=element_line(color="grey"),
              axis.ticks.x=element_line(color="grey"))
    }else{
      Plot_DWLE <- Plot_DWLE +
        labs(x="cal BP")+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_blank(), panel.grid=element_blank(),
              axis.ticks.y=element_blank(), axis.title.y=element_blank(),
              axis.text.y=element_blank(), axis.line.x=element_line(color="grey"),
              axis.ticks.x=element_line(color="grey"))
    }
    
    if(w==1){
      # On first outlier model iteration, initialize plot column
      plotmaster <- Plot_sdNLE + Plot_sdWLE + Plot_DNLE + 
        Plot_DWLE + plot_layout(ncol=1)
      
    }else{
      # On later outlier model iterations, add columns to the plot
      plotcolumn <- Plot_sdNLE + Plot_sdWLE + Plot_DNLE + 
        Plot_DWLE + plot_layout(ncol=1)
      plotmaster <- plotmaster | plotcolumn
    }
    
  }
  
  # Assign master simulation plot to list of events
  List_plots_simulations[[q]] <- plotmaster
  
  # Remove extraneous variables from environment
  rm(plotcolumn, plotmaster, ymaxsd, ymaxD, Years, xbreaks, w, q, 
     obs_sd, obs_D, summaryDFsdWLE, summaryDFsdNLE, summaryDFDWLE,
     summaryDFDNLE, simMDNLE, simMDWLE, simMsdNLE, simMsdWLE, Plot_DNLE,
     Plot_DWLE, Plot_sdNLE, Plot_sdWLE, lfDFDWLE, lfDFDNLE, lfDFsdNLE,
     lfDFsdWLE, event_data)
}

# If the user simulated both YDB and Laacher See, plot both 
# simulations in one figure.
if(("YDB" %in% events) & ("Laacher See" %in% events)){
  
  ###################################################################
  #########----- PLOTS OF PROPORTIONS OF SIMULATED RESULTS ---######
  ##########---- THAT EXCEED OBSERVED VALUES ----##################
  ################################################################
  
  # Create data frame of z scores across both events and all 
  # old wood models for plotting
  zscoreDF <- lapply(1:(oml +1), function(x){
    # Apply to each event
    modelx <- lapply(1:2, function(z){
      # Create empty data frames
      elen <- length(List_years_master[[z]])
      oVals <- List_stats_master[[z]]$Stats
      sdn <- sdw <- dn <- dw <- data.frame(event=rep(events[z], elen),
                                           BP=List_years_master[[z]],
                                           inHypRange=c(rep("no", ybuff),
                                                        rep("yes", elen-2*ybuff), 
                                                        rep("no", ybuff)), 
                                           outlierModel=rep(x, elen))
      sdn$obsStat <- sdw$obsStat<- rep("sd", elen)
      dn$obsStat <- dw$obsStat<- rep("D", elen)
      sdn$labError <- dn$labError <- rep("no", elen)
      sdw$labError <- dw$labError <- rep("yes", elen)
      # Get logit and log z scores
      sdn$values <- sapply(1:elen, function(j){
        rowj <- log(oVals$Uncalibrated_mu_sd$No_lab_error[[x]][j,])
        zscore <- (log(List_observed_master[[z]][1])-mean(rowj))/sd(rowj)
        return(zscore)
      })
      sdw$values <- sapply(1:elen, function(j){
        rowj <- log(oVals$Uncalibrated_mu_sd$With_lab_error[[x]][j,])
        zscore <- (log(List_observed_master[[z]][1])-mean(rowj))/sd(rowj)
        return(zscore)
      })
      dn$values <- sapply(1:elen, function(j){
        rowj <- logit(oVals$Calibrated_mu_dist$No_lab_error[[x]][j,])
        zscore <- (logit(List_observed_master[[z]][2])-mean(rowj))/sd(rowj)
        return(zscore)
      })
      dw$values <- sapply(1:elen, function(j){
        rowj <- logit(oVals$Calibrated_mu_dist$With_lab_error[[x]][j,])
        zscore <- (logit(List_observed_master[[z]][2])-mean(rowj))/sd(rowj)
        return(zscore)
      })
      
      return(rbind(sdn, sdw, dn, dw))
    }) 
    
    return(do.call("rbind", modelx))
  })
  zscoreDF <- do.call("rbind", zscoreDF)
  zscoreDF$grp <- paste(zscoreDF$event, zscoreDF$obsStat, sep="")
  
  ###################################################################
  #########------------ PLOTS OF SIMULATED RESULTS ------------#####
  #################################################################
  
  # Extract results for Laacher See and YDB
  event_dataLS <- List_stats_master[[which(events=="Laacher See")]]
  event_dataYDB <- List_stats_master[[which(events=="YDB")]]
  # Extract sequence of years
  YearsYDB <- List_years_master[[which(events=="YDB")]]
  YearsLS <- List_years_master[[which(events=="Laacher See")]]
  Years <- seq(max(c(YearsYDB, YearsLS)), min(c(YearsYDB, YearsLS)), -1)
  # Extract observed sd and D statistics for Laacher See and YDB
  obs_sdLS <- List_observed_master[[which(events=="Laacher See")]][1]
  obs_DLS <- List_observed_master[[which(events=="Laacher See")]][2]
  obs_sdYDB <- List_observed_master[[which(events=="YDB")]][1]
  obs_DYDB <- List_observed_master[[which(events=="YDB")]][2]
  
  # Sampling distribution bands. Final vector value must be 0.5.
  SimIntervals <- c(0.975, 0.025, 0.75, 0.25, 0.5)
  
  # Set colors for each event
  ls_col <- "blue"
  ydb_col <- "red"
  
  # Set uniform y-axis limits for each statistic between Laacher See and YDB
  y_max <- sapply(1:(oml+1), function(w){
    
    # Extract matrices of results from event_dataLS and event_dataYDB
    simMsdNLEls <- event_dataLS$Stats$Uncalibrated_mu_sd$No_lab_error[[w]]
    simMsdWLEls <- event_dataLS$Stats$Uncalibrated_mu_sd$With_lab_error[[w]]
    simMDNLEls <- event_dataLS$Stats$Calibrated_mu_dist$No_lab_error[[w]]
    simMDWLEls <- event_dataLS$Stats$Calibrated_mu_dist$With_lab_error[[w]]
    simMsdNLEydb <- event_dataYDB$Stats$Uncalibrated_mu_sd$No_lab_error[[w]]
    simMsdWLEydb <- event_dataYDB$Stats$Uncalibrated_mu_sd$With_lab_error[[w]]
    simMDNLEydb <- event_dataYDB$Stats$Calibrated_mu_dist$No_lab_error[[w]]
    simMDWLEydb <- event_dataYDB$Stats$Calibrated_mu_dist$With_lab_error[[w]]
    
    # Get max 97.5% quantile for each simulation
    sdmax <- c(obs_sdLS, obs_sdYDB, rep(NA, 4))
    dmax <- c(obs_DLS, obs_DYDB, rep(NA, 4))
    sdmax[3] <- max(sapply(1:length(YearsLS), function(x){
      quantile(simMsdNLEls[x,], probs=0.975)}))
    sdmax[4] <- max(sapply(1:length(YearsLS), function(x){
      quantile(simMsdWLEls[x,], probs=0.975)}))
    dmax[3] <- max(sapply(1:length(YearsLS), function(x){
      quantile(simMDNLEls[x,], probs=0.975)}))
    dmax[4] <- max(sapply(1:length(YearsLS), function(x){
      quantile(simMDWLEls[x,], probs=0.975)}))
    sdmax[5] <- max(sapply(1:length(YearsYDB), function(x){
      quantile(simMsdNLEydb[x,], probs=0.975)}))
    sdmax[6] <- max(sapply(1:length(YearsYDB), function(x){
      quantile(simMsdWLEydb[x,], probs=0.975)}))
    dmax[5] <- max(sapply(1:length(YearsYDB), function(x){
      quantile(simMDNLEydb[x,], probs=0.975)}))
    dmax[6] <- max(sapply(1:length(YearsYDB), function(x){
      quantile(simMDWLEydb[x,], probs=0.975)}))
    
    return(c(max(sdmax), max(dmax)))
    
  })
  # Assign uniform y-axis limits to each statistic
  ymaxsd <- max(y_max[1,])
  ymaxD <- max(y_max[2,])
  
  # Set x-axis breaks based on number of calendar years in Years.
  if(length(Years) < 80){
    xbreaks <- 20} else if(length(Years) < 160){
      xbreaks <- 40} else if(length(Years) < 400){
        xbreaks <- 100} else if(length(Years) < 800){
          xbreaks <- 200} else {xbreaks <- 200}
  xbreaks <- seq(from=0, to=20000, by=xbreaks)
  
  
  # Set y-axis breaks for z-score plots
  breaksnle <- seq(0, 80, 20)
  breakswle <- seq(0, 12, 3)

  # For each z-score plot, obtain associated y-axis probabilities
  # if they can be calculated, otherwise return asymptotic probabilities.
  mAxnle <- sapply(breaksnle, function(x){
    ifelse(1-pnorm(x) < 1e-15, return("<1e-15"),
           return(as.character(format(1-round(pnorm(x),16), 
                                      scientific=TRUE, digits=2))))
  })
  mAxwle <- sapply(breakswle, function(x){
    ifelse(1-pnorm(x) < 1e-15, return("<1e-15"),
           return(as.character(format(1-round(pnorm(x),15), 
                                      scientific=TRUE, digits=2))))
  })

  
  # For each outlier model output, generate four plots for event q.
  for(w in 1:(1 + oml)){
    
    ###################################################################
    #########----- PLOTS OF PROPORTIONS OF SIMULATED RESULTS ---######
    ##########---- THAT EXCEED OBSERVED VALUES ----##################
    ################################################################
    
    # Subset zscoreDF for outlier model w for data generated from simulation
    # both with and without lhe lab measurement bias and repeatability model.
    zsNE <- zscoreDF[which(zscoreDF$outlierModel==w & zscoreDF$labError=="no"),]
    zsWE <- zscoreDF[which(zscoreDF$outlierModel==w & zscoreDF$labError=="yes"),]
    
    # Plot witn no lab measurement bias and repeatability model.
    Plot_zsNLE <- ggplot(data=zsNE, aes(x=BP, y=values, group=grp))+
      annotate("rect", ymin=rep(0, 2), ymax=rep(Inf, 2), 
               xmin=Range_years$end, xmax=Range_years$start, fill="grey", alpha=0.3)+
      annotate("rect", ymin=0, ymax=qnorm(0.95), xmin=-Inf, xmax=Inf, fill="orange", alpha=0.5)+
      geom_line(aes(color=event, alpha=obsStat), size=0.8)+
      annotate("text", label=c(paste0(LETTERS[w], 1), paste0(LETTERS[w], 1)),
               x=(Range_years$start+Range_years$end)/2, 
               y=rep(max(zscoreDF$values[which(zscoreDF$labError=="no")])*1.1,2),
               color=c("red", "blue"), vjust=1)+
      scale_color_manual(values=c("YDB"="red", "Laacher See"="blue"))+
      scale_alpha_manual(values=c("sd"=0.3, "D"=1))+
      scale_x_reverse(expand=c(0.015,0.015), limits=c(max(zscoreDF$BP), min(zscoreDF$BP)))+
      scale_y_continuous(limits=c(0, 
                        max(zscoreDF$values[which(zscoreDF$labError=="no")])*1.1),
                        breaks=breaksnle)
    if(w==1){
      Plot_zsNLE <- Plot_zsNLE +
        labs(y="z-score", title="No OWM")+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_line(color="grey"), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.line.x=element_blank(),
              axis.ticks.x=element_line(color="grey"), axis.title.x=element_blank(),
              axis.text.x=element_blank(), legend.position="none",
              plot.title=element_text(hjust=0.5, size=11))
    }else{
      Plot_zsNLE <- Plot_zsNLE +
        labs(title=bquote(lambda==.(OM_lambda[w-1])))+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_blank(), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.title.y.left=element_blank(),
              axis.text.y.left=element_blank(), axis.line.x=element_blank(),
              axis.ticks.x=element_line(color="grey"), axis.title.x=element_blank(),
              axis.text.x=element_blank(), legend.position="none",
              plot.title=element_text(hjust=0.5, size=11))
    }
    
    # Plot with lab measurement bias and repeatability model.
    Plot_zsWLE <- ggplot(data=zsWE, aes(x=BP, y=values, group=grp))+
      annotate("rect", ymin=rep(0, 2), ymax=rep(Inf, 2), 
               xmin=Range_years$end, xmax=Range_years$start, fill="grey", alpha=0.3)+
      annotate("rect", ymin=0, ymax=qnorm(0.95), xmin=-Inf, xmax=Inf, fill="orange", alpha=0.5)+
      geom_line(aes(color=event, alpha=obsStat), size=0.8)+
      annotate("text", label=c(paste0(LETTERS[w], 2), paste0(LETTERS[w], 2)),
               x=(Range_years$start+Range_years$end)/2, 
               y=rep(max(zscoreDF$values[which(zscoreDF$labError=="yes")]),2)*1.1,
               color=c("red", "blue"), vjust=1)+
      scale_color_manual(values=c("YDB"="red", "Laacher See"="blue"))+
      scale_alpha_manual(values=c("sd"=0.3, "D"=1))+
      scale_x_reverse(expand=c(0.015,0.015), limits=c(max(zscoreDF$BP), min(zscoreDF$BP)))+
      scale_y_continuous(limits=c(0, 
                        max(zscoreDF$values[which(zscoreDF$labError=="yes")])*1.1),
                        breaks=breakswle)+
      labs(x="cal BP")

    if(w==1){
      Plot_zsWLE <- Plot_zsWLE + labs(y="z-score") +
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_line(color="grey"), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.line.x=element_blank(),
              axis.ticks.x=element_line(color="grey"), legend.position="none")
    }else{
      Plot_zsWLE <- Plot_zsWLE + 
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_line(color="grey"), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.line.x=element_blank(),
              axis.ticks.x=element_line(color="grey"), axis.text.y.left=element_blank(),
              axis.title.y.left=element_blank(), legend.position="none")
    }
    
    ###################################################################
    #########------------ PLOTS OF SIMULATED RESULTS ------------#####
    #################################################################
    
    # Extract matrices of results from event_dataLS and event_dataYDB
    simMsdNLEls <- event_dataLS$Stats$Uncalibrated_mu_sd$No_lab_error[[w]]
    simMsdWLEls <- event_dataLS$Stats$Uncalibrated_mu_sd$With_lab_error[[w]]
    simMDNLEls <- event_dataLS$Stats$Calibrated_mu_dist$No_lab_error[[w]]
    simMDWLEls <- event_dataLS$Stats$Calibrated_mu_dist$With_lab_error[[w]]
    simMsdNLEydb <- event_dataYDB$Stats$Uncalibrated_mu_sd$No_lab_error[[w]]
    simMsdWLEydb <- event_dataYDB$Stats$Uncalibrated_mu_sd$With_lab_error[[w]]
    simMDNLEydb <- event_dataYDB$Stats$Calibrated_mu_dist$No_lab_error[[w]]
    simMDWLEydb <- event_dataYDB$Stats$Calibrated_mu_dist$With_lab_error[[w]]
    
    # Pull summaries of the number of distance and sd values that exceeded
    # the observed for each simulation and store these in a data frame.
    Iter_summary <- data.frame(Event=rep(c(1,2), 4),
       Iter=c(length(which(simMsdNLEls[(ybuff+1):(length(YearsLS)-ybuff),]>=
                    List_observed_master[[which(events=="Laacher See")]][1])),
              length(which(simMsdNLEydb[(ybuff+1):(length(YearsYDB)-ybuff),]>=
                    List_observed_master[[which(events=="YDB")]][1])),
              length(which(simMsdWLEls[(ybuff+1):(length(YearsLS)-ybuff),]>=
                    List_observed_master[[which(events=="Laacher See")]][1])),
              length(which(simMsdWLEydb[(ybuff+1):(length(YearsYDB)-ybuff),]>=
                    List_observed_master[[which(events=="YDB")]][1])),
              length(which(simMDNLEls[(ybuff+1):(length(YearsLS)-ybuff),]>=
                    List_observed_master[[which(events=="Laacher See")]][2])),
              length(which(simMDNLEydb[(ybuff+1):(length(YearsYDB)-ybuff),]>=
                    List_observed_master[[which(events=="YDB")]][2])),
              length(which(simMDWLEls[(ybuff+1):(length(YearsLS)-ybuff),]>=
                    List_observed_master[[which(events=="Laacher See")]][2])),
              length(which(simMDWLEydb[(ybuff+1):(length(YearsYDB)-ybuff),]>=
                    List_observed_master[[which(events=="YDB")]][2]))))
    # Add a variables that expresses the number of iterations as a %
    Iter_summary$Perc <- sapply(Iter_summary$Iter, function(x){
      p <- 100*x/(setsn*(length(YearsYDB)-2*ybuff))
      formattednum <- sprintf(paste0("%4.", nchar(setsn)-1, "f"), p)
      return(paste0(formattednum, "%"))
      })
    
    # Initialize summary data frames
    summaryDFsdNLEls <- data.frame(BP=YearsLS)
    summaryDFsdWLEls <- data.frame(BP=YearsLS)
    summaryDFDNLEls <- data.frame(BP=YearsLS)
    summaryDFDWLEls <- data.frame(BP=YearsLS)
    summaryDFsdNLEydb <- data.frame(BP=YearsYDB)
    summaryDFsdWLEydb <- data.frame(BP=YearsYDB)
    summaryDFDNLEydb <- data.frame(BP=YearsYDB)
    summaryDFDWLEydb <- data.frame(BP=YearsYDB)
    
    # Get 95% and 50% quantiles for each statistic
    summaryDFsdNLEls <- cbind(summaryDFsdNLEls, 
                              t(sapply(1:nrow(simMsdNLEls), function(x){
                                quantile(simMsdNLEls[x,], probs=SimIntervals)})))
    
    summaryDFsdWLEls <- cbind(summaryDFsdWLEls, 
                              t(sapply(1:nrow(simMsdWLEls), function(x){
                                quantile(simMsdWLEls[x,], probs=SimIntervals)})))
    
    summaryDFDNLEls <- cbind(summaryDFDNLEls, 
                             t(sapply(1:nrow(simMDNLEls), function(x){
                               quantile(simMDNLEls[x,], probs=SimIntervals)})))
    
    summaryDFDWLEls <- cbind(summaryDFDWLEls, 
                             t(sapply(1:nrow(simMDWLEls), function(x){
                               quantile(simMDWLEls[x,], probs=SimIntervals)})))
    
    summaryDFsdNLEydb <- cbind(summaryDFsdNLEydb, 
                               t(sapply(1:nrow(simMsdNLEydb), function(x){
                                 quantile(simMsdNLEydb[x,], probs=SimIntervals)})))
    
    summaryDFsdWLEydb <- cbind(summaryDFsdWLEydb, 
                               t(sapply(1:nrow(simMsdWLEydb), function(x){
                                 quantile(simMsdWLEydb[x,], probs=SimIntervals)})))
    
    summaryDFDNLEydb <- cbind(summaryDFDNLEydb, 
                              t(sapply(1:nrow(simMDNLEydb), function(x){
                                quantile(simMDNLEydb[x,], probs=SimIntervals)})))
    
    summaryDFDWLEydb <- cbind(summaryDFDWLEydb, 
                              t(sapply(1:nrow(simMDWLEydb), function(x){
                                quantile(simMDWLEydb[x,], probs=SimIntervals)})))
    
    # Rename quantile columns
    colnames(summaryDFsdNLEls)[2:ncol(summaryDFsdNLEls)] <- 
      colnames(summaryDFsdWLEls)[2:ncol(summaryDFsdWLEls)] <-
      colnames(summaryDFDNLEls)[2:ncol(summaryDFDNLEls)] <-
      colnames(summaryDFDWLEls)[2:ncol(summaryDFDWLEls)] <-
      colnames(summaryDFsdNLEydb)[2:ncol(summaryDFsdNLEydb)] <- 
      colnames(summaryDFsdWLEydb)[2:ncol(summaryDFsdWLEydb)] <-
      colnames(summaryDFDNLEydb)[2:ncol(summaryDFDNLEydb)] <-
      colnames(summaryDFDWLEydb)[2:ncol(summaryDFDWLEydb)] <-
      sapply(colnames(summaryDFsdNLEls)[2:ncol(summaryDFsdNLEls)],
             function(x){paste("b", substr(x, 1, nchar(x)-1), sep="")})
    
    # Create long form quantile lists
    lfDFsdNLEls <- lapply(1:2, function(x){
      return(data.frame(BP=YearsLS, upper=summaryDFsdNLEls[,x*2], 
                        lower=summaryDFsdNLEls[,x*2+1],
                        QuantileID=rep(c("a","b")[x], 
                                       nrow(summaryDFsdNLEls))))
    })
    
    lfDFsdWLEls <- lapply(1:2, function(x){
      return(data.frame(BP=YearsLS, upper=summaryDFsdWLEls[,x*2], 
                        lower=summaryDFsdWLEls[,x*2+1],
                        QuantileID=rep(c("a","b")[x], 
                                       nrow(summaryDFsdWLEls))))
    })
    
    lfDFDNLEls <- lapply(1:2, function(x){
      return(data.frame(BP=YearsLS, upper=summaryDFDNLEls[,x*2], 
                        lower=summaryDFDNLEls[,x*2+1],
                        QuantileID=rep(c("a","b")[x], 
                                       nrow(summaryDFDNLEls))))
    })
    
    lfDFDWLEls <- lapply(1:2, function(x){
      return(data.frame(BP=YearsLS, upper=summaryDFDWLEls[,x*2], 
                        lower=summaryDFDWLEls[,x*2+1],
                        QuantileID=rep(c("a","b")[x], 
                                       nrow(summaryDFDWLEls))))
    })
    
    lfDFsdNLEydb <- lapply(1:2, function(x){
      return(data.frame(BP=YearsYDB, upper=summaryDFsdNLEydb[,x*2], 
                        lower=summaryDFsdNLEydb[,x*2+1],
                        QuantileID=rep(c("a","b")[x], 
                                       nrow(summaryDFsdNLEydb))))
    })
    
    lfDFsdWLEydb <- lapply(1:2, function(x){
      return(data.frame(BP=YearsYDB, upper=summaryDFsdWLEydb[,x*2], 
                        lower=summaryDFsdWLEydb[,x*2+1],
                        QuantileID=rep(c("a","b")[x], 
                                       nrow(summaryDFsdWLEydb))))
    })
    
    lfDFDNLEydb <- lapply(1:2, function(x){
      return(data.frame(BP=YearsYDB, upper=summaryDFDNLEydb[,x*2], 
                        lower=summaryDFDNLEydb[,x*2+1],
                        QuantileID=rep(c("a","b")[x], 
                                       nrow(summaryDFDNLEydb))))
    })
    
    lfDFDWLEydb <- lapply(1:2, function(x){
      return(data.frame(BP=YearsYDB, upper=summaryDFDWLEydb[,x*2], 
                        lower=summaryDFDWLEydb[,x*2+1],
                        QuantileID=rep(c("a","b")[x], 
                                       nrow(summaryDFDWLEydb))))
    })
    
    # Collapse long form quantile lists into data frames
    lfDFsdNLEls <- do.call("rbind", lfDFsdNLEls)
    lfDFsdWLEls <- do.call("rbind", lfDFsdWLEls)
    lfDFDNLEls <- do.call("rbind", lfDFDNLEls)
    lfDFDWLEls <- do.call("rbind", lfDFDWLEls)
    lfDFsdNLEydb <- do.call("rbind", lfDFsdNLEydb)
    lfDFsdWLEydb <- do.call("rbind", lfDFsdWLEydb)
    lfDFDNLEydb <- do.call("rbind", lfDFDNLEydb)
    lfDFDWLEydb <- do.call("rbind", lfDFDWLEydb)
    
    # Plot for standard deviation without lab error
    Plot_sdNLE <- ggplot(data=NULL, aes(x=BP, ymin=lower, ymax=upper))+
      annotate("segment", x=Range_years$start[which(events=="Laacher See")], 
               xend=Range_years$end[which(events=="Laacher See")], y=obs_sdLS, 
               yend=obs_sdLS, color=ls_col, size=0.7)+
      geom_ribbon(data=lfDFsdNLEls, aes(group=QuantileID), fill=ls_col, alpha=0.2)+
      annotate("line", x=summaryDFsdNLEls$BP, y=summaryDFsdNLEls$b50, color=ls_col, alpha=0.4)+
      annotate("segment", x=Range_years$start[which(events=="YDB")], 
               xend=Range_years$end[which(events=="YDB")], y=obs_sdYDB, 
               yend=obs_sdYDB, color=ydb_col, size=0.7)+
      geom_ribbon(data=lfDFsdNLEydb, aes(group=QuantileID), fill=ydb_col, alpha=0.2)+
      annotate("line", x=summaryDFsdNLEydb$BP, y=summaryDFsdNLEydb$b50, color=ydb_col, alpha=0.4)+
      annotate("text", label=paste(LETTERS[w],1,sep=""), size=5, color="grey",
               y=0-ymaxsd*.28, x=max(Years), vjust=0, hjust=0)+
      annotate("text", label=("'Expected'>=~'Observed values'"), x=mean(Years),
               y=0, vjust=1, parse=TRUE)+
      annotate("text", label=Iter_summary$Iter[1:2], x=c(mean(YearsLS),
               mean(YearsYDB)), color=c("blue", "red"), y=rep(0-ymaxsd*0.14, 2))+
      annotate("text", label=Iter_summary$Perc[1:2], x=c(mean(YearsLS),
               mean(YearsYDB)), color=c("blue", "red"), vjust=0,
               y=rep(0-ymaxsd*0.28, 2))+
      scale_x_reverse(expand=c(0.015,0.015), limits=c(max(Years), min(Years)))+
      scale_y_continuous(limits=c(0-ymaxsd*0.28, ymaxsd), breaks=seq(0,1e3,100))
    if(w==1){
      Plot_sdNLE <- Plot_sdNLE +
        labs(y=expression(sigma[mu][""^14*"C"]), title="No OWM")+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_line(color="grey"), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.line.x=element_blank(),
              axis.ticks.x=element_line(color="grey"), axis.title.x=element_blank(),
              axis.text.x=element_blank(), plot.title=element_text(hjust=0.5, size=11))
    }else{
      Plot_sdNLE <- Plot_sdNLE +
        labs(title=bquote(lambda==.(OM_lambda[w-1])))+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_blank(), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.title.y=element_blank(),
              axis.text.y=element_blank(), axis.line.x=element_blank(),
              axis.ticks.x=element_line(color="grey"), axis.title.x=element_blank(),
              axis.text.x=element_blank(), plot.title=element_text(hjust=0.5, size=11))
    }
    
    # Plot for standard deviation with lab error
    Plot_sdWLE <- ggplot(data=NULL, aes(x=BP, ymin=lower, ymax=upper))+
      annotate("segment", x=Range_years$start[which(events=="Laacher See")], 
               xend=Range_years$end[which(events=="Laacher See")], y=obs_sdLS, 
               yend=obs_sdLS, color=ls_col, size=0.7)+
      geom_ribbon(data=lfDFsdWLEls, aes(group=QuantileID), fill=ls_col, alpha=0.2)+
      annotate("line", x=summaryDFsdWLEls$BP, y=summaryDFsdWLEls$b50, color=ls_col, alpha=0.4)+
      annotate("segment", x=Range_years$start[which(events=="YDB")], 
               xend=Range_years$end[which(events=="YDB")], y=obs_sdYDB, 
               yend=obs_sdYDB, color=ydb_col, size=0.7)+
      geom_ribbon(data=lfDFsdWLEydb, aes(group=QuantileID), fill=ydb_col, alpha=0.2)+
      annotate("line", x=summaryDFsdWLEydb$BP, y=summaryDFsdWLEydb$b50, color=ydb_col, alpha=0.4)+
      annotate("text", label=paste(LETTERS[w],2,sep=""), size=5, color="grey",
               y=0-ymaxsd*.28, x=max(Years), vjust=0, hjust=0)+
      annotate("text", label=("'Expected'>=~'Observed values'"), x=mean(Years),
               y=0, vjust=1, parse=TRUE)+
      annotate("text", label=Iter_summary$Iter[3:4], x=c(mean(YearsLS),
               mean(YearsYDB)), color=c("blue", "red"), y=rep(0-ymaxsd*0.14, 2))+
      annotate("text", label=Iter_summary$Perc[3:4], x=c(mean(YearsLS),
               mean(YearsYDB)), color=c("blue", "red"), vjust=0,
               y=rep(0-ymaxsd*0.28, 2))+
      scale_x_reverse(expand=c(0.015,0.015), breaks=xbreaks, limits=c(max(Years), min(Years)))+
      scale_y_continuous(limits=c(0-ymaxsd*0.28, ymaxsd), breaks=seq(0,1e3,100))
    if(w==1){
      Plot_sdWLE <- Plot_sdWLE +
        labs(y=expression(sigma[mu][""^14*"C"]), x="cal BP")+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_line(color="grey"), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.line.x=element_blank(),
              axis.ticks.x=element_line(color="grey"))
    }else{
      Plot_sdWLE <- Plot_sdWLE +
        labs(x="cal BP")+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_blank(), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.title.y=element_blank(),
              axis.text.y=element_blank(), axis.line.x=element_blank(),
              axis.ticks.x=element_line(color="grey"))
    }
    
    # Plot for distance without lab error
    Plot_DNLE <- ggplot(data=NULL, aes(x=BP, ymin=lower, ymax=upper))+
      annotate("segment", x=Range_years$start[which(events=="Laacher See")], 
               xend=Range_years$end[which(events=="Laacher See")], y=obs_DLS, 
               yend=obs_DLS, color=ls_col, size=0.7)+
      geom_ribbon(data=lfDFDNLEls, aes(group=QuantileID), fill=ls_col, alpha=0.2)+
      annotate("line", x=summaryDFDNLEls$BP, y=summaryDFDNLEls$b50, color=ls_col, alpha=0.4)+
      annotate("segment", x=Range_years$start[which(events=="YDB")], 
               xend=Range_years$end[which(events=="YDB")], y=obs_DYDB, 
               yend=obs_DYDB, color=ydb_col, size=0.7)+
      geom_ribbon(data=lfDFDNLEydb, aes(group=QuantileID), fill=ydb_col, alpha=0.2)+
      annotate("line", x=summaryDFDNLEydb$BP, y=summaryDFDNLEydb$b50, color=ydb_col, alpha=0.4)+
      annotate("text", label=paste(LETTERS[w],1,sep=""), size=5, color="grey",
               y=0-ymaxD*.28, x=max(Years), vjust=0, hjust=0)+
      annotate("text", label=("'Expected'>=~'Observed values'"), x=mean(Years),
               y=0, vjust=1, parse=TRUE)+
      annotate("text", label=Iter_summary$Iter[5:6], x=c(mean(YearsLS),
              mean(YearsYDB)), color=c("blue", "red"), y=rep(0-ymaxD*0.14, 2))+
      annotate("text", label=Iter_summary$Perc[5:6], x=c(mean(YearsLS),
               mean(YearsYDB)), color=c("blue", "red"), vjust=0,
               y=rep(0-ymaxD*0.28, 2))+
      scale_x_reverse(expand=c(0.015,0.015), limits=c(max(Years), min(Years)))+
      scale_y_continuous(limits=c(0-ymaxD*0.28, ymaxD), breaks=seq(0, 1, 0.25))
    if(w==1){
      Plot_DNLE <- Plot_DNLE +
        labs(y="MPMD", title="No OWM")+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_line(color="grey"), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.line.x=element_blank(),
              axis.ticks.x=element_line(color="grey"), axis.title.x=element_blank(),
              axis.text.x=element_blank(), plot.title=element_text(hjust=0.5, size=11))
    }else{
      Plot_DNLE <- Plot_DNLE +
        labs(title=bquote(lambda==.(OM_lambda[w-1])))+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_blank(), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.title.y=element_blank(),
              axis.text.y=element_blank(), axis.line.x=element_blank(),
              axis.ticks.x=element_line(color="grey"), axis.title.x=element_blank(),
              axis.text.x=element_blank(), plot.title=element_text(hjust=0.5, size=11))
    }
    
    # Plot for distance with lab error
    Plot_DWLE <- ggplot(data=NULL, aes(x=BP, ymin=lower, ymax=upper))+
      annotate("segment", x=Range_years$start[which(events=="Laacher See")], 
               xend=Range_years$end[which(events=="Laacher See")], y=obs_DLS, 
               yend=obs_DLS, color=ls_col, size=0.7)+
      geom_ribbon(data=lfDFDWLEls, aes(group=QuantileID), fill=ls_col, alpha=0.2)+
      annotate("line", x=summaryDFDWLEls$BP, y=summaryDFDWLEls$b50, color=ls_col, alpha=0.4)+
      annotate("segment", x=Range_years$start[which(events=="YDB")], 
               xend=Range_years$end[which(events=="YDB")], y=obs_DYDB, 
               yend=obs_DYDB, color=ydb_col, size=0.7)+
      geom_ribbon(data=lfDFDWLEydb, aes(group=QuantileID), fill=ydb_col, alpha=0.2)+
      annotate("line", x=summaryDFDWLEydb$BP, y=summaryDFDWLEydb$b50, color=ydb_col, alpha=0.4)+
      annotate("text", label=paste(LETTERS[w],2,sep=""), size=5, color="grey",
               y=0-ymaxD*.28, x=max(Years), vjust=0, hjust=0)+
      annotate("text", label=("'Expected'>=~'Observed values'"), x=mean(Years),
               y=0, vjust=1, parse=TRUE)+
      annotate("text", label=Iter_summary$Iter[7:8], x=c(mean(YearsLS),
               mean(YearsYDB)), color=c("blue", "red"), y=rep(0-ymaxD*0.14, 2))+
      annotate("text", label=Iter_summary$Perc[7:8], x=c(mean(YearsLS),
               mean(YearsYDB)), color=c("blue", "red"), vjust=0,
               y=rep(0-ymaxD*0.28, 2))+
      scale_x_reverse(expand=c(0.015,0.015), breaks=xbreaks, limits=c(max(Years), min(Years)))+
      scale_y_continuous(limits=c(0-ymaxD*0.28, ymaxD), breaks=seq(0, 1, 0.25))

    if(w==1){
      Plot_DWLE <- Plot_DWLE +
        labs(x="cal BP", y="MPMD")+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_line(color="grey"), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.line.x=element_line(color="grey"),
              axis.ticks.x=element_line(color="grey"))
      
    }else{
      Plot_DWLE <- Plot_DWLE +
        labs(x="cal BP")+
        theme(panel.background=element_rect(fill=NA, color="grey"), 
              axis.line.y=element_blank(), panel.grid=element_blank(),
              axis.ticks.y=element_line(color="grey"), axis.title.y=element_blank(),
              axis.text.y=element_blank(), axis.line.x=element_line(color="grey"),
              axis.ticks.x=element_line(color="grey"))
    }
    
    if(w==1){
      # On first outlier model iteration, initialize plot columns
      Plot_YDB_LS_sim_sd <- Plot_sdNLE + Plot_sdWLE + plot_layout(ncol=1)
      Plot_YDB_LS_sim_D <- Plot_DNLE + Plot_DWLE + plot_layout(ncol=1)
      Plot_YDB_LS_z <- Plot_zsNLE + Plot_zsWLE + plot_layout(ncol=1)
      
    }else{
      # On later outlier model iterations, add columns to the plot
      plotcolumn <- Plot_sdNLE + Plot_sdWLE + plot_layout(ncol=1)
      Plot_YDB_LS_sim_sd <- Plot_YDB_LS_sim_sd | plotcolumn
      plotcolumn <- Plot_DNLE + Plot_DWLE + plot_layout(ncol=1)
      Plot_YDB_LS_sim_D <- Plot_YDB_LS_sim_D | plotcolumn
      plotcolumn <- Plot_zsNLE + Plot_zsWLE + plot_layout(ncol=1)
      Plot_YDB_LS_z <- Plot_YDB_LS_z | plotcolumn
    }
    
  }
  
  # Create plot for lab measurement bias and repeatability model titles
  modeltitleplot <- ggplot()+
    annotate("text", label="With LBM", y=10, x=1, angle=270)+
    scale_x_continuous(limits=c(0,2))+
    scale_y_continuous(limits=c(0,50))+
    theme_void()
  
  Plot_YDB_LS_sim_sd <- Plot_YDB_LS_sim_sd + modeltitleplot +
    plot_layout(widths=c(rep(8, oml+1),1))
  
  Plot_YDB_LS_sim_D <- Plot_YDB_LS_sim_D + modeltitleplot +
    plot_layout(widths=c(rep(8, oml+1),1))
  
  Plot_YDB_LS_z <- Plot_YDB_LS_z + modeltitleplot +
    plot_layout(widths=c(rep(8, oml+1),1))
  
  # Remove extraneous variables from environment
  rm(plotcolumn, ymaxsd, ymaxD, xbreaks, w, obs_sdLS, obs_sdYDB, obs_DLS,
     obs_DYDB, summaryDFsdWLEls, summaryDFsdNLEls, summaryDFDWLEls, 
     summaryDFDNLEls, simMDNLEls, simMDWLEls, simMsdNLEls, simMsdWLEls, 
     Plot_DNLE, Plot_DWLE, Plot_sdNLE, Plot_sdWLE, lfDFDWLEls, lfDFDNLEls, 
     lfDFsdNLEls, lfDFsdWLEls, event_dataLS, event_dataYDB, summaryDFsdWLEydb, 
     summaryDFsdNLEydb, summaryDFDWLEydb, summaryDFDNLEydb, simMDNLEydb, 
     simMDWLEydb, simMsdNLEydb, simMsdWLEydb, lfDFDWLEydb, lfDFDNLEydb, 
     lfDFsdNLEydb,lfDFsdWLEydb, SimIntervals, ls_col, ydb_col, YearsLS, 
     YearsYDB, modeltitleplot, Plot_zsNLE, Plot_zsWLE, zsNE, zsWE, zscoreDF,
     breaksnle, breakswle, mAxnle, mAxwle, y_max, Years, Iter_summary)
  
}

# Remove extraneous variables from environment
rm(u, y, SimScript, sc, SBatchCmd, filename, f, params, paramsl, RCoffsets)

# Save variables from environment
save.image(paste("SimDat", setsn, ".RData", sep=""))
