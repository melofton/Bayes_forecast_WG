##First cut at Sobol analysis of hindcast variance
##Author: Mary Lofton
##Date: 16JUL21

#clear environment
rm(list = ls())

#load packages
pacman::p_load(tidyverse, readxl, rjags, runjags, moments, coda)

#set a directory to use as a local file repository for plots if desire to write to file
my_directory <- "C:/Users/Mary Lofton/Documents/RProjects/Bayes_forecast_WG/11_Sobol_analysis/Calibration_output"
write_plots <- FALSE

##List of steps

####1. Calibrate model####

#a. Source helper functions ---------------------------------------------------------
source('11_Sobol_analysis/Function_library/model_calibration_plug_n_play.R')
source('11_Sobol_analysis/Function_library/model_calibration_get_data.R')
source('11_Sobol_analysis/Function_library/model_calibration_plots.R')

#b. Model options => pick model -----------------------------------------------------

model_name = "temp_and_wind"
model=paste0("11_Sobol_analysis/Models/",model_name, '.R') #Do not edit


#c. Read in data for model ------------------------------------------------------------------------------------------------------------

#see 0_Function_library/model_calibration_get_data.R for this function
cal_data <- get_calibration_data(model_name)


#d. JAGS Plug-Ins => initial conditions, priors, data, etc. --------------------------------------------------------------------------------------

#see 0_Function_library/model_calibration_plug_n_play.R for this function
jags_plug_ins <- jags_plug_ins(model_name = model_name)


#e. Run model  -------------------------------------------------------------
j.model   <- jags.model (file = model,
                         data = jags_plug_ins$data.model,
                         inits = jags_plug_ins$init.model,
                         n.chains = 3)

jags.out <- run.jags(model = model,
                     data = jags_plug_ins$data.model,
                     adapt =  5000,
                     burnin =  10000,
                     sample = 50000,
                     n.chains = 3,
                     inits=jags_plug_ins$init.model,
                     monitor = jags_plug_ins$variable.namesout.model)

#convert to an MCMC list to calculate cross-correlation later
jags.out.mcmc <- as.mcmc.list(jags.out)


#f. Save output for calibration assessment -------------------------------

#plot parameters
plot_parameters(params = jags_plug_ins$params.model,
                write_plots = write_plots,
                my_directory = my_directory)

#calculate parameter summaries, effective sample size, and cross-correlations
sum <- summary(jags.out, vars = jags_plug_ins$variable.names.model)
crosscorr <- crosscorr(jags.out.mcmc[,c(jags_plug_ins$params.model)])

#save results
sink(file = file.path("./11_Sobol_analysis/Calibration_output/",paste0(model_name,'_param_summary.txt')))
print("Parameter summary")
print(sum)
print("Parameter cross-correlations")
print(crosscorr)
sink()

#save parameter posteriors
Nmc = 7500
out <- as.matrix(jags.out.mcmc)
rows_A <- sample.int(nrow(out),Nmc,replace=TRUE)
rows_B <- sample.int(nrow(out),Nmc,replace=TRUE)

####2. Grab what you need for start-of-season hindcast####

#a. Parameter estimates (including process error) -----------

params_A <- out[rows_A,jags_plug_ins$params.model]
params_B <- out[rows_B,jags_plug_ins$params.model]


params.A <- list(sd_obs = 0, sd_proc = 1/sqrt(params_A[,"tau_proc"]), beta1 = params_A[,"beta1"],
                 beta2 = params_A[,"beta2"], beta3 = params_A[,"beta3"], beta4 = params_A[,"beta4"],
                 sd_C1 = 1/sqrt(params_A[,"tau_C1_proc"]), sd_C1 = 1/sqrt(params_A[,"tau_C2_proc"]))

params.B <- list(sd_obs = 0, sd_proc = 1/sqrt(params_B[,"tau_proc"]), beta1 = params_B[,"beta1"],
                 beta2 = params_B[,"beta2"], beta3 = params_B[,"beta3"], beta4 = params_B[,"beta4"],
                 sd_C1 = 1/sqrt(params_B[,"tau_C1_proc"]), sd_C1 = 1/sqrt(params_B[,"tau_C2_proc"]))


#b. Initial conditions from informed prior -------------

x_ic=-5
tau_ic = 100

ic = rnorm(150000,x_ic,tau_ic)

ic_A = ic[rows_A]
ic_B = ic[rows_B]

#c. Z-scores ---------------

z = rnorm(Nmc,0,1)

#d. Driver hindcasts ---------------

#NEED TWO HINDCASTS FOR A AND B ENSEMBLES

#set temporary counters for yr and wks till get for-loop up and running
yrs <- c(2013:2016)
year_no = as.numeric(as.factor(yrs))
wks <- c(1:20)
j = 1
k = 1

#source helper functions
source('11_Sobol_analysis/Function_library/hindcasting_get_data.R') #check
source('11_Sobol_analysis/Function_library/hindcasting_get_covar_hindcasts.R') #check

###grab covar data
  #see 0_Function_library/hindcasting_get_data.R for this function
  hindcast_data <- get_hindcast_data(model_name = model_name,
                                     year = yrs[j],
                                     season_week = wks[k])

###gap-fill missing covariate values using latent states from calibrated model

  #NOTES: need to change calibration to also monitor latent states of drivers
  #then need to edit covar hindcasts to only draw from 2009-2012
  covar_ls1 <- out[,grep("covar1", colnames(out))]
  missing1 <- which(is.na(hindcast_data$covar1_hindcast))

  covar_ls2 <- out[,grep("covar2", colnames(out))]
  missing2 <- which(is.na(hindcast_data$covar2_hindcast))

  #this will bonk as soon as start hindcasting for 2014 b/c won't have imputed
  #driver values for 2013 since calibration does not currently include 2013
  for (m in 1:length(missing1)){
    hindcast_data$covar1_hindcast[missing1[m]] <- mean(covar_ls1[,missing1[m]],na.rm = TRUE)
  }
  for (m in 1:length(missing2)){
    hindcast_data$covar2_hindcast[missing2[m]] <- mean(covar_ls2[,missing2[m]],na.rm = TRUE)
  }


#set up sampling of covariate hindcasting ensemble for models with covariates

  #set up sampling for hindcasted covariates
  if(yrs[j] == 2013){year_no <- c(1:4)}
  if(yrs[j] == 2014){year_no <- c(1:5)}
  if(yrs[j] == 2015){year_no <- c(1:6)}
  if(yrs[j] == 2016){year_no <- c(1:7)}
  #set up draws from covariate matrix to account for intraannual correlation
  yrsamp.A <- sample(year_no, Nmc, replace = TRUE)
  yrsamp.B <- sample(year_no, Nmc, replace = TRUE)


#get hindcasted covariates
covar.hindcast.A <- get_covar_hindcasts(model_name = model_name,
                                             wk = wks[k],
                                             yrsamp = yrsamp.A,
                                             Nmc = Nmc,
                                             covar_ensemble = list(covar1 = hindcast_data$covar1_hindcast, covar2 = hindcast_data$covar2_hindcast))

covar.hindcast.B <- get_covar_hindcasts(model_name = model_name,
                                        wk = wks[k],
                                        yrsamp = yrsamp.B,
                                        Nmc = Nmc,
                                        covar_ensemble = list(covar1 = hindcast_data$covar1_hindcast, covar2 = hindcast_data$covar2_hindcast))



####3. Run start-of-season hindcast (4 weeks)####

#use existing hindcast function code as starting point for this; will need to
#modify process error so that function is deterministic
#REMEMBER TO RUN TWICE! (ENSEMBLE A AND B)
source('11_Sobol_analysis/Function_library/hindcasting_run_hindcast.R') #check


# Each team has a forecast function in R: forecast(p,i,d,z)
# Run the following 12 ensemble forecast experiments [exp]

#a. Ensemble A -----
hindcast.A <- run_hindcast(model_name = model_name,
                           params = params.A,
                           proc = params.A$sd_proc,
                           Nmc = Nmc,
                           IC = ic_A,
                           wk = wks[k],
                           covar_hindcast = covar.hindcast.A,
                           z = z)
#b. Ensemble B ----
hindcast.B <- run_hindcast(model_name = model_name,
                           params = params.B,
                           proc = params.B$sd_proc,
                           Nmc = Nmc,
                           IC = ic_B,
                           wk = wks[k],
                           covar_hindcast = covar.hindcast.B,
                           z = z)

#c. ABp (A with B’s parameter samples) ----
hindcast.ABp <- run_hindcast(model_name = model_name,
                           params = params.B,
                           proc = params.A$sd_proc,
                           Nmc = Nmc,
                           IC = ic_A,
                           wk = wks[k],
                           covar_hindcast = covar.hindcast.A,
                           z = z)


#d. ABi (A with B’s initial condition samples) ----
hindcast.ABi <- run_hindcast(model_name = model_name,
                           params = params.A,
                           proc = params.A$sd_proc,
                           Nmc = Nmc,
                           IC = ic_B,
                           wk = wks[k],
                           covar_hindcast = covar.hindcast.A,
                           z = z)

#e. ABd (A with B’s driver samples) ----
hindcast.ABd <- run_hindcast(model_name = model_name,
                           params = params.A,
                           proc = params.A$sd_proc,
                           Nmc = Nmc,
                           IC = ic_A,
                           wk = wks[k],
                           covar_hindcast = covar.hindcast.B,
                           z = z)

#f. ABz (A with B’s z-score samples) ----
hindcast.ABz <- run_hindcast(model_name = model_name,
                           params = params.A,
                           proc = params.B$sd_proc,
                           Nmc = Nmc,
                           IC = ic_A,
                           wk = wks[k],
                           covar_hindcast = covar.hindcast.A,
                           z = z)

#g. ABpi (A with B’s parameters & initial conditions) ----
hindcast.ABpi <- run_hindcast(model_name = model_name,
                           params = params.B,
                           proc = params.A$sd_proc,
                           Nmc = Nmc,
                           IC = ic_B,
                           wk = wks[k],
                           covar_hindcast = covar.hindcast.A,
                           z = z)

#h. ABpd ----
hindcast.ABpd <- run_hindcast(model_name = model_name,
                           params = params.B,
                           proc = params.A$sd_proc,
                           Nmc = Nmc,
                           IC = ic_A,
                           wk = wks[k],
                           covar_hindcast = covar.hindcast.B,
                           z = z)

#i. ABpz ----
hindcast.ABpz <- run_hindcast(model_name = model_name,
                           params = params.B,
                           proc = params.B$sd_proc,
                           Nmc = Nmc,
                           IC = ic_A,
                           wk = wks[k],
                           covar_hindcast = covar.hindcast.A,
                           z = z)

#j. ABid ----
hindcast.ABid <- run_hindcast(model_name = model_name,
                           params = params.A,
                           proc = params.A$sd_proc,
                           Nmc = Nmc,
                           IC = ic_B,
                           wk = wks[k],
                           covar_hindcast = covar.hindcast.B,
                           z = z)

#k. ABiz ----
hindcast.ABiz <- run_hindcast(model_name = model_name,
                           params = params.A,
                           proc = params.B$sd_proc,
                           Nmc = Nmc,
                           IC = ic_B,
                           wk = wks[k],
                           covar_hindcast = covar.hindcast.A,
                           z = z)

#l. ABdz ----
hindcast.ABdz <- run_hindcast(model_name = model_name,
                           params = params.A,
                           proc = params.B$sd_proc,
                           Nmc = Nmc,
                           IC = ic_A,
                           wk = wks[k],
                           covar_hindcast = covar.hindcast.B,
                           z = z)


####stopped developing script here - analysis code below is non-functional#####
##4. Analysis

# When you encounter data, run the Analysis on BOTH forecasts A and forecast B,
# but NOT any of the other experiments

#set temporary counter for wks till get for-loop up and running
k = 2

#a. Retrieve next week of data for both Gloeo and drivers
source('11_Sobol_analysis/Function_library/model_analysis_get_data.R') #check

analysis.data <- get_analysis_data(model_name = model_name,
                                   year_no = year_no[j],
                                   wk = wks[k])

#b. Name priors and data

#see 0_Function_library/model_calibration_plug_n_play.R for this function
#Linear_2var
data.Linear_2var <- list(y=analysis.data$y, covar1=analysis.data$covar1, covar2=analysis.data$covar2, week_avg1=analysis.data$week_avg1,week_avg2=analysis.data$week_avg2, beta.m1=0,  beta.m2=0,beta.m3=0,beta.m4=0, beta.v1=0.001, beta.v2=0.001,beta.v3=0.001,beta.v4=0.001, x_ic=-5,tau_ic = 100,a_proc = 0.001,r_proc = 0.001, a_obs = 15.37, r_obs = 7.84)
variable.names.Linear_2var <- c("tau_proc", "beta1","beta2", "beta3","beta4", "tau_obs","tau_C1_proc", "tau_C2_proc")
variable.namesout.Linear_2var <- c("tau_proc", "beta1", "beta2","beta3","beta4",  "mu", "tau_obs", "tau_C1_proc", "tau_C2_proc", "covar1","covar2")
init.Linear_2var <- list(list(tau_proc=0.001, tau_obs = 0.1,  tau_C1_proc = 0.01,tau_C2_proc = 0.01, beta1=-0.5, beta2=-0.5, beta3=-0.5, beta4=-0.5), list(tau_proc=0.1,  tau_obs = 1,tau_C1_proc = 0.1,tau_C2_proc = 0.1, beta1=0, beta2=0, beta3=0, beta4=0), list(tau_proc=1, tau_obs = 5,tau_C1_proc = 1,tau_C2_proc = 1, beta1=0.5,beta2=0.5, beta3=0.5, beta4=0.5))
params.Linear_2var <- c("tau_proc","beta1", "beta2", "beta3","beta4","tau_obs","tau_C1_proc", "tau_C2_proc")

j.model   <- jags.model (file = model,
                         data = jags_plug_ins$data.model,
                         inits = jags_plug_ins$init.model,
                         n.chains = 3)

#c. Call to JAGS



jags.out <- run.jags(model = model,
                     data = jags_plug_ins$data.model,
                     adapt =  5000,
                     burnin =  10000,
                     sample = 50000,
                     n.chains = 3,
                     inits=jags_plug_ins$init.model,
                     monitor = jags_plug_ins$variable.namesout.model)

  #i. parameter priors are posteriors of calibrated parameters
  #ii. initial conditions of Gloeo are hindcasted Gloeo for that week
  #iii. run with single week of Gloeo/driver data and save estimated latent state

#AGAIN WILL NEED TO DO TWICE WITH A AND B HINDCASTS

##5. Grab what you need for week 2 hindcast

# Same as above, but the initial conditions come from the previous A and B
# forecasts (or Analysis posterior where appropriate).
# You do NOT grab IC from any of the other 10 experiments

#a. updated latent state for A and B from JAGS or previous forecast
#b. everything else stays the same

