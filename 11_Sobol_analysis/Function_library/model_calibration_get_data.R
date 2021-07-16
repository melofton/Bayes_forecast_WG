# Title: Pull appropriate data files for model calibration runs
# History:
# created MEL 27MAR20

#install to load and install other packages as needed
#install.packages('pacman')

#load packages
pacman::p_load(tidyverse)

get_calibration_data <- function(model_name){

#set calibration years and weeks of season - do not edit
years <- c(2009:2012)
year_no = as.numeric(as.factor(years))
season_weeks = c(1:20)

#read in Gloeo data
y <- as.matrix(read_csv("./11_Sobol_analysis/Data/Gechinulata_Site1.csv"))
#remove 2015-2016 data
y <- y[-c(5:8),]

###############################TWO COVARIATE MODELS#####################################

#for temp_and_wind
if(model_name == "temp_and_wind"){

  #read in covariate 1 data
  Temp <- as.matrix(read_csv("./11_Sobol_analysis/Data/wtrtemp_min_Site1.csv"))
  #remove 2015-2016 data
  Temp <- Temp[-c(5:8),]
  #center covariate data
  Temp <- (Temp - mean(Temp, na.rm = TRUE))/sd(Temp, na.rm = TRUE)

  #read in data from Site 2 for data gap-filling
  Temp_prior <- as.matrix(read_csv("./11_Sobol_analysis/Data/wtrtemp_min_Site2.csv"))
  #remove 2015-2016 data
  Temp_prior <- Temp_prior[-c(5:8),]
  #center water temp data
  Temp_prior <- (Temp_prior - mean(Temp_prior, na.rm = TRUE))/sd(Temp_prior, na.rm = TRUE)

  #calculate weekly average of covariate from past years for gap filling
  week_avg1 = colMeans(Temp_prior, na.rm = TRUE)
  #use weekly average from last sampled week (18) to serve as prior for weeks 19 & 20
  week_avg1[is.na(week_avg1)] <- week_avg1[19]

  #read in covariate 2 data
  Wind <- as.matrix(read_csv("./11_Sobol_analysis/Data/wnd_dir_mean_2daylag.csv"))
  #remove 2015-2016 data
  Wind <- Wind[-c(5:8),]
  #center covariate data
  Wind <- (Wind - mean(Wind, na.rm = TRUE))/sd(Wind, na.rm = TRUE)

  #calculate weekly average of covariate from past years for gap filling
  week_avg2 = colMeans(Wind, na.rm = TRUE)

  return(list(year_no = year_no, season_weeks = season_weeks, y = y, covar1 = Temp, covar2 = Wind, week_avg1 = week_avg1, week_avg2 = week_avg2))

}

}

