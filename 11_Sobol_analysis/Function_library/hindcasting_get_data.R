# Title: Pull appropriate data files for model calibration runs
# History:
# created MEL 27MAR20

#install to load and install other packages as needed
#install.packages('pacman')

#load packages
pacman::p_load(tidyverse)

get_hindcast_data <- function(model_name, year, season_week){

###############################FOR ALL MODELS#####################################

#set calibration years and weeks of season - do not edit
if(year == 2015){years <- c(2009:2015)} else {years <- c(2009:2016)}
year_no = as.numeric(as.factor(years))
season_weeks = c(1:20)

#read in Gloeo data
y0 <- as.matrix(read_csv("./11_Sobol_analysis/Data/Gechinulata_Site1.csv"))

#subset data depending on year and season_week
if(year == 2013){
  y <- y0[c(1:4),]

}
if(year == 2014){
  y <- y0[c(1:5),]

}
if(year == 2015){
  y <- y0[c(1:6),]

}

if(year == 2016){
  y <- y0[c(1:7),]

}

if (model_name == "temp_and_wind"){
  #read in covar1 data
  covar01 <- as.matrix(read_csv("./11_Sobol_analysis/Data/wtrtemp_min_Site1.csv"))
  #read in prior data
  prior1 <- as.matrix(read_csv("./11_Sobol_analysis/Data/wtrtemp_min_Site2.csv"))

  #read in covar2 data
  covar02 <- as.matrix(read_csv("./11_Sobol_analysis/Data/wnd_dir_mean_2daylag.csv"))
  #read in prior2 data
  prior2 <- as.matrix(read_csv("./11_Sobol_analysis/Data/wnd_dir_mean_2daylag.csv"))
}


#subset covar data depending on year and season_week
if(model_name == "temp_and_wind"){
  if(year == 2013){

    #create covar1 timeseries
    covar1 <- covar01[c(1:4),]

    #create covar1 hindcast
    covar1_hindcast <- covar01[c(1:4),]
    #standardize covar hindcast
    covar1_hindcast <- (covar1_hindcast - mean(covar1_hindcast, na.rm = TRUE))/sd(covar1_hindcast, na.rm = TRUE)

    #standardize covar1 gap-filling dataset
    prior1 <- (prior1 - mean(prior1, na.rm = TRUE))/sd(prior1, na.rm = TRUE)
    #create gap-filling weekly avg
    week_avg1 = colMeans(prior1[c(1:4),], na.rm = TRUE)
    week_avg1[is.na(week_avg1)] <- week_avg1[19]

    #create covar2 timeseries
    covar2 <- covar02[c(1:4),]

    #create covar1 hindcast
    covar2_hindcast <- covar02[c(1:4),]
    #standardize covar hindcast
    covar2_hindcast <- (covar2_hindcast - mean(covar2_hindcast, na.rm = TRUE))/sd(covar2_hindcast, na.rm = TRUE)

    #standardize covar1 gap-filling dataset
    prior2 <- (prior2 - mean(prior2, na.rm = TRUE))/sd(prior2, na.rm = TRUE)
    #create gap-filling weekly avg
    week_avg2 = colMeans(prior2[c(1:4),], na.rm = TRUE)

  }

  if(year == 2014){

    #create covar1 timeseries
    covar1 <- covar01[c(1:5),]

    #create covar1 hindcast
    covar1_hindcast <- covar01[c(1:5),]
    #standardize covar hindcast
    covar1_hindcast <- (covar1_hindcast - mean(covar1_hindcast, na.rm = TRUE))/sd(covar1_hindcast, na.rm = TRUE)

    #standardize covar1 gap-filling dataset
    prior1 <- (prior1 - mean(prior1, na.rm = TRUE))/sd(prior1, na.rm = TRUE)
    #create gap-filling weekly avg
    week_avg1 = colMeans(prior1[c(1:5),], na.rm = TRUE)
    week_avg1[is.na(week_avg1)] <- week_avg1[19]

    #create covar2 timeseries
    covar2 <- covar02[c(1:5),]

    #create covar1 hindcast
    covar2_hindcast <- covar02[c(1:5),]
    #standardize covar hindcast
    covar2_hindcast <- (covar2_hindcast - mean(covar2_hindcast, na.rm = TRUE))/sd(covar2_hindcast, na.rm = TRUE)

    #standardize covar1 gap-filling dataset
    prior2 <- (prior2 - mean(prior2, na.rm = TRUE))/sd(prior2, na.rm = TRUE)
    #create gap-filling weekly avg
    week_avg2 = colMeans(prior2[c(1:5),], na.rm = TRUE)

  }

  if(year == 2015){

    #create covar1 timeseries
    covar1 <- covar01[c(1:6),]

    #create covar1 hindcast
    covar1_hindcast <- covar01[c(1:6),]
    #standardize covar hindcast
    covar1_hindcast <- (covar1_hindcast - mean(covar1_hindcast, na.rm = TRUE))/sd(covar1_hindcast, na.rm = TRUE)

    #standardize covar1 gap-filling dataset
    prior1 <- (prior1 - mean(prior1, na.rm = TRUE))/sd(prior1, na.rm = TRUE)
    #create gap-filling weekly avg
    week_avg1 = colMeans(prior1[c(1:6),], na.rm = TRUE)
    week_avg1[is.na(week_avg1)] <- week_avg1[19]

    #create covar2 timeseries
    covar2 <- covar02[c(1:6),]

    #create covar1 hindcast
    covar2_hindcast <- covar02[c(1:6),]
    #standardize covar hindcast
    covar2_hindcast <- (covar2_hindcast - mean(covar2_hindcast, na.rm = TRUE))/sd(covar2_hindcast, na.rm = TRUE)

    #standardize covar1 gap-filling dataset
    prior2 <- (prior2 - mean(prior2, na.rm = TRUE))/sd(prior2, na.rm = TRUE)
    #create gap-filling weekly avg
    week_avg2 = colMeans(prior2[c(1:6),], na.rm = TRUE)

  }

  if(year == 2016){

    #create covar1 timeseries
    covar1 <- covar01[c(1:7),]

    #create covar1 hindcast
    covar1_hindcast <- covar01[c(1:7),]
    #standardize covar hindcast
    covar1_hindcast <- (covar1_hindcast - mean(covar1_hindcast, na.rm = TRUE))/sd(covar1_hindcast, na.rm = TRUE)

    #standardize covar1 gap-filling dataset
    prior1 <- (prior1 - mean(prior1, na.rm = TRUE))/sd(prior1, na.rm = TRUE)
    #create gap-filling weekly avg
    week_avg1 = colMeans(prior1[c(1:7),], na.rm = TRUE)
    week_avg1[is.na(week_avg1)] <- week_avg1[19]

    #create covar2 timeseries
    covar2 <- covar02[c(1:7),]

    #create covar1 hindcast
    covar2_hindcast <- covar02[c(1:7),]
    #standardize covar hindcast
    covar2_hindcast <- (covar2_hindcast - mean(covar2_hindcast, na.rm = TRUE))/sd(covar2_hindcast, na.rm = TRUE)

    #standardize covar1 gap-filling dataset
    prior2 <- (prior2 - mean(prior2, na.rm = TRUE))/sd(prior2, na.rm = TRUE)
    #create gap-filling weekly avg
    week_avg2 = colMeans(prior2[c(1:7),], na.rm = TRUE)

  }

  #standardize covar timeseries
  covar1 <- (covar1 - mean(covar1, na.rm = TRUE))/sd(covar1, na.rm = TRUE)
  covar2 <- (covar2 - mean(covar2, na.rm = TRUE))/sd(covar2, na.rm = TRUE)

  return(list(year_no = year_no, season_weeks = season_weeks, y = y, covar1 = covar1, covar2 = covar2, covar1_hindcast = covar1_hindcast, covar2_hindcast = covar2_hindcast, week_avg1 = week_avg1, week_avg2 = week_avg2))

}


}
