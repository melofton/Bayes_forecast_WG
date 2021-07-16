# Title: Create hindcasted covariates for various forms of hindcast
# History:
# created MEL 17APR20

#install to load and install other packages as needed
#install.packages('pacman')

#load packages
pacman::p_load(tidyverse)

#function
get_covar_hindcasts <- function(model_name, wk, yrsamp, Nmc, covar_ensemble){

  #assign same model name for models with the same structure
  if(model_name %in% c("wtrtemp_min","wtrtemp_min_lag","wtrtemp_MA7","wnd_dir_2day_lag","GDD","schmidt_max_lag")){
    model_type <- "1var"
  }
  if(model_name %in% c("schmidt_and_wind","temp_and_wind","wind_and_GDD")){
    model_type <- "2var"
  }

  #make vector of columns corresponding to week of hindcast
  if(wk %in% c(1:17)){
    hindcast_wks <- c(wk,wk+1,wk+2,wk+3)
    } else if (wk == 18){
    hindcast_wks <- c(wk,wk+1,wk+2)
    } else if (wk == 19){
    hindcast_wks <- c(wk,wk+1)
    } else {hindcast_wks <- c(wk)}

  if(model_type == "2var"){
    #set up output matrices
    covar1 <- matrix(NA, Nmc, length(hindcast_wks))
    covar2 <- matrix(NA, Nmc, length(hindcast_wks))


      for(i in 1:length(hindcast_wks)){
        covar1[,i] <- covar_ensemble$covar1[yrsamp,hindcast_wks[i]]
        covar2[,i] <- covar_ensemble$covar2[yrsamp,hindcast_wks[i]]
      }


    return(list(covar1 = covar1, covar2 = covar2))
  }

}


