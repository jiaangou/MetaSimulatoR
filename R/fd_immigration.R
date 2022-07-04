#' @title Frequency-dependent immigration

#' @description Applies weights to each neighborhood cell and distributes emigrants to each accordingly

#' @param neighborhood_df neighborhood data with, patch ID, coordinates (x,y) and species densities
#' @param disp_inds number immigrants to distribute (for each species)
#' @param delta coefficients determining strength of frequency-dependence
#' @param discrete whether or not immigrants should be distributed as discrete individuals

#' @return a dataframe with neighborhood coordinates and immigrants densities (NOT updated densities!)

#' @import dplyr
#' @import plyr
#' @import tidyr
#' @export



#FUNCTION: Returns neighborhood with new immigrants ----
fd_immigration <- function(neighborhood_df, disp_inds, delta, discrete = FALSE, no_spp = 2){

  # require(dplyr)
  # require(plyr)
  # require(tidyr)

  #Species names for subsetting and generalization to >2spp
  sp_names <- paste0("N", 1:no_spp)


  #If dispersing individuals are 0, then return neighborhood_df with 0s
  if(sum(disp_inds) == 0){

    neighborhood_df[,sp_names] <- 0
    immi <- neighborhood_df


  }else{



    ######### Calculate weights #######################
    #Function: weighted frequencies (given delta)-------
    weighted_freq <- function(N, delta){

      w <- exp(delta*(N+1))

      w_freq <- w/sum(w)

      return(w_freq)
    }

    #Calculate weights for each species in the neighborhood
    weights <- neighborhood_df%>%
      select(sp_names)%>%
      apply(., 2, function(x)weighted_freq(x, delta = delta))


    ##############################################
    #Sample with replacement if discrete = TRUE
    ##############################################
    if(discrete == TRUE){

      ########### BUG #########################
      ########### BUG: number of inds  #########################
      #Iterates disp_inds and weights simultaneously

      immi <- purrr::map2(.x = disp_inds,
                          .y = weights%>%split(col(.)),
                          .f = function(x,y)sample(neighborhood_df$ID, size = x, prob = y, replace = TRUE))%>%
        lapply(function(x)count(x))%>%
        bind_rows(.id = 'Species')%>%
        dplyr::rename(`ID` = `x`)%>%
        mutate(ID = factor(ID, levels = neighborhood_df$ID))%>%
        mutate(Species = paste0('N', Species))%>%
        mutate(Species = factor(Species, levels = sp_names))%>%
        complete(Species, ID, fill = list(freq = 0))%>%
        pivot_wider(names_from = 'Species', values_from = freq, values_fill = 0)%>%
        mutate(ID = ID%>%
                 as.character()%>% #so that values do not change
                 as.integer())%>% #to match original class for combining
        left_join(neighborhood_df%>%select(-sp_names), ., by = 'ID')%>%
        select(colnames(neighborhood_df)) #rearranges order to be same as input




    }

    ##############################################
    #If detereministic and continous if discrete = FALSE
    ##############################################
    else{
      #Multiple frequencies with disp inds
      immi <- weights%*%diag(disp_inds)%>%
        as.data.frame()%>%
        setNames(sp_names)%>%
        cbind(neighborhood_df%>%select(-sp_names),.)%>%
        select(colnames(neighborhood_df)) #rearranges order to be same as input


    }



  }




  return(immi)



}

