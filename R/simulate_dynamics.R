#'@title Simulate metacommunity dynamics

#' @param initial_df dataframe specifying the initial conditions with ID, xy coordinates, speices denisities
#' @param r growth rate of species
#' @param alpha length(r)^2 competition matrix
#' @param delta aggregation parameter where 0 is uniform weighting, >0 high densities are weighted more, <0 less
#' @param timesteps length of time to simulate dynamics
#' @param disp_rate dispersal rate which can range from 0 (no dispersal) to 1 (maximum dispersal)
#' @param nh_size Manhattan distance specifying size of neighborhood
#' @param nh Type of neighborhood for dispersal to occur. Can be either 'Von Neummann' or 'Moore'
#' @param dd_emi If TRUE, emigration is density-dependent and follows type-II functional response
#' @param torus If TRUE, the lattice is wrapped around to remove edge effects
#' @param stochastic If TRUE, denisities are discretized and the process is stochastic


#' @import dplyr
#' @import tidyr
#' @export

simulate_dynamics <- function(initial_df, r, alpha, delta, timesteps, disp_rate = 0.1, nh_size = 1,
                              nh = "vonNeumann", stochastic = FALSE, dd_emi = FALSE,  torus = TRUE){

  # require(dplyr)
  # options(dplyr.summarise.inform = FALSE)
  # require(tidyr)

  #Species names
  sp_names <- paste0('N', 1:length(r))

  #Empty list to store iterations -----
  out <- vector(mode = 'list', length = timesteps+1) #spp densities


  #Set initial data as first element of list
  out[[1]] <- initial_df

  #Iterate dynamics ----------
  for(i in 2:(timesteps+1)){

    #=================================
    #One iteration of the simulation -------
    #=================================

    # 1. Growth --------------------

    #Species densities at time t (coordinates removed)
    # sp_densities <- out[[i-1]]%>%
    #   select(sp_names)
    #
    # #Spatial ID and coordinates
    # coords <- out[[i-1]]%>%
    #   select(-sp_names)
    #

    #apply bh_competition to each cell (single iteration)
    comp <- out[[i-1]]%>%
      select(sp_names)%>%
      apply(., 1, function(x)bh_competition(r = r,  #applies function to each row (i.e. cell)
                                        alpha = alpha,
                                        timesteps = 1,
                                        discrete = stochastic, #if stochastic is TRUE, then discrete individuals are used (w/ Poisson)
                                        initial_N = x, final_N = TRUE))%>%
      t()%>%
      as.data.frame()%>%
      setNames(., sp_names)%>% #put species names back in dataframe
      cbind(out[[i-1]]%>%select(-sp_names), .) #put back ID and coordinates


    # 2. Emigrate --------------------

    #note: spp have same emigration rates
    emi <- out[[i-1]]%>% #this DOES NOT implement order into simulation: is emigration(n(t)) and NOT emigration(n'(t))
      rowwise()%>%
      mutate_at(vars(-ID, -x, -y), .funs = function(x)dd_emigration(N = x, dd = dd_emi, discrete = stochastic, max = disp_rate))%>%
      ungroup()



    #3. Immigrate --------------------
    #Calculate neighborhood just once (doing it every time slows the process)
    immi <- purrr::map2(emi$x, emi$y, .f = function(x,y)neighborhood(x = x, y = y, df = out[[i-1]],
                                                                     nh_size = nh_size,
                                                                     neighborhood = nh,
                                                                     torus = torus))%>% #Get neighborhoods of each cell
      purrr::map2(.,
                  emi%>%select(sp_names)%>%as.matrix()%>%split(1:nrow(.)),
                  .f = function(x,y)fd_immigration(neighborhood_df = x,
                                                   disp_inds = y,
                                                   delta = delta,
                                                   discrete = stochastic))%>% #Distribute emigrants across neighborhood patches
      bind_rows()%>%
      group_by(ID)%>% #take the sum of all immigrants coming in to the patch (ID)
      dplyr::summarise(N1 = sum(N1), N2 = sum(N2)) #explicit from dplyr so that there is no conflict with plyr





    # 4. Combine steps: growth - emigration + immigration --------------------
    out[[i]] <- out[[i-1]]%>%
      select(-sp_names)%>%
      cbind(comp%>%select(sp_names) - emi%>%select(sp_names) + immi%>%select(sp_names)) #Add components together (C - E + I)


  }



  #Export data
  out <- out%>%
    bind_rows(.id = 'time')%>%
    mutate(time = as.integer(time))
  #mutate(delta = delta)



  return(out)

}
