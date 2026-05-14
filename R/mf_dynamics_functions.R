#' @title
#' rotate matrix
#' @description
#' function to rotate matrix, used in mf function.
#'
#' @param x matrix to rotate.
#'
#' @returns rotated matrix (?)
rotate <- function(x) {
  x_rotated <- t(apply(x, 2, rev))

  return(x_rotated)
}


#HERE IS WHERE MY CHANGES ARE
#' @title
#' change.micro()
#' @description
#' To evaluate/change? the number of MF in one MF age class for all people. Replaces change.micro and derivmf.one/rest
#' 
#' @param give.treat  value (1/0) indicating whether treatment occurs 
#' @param mu.rates.mf mf mortality rate (modelled as a weibull distribution, see other param defn)
#' @param mf.cpt age class under consideration.
#' @param num.comps number of age classes in adult worms
#' @param ws  column in main matrix where worm age compartments start for each individual 
#' @param ep.in fecundity of female worms (based on fec.rates in adult female worms (m(a) in Hamley et al. 2019))
#' @param mf.st (initialised at 7)  column (first age compartment) in matrix where mf begin, in order to calculate which column to use depending on compartment
#' @param dat main matrix containing different worm compartments for age-classes; can i just use dat? or does this skip that param defn?
#' @param treat.vec vector of treatment values for each individuals, to specify human population length with/without treatment 
# mf.mort mf.mu rate specified for first age class in each human
# mf.move determines mf aging (moving rate to the next age class)
#' @param DT timestep
#' @param time.each.comp Duration of each age class for MF (q_M in Hamley et al. 2019 supp)
#'
#' @returns vector/matrix? of values for MF count per individual

change.micro <- function(mf.st = 7, give.treat= give.treat, mu.rates.mf= mort.rates.mf,
		 iteration = i,  mf.cpt=mf.c, num.comps=num.comps.worm, ws=worms.start,
                 dat= all.mats.cur,
		 treat.vec=treat.vec.in,
                 DT=DT, time.each.comp=time.each.comp.worms,
		 aging_in = rep(0,N),
		 num.mf.comps = num.mf.comps, fec.rates = fec.rates.worms, mf.move.rate = mf.move.rate,
		 up = up, kap = kap, treat.start = treat.start){
  
  #for both scenarios (treatment & no treatment)
   ep.in <- fec.rates
   mf.move <- mf.move.rate
   N <- length(treat.vec)
   mf.mu <- rep(mu.rates.mf[mf.cpt], N)

#there's a number of repetitions, could make this more concise   

  #######FOR TREATMENT SCENARIO#####

   if(give.treat == 1 & iteration >= treat.start){
    #taken from mf_dynamics_function directly
    tao <- ((iteration - 1) * DT) - treat.vec # tao is zero if treatment has been given at this timestep

    mu.mf.prime <- ((tao + up) ^ (- kap)) # additional mortality due to ivermectin treatment

    mu.mf.prime[which(is.na(mu.mf.prime) == TRUE)] <- 0

    mf.mu <- mf.mu + mu.mf.prime #this is the big change

  # indexes for fertile worms (to use in production of mf)
    fert.worms.start <-  ws + num.comps * 2
    fert.worms.end <-  (ws - 1) + num.comps * 3

  # indexes to check if there are males (males start is just 'ws')
  # there must be >= 1 male worm for females to produce microfilariae
    mal.worms.end <- (ws - 1) + num.comps

    fert.worms <- dat[, fert.worms.start:fert.worms.end] #number of fertile females worms
    #fertile worms + male worm -> new birthed mf
    new.in <- (rotate(fert.worms) * ep.in) #need to rotate matrix to each column is multiplied by respective fecundity rate, not each row
    new.in <- rotate(rotate(rotate(new.in)))
    new.in <- rowSums(new.in)
#indexing  compartment location    
    compartments_ind <- (mf.st -1) + mf.cpt
    mf.birthed <-rpois(N, new.in*DT ) # check, might be wrong scale, also check if Poisson?


    mf.cur <- dat[,compartments_ind] #current MF load in compartment
    #mortality -- check if correct
    mf.mu <- rep(1 - exp(-mu.rates.mf[mf.cpt] * DT), N)
    mf.mort <- mf.mu #unnecessary, and sloppy notation
    #making death stochastic using rbinom()
    mf.die <- rbinom(N, pmax(0L, mf.cur), mf.mort) 
    #aging out using rbinom()
    mf.loss.aged <- rbinom(N, pmax(0L, mf.cur - mf.die), (DT/time.each.comp))

   #influx/movement param
    in.param <- mf.move
    in.param <- rep(1 - exp(-mf.move * DT), N) #hopefully will fix NAs

#change by compartment
    if(mf.cpt == 1){
      mf.out <- mf.cur + mf.birthed - mf.die - mf.loss.aged
     }

    if(mf.cpt > 1){
      mf.out <- mf.cur + aging_in - mf.die - mf.loss.aged  # aging_in comes from previous compartment
     }

   mf.out <- pmax(0L, mf.out)  # safety floor
  
  
}  
  
  ###########NON-TREATMENT SCENARIO##############
   
   
   #######ATM, also for the non-treatment period of treatment scenario
   #give.treat == 1 & iteration >= treat.start
 if(give.treat == 0 | give.treat == 1 & iteration < treat.start | give.treat == 1 & iteration > treat.end){
      
  #mortality      
  mf.mort <- 1 - exp(-mu.rates.mf[mf.cpt] * DT)
  mf.mu <- rep(mf.mort, N)
  mf.mort <- mf.mu #unnecessary, and sloppy notation, but i'll fix later

  #checking for fertile female worms to birth new MF    
  # indexes for fertile worms (to use in production of mf)
  fert.worms.start <-  ws + num.comps * 2
  fert.worms.end <-  (ws - 1) + num.comps * 3

  # indexes to check if there are males (males start is just 'ws')
  # there must be >= 1 male worm for females to produce microfilariae
  mal.worms.end <- (ws - 1) + num.comps

  fert.worms <- dat[, fert.worms.start:fert.worms.end] #number of fertile females worms
  
  #fertile worms + male worm -> new birthed mf
  new.in <- (rotate(fert.worms) * ep.in) #need to rotate matrix to each column is multiplied by respective fecundity rate, not each ro
  new.in <- rotate(rotate(rotate(new.in)))
  new.in <- rowSums(new.in)
    

  compartment <- mf.cpt

  compartments_ind <- (mf.st -1) + compartment
  mf.birthed <-rpois(N, new.in*DT)# CHECK
  #current load   
  mf.cur <- dat[,compartments_ind]

  N <- length(treat.vec)
    
  mf.die <- rbinom(N, pmax(0L, mf.cur), mf.mort)    


  mf.mu <- rep(1 - exp(-mu.rates.mf[mf.cpt] * DT), N)


  mf.loss.aged <- rbinom(N, pmax(0L, mf.cur - mf.die), (DT/time.each.comp))

  in.param <- mf.move
  in.param <- rep(1 - exp(-mf.move * DT), N)   
   


  if(mf.cpt == 1){
    mf.out <- mf.cur + mf.birthed - mf.die - mf.loss.aged
  }

  if(mf.cpt > 1){
    mf.out <- mf.cur + aging_in - mf.die - mf.loss.aged  # aging_in comes from previous compartment
  }

 mf.out <- pmax(0L, mf.out)  # safety floor

}
    
    #output, depending on scenario
return(list(mf_out = mf.out, loss_aged = mf.loss.aged))
  }






#' @title
#' calculate mf per skin snip
#' @description
#' function calculates number of mf in skin snip for all people.
#' people are tested for the presence of mf using a skin snip, we assume mf are overdispersed in the skin.
#' @param ss.wt weight of the skin snip
#' @param num.ss number of skin snips taken (default set to 2)
#' @param slope.kmf slope value governing linear relationship between degree of mf overdispersion and adult female worms
#' @param int.kMf initial value governing linear relationship between degree of mf overdispersion and adult female worms
#' @param data data is the matrix tracking age compartments of mf and W per individual
#' @param nfw.start column (first age compartment) in matrix where non-fertile female worms begin
#' @param fw.end column (last age compartment) in matrix where fertile female worms end
#' @param mf.start column (first age compartment) in matrix where mf begin
#' @param mf.end column (last age compartment) in matrix where mf ends
#' @param pop.size human population size
#' @param kM.const.toggle if set to YES then kM is a constant (default = 15)
#'
#' @returns element (1) in list is mean of mf per skin snip; element (2) contains all mf per skin snip for each individual
mf.per.skin.snip <- function(ss.wt, num.ss, slope.kmf, int.kMf, data, nfw.start, fw.end,  ###check vectorization
                             mf.start, mf.end, pop.size, kM.const.toggle)

{

  all.mfobs <- c()

  if(isTRUE(kM.const.toggle)){
    kmf <- 0 * (rowSums(data[,nfw.start:fw.end])) + 15}
  else {
     kmf <- slope.kmf * (rowSums(data[,nfw.start:fw.end])) + int.kMf #rowSums(da... sums up adult worms for all individuals giving a vector of kmfs
  }

  mfobs <- rnbinom(pop.size, size = kmf, mu = ss.wt * (rowSums(data[,mf.start:mf.end])))

  nans <- which(mfobs == 'NaN'); mfobs[nans] <- 0

  if(num.ss > 1)

  {

    tot.ss.mf <- matrix(, nrow = length(data[,1]), ncol = num.ss) # error?
    tot.ss.mf[,1] <- mfobs

    for(j in 2 : (num.ss)) #could be vectorized

    {

      temp <- rnbinom(pop.size, size = kmf, mu = ss.wt * (rowSums(data[,mf.start:mf.end])))

      nans <- which(temp == 'NaN'); temp[nans] <- 0

      tot.ss.mf[,j] <- temp

    }

    mfobs <- rowSums(tot.ss.mf)

  }

  mfobs <- mfobs / (ss.wt * num.ss)

  list(mean(mfobs), mfobs)

}

#' @title
#' mf prevalence by strata
#' @description
#' calculates mf prevalence in people based on a skin snip for different population strata (age strata, lower & upper ages)
#' @param ss.in takes mf per skin snip count object for each individual to convert to binary variable for prevalence
#' @param main.dat main matrix contain age of each individual (to ensure only calculate prevalence based o age from which skin snips are taken)
#' @param lwr.age lower age to measure prevalence from
#' @param upr.age lower age to measure prevalence to
#'
#' @returns value for prevalence
prevalence.for.age <- function(ss.in, main.dat, lwr.age=5, upr.age=81)

{
  inds <- which(main.dat[,2] >= lwr.age & main.dat[,2] < upr.age)

  out <- length(which(ss.in[[2]][inds] > 0)) / length(inds)

  return(out)
}



#' @title
#' mf prevalence by strata
#' @description
#' calculates mf prevalence in people based on a skin snip for different population strata (age and sex)
#' @param ss.in takes mf per skin snip count object for each individual to convert to binary variable for prevalence
#' @param main.dat main matrix contain age of each individual (to ensure only calculate prevalence based o age from which skin snips are taken)
#' @param lwr_age lower age to measure prevalence from
#' @param upr_age lower age to measure prevalence to
#' @param sex  sex of strata to measure prevalence in
#'
#' @returns value for prevalence
prevalence.for.age_sex.strata <- function(ss.in, main.dat, lwr_age, upr_age, sex)

{
  if(sex == "male"){
    sex_ind <- 1
  } else {
    sex_ind <- 0
  }

  inds <- which(main.dat[,2] >= lwr_age & main.dat[,2] < upr_age & main.dat[,3] == sex_ind)

  out <- length(which(ss.in[[2]][inds] > 0)) / length(inds)

  return(out)
}


#' @title
#' mf prevalence by strata (including compliance)
#' @description
#' calculates mf prevalence in people based on a skin snip for different population strata (age, sex and compliance)
#' @param ss.in takes mf per skin snip count object for each individual to convert to binary variable for prevalence
#' @param main.dat main matrix contain age of each individual (to ensure only calculate prevalence based o age from which skin snips are taken)
#' @param lwr_age lower age to measure prevalence from
#' @param upr_age lower age to measure prevalence to
#' @param sex  sex of strata to measure prevalence in
#'
#' @returns value for prevalence
prevalence.for.age_sex_compl.strata <- function(ss.in, main.dat, lwr_age, upr_age, sex, compliance)

{
  if(sex == "male"){
    sex_ind <- 1
  } else {
    sex_ind <- 0
  }

  inds <- which(main.dat[,2] >= lwr_age & main.dat[,2] < upr_age)

  out <- length(which(ss.in[[2]][inds] > 0)) / length(inds)

  return(out)
}

calculate_mf_stats_across_age_groups <- function(stat_type, temp_mf, main_dat, age_groups) {
  output_mf_data <- rep(NA, length(age_groups))
  if (stat_type == "prevalence") {
    for (age_group_index in 1:length(age_groups)) {
      age_group <- age_groups[[age_group_index]]
      output_mf_data[age_group_index] <- prevalence.for.age(
        ss.in = temp_mf,
        main.dat = main_dat,
        lwr.age = age_group[1], upr.age = age_group[2]
      )
    }
  } else if (stat_type == "intensity") {
    for (age_group_index in 1:length(age_groups)) {
      age_group <- age_groups[[age_group_index]]
      output_mf_data[age_group_index] <- mean(
        temp_mf[[2]][which(
          main_dat[, 2] >= age_group[1] & main_dat[, 2] < age_group[2]
        )]
      )
    }
  }
  return(output_mf_data)
}
