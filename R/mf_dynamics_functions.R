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

#change.micro.stochastic, if change rest of functinos; this is the stochastic version
change.micro <- function(mf.st = 7, give.treat= give.treat, mu.rates.mf= mort.rates.mf,
		 iteration = i,  mf.cpt=mf.c, num.comps=num.comps.worm, ws=worms.start,
                 dat= all.mats.cur,
		 treat.vec=treat.vec.in,
                 DT=DT, time.each.comp=time.each.comp.worms,
		 aging_in = rep(0,N),
		 num.mf.comps = num.mf.comps, fec.rates = fec.rates.worms, mf.move.rate = mf.move.rate,
		 up = up, kap = kap, treat.start = treat.start){
  




  #for both
  #ignore #just notation, keeping it consistent with model_wrappers, but also with object in this function
   ep.in <- fec.rates
   mf.move <- mf.move.rate



 #mp is just length of human population, so to account for treatment we can just use human population definition used for adult worms
   N <- length(treat.vec)

   mf.mu <- rep(mu.rates.mf[mf.cpt], N)

  
  #######FOR TREATMENT SCENARIO#####

   if(give.treat == 1 & iteration >= treat.start){
    #taken from mf_dynamics_function directly
    tao <- ((iteration - 1) * DT) - treat.vec # tao is zero if treatment has been given at this timestep

    mu.mf.prime <- ((tao + up) ^ (- kap)) # additional mortality due to ivermectin treatment

    mu.mf.prime[which(is.na(mu.mf.prime) == TRUE)] <- 0

    mf.mu <- mf.mu + mu.mf.prime #this is the big change

    ####DOES TERATMENT AFFECT THE NUMBER OF FERTILE WORMS??? OR IS IT MICROFILARICIDAL ONLY??
  
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
    
    compartments_ind <- (mf.st -1) + mf.cpt
    #there's mf.st to mf.end
    mf.birthed <-rpois(N, new.in) # new.in updating, maybe change, just trying to make right shape

    
#debug
cat("NAs in new.in / mf.birthed:", sum(is.na(new.in)), "\n")
cat("Length of new.in:", length(new.in), "\n")
cat("new.in range:", range(new.in, na.rm = TRUE), "\n")






    mf.cur <- dat[,compartments_ind]
    #this might be where the big error is^^ check
 

#DEBUG CHECKS
# --- DIAGNOSTICS (remove once fixed) ---
cat("mf.cpt:", mf.cpt, "\n")
cat("compartments_ind:", compartments_ind, "\n")
cat("mf.cur range:", range(mf.cur, na.rm = TRUE), "\n")
cat("NAs in mf.cur:", sum(is.na(mf.cur)), "\n")
cat("Negatives in mf.cur:", sum(mf.cur < 0, na.rm = TRUE), "\n")
cat("mf.mort range:", range(mf.mort, na.rm = TRUE), "\n")
# ---------------------------------------






   
#    mf.mort <- mf.mu
mf.mort <- 1 - exp(-mu.rates.mf[mf.cpt] * DT)
mf.mu <- rep(1 - exp(-mu.rates.mf[mf.cpt] * DT), N)

#was rates not probability, causing errors
mf.mort <- mf.mu #unnecessary, and sloppy notation, but i'll fix later




   
    #mf.mu versus mf.mort?
#    mf.die <- rbinom(N, mf.cur , mf.mort) #death from compartment
#getting errors, need to curb nonneg
mf.die <- rbinom(N, pmax(0L, mf.cur), mf.mort) 


   
    #aging out at rate mf.move
  #  mf.loss.aged <- rbinom(N, (mf.cur - mf.die),rep((DT/time.each.comp), mf.move))
#   remaining <- pmax(0L, mf.cur - mf.die)
 # mf.loss.aged <- rbinom(N, remaining, rep((DT/time.each.comp), mf.move)) 
 
mf.loss.aged <- rbinom(N, pmax(0L, mf.cur - mf.die), (DT/time.each.comp))


   
    ##should i add some checks for <0 values?
#^did


    #now to make mf.age.in stochastic, we have mf.birthed <- new.in, so the age 1 is fine, just need for future
#can use from det defn mf.in = dat[, 6 + mf.cpt], should mimic mf.loss.aged
#should define it as the aging in compartment that still has MF after death and aging
    

 #   mf.age.in <- rep(0,N)
#    temp.in <- which((mf.cur - mf.die - mf.loss.aged) >0)
    #checking for MF in age compartment after death and aging, will deal with age class 1 in a sec
#gotta define a in.param var, check if there is one for MF aging [CHECK] 
#just use mf.move.rate? mf.move
 


in.param <- mf.move
in.param <- rep(1 - exp(-mf.move * DT), N) #hopefully will fix NAs

# mf.age.in[1:N] <- mf.birthed #CHECK WHAT mf.birthed outputs, if it's a vector for all individuals or just one?



#   if(length(temp.in) >0) {
 #      mf.age.in[temp.in] <- rbinom(length(temp.in),
  #            (mf.cur[temp.in] - mf.die[temp.in] - mf.loss.aged[temp.in]),
   #            in.param[temp.in])


#}


    
    #now the changes, by compartment
#    if(mf.cpt  ==1){
 #   mf.out <- mf.cur + mf.birthed - mf.die - mf.loss.aged
    
#debug
#cat("NAs in mf.out before return:", sum(is.na(mf.out)), "\n")


#    trans.mf.out <- which(mf.out < 0)
  #  if(length(trans.mf.out)>0){print('MF negative output3')}
    
  #  }
    
  #  if(mf.cpt >1){
 # 
#
#mf.age.in <- rep(0,N)
#temp.in <- which((mf.cur - mf.die - mf.loss.aged) >0)
#if(length(temp.in) > 0) {
 #       mf.age.in[temp.in] <- rbinom(length(temp.in),
  #          (mf.cur[temp.in] - mf.die[temp.in] - mf.loss.aged[temp.in]),
   #         in.param[temp.in])
   # }


    
    #  mf.out <- mf.cur + mf.age.in - mf.die - mf.loss.aged

#debug
#cat("NAs in mf.out before return:", sum(is.na(mf.out)), "\n")
 #     
  #    trans.mf.out <- which(mf.out < 0)
   #   if(length(trans.mf.out)>0){print('MF negative output4')}
    #  
   # }
    #  }
  
  mf.loss.aged <- rbinom(N, pmax(0L, mf.cur - mf.die), (DT/time.each.comp))

if(mf.cpt == 1){
  mf.out <- mf.cur + mf.birthed - mf.die - mf.loss.aged
}

if(mf.cpt > 1){
  mf.out <- mf.cur + aging_in - mf.die - mf.loss.aged  # aging_in comes from previous compartment
}

mf.out <- pmax(0L, mf.out)  # safety floor
  
  
}  
  
  
  
  
  
  
  
  
  
  
  
  ###########NON-TREATMENT SCENARIO##############
    if(give.treat == 0){
      
      
#  mf.mort <- mf.mu
 
#    mf.mort <- mf.mu
mf.mort <- 1 - exp(-mu.rates.mf[mf.cpt] * DT)
mf.mu <- rep(mf.mort, N)
#was rates not probability, causing errors
mf.mort <- mf.mu #unnecessary, and sloppy notation, but i'll fix later
    





 
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
    

    compartment <- mf.cpt

    compartments_ind <- (mf.st -1) + compartment
    #there's mf.st to mf.end
    mf.birthed <-rpois(N, new.in)# new.in changing to make right shape? check
    

#debug
cat("NAs in new.in / mf.birthed:", sum(is.na(new.in)), "\n")
cat("Length of new.in:", length(new.in), "\n")
cat("new.in range:", range(new.in, na.rm = TRUE), "\n")



    mf.cur <- dat[,compartments_ind]




####FURTHER DEBUGGING
  cat("ncol(dat):", ncol(dat), "\n")
cat("mf.st:", mf.st, "\n")
cat("Column calculation: (mf.st - 1) + mf.cpt =", (mf.st - 1) + mf.cpt, "\n")
cat("First few non-NA columns:", which(colSums(!is.na(dat)) > 0), "\n")




#DEBUG CHECKS
# --- DIAGNOSTICS (remove once fixed) ---
cat("mf.cpt:", mf.cpt, "\n")
cat("compartments_ind:", compartments_ind, "\n")
cat("mf.cur range:", range(mf.cur, na.rm = TRUE), "\n")
cat("NAs in mf.cur:", sum(is.na(mf.cur)), "\n")
cat("Negatives in mf.cur:", sum(mf.cur < 0, na.rm = TRUE), "\n")
cat("mf.mort range:", range(mf.mort, na.rm = TRUE), "\n")
# ---------------------------------------    
    
    #mp is just length of human population, so to account for treatment we can just use human population definition used for adult worms
    N <- length(treat.vec)
    
    #mf.mu versus mf.mort? <- check
   # mf.die <- rbinom(N, mf.cur , mf.mort) #death from compartment
#mf.die getting errors so need to curb it to ensure size nonneg
mf.die <- rbinom(N, pmax(0L, mf.cur), mf.mort)    


mf.mu <- rep(1 - exp(-mu.rates.mf[mf.cpt] * DT), N)




    #aging out at rate mf.move
#    mf.loss.aged <- rbinom(N, (mf.cur - mf.die),rep((DT/time.each.comp), mf.move))
# remaining <- pmax(0L, mf.cur - mf.die)
# mf.loss.aged <- rbinom(N, remaining, rep((DT/time.each.comp), mf.move))  
mf.loss.aged <- rbinom(N, pmax(0L, mf.cur - mf.die), (DT/time.each.comp))

    #now to make mf.age.in stochastic, we have mf.birthed <- new.in, so the age 1 is fine, just need f$
#can use from det defn mf.in = dat[, 6 + mf.cpt], should mimic mf.loss.aged
#should define it as the aging in compartment that still has MF after death and aging
    
    
#    mf.age.in <- rep(0,N)
 #   temp.in <- which((mf.cur - mf.die - mf.loss.aged) >0)
    #checking for MF in age compartment after death and aging, will deal with age class 1 in a sec
#gotta define a in.param var, check if there is one for MF aging [CHECK]
#just use mf.move.rate? mf.move
 

 in.param <- mf.move
 in.param <- rep(1 - exp(-mf.move * DT), N)   
   
#   if(length(temp.in) >0) {
 #      mf.age.in[temp.in] <- rbinom(length(temp.in),
 #             (mf.cur[temp.in] - mf.die[temp.in] - mf.loss.aged[temp.in]),
 #              in.param[temp.in])
       
 
# }   

# mf.age.in[1:N] <- mf.birthed 



 
#    #now the changes, by compartment
 #   if(mf.cpt ==1){


  #  mf.out <- mf.cur + mf.birthed - mf.die - mf.loss.aged

#debug
#cat("NAs in mf.out before return:", sum(is.na(mf.out)), "\n")

 #   trans.mf.out <- which(mf.out < 0)
  #  if(length(trans.mf.out)>0){print('MF negative output1')}
  #  
  #  }
    
  #  if(mf.cpt >1){
  #   mf.age.in <- rep(0, N)
  #  temp.in <- which((mf.cur - mf.die - mf.loss.aged) > 0)
  #  if(length(temp.in) > 0) {
  #      mf.age.in[temp.in] <- rbinom(length(temp.in),
  #          (mf.cur[temp.in] - mf.die[temp.in] - mf.loss.aged[temp.in]),
  #          in.param[temp.in])
  #  }#



   # mf.out <- mf.cur + mf.age.in - mf.die - mf.loss.aged




#debug
#cat("NAs in mf.out before return:", sum(is.na(mf.out)), "\n")


   
 #   trans.mf.out <- which(mf.out < 0)
  #  if(length(trans.mf.out)>0){print('MF negative output2')}
   #   
      
   # }}
    


mf.loss.aged <- rbinom(N, pmax(0L, mf.cur - mf.die), (DT/time.each.comp))

if(mf.cpt == 1){
  mf.out <- mf.cur + mf.birthed - mf.die - mf.loss.aged
}

if(mf.cpt > 1){
  mf.out <- mf.cur + aging_in - mf.die - mf.loss.aged  # aging_in comes from previous compartment
}

mf.out <- pmax(0L, mf.out)  # safety floor


}
    
    #output, depending on scenario
 #

return(list(mf_out = mf.out, loss_aged = mf.loss.aged))

 #  return(mf.out)
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
