#EPIONCHO-IBM
#30/04/2020
#Jonathan Hamley


#' @title
#' Run EPIONCHO-IBM epidemiological model with or without interventions
#' @description
#' Runs the individual-based Onchocerciasis transmission model (EPIONCHO-IBM), based on Hamley et al. 2019
#' ep.equi.sim will run one repeat of the model, typically the mean of 500 repeats is required
#' model must be run to equilibrium (100 years), before treatment can begin
#' treatment is parameterised based on ivermectin
#' if the mf prevalence is zero 50 years after the final treatment, we assume elimination has occured
#' code is available which saves the equilibrium and receives it as an input
#'
#' @param time.its number of iterations (this input is a single value)
#' @param ABR annual biting rate (this input is a single value)
#' @param N.in human population size
#' @param DT timestep (must be one day e.g., 1/366)
#' @param treat.int treatment interval in years e.g., 1 is every 1 year, 0.5 is every 6 months (this input is a single value)
#' @param treat.timing specific timing of treatment rounds (default is NA)
#' @param treat.prob total population coverage (this input is a single value between 0 - 1)
#' @param treat.prob.variable variable total population coverage (value between 0 - 1) specified as a vector of coverages for each treatment round
#' @param give.treat takes 1 (MDA) or 0 (no MDA)
#' @param treat.start iteration where treatment starts
#' @param treat.stop iteration where treatment stops
#' @param pnc proportion of population which never receive treatment (single input value between 0 - 1)
#' @param min.mont.age minimum age for giving a skin snip (single input value, default is 5)
#' @param vector.control.strt start year for vector control
#' @param vector.control.duration duration in years for vector control
#' @param vector.control.efficacy efficacy of vector (proportion of original ABR value)
#' @param delta.hz.in (if gam.dis.in is 0.2, 0.3, or 0.4, this will not be used) this is a new user-input for density-dependence in humans (proportion of L3 larvae establishing/ developing to adult stage within human host, per bit, when ATP tends to 0)
#' @param delta.hinf.in (if gam.dis.in is 0.2, 0.3, or 0.4, this will not be used) this is a new user-input for density-dependence in humans (proportion of L3 larvae establishing/ developing to adult stage within human host, per bit, when ATP tends to infinity)
#' @param c.h.in (if gam.dis.in is 0.2, 0.3, or 0.4, this will not be used) this is a new user-input for density-dependence in humans (severity of transmission intensity - dependent parasite establishment within humans)
#' @param gam.dis.in this is a new user-input for individual-level exposure in humans (overdispersion parameter k_E determining degree of individual exposure heterogeneity; default value is 0.3)
#' @param Q This is the parameter that controls the Ratio of Male Exposure vs Female Exposure. Default is 1.2
#' @param kM.const.toggle specifies whether kM set to constant (if yes, then overdispersion in mf in skin set to 15)
#' @param output_age_groups age groups to output prevalences at
#' @param run_equilibrium specify whether to run to equilibrium first
#' @param equilibrium equilibrium input given to continue model
#' @param print_progess print the current status of the model run
#' @param morbidity_module whether to use and output morbidity in the model. Default is "YES"
#' @param morbidities the morbidities to track, options from c("SevereItch", "RSD", "Atrophy", "HG", "Depig", "Blindness", "VI", "OAE"). By default all are enabled
#' @param correlated_compliance correlated compliance structure specified, to toggle on, set to "YES". Default is "NO".
#' @param comp.correlation this is the probability associated with the compliance correlation (rho). Default is 0. correlated_compliance parameter must be "YES"
#' @param treat.switch An array designating which drug (IVM or MOX) to be used. Must be a list of the same size as treatment time. Default is NA. ex: c("IVM", "MOX", "IVM")
#' @param treat.type If you are not using treat.switch, this should be the drug to be used for treatment. Can be "MOX" or "IVM". Default is "IVM"
#' @param ov16_store_times The timesteps at which you want to store the ov16 status for each individual. Provide in terms of actual time (i.e year 10). Default is c(). ex: c(50, 100)
#' @param ov16_diagnostic_adjustment The diagnostic adjustment to provide for the adjusted ov16 seroprevalence output. Must be a vector of exactly 2 values. Default is c(0.8, 0.99)
#'
#' @export

ep.equi.sim <- function(time.its,
                        DT.in = 1 / 366,
                        ABR,
                        N.in=440,
                        treat.int,
                        treat.timing,
                        treat.prob,
                        treat.prob.variable = NA,
                        give.treat,
                        treat.start,
                        treat.stop,
                        pnc,
                        min.mont.age,
                        vector.control.strt,
                        vector.control.duration,
                        vector.control.efficacy,
                        delta.hz.in=0.186, # these inputs are new (matt) for testing DD
                        delta.hinf.in=0.003,
                        c.h.in=0.005,
                        gam.dis.in=0.3,
                        Q = 1.2,
                        k0_in = 0.0054, # 95% CI = (0.0012, 0.0097)
                        k1_in = 0.1459, # 95% CI = (0.0548, 0.237)
                        kM.const.toggle = FALSE,
                        output_age_groups = list(c(0, 5), c(5, 10), c(10, 15), c(15, 20), c(20, 30), c(30, 50), c(50, 81)),
                        run_equilibrium,
                        equilibrium,
                        print_progress = TRUE,
                        morbidity_module = "YES",
                        morbidities = c("SevereItch", "RSD", "Atrophy", "HG", "Depig", "Blindness", "VI", "OAE"),
                        correlated_compliance = "NO",
                        comp.correlation = 0,
                        treat.switch = NA,
                        treat.type = "IVM",
                        ov16_store_times = c(),
                        ov16_diagnostic_adjustment = c(0.8, 0.99)) {

  output_age_groups_as_strings <- rep("", length(output_age_groups))
  for (output_age_group_index in 1:length(output_age_groups)) {
    output_age_groups_as_strings[output_age_group_index] <- paste(
      output_age_groups[[output_age_group_index]], collapse = "_"
    )
  }
  # ====================== #
  # Set-up time parameters #

  DT <- DT.in # default timestep
  if ((treat.type == "MOX") || (all(!is.na(treat.switch)) && "MOX" %in% treat.switch)) {
    # half-day timestep required for MOX dynamics
    DT <- min(DT.in / 2, (1 / 366) / 2)
  }

  time.its <- round(time.its / (DT))
  year_its <- seq(0, time.its, 1/DT)

  # if(give.treat == 1) #calculate timesteps at which treatment is given
  # {times.of.treat.in <- seq(treat.start, treat.stop - (treat.int / DT), treat.int / DT)}
  # else {times.of.treat.in <- 0}

  if (gam.dis.in == 0.2) {
    delta.hz.in =  0.385
    delta.hinf.in = 0.003
    c.h.in = 0.008
  } else if (gam.dis.in == 0.3) {
    delta.hz.in =  0.186
    delta.hinf.in = 0.003
    c.h.in = 0.005
  } else if (gam.dis.in == 0.4) {
    delta.hz.in = 0.118
    delta.hinf.in = 0.002
    c.h.in = 0.004
  }

  if(give.treat == 1)
  {
    treat.stop <- round(treat.stop / (DT))
    if(treat.start >= 1) {treat.start <-  round( (treat.start) / (DT)) + 1}
    if(treat.start == 0) {treat.start <-  1}
  }

  if(length(ov16_store_times) > 0) {
    ov16_store_times <- round(ov16_store_times / (DT))
  } else {
    if (give.treat == 1) {
      ov16_store_times <- c(treat.start-1, treat.stop, treat.stop+(0.5/DT), treat.stop+(1/DT))
    } else {
      ov16_store_times <- c(time.its-1)
    }
  }

  print("Ov16 full population store times")
  print(ov16_store_times)

  # vector control #
  if(!is.na(vector.control.strt)){

    vc.iter.strt <- round(vector.control.strt / (DT))
    vc.iter.stp <- vc.iter.strt + round(vector.control.duration / (DT))

  } else {
    vc.iter.strt <- NA
    vc.iter.stp <- NA
    vector.control.duration <- NA
  }


  # ================ #
  # hard coded parms #

  # density dep pars (worm establishment in humans)
  delta.hz <- delta.hz.in
  delta.hinf <- delta.hinf.in
  c.h <- c.h.in

  m = ABR * ((1/104) / 0.63) # matt: m = vector to host ratio (V/H) ?; ABR = beta * V/H, where V/H = ABR/beta or V/H = ABR/(h/g), or V/H = ABR * (g/h) which you see here (note, V/H is inferred from the ABR - KEY INPUT for adjusting endemicity level of EPIONCHO-IBM sims)
  beta = 0.63 / (1/104) # matt: beta = per blackfly biting rate on humans (h/g; where h is the human blood index, g is duration of the gonotrophic cycle)
  mu.v = 26
  int.mf = 0
  sigma.L0 = 52
  a.H = 0.8
  g = 0.0096
  a.v = 0.39
  real.max.age = 80 #no humans live longer than 80 years
  N = N.in #human population size
  mean.age = 50 #mean human age (before truncation)
  int.L3 = 0.03; int.L2 = 0.03; int.L1 = 0.03
  lambda.zero = 0.33 # (matt:) per-capita rate that female worms lose their fertility (W_FF) & return to non-fertile state (W_FN)
  omeg = 0.59 # (matt:) per-capita rate that female worms progress from non-fertile (W_FN) to fertile (W_FF)
  delta.vo = 0.0166 # matt : delta V0, density dependence when microfilarae tend to 0
  c.v = 0.0205 # matt: severity of constraining density-dep larval development (per dermal mf) : Table F Supp
  num.mf.comps = 21; num.comps.worm = 21 #number of age classes for mf and adult worms (matt: c_max ?)
  time.each.comp.worms = 1; time.each.comp.mf = 0.125; mf.move.rate = 8.133333 #for aging in parasites (matt: time.each.comp.worms = q_W & time.each.comp.mf = q_M in supp table E )
  int.worms=1 #initial number of worms in each worm age compartment
  ss.wt = 2; num.ss = 2 #skin snip parameters (matt weight and number)
  slope.kmf = 0.0478 # matt: parameter associated with decreasing degree of aggregation of skin mf with increasing no. of adult female worms (slope in linear model -> k_M = 0.0478 * W_F + 0.313)
  int.kMf = 0.313 # matt: parameter associated with decreasing degree of aggregation of skin mf with increasing no. of adult female worms (inital y value in linear model -> k_M = 0.0478 * W_F + 0.313)

  sex.rat = 0.5 #sex ratio (matt: inidividual assigned a sex randomly - equal probability psi_F = 0.5, psi_m = 0.5)

  nuone = 201.6189; nutwo = 207.7384 #movement of fly parasite life stages

  mu.w1 = 0.09953; mu.w2 = 6.00569 #parameters controlling age-dependent mortality in adult worms (matt: these are y_l = y_w and d_l = d_w in equation S6/S7 & Table E)
  mu.mf1 = 1.089; mu.mf2 = 1.428 #parameters controlling age-dependent mortality in mf (matt: these are y_l = y_m and d_l = d_m in equation S6/S7 & Table E)
  fec.w.1 = 70; fec.w.2 = 0.72 #parameters controlling age-dependent fecundity in adult worms (matt: fec.w.1 = F and fec.w.2 = G in Supp table E)
  l3.delay = 10; dt.days = DT * 366 #delay in worms entering humans and joining the first adult worm age class (dt.days = DT.in*366)

  # Morbidity Probabilities
  SI_prob = 0.1636701; RSD_prob = 0.04163095
  Atrp_prob = 0.036; HG_prob = 0.018
  Depigm_prob = 0.059
  avg_age_dpm_atr_hg = 22.876


  # Treatment parameters #

  if(all(is.na(treat.switch)) && (treat.type == "IVM")){
    lam.m = 32.4; phi = 19.6 #effects of ivermectin (matt: embryostatic effect - lam.m is the max rate of treatment-induced sterility; phi is the rate of decay of this effect - Table G in Supp)
    cum.infer = 0.345 # permanent infertility in worms due to ivermectin (irreversible sterlising effect- "global") - this could be changed as a macrofilaricidal to 0.9 (90%)
    up = 0.0096; kap = 1.25 #effects of ivermectin (matt: parameters u (up) and k (kap) define the microfilaricidal effect curve, u = finite effect follwoed by decline (rebound) = k - table G in Supp)

  }

  if(all(is.na(treat.switch)) && treat.type == "MOX"){
    #print("default pars")

    lam.m = 462; phi = 4.83 #effects of moxidectin (matt: embryostatic effect - lam.m is the max rate of treatment-induced sterility; phi is the rate of decay of this effect - Table G in Supp)
    cum.infer = 0.345 # permanent infertility in worms due to ivermectin (irreversible sterlising effect- "global") - this could be changed as a macrofilaricidal to 0.9 (90%)
    up = 0.04; kap = 1.82 #effects of Moxidectin from Kura et al. 2023(matt: parameters u (up) and k (kap) define the microfilaricidal effect curve

  }

  # Exposure parameters #

  # gam.dis = 0.3 #individual level exposure heterogeneity (matt: shape par in gamma dist, K_E)
  gam.dis <- gam.dis.in # when specifying user input (K_E)
  E0 = 0; q = 0; m.exp = 1.08; f.exp = 0.9; age.exp.m = 0.007; age.exp.f = -0.023 #age-dependent exposure to fly bites age.exp.m or .f = alpha_m or alpha_f)


  # print error messages when incorrect inputs / combination of inputs specified #

  if(give.treat == 1)
  {
    if(all(is.na(treat.prob.variable))){ if(treat.prob > 1) stop('treatment probability must be between 0 & 1') }
    else {if(any(treat.prob.variable > 1)) stop('treatment probability must be between 0 & 1')}

    if(treat.stop > time.its) stop('not enough time for requested MDA duration')

    #times.of.treat.in <- seq(treat.start, treat.stop - (treat.int / DT), treat.int / DT)

    if(all(!is.na(treat.timing))) {treat.timing <- treat.timing + ((treat.start - (1/DT)+round((1/DT)/366)) * DT)} # # 1 day dt
    if(all(is.na(treat.timing)))
      {times.of.treat.in <- seq(treat.start, treat.stop, treat.int / DT)}
    else {times.of.treat.in <- round((treat.timing) / (DT)) + 1}

    print(paste(length(times.of.treat.in), 'MDA rounds to be given', sep = ' '))

    print('MDA given at')
    print(paste(round(times.of.treat.in * DT, digits = 2), 'yrs', sep = ''))

    print(times.of.treat.in)

    print('Target coverage at each round')
    if(all(!is.na(treat.prob.variable))) {print(paste(treat.prob.variable*100, "%", sep = ''))}
    else{print(paste(treat.prob*100, "%", sep = ''))}

    print('Treatment to be given at each round')
    if(all(is.na(treat.switch)) && (treat.type == "IVM")){
      print("IVM only")
    }
    if(all(is.na(treat.switch)) && treat.type == "MOX"){
      print("MOX only")
    }
    if(all(!is.na(treat.switch))){
      print(treat.switch)
    }

    print('ABR is')
    print(paste(ABR, 'bites person-1 yr-1 at endemic equilibria'))
    if(!is.na(vector.control.strt)){
      print('ABR changes to the following due to vector control')
      print(paste((ABR - (ABR * vector.control.efficacy)), 'bites person-1 yr-1 at', vector.control.strt, 'yrs'))}

  }
  else {times.of.treat.in <- 0; print('no MDA to be simulated')}

  if(is.logical(run_equilibrium) == FALSE) stop('no input indicating if model needs to be run to equilibrium')


  #columns to set to zero when an individual dies
  cols.to.zero <- seq(from = 1, to = (6 + num.mf.comps + 3*num.comps.worm))
  cols.to.zero <- cols.to.zero[-c(1,5, 6)] #compliance, L2 and L3 do not become zero when an individual dies

  #columns, used to perform operations on different worm and mf compartments
  tot.worms <- num.comps.worm*3
  num.cols <- 6 + num.mf.comps + tot.worms
  worms.start <- 7 + num.mf.comps

  nfw.start <- 7 + num.mf.comps + num.comps.worm # start of infertile worms
  fw.start <- nfw.start + num.comps.worm # start of fertile worms
  fw.end <- num.cols # end of fertile worms
  mf.start <- 7
  mf.end <- 6 + num.mf.comps

  #age-dependent mortality and fecundity rates of parasite life stages

  age.cats <- seq(0, 20, length = num.comps.worm) #up to 20 years old (assume all worms die after age 20 years)

  mort.rates.worms <- weibull.mortality(DT = DT, par1 = mu.w1, par2 = mu.w2, age.cats = age.cats)

  fec.rates.worms <- 1.158305 * fec.w.1 / (fec.w.1 + (fec.w.2 ^ -age.cats) - 1) #no DT - Rk4

  age.cats.mf <- seq(0, 2.5, length = num.mf.comps) #up to 2.5 years old (assume mf die after age 2.5 years)

  #DT not relevent here because RK4 is used to calculate change in mf
  mort.rates.mf <- weibull.mortality(DT = 1, par1 = mu.mf1, par2 = mu.mf2, age.cats = age.cats.mf)


  # create output variables
  mf_prevalence_outputs <- matrix(NA, nrow = time.its - 1, ncol = length(output_age_groups) + 1)
  colnames(mf_prevalence_outputs) <- c("prev", paste0("prev_", output_age_groups_as_strings))

  never_treated_values <- rep(NA, time.its - 1)
  has_been_treated <- rep(FALSE, N)
  mfp_recorded_year_tracker <- rep(NA, time.its - 1)

  mf_intensity_outputs <- matrix(NA, nrow = time.its - 1, ncol = length(output_age_groups) + 1)
  colnames(mf_intensity_outputs) <- c(
    "mean.mf.per.snip",
    paste0("mean.mf.per.snip", output_age_groups_as_strings)
  )

  worm_burden_outputs <- matrix(NA, nrow = time.its - 1, ncol = (length(output_age_groups) + 1) * 3)
  colnames(worm_burden_outputs) <- c(
    "male_worm_burden",
    paste0("male_worm_burden", output_age_groups_as_strings),
    "fertile_female_worm_burden",
    paste0("fertile_female_worm_burden", output_age_groups_as_strings),
    "infertile_female_worm_burden",
    paste0("infertile_female_worm_burden", output_age_groups_as_strings)
  )

  L3_vec <- rep(NA, time.its - 1)
  l3_prev_blackflies <- rep(NA, time.its - 1)
  ABR_recorded <- rep(NA, time.its - 1)
  coverage.recorded <- rep(NA, time.its - 1)

  # Ov16 Variables
  ov16_timetrend_outputs <-  matrix(NA, nrow = time.its - 1, ncol = (length(output_age_groups) + 1) * 2)
  colnames(ov16_timetrend_outputs) <- c(
    "ov16_seroprevalence_no_seroreversion",
    paste0("ov16_seroprevalence_no_seroreversion", output_age_groups_as_strings),
    "ov16_seroprevalence_finite_seroreversion",
    paste0("ov16_seroprevalence_finite_seroreversion", output_age_groups_as_strings)
  )
  ov16_timetrend_outputs_adj <-  matrix(NA, nrow = time.its - 1, ncol = (length(output_age_groups) + 1) * 2)
  colnames(ov16_timetrend_outputs_adj) <- c(
    "ov16_seroprevalence_no_seroreversion_adj",
    paste0("ov16_seroprevalence_no_seroreversion_adj", output_age_groups_as_strings),
    "ov16_seroprevalence_finite_seroreversion_adj",
    paste0("ov16_seroprevalence_finite_seroreversion_adj", output_age_groups_as_strings)
  )

  ov16_seropositive_no_seroreversion <- rep(0, N)
  ov16_seropositive_finite_serorevert <- rep(0, N)
  ov16_indiv_matrix <- matrix(0, nrow = N, ncol = length(ov16_store_times) * 5)
  matrix_index <- 1

  if (morbidity_module == "YES") {
    oae_outputs_dt <- matrix(NA, nrow = time.its - 1, ncol = 7)
    colnames(oae_outputs_dt) <- c(
      "prev_OAE", "OAE_incidence", "OAE_incidence_under_5",
      "OAE_incidence-5_10", "OAE_incidence-11_15",
      "OAE_incidence_M", "OAE_incidence_F"
    )
    morbidity_prevalence_outputs <- matrix(
      NA, nrow = time.its - 1, ncol = length(morbidities) * (length(output_age_groups) + 1)
    )
    morbidity_column_names <- c()
    for (morbidity in morbidities) {
      morbidity_column_names <- c(
        morbidity_column_names,
        paste0(morbidity, "_prev"),
        paste0(morbidity, "_prev", output_age_groups_as_strings)
      )
    }
    colnames(morbidity_prevalence_outputs) <- morbidity_column_names

    # extract probabilities for each condition
    eye_disease_file_path <- file.path("inst", "extdata", "eye_disease_probabilties_updated.rds")
    if (!file.exists(eye_disease_file_path)) {
      eye_disease_file_path <- system.file("extdata", "eye_disease_probabilties_updated.rds", package = "EPIONCHO.IBM")
    }
    eye.disease.probs <- readRDS(eye_disease_file_path) # estimated from Little et al. 2004

    eye.dis.probs <- eye.disease.probs$fit
  }

  # ========================================================================================================== #
  # create inital age distribution and simulate stable age distribution (where equilibrium input not provided) #
  if (isTRUE(run_equilibrium)) {
    cur.age <- rep(0, N)

    #(the approach below must be used, drawing human lifespans from an exponential distribution eventually leads to a non-exponential distribution)
    for(i in 1 : 75000) #if at equilibrium you saved the age at which inds die and simulated further, you should get an exponential distribution
    {
      cur.age <- cur.age + DT

      death.vec <- rbinom(N, 1, (1/mean.age) * DT) # Matt: human mortality (constant with age) - no. of deaths at time step t is a random variable drawn from binomial distribution (N = human pop size?)

      cur.age[which(death.vec == 1)] <- 0 #set individuals which die to age 0
      cur.age[which(cur.age >= real.max.age)] <- 0 #all individuals >= maximum imposed age die (matt: distribution truncated to prevent excessively long life spans - a_max)
    }


    ex.vec <- rgamma(N, gam.dis, gam.dis) #individual level exposure to fly bites (matt: individual-specific exposure factor assigned at birth - drawn from gamma dist, with shape par (K_E = gam.dis, and rate par set to this))

    ###############################################
    #matrix for delay in L3 establishment in humans
    num.delay.cols <- l3.delay * (28 / dt.days)
    l.extras <- matrix(0, ncol= num.delay.cols, nrow= N)
    inds.l.mat <- seq(2,(length(l.extras[1,]))) #for moving columns along with time

    ################################################
    #L1 delay in flies
    l1.delay <- rep(int.L1, N)

    ###############################################
    #matrix for tracking mf for L1 delay
    num.mfd.cols <- 4 / dt.days
    mf.delay <- matrix(int.mf, ncol= num.mfd.cols, nrow= N)
    inds.mfd.mats <- seq(2,(length(mf.delay[1,])))

    ###############################################
    #matrix for exposure (to fly bites) for L1 delay
    num.exp.cols <- 4 / dt.days
    exposure.delay <- matrix(ex.vec, ncol= num.exp.cols, nrow= N)
    inds.exp.mats <- seq(2,(length(exposure.delay[1,])))

    #matrix for first timestep, contains all parasite values, human age, sex and compliance

    all.mats.temp <- matrix(nrow=N, ncol=num.cols)

    all.mats.temp[,  (worms.start) : num.cols] <- int.worms

    all.mats.temp[, 4] <- int.L1

    all.mats.temp[, 5] <- int.L2

    all.mats.temp[, 6] <- int.L3

    all.mats.temp[, 7 : (7 + (num.mf.comps-1))] <- int.mf

    all.mats.temp[,1] <- rep(0, N) #column used during treatment
    all.mats.temp[,2] <- cur.age

    #assign sex to humans

    sex <- rbinom(N, 1, sex.rat) # matt: randomly assigned (binomal dist) with equal probability (e.g., psi_F = 0.5, psi_m = 0.5)

    all.mats.temp[,3] <- sex

    #non-compliant people
    non.comp <- ceiling(N * pnc)
    out.comp <- rep(0, N)
    s.comp <- sample(N, non.comp)
    out.comp[s.comp] <- 1
    all.mats.temp[,1] <- out.comp

    treat.vec.in <- rep(NA, N) #for time since treatment calculations

    if (morbidity_module == "YES") {

      # OAE set-up #
      infected_at_all <- rep(1, N)
      age_to_samp <- sample(seq(3, 15, DT), size = N, replace = TRUE)
      OAE <- rep(1, N)
      tested_OAE <- rep(0, N)
      check_ind <- c()

      # data and function to obtain OAE probability for a given mf count
      Chesnais_dat <- data.frame(prob = c(0.0061, 0.0439, 0.0720, 0.0849, 0.1341, 0.1538, 0.20),
                                  mf = c(0, 3, 13, 36, 76, 151, 200))

      OAE_probs <- OAE_mfcount_prob_func(dat = Chesnais_dat)

      # ========================= #
      # Skin disease (OSD) set-up #

      num.cols.morb2 <- 53 # 4 cols for age, sex, compliance, mf count; 8 cols for whether undergo trial for each condition; 8 cols disease state 1/0

      # define Morbidity dataframe #
      all.morb.temp <- data.frame(
        Age = cur.age,
        Sex = sex,
        TrueMFCount = 0,
        ObservedMFCount = 0,
        AtrophySampleAges = 0,
        HangingGroinSampleAges = 0,
        DepigSampleAges = 0,
        Day3SevereItchStatus = 0,
        Day3RSDStatus = 0,
        ToTestSevereItch = 0,
        ToTestRSD = 0,
        ToTestAtrophy = 0,
        ToTestHG = 0,
        ToTestDepig = 0,
        SevereItchStatus = 0,
        RSDStatus = 0,
        AtrophyStatus = 0,
        HGStatus = 0,
        DepigStatus = 0
      )

      age_to_samp_vec_reversible <- seq(0+1/366, 79+1/366, 1) # between 5 and 80, sample once year year of age
      age_to_samp_vec_nonreversible <- seq(0+1/366, 79+1/366, 1) # between 5 and 80, sample once year year of age


      # set up reversible sequela (SI and RSD) 3-day delay matrix #
      num.days.cols <- 3 + 1 # day 0 to day 3 (so 4)
      sequela.postive.mat <- matrix(0, ncol = num.days.cols, nrow = N)
      inds.sequela.mat <- seq(2,(length(sequela.postive.mat[1,]))) # for moving columns along with time

      sequela.postive.mat1 <- sequela.postive.mat # for severe itch
      sequela.postive.mat2 <- sequela.postive.mat # for RSD

      # ========================= #
      # Eye disease set-up        #

      num.cols.blindness <- 12 # 9 + 2 for lagged age + 1 for updated blindness status

      # define Morbidity matrix #
      # set other cols to 0 at i = 1
      all.blind.temp <- data.frame(
        Age = cur.age,
        Sex = sex,
        TrueMFCount = 0,
        ObservedMFCount = 0,
        BlindnessProb = 0,
        BlindnessStatus = 0,
        BlindnessPending = 0,
        BlindnessCountdown = 730
      )
      # TODO? #
      #   lagged vector for blindness / visual impairment (start filling from 2 years)  #
    }
  }


  # =============================================================================================#
  # set-up main structures for tracking human/parasite stages if equilibrium dataframe given     #
  if(isFALSE(run_equilibrium))

  {
    if(is.list(equilibrium) == FALSE) stop('equilibrium condition not in correct format')
    ex.vec <- equilibrium[[2]] #exposure

    ###############################################
    #matrix for delay in l3 establishment in humans
    num.delay.cols <- l3.delay * (28 / dt.days)
    l.extras <- equilibrium[[4]]
    inds.l.mat <- seq(2,(length(l.extras[1,]))) #for moving columns along with time

    ################################################
    l1.delay <- equilibrium[[6]]

    ###############################################
    #matrix for tracking mf for l1 delay
    num.mfd.cols <- 4 / dt.days
    mf.delay <- equilibrium[[5]]
    inds.mfd.mats <- seq(2,(length(mf.delay[1,])))

    ###############################################
    #matrix for exposure for L1 delay
    num.exp.cols <- 4 / dt.days

    #matrix for first time step
    all.mats.temp <- equilibrium[[1]]


    exposure.delay <- equilibrium[[8]]
    inds.exp.mats <- seq(2,(length(exposure.delay[1,])))

    temp <- mf.per.skin.snip(ss.wt = 2, num.ss = 2, slope.kmf = 0.0478, int.kMf = 0.313, data = all.mats.temp, nfw.start, fw.end,
                             mf.start, mf.end, pop.size = N, kM.const.toggle)

    treat.vec.in <- equilibrium[[3]]
    has_been_treated <- equilibrium[[9]]

    if(morbidity_module == "YES"){
      OAE <- equilibrium$morbidity_equilibrium_outputs$all_OAE_equilibrium_outputs$OAE
      age_to_samp <- equilibrium$morbidity_equilibrium_outputs$all_OAE_equilibrium_outputs$age_to_samp
      tested_OAE <- equilibrium$morbidity_equilibrium_outputs$all_OAE_equilibrium_outputs$tested_OAE
      infected_at_all <- equilibrium$morbidity_equilibrium_outputs$all_OAE_equilibrium_outputs$infected_at_all
      check_ind <- equilibrium$morbidity_equilibrium_outputs$all_OAE_equilibrium_outputs$check_ind
      tot_ind_ep_samp <- equilibrium$morbidity_equilibrium_outputs$all_OAE_equilibrium_outputs$tot_ind_ep_samp
      OAE_probs <- equilibrium$morbidity_equilibrium_outputs$all_OAE_equilibrium_outputs$OAE_probs

      #==================#
      all.morb.temp <- equilibrium$morbidity_equilibrium_outputs$other_morbidity_equilibrium_outputs$all.morb.updated

      all.blind.temp <- equilibrium$morbidity_equilibrium_outputs$other_morbidity_equilibrium_outputs$all.blind.updated

      age_to_samp_vec_reversible <- seq(0+1/366, 79+1/366, 1) # between 5 and 80, sample once year year of age
      age_to_samp_vec_nonreversible <- seq(0+1/366, 79+1/366, 1) # between 5 and 80, sample once year year of age

      sequela.postive.mat1 <- equilibrium$morbidity_equilibrium_outputs$other_morbidity_equilibrium_outputs$sequela_pos_mat_1 # 3-day delay sequelae matrix (SI)
      sequela.postive.mat2 <- equilibrium$morbidity_equilibrium_outputs$other_morbidity_equilibrium_outputs$sequela_pos_mat_2 # 3-day delay sequelae matrix (RSD)

      inds.sequela.mat <- seq(2,(length(sequela.postive.mat1[1,]))) # for moving columns along with time
    }

    # Ov16 Variables
    ov16_seropositive_no_seroreversion <- equilibrium$ov16_equilibrium$ov16_seropositive_no_seroreversion
    ov16_seropositive_finite_serorevert <- equilibrium$ov16_equilibrium$ov16_seropositive_finite_serorevert
  }

  i <- 1

  # ================================================================================================= #
  #               Main loop for running through processes in model with i > 1                         #
  while(i < time.its) #over time

  {
    #print(paste(round(i * DT, digits = 2), 'yrs', sep = ' '))

    # if(isTRUE(print_progress)) {print(paste(round(i * DT, digits = 2), 'yrs', sep=' '))}

    # if(isTRUE(print_progress) & (any(i == year_its))) {print(paste(round(i * DT, digits = 2), 'yrs', sep=' '))}
    # if(isTRUE(print_progress) & (any(i == year_its))) {print(paste(round(i/time.its * 100, digits = 1), '%', sep=' '))}

    if(isTRUE(print_progress) & (any(i == year_its))) {print(paste(round(i * DT, digits = 2), 'yrs;',
                                                                   (paste(round(i/time.its * 100, digits = 1), '%', sep=' '))))}

    #stores mean L3 and adult worms from previous timesteps

    all.mats.cur <- all.mats.temp

    if(any(i == times.of.treat.in)) {
      index.iter.treat <- match(i, times.of.treat.in) # find element where iteration number matches a time in times.of.treat vector

      # extract element from treat.prob.variable depending on specific treatment round (iteration in times.of.treat.in)
      if(all(!is.na(treat.prob.variable))){
        treat.prob <- treat.prob.variable[index.iter.treat] # index prob value from treat.prob.variable vector
      }

      # if IVM/MOX switch included, specify which treatment parameters to use if treatment iteration
      if(all(!is.na(treat.switch))){
        treat.type <- treat.switch[index.iter.treat] # index IVM or MOX value from treat.switch.in vector
        if(treat.type == "IVM") {
          lam.m = 32.4; phi = 19.6 # treatment induced embryostatic parameters
          cum.infer= 0.345 # permanent infertility in worms
          up = 0.0096; kap = 1.25 # microfilaricidal effect curve parameters
          print("IVM parameters updated")
        } else if (treat.type == "MOX") {
          lam.m = 462; phi = 4.83 # treatment induced embryostatic parameters
          cum.infer= 0.345 # permanent infertility in worms
          up = 0.04; kap = 1.82 # microfilaricidal effect curve parameters
          print("MOX parameters updated")
        } else {
          print("ERROR - either IVM or MOX not specified in treatment switch vector at this iteration")
        }
      }
    }


    # to track variable coverage
    if(i >= treat.start & i <= treat.stop & give.treat == 1) {
      coverage.upd <- treat.prob
    } else {
      coverage.upd <- 0
    }

    # which individuals will be treated if treatment is given (old compliance approach)
    if(i >= treat.start & give.treat == 1 & correlated_compliance != "YES") {
      cov.in <- os.cov(all.dt = all.mats.cur, pncomp = pnc, covrg = treat.prob, N = N)
    }


    # ============================= #
    # for new compliance structure;
    # initialize probability of treatment values (pTreat) for each individual if first round
    # subsequent rounds: check to see if coverage or correlation par values have changed since last treatment,
    # if so, need to edit pTreat values
    # always check for zero values in pTreat for subsequent rounds
    if(correlated_compliance == "YES" & any(i == times.of.treat.in)){
      probneverTreat <- pnc
      cov <- treat.prob

        # first MDA round
        if(i == times.of.treat.in[1]){

          # specify neever treat individuals
          compliance.mat <- matrix(nrow=N, ncol=6) # col 1 = age, col 2 = never_treat,
                                                   # col 3 = probability of treatment, col 4 = to be treated in this round
          compliance.mat[,2] = generateNeverTreat(N = N, probneverTreat) # never treat col (mat[,1])

          # individual probability of treatment values
          cov = treat.prob # whatever the coverage of this MDA is
          rho = comp.correlation # whatever the correlation of this MDA is
          compliance.mat[,3] = initializePTreat(N = N, cov, rho) # initialize pTreat - correlation for each individual (mat[,2])

          # record this value for previous coverage #
          prevCov = cov # set prevCov to coverage value used
          prevRho = rho # set prevRho to correlation value used
        }

        # subsequent MDA rounds

        if (i %in% times.of.treat.in[-1]) {
          if((prevCov != treat.prob) | (prevRho != comp.correlation)){
            # 1) check and update/redraw any zero values introduced in pTreat for individuals since last MDA round
            #compliance.mat[,3] = checkForZeroPTreat(pTreat = compliance.mat[,2], prevCov, prevRho)
            compliance.mat[,3] = checkForZeroPTreat(pTreat = compliance.mat[,3], prevCov, prevRho)

            # 2) assign everyone a new/updated pTreat value for the next MDA round if cov and/or rho have changed
            cov = treat.prob
            rho = comp.correlation

            compliance.mat[,3] = editPTreat(pTreat = compliance.mat[,3], cov, rho)

            prevCov = cov
            prevRho = rho

          }

          # check for zero pTreat values since last MDA regardless of whether new cov/rho values
          compliance.mat[,3] = checkForZeroPTreat(pTreat = compliance.mat[,3], prevCov, prevRho)
        }

      # specify if individuals are to be treated in this round in compliance.mat (column 6)
      eligible_out <- check_eligibility(comp.mat = compliance.mat, all.dt = all.mats.cur, minAgeMDA = 5, maxAgeMDA = 80)

      compliance.mat <- eligible_out[[1]] # extract updated compliance matrix

      cov.in <- compliance.mat[,6] # this is vector of individuals to be treated from compliance mat, to feed into change.worm.per.ind.treat

      # Count the number of treated hosts
      hostsEligibleAge <- compliance.mat[,4]
      eligible_hosts <- eligible_out[[2]]
      hostsTreated <- length(eligible_hosts)
      CovEligibles = hostsTreated / hostsEligibleAge * 100
      CovTotal = hostsTreated / N * 100

    }
    if (i >= treat.start & give.treat == 1 ) {
      has_been_treated <- has_been_treated | (cov.in == 1)
    }
    #         TODO: Update print output here with new coverage vals from compliance structure                       #

    # # update when no MDA round
    # if(isTRUE(print_progress) & (any(i == year_its + 1)) & any(i == times.of.treat.in))
    # {print(paste(round(CovEligibles, digits = 3), '% coverage eligibles', paste(round(CovTotal, digits = 3), '% coverage total', sep=' ')))}

    #sex and age dependent exposure, mean exposure must be 1, so ABR is meaningful

    mls <- which(all.mats.cur[,3] == 1) # index males
    fmls <- which(all.mats.cur[,3] == 0) # index females

    e.f <- 1/(sex.rat*(Q-1) + 1) # 0.9
    e.m <- Q * e.f # 1.08

    s.a.exp <- rep(0, N)

    alpha.m <- exp(-age.exp.m * (all.mats.cur[mls, 2]))
    gamma.m <- 1/mean(alpha.m)
    s.a.exp[mls] <- alpha.m * e.m * gamma.m

    alpha.f <- exp(-age.exp.f * (all.mats.cur[fmls, 2]))
    gamma.f <- 1/mean(alpha.f)
    s.a.exp[fmls] <- alpha.f * e.f * gamma.f

    norm.ex.vec <- ex.vec * (1 / mean(ex.vec)) #normalize so mean = 1 (matt: normalising the indvidual-specific exposure from line 565)

    tot.ex.ai <- s.a.exp * norm.ex.vec # matt: combine sex/age specific exposure + individual specific exposure (total exposure to blackfly bites)
    tot.ex.ai <- tot.ex.ai * (1 / mean(tot.ex.ai)) #normalize so mean = 1

    #increase age (for next time step)

    all.mats.temp[,2] <- (all.mats.cur[,2]) + DT #increase age for all individuals

    death.vec <- rbinom(N, 1, (1/mean.age) * DT) #select individuals to die

    to.die <- which(death.vec == 1)

    at.ab.max <- which(all.mats.temp[,2] >= real.max.age)

    to.die <- c(to.die, at.ab.max)

    to.die <- unique(to.die) #may have repeated indivudals i.e selected by binom and >80

    ##################
    #delay calculations
    ##################

    #there is a delay in new parasites entering humans (from fly bites) and entering the first adult worm age class

    new.worms.m <- c()
    new.worms.nf <- c()

    new.worms.m <- rbinom(N, size = l.extras[,length(l.extras[1,])], prob = 0.5) #draw males and females from last column of delay matrix
    new.worms.nf <- l.extras[,length(l.extras[1,])] - new.worms.m

    #move individuals along
    l.extras[,inds.l.mat] <- l.extras[,(inds.l.mat-1)]

    #mean number of L3 in fly population
    L3.in <- mean(all.mats.cur[, 6])


    # change m based on ABR change due to vector control if called (during vector control duration iteration period)
    if(!is.na(vector.control.strt)){

      if (i >= vc.iter.strt && i < vc.iter.stp) {

      ABR_updated <- ABR - (ABR * vector.control.efficacy) # proportional reduction in ABR (x efficacy) during VC

      m = ABR_updated * ((1/104) / 0.63) # update m
      }
    }

    # change m back to original after VC finishes (can change when this occurs i.e, one year after VC ends)
    if(!is.na(vector.control.strt)){

      if (i >= vc.iter.stp) {

        m = ABR * ((1/104) / 0.63) # update m
      }
    }


    # to track #
    if (!is.na(vector.control.strt) && i >= vc.iter.strt && i < vc.iter.stp) {

      ABR_upd <- ABR_updated

      } else {

      ABR_upd <- ABR
    }


    #rate of infections in humans
    #delta.hz, delta.hinf, c.h are density dependence parameters, expos is the exposure of each person to bites
    nw.rate <- Wplus1.rate(delta.hz, delta.hinf, c.h, L3 = L3.in, m ,
                           beta, expos = tot.ex.ai, DT)


    new.worms <- rpois(N, nw.rate) #total new establishing L3 for each individual

    l.extras[,1] <- new.worms


    for(k in 1 : num.comps.worm) #go through each adult worm compartment

    {

      if(k == 1) {from.last <- rep(0, N)} #create vector for worms coming from previous compartment (needs to be 0 when k ==1)


       res.w1 <- change.worm.per.ind1(treat.vec = treat.vec.in, lambda.zero = lambda.zero, DT = DT, omeg = omeg,
                                      ws = worms.start, compartment = k, total.dat = all.mats.cur, mort.rates = mort.rates.worms,
                                      time.each.comp = time.each.comp.worms, new.worms.m = new.worms.m, w.f.l.c = from.last,
                                      num.comps = num.comps.worm)

       res.w.treat <- change.worm.per.ind.treat(give.treat = give.treat, iteration = i, treat.start = treat.start,
                                                times.of.treat = times.of.treat.in, treat.stop = treat.stop,
                                                onchosim.cov = cov.in, treat.vec = treat.vec.in, DT = DT, cum.infer = cum.infer,
                                                lam.m = lam.m, phi = phi, N = res.w1[[3]], mort.fems = res.w1[[2]], lambda.zero.in = res.w1[[1]])

       res.w2 <- change.worm.per.ind2(DT = DT, time.each.comp = time.each.comp.worms, compartment = k, new.worms.nf.fo = new.worms.nf, w.f.l.c = from.last,
                                      N = res.w1[[3]], cur.Wm.nf = res.w1[[4]], mort.fems = res.w.treat[[3]], cur.Wm.f = res.w1[[5]], omeg = res.w1[[7]],
                                      male.tot.worms = res.w1[[8]], worm.loss.males = res.w1[[9]],
                                      lambda.zero.in = res.w.treat[[1]], treat.vec = res.w.treat[[2]])


      res <- res.w2 # (matt: re-label the final result output to res so do not have to change res below)

      from.last <- res # (matt: re-label the final result output to res so do not have to change res below)

      # update male worms in matrix for compartment k

      all.mats.temp[, (6 + num.mf.comps + k)] <- res[[1]]

      # update females worms in matrix

      all.mats.temp[, (6 + num.mf.comps + num.comps.worm + k)] <- res[[3]] # infertile, num.comps.worm skips over males
      all.mats.temp[, (6 + num.mf.comps + 2*num.comps.worm + k)] <- res[[4]] # fertile, num.comps.worm skips over males and infertile females


    }

    if(give.treat == 1 & i >= treat.start) {treat.vec.in <- res[[7]]} #treated individuals




#######CHANGING THIS ONE CALLING FUNCTION####
  #  for(mf.c in 1 : num.mf.comps)

 #   {
#
#      res.mf <- change.micro(dat = all.mats.cur, num.comps =num.comps.worm, mf.cpt = mf.c,
#                             num.mf.comps = num.mf.comps, ws=worms.start, DT=DT, time.each.comp = time.each.comp.mf,
#                             mu.rates.mf = mort.rates.mf, fec.rates = fec.rates.worms, mf.move.rate = mf.move.rate, up = up, kap = kap, iteration = i,
#                             treat.vec = treat.vec.in, give.treat = give.treat, treat.start = treat.start)#
#
 #     all.mats.temp[, 6 + mf.c] <- res.mf
 #   }
#replacing with :
aging_in_from_prev <- rep(0, N)  # compartment 1 receives no aging-in
for(mf.c in 1 : num.mf.comps) {
  res.mf <- change.micro(dat = all.mats.cur, num.comps =num.comps.worm, mf.cpt = mf.c,
                             num.mf.comps = num.mf.comps, ws=worms.start, DT=DT, time.each.comp = time.each.comp.mf,
                             mu.rates.mf = mort.rates.mf, fec.rates = fec.rates.worms, mf.move.rate = mf.move.rate, up = up, kap = kap, iteration = i,
                             treat.vec = treat.vec.in, give.treat = give.treat, treat.start = treat.start,
 aging_in = aging_in_from_prev)
  aging_in_from_prev <- res.mf$loss_aged   # pass this compartment's aged-out mf to next
  all.mats.temp[, 6 + mf.c] <- res.mf$mf_out
}











    #inputs for delay in L1
    exp.delay.temp <- exposure.delay[, length(exposure.delay[1,])]
    mf.delay.temp <- mf.delay[, length(mf.delay[1,])]
    l1.delay.temp <- l1.delay #L1 from previous timestep

    #move values along
    exposure.delay[, inds.exp.mats] <- exposure.delay[, (inds.exp.mats -1)]
    mf.delay[, inds.mfd.mats] <- mf.delay[, (inds.mfd.mats - 1)]

    #update L1, L2 and L3

    #total number of mf in each person
    mf.temp <- rowSums(all.mats.cur[, 7 : (6 + num.mf.comps)]) #sum mf over compartments, mf start in column 7

    all.mats.temp[, 4] <- calc.L1(beta, mf = mf.temp, mf.delay.in = mf.delay.temp, expos = tot.ex.ai, delta.vo, c.v, nuone, mu.v, a.v, expos.delay = exp.delay.temp)
    all.mats.temp[, 5] <- calc.L2(nuone, L1.in = l1.delay.temp, mu.v, nutwo, mf = mf.delay.temp, a.v, expos = exp.delay.temp)
    all.mats.temp[, 6] <- calc.L3(nutwo, L2.in = all.mats.cur[, 5], a.H, g, mu.v, sigma.L0)

    #new values for delay parts
    l1.delay <- all.mats.temp[, 4]
    mf.delay[, 1] <- rowSums(all.mats.cur[, 7 : (6 + num.mf.comps)])
    exposure.delay[, 1] <- tot.ex.ai

    # ======================================== #
    # Morbidities                              #
    # ======================================== #

    if(morbidity_module == "YES"){
      #===========================#
      #     OAE module funcs      #

      OAE_out1 <- find_indiv_OAE_func(dat = all.mats.temp, mf.start = mf.start, mf.end = mf.end, worms.start = worms.start, nfw.start = nfw.start, fw.end = fw.end,
                  infected_at_all = infected_at_all, age_to_samp = age_to_samp, OAE = OAE, tested_OAE = tested_OAE, check_ind = check_ind) # step 1

      infected_at_all = OAE_out1[[2]] # updated (when i = 1)
      check_ind = OAE_out1[[3]] # updated (when i = 1)

      temp.mf <- mf.per.skin.snip(ss.wt = 2, num.ss = 2, slope.kmf = 0.0478, int.kMf = 0.313, data = all.mats.temp, nfw.start, fw.end,
                              mf.start, mf.end, pop.size = N, kM.const.toggle)

      OAE_out2 <- new_OAE_cases_func(temp.mf = temp.mf, tot_ind_ep_samp = OAE_out1[[1]], OAE_probs = OAE_probs, dat = all.mats.temp,
                                      OAE = OAE, tested_OAE = tested_OAE, age_groups = output_age_groups) # step 2

      oae_outputs_dt[i,] <- OAE_out2$output_vals
      OAE = OAE_out2[[3]]
      tested_OAE = OAE_out2[[4]]


      # ============ #
      # Skin disease and blindness #

      if(i == 1) {
        all.morb.updated <- all.morb.temp
        all.blind.updated <- all.blind.temp
      } else if (i > 1) {

        # =================== #
        #     Skin disease    #

        temp.mf <- mf.per.skin.snip(ss.wt = 2, num.ss = 2, slope.kmf = 0.0478, int.kMf = 0.313, data = all.mats.temp, nfw.start, fw.end,
                                    mf.start, mf.end, pop.size = N, kM.const.toggle)

        all.morb.updated <- find_indiv_totest_func(dat = all.mats.temp, mf.start = mf.start, mf.end = mf.end,
                                                   morb.mat.tmp = all.morb.updated, temp_mf = temp.mf,
                                                   age_to_samp_vec_reversible = age_to_samp_vec_reversible,
                                                   age_to_samp_vec_nonreversible = age_to_samp_vec_nonreversible)

        all.morb.updated <- new_cases_morbidity_func(morb.mat.tmp = all.morb.updated,
                                                      SI_probs = SI_prob, RSD_probs = RSD_prob,
                                                      Atrp_probs = Atrp_prob,
                                                      Hg_probs = HG_prob, Depigm_probs = Depigm_prob,
                                                      avg_age = avg_age_dpm_atr_hg)


        # un hash below if want to shift to probabilities based on mf count #

        # all.morb.updated <- new_cases_morbidity_func2(morb.mat.tmp = all.morb.updated, temp_mf = temp.mf,
        #                                               SI_probs = SI_prob, RSD_probs = RSD_prob,
        #                                               Atrp_probs = Atrp_prob,
        #                                               Hg_probs = Hg_prob, Depigm_probs = Depigm_prob)

        reversible.morb.updated <- update_reversible_sequela_func(sequela.postive.mat1 = sequela.postive.mat1,
                                       sequela.postive.mat2 = sequela.postive.mat2,
                                       inds.sequela.mat = inds.sequela.mat, morb.mat.tmp = all.morb.updated)

        sequela.postive.mat1 <- reversible.morb.updated[[1]] # updated SI 3-day delay matrix
        sequela.postive.mat2 <- reversible.morb.updated[[2]] # updated RSD 3-day delay matrix
        all.morb.updated <- reversible.morb.updated[[3]] # updated morb.mat

        # ============ #
        # eye disease  #
        all.blind.updated <- find_indiv_totest_func2(dat = all.mats.temp, mf.start = mf.start, mf.end = mf.end,
                                                     morb.mat.tmp = all.blind.updated,
                                                     age_to_samp_vec_nonreversible = age_to_samp_vec_nonreversible)

        all.blind.updated <- new_cases_morbidity_func2(morb.mat.tmp = all.blind.updated, temp.mf = temp.mf,
                                                         blind.probs = eye.dis.probs)
      }

      morbidity_prevalence_outputs[i,] <- morbidity_prev_func(
        non_blindness_morb_dat = all.morb.updated,
        blidness_dat = all.blind.updated,
        oae_status = OAE,
        N = N,
        age_groups = append(list(c(5, 81)), output_age_groups)
      )
    }

    #save prevalence at current timestep
    temp.mf <- mf.per.skin.snip(ss.wt = 2, num.ss = 2, slope.kmf = 0.0478, int.kMf = 0.313, data = all.mats.temp, nfw.start, fw.end,
                                mf.start, mf.end, pop.size = N, kM.const.toggle)

    mfp_recorded_year_tracker[i] <- i * DT

    mf_prevalence_outputs[i, ] <- calculate_mf_stats_across_age_groups(
      "prevalence",
      temp_mf = temp.mf,
      main_dat = all.mats.temp,
      age_groups=append(list(c(min.mont.age, 81)), output_age_groups)
    )

    never_treated_values[i] <- 1 - mean(has_been_treated, na.rm=TRUE)

    mf_intensity_outputs[i, ] <- calculate_mf_stats_across_age_groups(
      "intensity",
      temp_mf = temp.mf,
      main_dat = all.mats.temp,
      age_groups=append(list(c(min.mont.age, 81)), output_age_groups)
    )

    worm_burden_outputs[i, ] <- calculate_worm_burden_across_age_groups(
      main_dat = all.mats.temp,
      age_groups = append(list(c(min.mont.age, 81)), output_age_groups),
      male_start = worms.start,
      infertile_female_start = nfw.start,
      fertile_female_start = fw.start,
      fertile_female_end = fw.end
    )

    mf.per.skin.snp.out <- temp.mf[[2]] #to extract out mf per skin snip for each individual?

    L3_vec[i] <- mean(all.mats.temp[, 6])

    l3_prev_blackflies[i] <- calculate_l3_prevalence_in_blackflies(
      l3_per_blackfly = mean(all.mats.temp[, 6]),
      k0 = k0_in,
      k1 = k1_in
    )

    ABR_recorded[i] <- ABR_upd # tracking changing ABR

    coverage.recorded[i] <- coverage.upd # track changing coverage if specified

    # ======================================================================================== #
    # new individual exposure for newborns, clear rows for new borns & update various matrices #

    # Ov16 exposure checks
    any_larvae <- (rowSums(l.extras) > 0)
    any_worms <- (rowSums(all.mats.temp[,worms.start:fw.end]) > 0)
    mating_worm <- ((rowSums(all.mats.temp[,worms.start:nfw.start])) > 0 & (rowSums(all.mats.temp[, fw.start:fw.end]) > 0))
    mating_worm_any_mf <- (mating_worm & (rowSums(all.mats.temp[,mf.start:mf.end]) > 0))

    ov16_seropositive_no_seroreversion <- determine_serostatus(
      exposure_array = mating_worm_any_mf,
      curr_array = ov16_seropositive_no_seroreversion
    )
    ov16_seropositive_finite_serorevert <- determine_serostatus(
      exposure_array = mating_worm_any_mf, curr_array = ov16_seropositive_finite_serorevert, do_serorevert=TRUE,
      seroreversion_arrays = list("any_larvae_arr" = any_larvae, "any_worms_arr" = any_worms)
    )

    ov16_timetrend_outputs[i, ] <- calculate_seroprevalence_across_age_groups(
      ov16_seropositive_no_seroreversion, ov16_seropositive_finite_serorevert,
      ages = all.mats.temp[,2], age_groups = append(list(c(min.mont.age, 81)), output_age_groups),
      diagnostic_adjustments = c(1, 1)
    )
    ov16_timetrend_outputs_adj[i, ] <- calculate_seroprevalence_across_age_groups(
      ov16_seropositive_no_seroreversion, ov16_seropositive_finite_serorevert,
      ages = all.mats.temp[,2], age_groups = append(list(c(min.mont.age, 81)), output_age_groups),
      diagnostic_adjustments = ov16_diagnostic_adjustment
    )

    # new individual exposure for newborns, clear rows for new borns
    if(length(to.die) > 0)
    {
      ex.vec[to.die] <- rgamma(length(to.die), gam.dis, gam.dis)

      l.extras[to.die, ] <- 0 #establishing adult worms

      mf.delay[to.die, 1] <- 0 #individual dies so no contribution to L1s at this timestep

      l1.delay[to.die] <- 0

      treat.vec.in[to.die] <- NA

      all.mats.temp[to.die, cols.to.zero] <- 0 #set age, sex and parasites to 0 (includes L1, but not L2 L3)
      all.mats.temp[to.die, 3] <- rbinom(length(to.die), 1, 0.5) #draw sex
      has_been_treated[to.die] <- FALSE

      # Ov16 reset dead individuals
      ov16_seropositive_no_seroreversion[to.die] <- 0
      ov16_seropositive_finite_serorevert[to.die] <- 0

      if(correlated_compliance == "YES" & any(i > times.of.treat.in)){
        compliance.mat[to.die, 3] <- 0
      } # if individual dies update pTreat to 0 in compliance matrix

      if(morbidity_module == "YES"){
        # OAE #
        infected_at_all[to.die] <- 0 # index those individuals to die as no longer ever infected

        age_to_samp[to.die] <- sample(seq(3, 15, DT), size = length(to.die), replace = TRUE) # for those individuals set to die, resample

        OAE[to.die] <-  0 # index those individuals to die as no longer with OAE

        tested_OAE[to.die] <-  0 # index those individuals to die as no longer tested

        # ==================== #
        #    Skin disease      #

        # set cols to zero for those individuals that die and replaced by newborn
        all.morb.updated[to.die, ] <- 0 #set all to 0

        # update sequela 3-day delay matrices if any individual dies
        sequela.postive.mat1[to.die, ] <- 0

        sequela.postive.mat2[to.die, ] <- 0

        # ==================== #
        #    Eye disease       #
        all.blind.updated[to.die, c("Age", "Sex", "TrueMFCount", "ObservedMFCount", "AgeToSampleEyeDist", "ToTestBlindness", "BlindnessProb", "BlindnessStatus", "BlindnessPending")] <- 0 # Reset dataframe to 0
        all.blind.updated[to.die, c("BlindnessCountdown")] <- 730
      }
    }

    # Ov16 add individual data to matrix
    if(!is.na(match(i, ov16_store_times))) {
      ov16_indiv_matrix[,5*matrix_index-4] <- all.mats.temp[,2]
      ov16_indiv_matrix[,5*matrix_index-3] <- all.mats.temp[,3]
      ov16_indiv_matrix[,5*matrix_index-2] <- as.integer(temp.mf[[2]] > 0)
      ov16_indiv_matrix[,5*matrix_index-1] <- ov16_seropositive_no_seroreversion
      ov16_indiv_matrix[,5*matrix_index] <- ov16_seropositive_finite_serorevert
      matrix_index <- matrix_index + 1
    }

    # ========================================================================================================== #
    # update sequela matrix if any individual has reached 3rd day with SI/RSD (4th column) - reversible OSD only #

    if(morbidity_module == "YES"){
      # SI 3-day matrix #
      current_day3_SI <- which(sequela.postive.mat1[,4] == 1)
      if(length(current_day3_SI) > 0) {
        sequela.postive.mat1[current_day3_SI, ] <- 0
      }

      # RSD 3-day matrix #
      current_day3_RSD <- which(sequela.postive.mat2[,4] == 1)
      if(length(current_day3_RSD) > 0) {
        sequela.postive.mat2[current_day3_RSD, ] <- 0
      }
    }

    i <- i + 1
  }

  general_outputs <- list(
    'mf_prev' = mf_prevalence_outputs[,'prev'], 'mf_intens' = mf_intensity_outputs[,'mean.mf.per.snip'],
    "ov16_seroprevalence_no_seroreversion" = ov16_timetrend_outputs[,"ov16_seroprevalence_no_seroreversion"],
    "ov16_seroprevalence_finite_seroreversion" = ov16_timetrend_outputs[,"ov16_seroprevalence_finite_seroreversion"],
    'L3' = L3_vec, 'ABR' = ABR, 'all_infection_burdens' = all.mats.temp, "Ke" = gam.dis,
    'blackfly_l3_prevalence' = l3_prev_blackflies, 'years' = mfp_recorded_year_tracker, 'all_mf_prevalence_age_grouped' = mf_prevalence_outputs,
    'all_mf_intensity_age_grouped' = mf_intensity_outputs, 'ov16_indiv_matrix' = ov16_indiv_matrix,
    "ov16_timetrend_outputs" = ov16_timetrend_outputs, 'ov16_timetrend_outputs_adj' = ov16_timetrend_outputs_adj,
    "worm_burden_outputs" = worm_burden_outputs, 'ABR_recorded' = ABR_recorded, 'coverage.recorded' = coverage.recorded,
    'percent_never_treated' = never_treated_values
  )

  if (morbidity_module == "YES"){
    morbidity_outputs <- list(
      'oae_incidence_outputs' = oae_outputs_dt,
      'all_morbidity_prevalence_outputs' = morbidity_prevalence_outputs
    )
    general_outputs <- c(general_outputs, morbidity_outputs)
  }

  # add outputs needed to restart sims
  if(isTRUE(run_equilibrium))
  {
    equilibrium_outputs <- list(all.mats.temp, ex.vec, treat.vec.in, l.extras, mf.delay, l1.delay, ABR, exposure.delay, has_been_treated)
    # Ov16 Equilibrium Values
    ov16_equib <- list(ov16_seropositive_no_seroreversion, ov16_seropositive_finite_serorevert)
    names(ov16_equib) <- c('ov16_seropositive_no_seroreversion', 'ov16_seropositive_finite_serorevert')
    equilibrium_outputs[["ov16_equilibrium_outputs"]] <- ov16_equib

    if (morbidity_module == "YES") {
      equilibrium_outputs[["morbidity_equilibrium_outputs"]] <- list(
        "all_OAE_equilibrium_outputs" = list(
          'OAE' = OAE, 'age_to_samp' = age_to_samp, 'tested_OAE' = tested_OAE, 'infected_at_all' = infected_at_all,
          'check_ind' = OAE_out1[[3]], 'tot_ind_ep_samp' = OAE_out1[[1]], 'OAE_probs' = OAE_probs,
          'oae_incidence_outputs' = oae_outputs_dt
        ),
        "other_morbidity_equilibrium_outputs" = list(
          'all.morb.updated' = all.morb.updated, 'all.blind.updated' = all.blind.updated,
          'sequela_pos_mat_1' = sequela.postive.mat1, 'sequela_pos_mat_2' = sequela.postive.mat2,
          'all_morbidity_prev_outputs' = morbidity_prevalence_outputs
        )
      )
    }
    general_outputs[["all_equilibrium_outputs"]] = equilibrium_outputs
  }
  return(general_outputs)
}


