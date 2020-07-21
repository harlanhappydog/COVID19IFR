# Sampler in the space of
# alpha = f1(eta, phi) = log(phi) + eta
# beta = f2(eta, phi) = log(phi) - eta
#
# where eta = cloglog_infectionrate[k] = cloglogr. (r = infectionrate for compactness)
# We sample beta holding alpha constant
# and provide -log(Det(Jacobian)) adjustment
sampler_infrate_phi <- nimbleFunction(
  name = 'sampler_infrate_phi',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    adaptive            <- if(!is.null(control$adaptive))            control$adaptive            else TRUE
    adaptInterval       <- if(!is.null(control$adaptInterval))       control$adaptInterval       else 20
    adaptFactorExponent <- if(!is.null(control$adaptFactorExponent)) control$adaptFactorExponent else 0.8
    scale               <- if(!is.null(control$scale))               control$scale               else 1
    ## node list generation
    cloglogrNode <- target[1]
    phiNode <- target[2]
    calcNodes <- model$getDependencies(target)
    ## This assumes there is no graph relation directly between cloglogr and phi
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    ## numeric value generation
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    scaleHistory  <- c(0, 0)   ## scaleHistory
    acceptanceHistory  <- c(0, 0)   ## scaleHistory
    if(nimbleOptions('MCMCsaveHistory')) {
      saveMCMChistory <- TRUE
    } else saveMCMChistory <- FALSE
    optimalAR     <- 0.44
    gamma1        <- 0
    ## checks
  },
  run = function() {
    current_cloglogr <- model[[cloglogrNode]] ## same as eta
    current_phi <- model[[phiNode]]
    current_logphi <- log(current_phi)
    current_alpha <- current_logphi + current_cloglogr
    current_beta <- current_logphi - current_cloglogr
    proposal_beta <- rnorm(1, current_beta, sd = scale)
    proposal_alpha <- current_alpha
    proposal_cloglogr <- 0.5 * (proposal_alpha - proposal_beta)
    proposal_logphi <- 0.5 * (proposal_alpha + proposal_beta)
    proposal_phi <- exp(proposal_logphi)
    model[[cloglogrNode]] <<- proposal_cloglogr
    model[[phiNode]] <<- proposal_phi
    logDetJacobianRatio <- proposal_logphi - current_logphi

    ## First we check for invalid priors to avoid further calculation
    logMHR <- calculateDiff(model, target)
    if(logMHR == -Inf) {
      nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
      jump <- FALSE
    } else {
      ## Next we calculate the full log Metropolis-Hastings ratio
      logMHR <- logMHR + calculateDiff(model, calcNodesNoSelf) + logDetJacobianRatio
      jump <- decide(logMHR)
      if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
      else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
    if(adaptive)     adaptiveProcedure(jump)
  },
  ## From here down is adaptation and reset code that is unmodified from sampler_RW except for removing reflective part
  methods = list(
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        if(saveMCMChistory) {
          setSize(scaleHistory, timesAdapted)                 ## scaleHistory
          scaleHistory[timesAdapted] <<- scale                ## scaleHistory
          setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
          acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
        }
        gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        ## If there are upper and lower bounds, enforce a maximum scale of
        ## 0.5 * (upper-lower).  This is arbitrary but reasonable.
        ## Otherwise, for a poorly-informed posterior,
        ## the scale could grow without bound to try to reduce
        ## acceptance probability.  This creates enormous cost of
        ## reflections.
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    getScaleHistory = function() {       ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(scaleHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },          
    getAcceptanceHistory = function() {  ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(acceptanceHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },          
    reset = function() {
      scale <<- scaleOriginal
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      if(saveMCMChistory) {
        scaleHistory  <<- c(0, 0)    ## scaleHistory
        acceptanceHistory  <<- c(0, 0)
      }
      gamma1 <<- 0
    }
  )
)

# Sampler in the space of
# alpha = f1(eta, gamma) = exp(eta) + exp(gamma)
# beta = f2(eta, phi) = exp(eta) - exp(gamma)
#
# where eta = cloglog_infectionrate[k] = cloglogr. (r = infectionrate for compactness)
# and   gamma = cloglog_IFR[k]
# We sample beta holding alpha constant
# and provide -log(Det(Jacobian)) adjustment
sampler_infrate_IFR <- nimbleFunction(
  name = 'sampler_infrate_IFR',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    adaptive            <- if(!is.null(control$adaptive))            control$adaptive            else TRUE
    adaptInterval       <- if(!is.null(control$adaptInterval))       control$adaptInterval       else 20
    adaptFactorExponent <- if(!is.null(control$adaptFactorExponent)) control$adaptFactorExponent else 0.8
    scale               <- if(!is.null(control$scale))               control$scale               else 1
    ## node list generation
    cloglogrNode <- target[1] ## eta
    cloglogIFRNode <- target[2] ## gamma
    calcNodes <- model$getDependencies(target)
    ## This assumes there is no graph relation directly between cloglogr and phi
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    ## numeric value generation
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    scaleHistory  <- c(0, 0)   ## scaleHistory
    acceptanceHistory  <- c(0, 0)   ## scaleHistory
    if(nimbleOptions('MCMCsaveHistory')) {
      saveMCMChistory <- TRUE
    } else saveMCMChistory <- FALSE
    optimalAR     <- 0.44
    gamma1        <- 0
    ## checks
  },
  run = function() {
    current_cloglogr <- model[[cloglogrNode]] ## same as eta
    current_cloglogIFR <- model[[cloglogIFRNode]]
    current_expeta <- exp(current_cloglogr)
    current_expgamma <- exp(current_cloglogIFR)
    current_alpha <- current_expeta + current_expgamma
    current_beta <- current_expeta - current_expgamma
    proposal_beta <- rnorm(1, current_beta, sd = scale)
    proposal_alpha <- current_alpha
    proposal_expeta <- 0.5 * (proposal_alpha + proposal_beta)
    proposal_expgamma <- 0.5 * (proposal_alpha - proposal_beta)
    proposal_cloglogr <- log(proposal_expeta)
    proposal_cloglogIFR <- log(proposal_expgamma)
    model[[cloglogrNode]] <<- proposal_cloglogr
    model[[cloglogIFRNode]] <<- proposal_cloglogIFR
    logDetJacobianRatio <- -(proposal_cloglogr + proposal_cloglogIFR) + (current_cloglogr + current_cloglogIFR)

    ## First we check for invalid priors to avoid further calculation
    logMHR <- calculateDiff(model, target)
    if(logMHR == -Inf) {
      nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
      jump <- FALSE
    } else {
      ## Next we calculate the full log Metropolis-Hastings ratio
      logMHR <- logMHR + calculateDiff(model, calcNodesNoSelf) + logDetJacobianRatio
      jump <- decide(logMHR)
      if(jump) nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
      else     nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
    if(adaptive)     adaptiveProcedure(jump)
  },
  ## From here down is adaptation and reset code that is unmodified from sampler_RW except for removing reflective part
  methods = list(
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        if(saveMCMChistory) {
          setSize(scaleHistory, timesAdapted)                 ## scaleHistory
          scaleHistory[timesAdapted] <<- scale                ## scaleHistory
          setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
          acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
        }
        gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        ## If there are upper and lower bounds, enforce a maximum scale of
        ## 0.5 * (upper-lower).  This is arbitrary but reasonable.
        ## Otherwise, for a poorly-informed posterior,
        ## the scale could grow without bound to try to reduce
        ## acceptance probability.  This creates enormous cost of
        ## reflections.
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    getScaleHistory = function() {       ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(scaleHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },          
    getAcceptanceHistory = function() {  ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(acceptanceHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },          
    reset = function() {
      scale <<- scaleOriginal
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      if(saveMCMChistory) {
        scaleHistory  <<- c(0, 0)    ## scaleHistory
        acceptanceHistory  <<- c(0, 0)
      }
      gamma1 <<- 0
    }
  )
)
