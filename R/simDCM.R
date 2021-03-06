
# 17.2 A general function to simulate data under the DCM model



simDCM <- function(nspec = 50, nsite = 100, nrep = 3, nyear = 10,
  mean.psi1 = 0.4, sig.lpsi1 = 1, mu.beta.lpsi1 = 0, sig.beta.lpsi1 = 0,
  range.mean.phi = c(0.8, 0.8), sig.lphi = 1,
  mu.beta.lphi = 0, sig.beta.lphi = 0,
  range.mean.gamma = c(0.2, 0.2), sig.lgamma = 1,
  mu.beta.lgamma = 0, sig.beta.lgamma = 0,
  range.mean.p = c(0.5, 0.5), sig.lp = 1,
  mu.beta.lp = 0, sig.beta.lp = 0,
  range.beta1.season = c(0, 0), range.beta2.season = c(0, 0),
  range.sd.site = c(0, 0), range.sd.survey = c(0, 0), show.plot = TRUE) {
#
# Written by Marc Kery, 28 Nov 2016
#
# Function is based on the dynocc function 'simDynocc.fn' (AHM2, Chap 16) and
# on the community occupancy function 'simComm' (AHM1, Chap 11).
#
# Function to simulate detection/nondetection data under a general
#     dynamic community (site-occ) model, including:
#   * annual variation in the probabilities of patch persistence, colonization
#     and detection is specified by the bounds of a uniform distribution.
#   * species heterogeneity around the means is specified by the SD of a normal
#     distribution and expressed on the logit scale
#   * one covariate is allowed to a parameter (site covariate for psi1,
#     site-year covariate for phi and gamma and site-year-rep for p).
#     Each covariate is allowed to differ among species again according
#     to a logit-normal model of heterogeneity.
#   * Additional detection heterogeneity at the site- or the survey level,
#     with the possibility of a temporal trend in this heterogeneity. E.g.,
#     an annual trend in detection heterogeneity at the site or
#     the survey level is specified by the value in the first and the last year.
#     Hence, range.sd.site = c(0, 1) will result in a linear trend in the
#     magnitude of site heterogeneity in detection from 0 in the first year to
#     1 in the last year.
#   * Additional detection heterogeneity that varies over the season (= occasion)
#     according to a quadratic effect of occasion number (to model phenology of
#     an insect species for instance).
#   * These last two types of detection heterogeneity are not (yet) allowed
#     to be species-specific.
#
# Function arguments:
# -------------------
# *** Sample size arguments ***
# nspec - Number of species (typically called N in AHM book)
# nsite - Number of sites (M)
# nrep - Number of replicate surveys within a year (= season) (J)
# nyear - Number of years (or 'seasons') (T)
#
# *** Arguments for mean parameters for the intercepts ***
# mean.psi1 - average occupancy probability in first year
# range.mean.p - bounds of uniform distribution from which annual p drawn
# range.mean.phi and range.mean.gamma - same for persistence and colonization prob.
# -------------------
# *** Arguments for mean parameters for the slopes ***
# mu.beta.lpsi1, mu.beta.lphi, mu.beta.lgamma, mu.beta.lp - coefficients of
#      covariates in probabilities of initial occupancy, persistence,
#      colonization and detection. These are the probability-scale means of
#      the normal distibutions, from which species-specific slopes are drawn
# -------------------
# *** Args. for species-specific heterogeneity in intercepts and slopes ***
# sig.lpsi1: sd of the normal distribution from which species-specific occupancy
#      intercepts are drawn (centered on logit(mean.psi1)), on logit scale
# sig.beta.lpsi1: sd of the normal distribution from which species-specific
#      slopes are drawn (centered on mu.beta.lpsi1)
# sig.lphi: sd of the normal distribution from which species-specific persistence
#      intercepts are drawn (centered on logit(mean.phi), which are year-specific),
#      on logit scale
# sig.beta.lphi: sd of the normal distribution from which species-specific
#      persistence slopes are drawn (centered on mu.beta.lphi)
# sig.lgamma: sd of the normal distribution from which species-specific
#      colonization intercepts are drawn (centered on logit(mean.gamma),
#      which are year-specific), on logit scale
# sig.beta.lgamma: sd of the normal distribution from which species-specific
#      colonization slopes are drawn (centered on mu.beta.lgamma)
# sig.lp: sd of the normal distribution from which species-specific
#      detection intercepts are drawn (centered on logit(mean.p),
#      which are year-specific), on logit scale
# sig.beta.lp: sd of the normal distribution from which species-specific
#      detection slopes are drawn (centered on mu.beta.lp)
# -------------------
# *** Args. for detection heterogeneity among sites and surveys
#      (this part of the model is NOT species-specific) ***
# range.sd.site: sd of normal distribution to model logit-normal noise in p
#      at the site level in the first and the last year of the simulation.
# range.sd.survey: sd of normal distribution to model logit-normal noise in p
#      at the site/year/rep = 'survey' level, in the first and the last year
# For the sd and error.rate arguments, if the two values in the range are the
#      same, a constant value is assumed over time, while if they are different,
#      a linear trend is assumed over time.
# -------------------
# *** Args. for detection heterogeneity among occasions within a season
#      (this part of model again NOT species-specific)
# range.beta1.season is the range of the annual variation in the linear effect
#     of season (i.e., of month 1-12) on the product of
#     availability and detection linear and quadratic effect of season
# range.beta2.season is the same for the quadratic effect of season
#
# show.plot: if TRUE, plots are produced. Usually set to FALSE when running sims.
#

# Set up arrays needed
spec <- 1:nspec                           # Species
site <- 1:nsite                           # Sites
year <- 1:nyear                           # Years
visit <- 1:nrep                           # Visit
month <- 1:nrep                           # Months (= seasons)
psi <- muZ <- z <- array(dim = c(nsite, nyear, nspec), dimnames =
   list(paste('Site', site, sep = ''), paste('Year', year, sep = ''),
   paste('Spec', spec, sep = ''))) # Occupancy, occurrence
phi <- gamma <- array(NA, dim = c(nsite, (nyear-1), nspec), dimnames =
   list(paste('Site', site, sep = ''), paste('Year', year[-nyear], sep = ''),
   paste('Spec', spec, sep = ''))) # Survival, colonisation
y <- p <- array(NA, dim = c(nsite, nrep, nyear, nspec), dimnames =
   list(paste('Site', site, sep = ''), paste('Visit', visit, sep = ''),
    paste('Year', year, sep = ''), paste('Spec', spec, sep = '')))# Det. hist and p

# Create covariates (same for all species)
Xpsi1 <- matrix(runif(nsite, -2, 2), ncol = 1, dimnames = list(paste('Site',site,sep=''), NULL))  # Site covariate for psi1
Xphi <- array(runif(nsite*nyear, -2, 2), dim = c(nsite,nyear), dimnames =
 list(paste('Site',site,sep=''), paste('Year',year,sep =''))) # Yearly-site cov
Xgamma <- array(runif(nsite*nyear, -2, 2), dim = c(nsite,nyear), dimnames =
 list(paste('Site',site,sep=''), paste('Year',year,sep =''))) # Yearly-site cov
Xp <- array(runif(nsite*nrep*nyear,-2,2),dim=c(nsite,nrep,nyear), dimnames =
   list(paste('Site', site, sep = ''), paste('Visit', visit, sep = ''),
    paste('Year', year, sep = '')))  # Observation cov.

# (1) Simulate all parameter values
# (a) State process parameters
# initial occupancy for all species
mu.lpsi1 <- ifelse(mean.psi1 == '1', 500, qlogis(mean.psi1))
beta0.lpsi <- rnorm(nspec, mu.lpsi1, sig.lpsi1)     # initial occupancy intercept
beta1.lpsi <- rnorm(nspec, mu.beta.lpsi1, sig.beta.lpsi1) # occ. slope on Xpsi1
for(s in 1:nspec){
  psi[,1,s] <- plogis(beta0.lpsi[s] + beta1.lpsi[s] * Xpsi1)    # psi1
}
# persistence and colonization for all species
beta0.lphi <- beta0.lgamma <- array(dim = c(nspec, nyear-1))
mean.phi <- runif(n = nyear-1, min = range.mean.phi[1], max = range.mean.phi[2])
mean.gamma <- runif(n = nyear-1, min = range.mean.gamma[1], max = range.mean.gamma[2])
mu.lphi <- ifelse(mean.phi == '1', 500, qlogis(mean.phi))
mu.lgamma <- ifelse(mean.gamma == '1', 500, qlogis(mean.gamma))
eps.lphi <- rnorm(nspec, 0, sig.lphi) # species effect in logit(phi) intercept
eps.lgamma <- rnorm(nspec, 0, sig.lgamma) # spec effect in logit(gam) intercept
for(t in 1:(nyear-1)){
  beta0.lphi[,t] <- mu.lphi[t] + eps.lphi       # logit(phi) intercept
  beta0.lgamma[,t] <- mu.lgamma[t] + eps.lgamma # logit(gamma) intercept
}
beta1.lphi <- rnorm(nspec, mu.beta.lphi, sig.beta.lphi) # slope of logit(phi) on Xphi
beta1.lgamma <- rnorm(nspec, mu.beta.lgamma, sig.beta.lgamma) # slope of logit(gamma) on Xphi
for(s in 1:nspec){
  for(t in 1:(nyear-1)){
    phi[,t, s] <- plogis(beta0.lphi[s, t] + beta1.lphi[s] * Xphi[,t])
    gamma[,t,s] <- plogis(beta0.lgamma[s, t] + beta1.lgamma[s] * Xgamma[,t])
  }
}


# (b) Observation process parameters
beta0.lp <- array(dim = c(nspec, nyear))
mean.p <- runif(n = nyear, min = range.mean.p[1], max = range.mean.p[2])
mu.lp <- ifelse(mean.p == '1', 500, qlogis(mean.p))
eps.lp <- rnorm(nspec, 0, sig.lp) # species effect in logit(p) intercept
for(t in 1:nyear){
  beta0.lp[,t] <- mu.lp[t] + eps.lp       # logit(p) intercept
}
beta1.lp <- rnorm(nspec, mu.beta.lp, sig.beta.lp) # slope of logit(p) on Xp
beta1 <- runif(n = nyear, min = range.beta1.season[1], max = range.beta1.season[2])
beta2 <- runif(n = nyear, min = min(range.beta2.season), max = max(range.beta2.season))
sd.site <- seq(from = range.sd.site[1], to = range.sd.site[2], length.out = nyear)
sd.survey <- seq(from = range.sd.survey[1], to = range.sd.survey[2], length.out = nyear)

# Create site and survey random effects
for(i in 1:nsite){
  for(t in 1:nyear){
    eps1 <- rnorm(n = nsite, mean = 0, sd = sd.site[t])  # Site random eff.
    eps2 <- rnorm(n = nrep, mean = 0, sd = sd.survey[t]) # Survey random eff.
  }
}
# Had this before
#for(s in 1:nspec){
#  for(i in 1:nsite){     # Sites
#    for(t in 1:nyear){   # Years
#      for(j in 1:nrep){ # Occasions interpreted as months
#        p[i,j,t,s] <- plogis(beta0.lp[s, t] + beta1.lp[s] * Xp[i,j,t] +
#        eps1[i] + eps2[j] + beta1[t] * (j - (nrep/2)) +
#        beta2[t] * (j - (nrep/2))^2)
#      }
#    }
#  }
#}

# This should work, too, and be faster
for(s in 1:nspec){
  for(t in 1:nyear){   # Years
    for(j in 1:nrep){ # Occasions interpreted as months
      p[,j,t,s] <- plogis(beta0.lp[s, t] + beta1.lp[s] * Xp[,j,t] +
      eps1 + eps2[j] + beta1[t] * (j - (nrep/2)) +
      beta2[t] * (j - (nrep/2))^2)
    }
  }
}

# (2) Simulate the true system dynamics (state process)
# First year
for(s in 1:nspec){
  z[,1, s] <- rbinom(nsite, 1, psi[,1,s])   # Initial occurrence state
}

# Years 2:nyear
# Had this before
#for(s in 1:nspec){                    # Loop over species
#  for(i in 1:nsite){                  # Loop over sites
#    for(t in 2:nyear){                # Loop over years
#      muZ[i,t,s] <- z[i,t-1,s] * phi[i,t-1,s] + (1-z[i,t-1,s]) * gamma[i,t-1,s]
#      z[i,t,s] <- rbinom(1, 1, muZ[i,t,s])
#    }
#  }
#}

# This should work, too, and be faster
for(s in 1:nspec){                    # Loop over species
  for(t in 2:nyear){                # Loop over years
    muZ[,t,s] <- z[,t-1,s] * phi[,t-1,s] + (1-z[,t-1,s]) * gamma[,t-1,s]
    z[,t,s] <- rbinom(nsite, 1, muZ[,t,s])
  }
}


# (3) Simulate observation process to get the observed data
# Had this before
#for(s in 1:nspec){                    # Loop over species
#  for(i in 1:nsite){                  # Loop over sites
#    for(t in 1:nyear){                # Loop over years
#      for(j in 1:nrep){               # Loop over replicates
#        prob <- z[i,t,s] * p[i,j,t,s] # zero out p for unoccupied sites
#        y[i,j,t,s] <- rbinom(1, 1, prob)
#      }
#    }
#  }
#}

# This should work, too, and be faster
for(s in 1:nspec){                    # Loop over species
  for(t in 1:nyear){                # Loop over years
    for(j in 1:nrep){               # Loop over replicates
      prob <- z[,t,s] * p[,j,t,s] # zero out p for unoccupied sites
      y[,j,t,s] <- rbinom(nsite, 1, prob)
    }
  }
}

# (4) Compute annual population occupancy
# Had this before
#for(s in 1:nspec){                    # Loop over species
#  for(i in 1:nsite){
#    for (t in 2:nyear){
#      psi[i,t,s] <- psi[i,t-1,s] * phi[i,t-1,s] + (1-psi[i,t-1,s]) * gamma[i,t-1,s]
#    }
#  }
#}

# This should work, too, and be faster
for(s in 1:nspec){                    # Loop over species
  for (t in 2:nyear){
    psi[,t,s] <- psi[,t-1,s] * phi[,t-1,s] + (1-psi[,t-1,s]) * gamma[,t-1,s]
  }
}

# Compute some derived stuff
n.occ <- apply(z, 2:3, sum)             # Number of occupied sites
psi.fs <- apply(z, 2:3, mean)           # Finite-sample occupancy proportion
mean.psi <- apply(psi, 2:3, mean)       # Average psi over sites
z.obs <- apply(y, c(1,3,4), max)        # Observed value of z matrix
n.occ.obs <- apply(z.obs, 2:3, sum)     # Observed number of occupied sites
psi.obs <- apply(z.obs, 2:3, mean)      # Observed occupancy (finite sample)

# Total number of species that occur in the sampled sites
tmp1 <- apply(z, 2:3, max)              # True presence per year and species
nyear.pres <- apply(tmp1, 2, sum)       # Number of years when species present
nspec.pres <- sum(nyear.pres > 0)       # Number of species ever present

# Total number of species that were detected anywhere in the sampled sites
tmp2 <- apply(z.obs, 2:3, max)          # Observed presence per year and species
nyear.det <- apply(tmp2, 2, sum)        # Number of years when species detected
nspec.det <- sum(nyear.det > 0)         # Number of species ever detected


# Print out number of occurring and detected species
cat(paste("\n *** Number of species ever occurring:", nspec.pres,
"\n *** Number of species ever detected:", nspec.det,
"\n *** Avg. number of years of occurrence:", round(mean(nyear.pres), 3),
"\n *** Avg. number of years with detection:", round(mean(nyear.det), 3), "\n\n"))



# Compute the average seasonal product of availability and detection
# (ignoring the other terms in the model for detection)
p.season <- array(NA, dim = c(nrep, nyear))
for(t in 1:nyear){   # Years
  p.season[,t] <- plogis(mean(beta0.lp[, t]) + beta1[t] * (month - (nrep/2)) + beta2[t] * (month - (nrep/2))^2)
}

# (5) Plots of stuff
if(show.plot){
  op <- par(mfrow = c(3, 2), mar = c(5,5,4,3), cex.lab = 1.2) ; on.exit(par(op))
  oldAsk <- devAskNewPage(ask = TRUE) ; on.exit(devAskNewPage(oldAsk), add=TRUE)

  # Get predicted covariate relationships and plot them in single graph
  pred.cov <- seq(-2, 2, length.out = 100)
  psi.pred <- phi.pred <- gamma.pred <- p.pred <- array(dim = c(length(pred.cov), nspec))
  for(s in 1:nspec){
    psi.pred[,s] <- plogis(beta0.lpsi[s] + beta1.lpsi[s] * pred.cov)
    phi.pred[,s] <- plogis(mean(beta0.lphi[s,]) + beta1.lphi[s] * pred.cov)
    gamma.pred[,s] <- plogis(mean(beta0.lgamma[s,]) + beta1.lgamma[s] * pred.cov)
    p.pred[,s] <- plogis(mean(beta0.lp[s,]) + beta1.lp[s] * pred.cov)
  }
  matplot(pred.cov, psi.pred, type = 'l', lty = 1, ylim = c(0,1), lwd = 2, main = paste('Occupancy (', nspec, ' species, ', nsite, ' sites)', sep = ''), xlab = 'Covariate', ylab = 'Initial occupancy prob.', las = 1, frame = FALSE)
  matplot(pred.cov, phi.pred, type = 'l', lty = 1, ylim = c(0,1), lwd = 2, main = paste('Persistence (averaged over years,\n', nspec, ' species, ', nsite, ' sites)', sep = ''), xlab = 'Covariate', ylab = 'Persistence prob.', las = 1, frame = FALSE)
  matplot(pred.cov, gamma.pred, type = 'l', lty = 1, ylim = c(0,1), lwd = 2, main = paste('Colonization (averaged over years,\n', nspec, ' species, ', nsite, ' sites)', sep = ''), xlab = 'Covariate', ylab = 'Colonization prob.', las = 1, frame = FALSE)
  matplot(pred.cov, p.pred, type = 'l', lty = 1, ylim = c(0,1), lwd = 2, main = paste('Detection (averaged over years,\n', nspec, ' species, ', nsite, ' sites)', sep = ''), xlab = 'Covariate', ylab = 'Detection prob.', las = 1, frame = FALSE)

  # Plot the average seasonal product of availability and detection
  # (ignoring the other terms in the model for detection)
  matplot(month, p.season, type = 'l', lty = 1, lwd = 2, main = 'Seasonal pattern in p over the years \n(only seasonal terms, same for all species)', xlab = 'Month', ylab = 'Detection probability', ylim = c(0,1))

  # Histo of detection
  hist(p, col = 'grey', breaks = 50, xlim = c(0,1), main = 'Detection probability p\n (all species, sites etc.)')


  # Annual (and species-specific) variation in persistence, colonisation, and detection
  matplot(t(plogis(beta0.lphi)), type = 'l', lty = 1, lwd = 2, ylim = c(0,1), xlab = 'Year', ylab = 'Persistence intercept', main = 'Average persistence per year and species', las = 1, frame = FALSE)
  matplot(t(plogis(beta0.lgamma)), type = 'l', lty = 1, lwd = 2, ylim = c(0,1), xlab = 'Year', ylab = 'Colonization intercept', main = 'Average colonization per year and species', las = 1, frame = FALSE)
  matplot(t(plogis(beta0.lp)), type = 'l', lty = 1, lwd = 2, ylim = c(0,1), xlab = 'Year', ylab = 'Detection intercept', main = 'Average detection per year and species', las = 1, frame = FALSE)

  # Histo of true mean occupancy probability (all species and years)
  hist(mean.psi, col = 'grey', breaks = 50, xlim = c(0,1), main = 'Mean occupancy probability psi1\n (all species and years)')

  # Plot realised and apparent proportion of occupied sites
  matplot(year, mean.psi, type = "l", lty = 1, xlab = "Year", ylab = "Occupancy prob.", xlim = c(0,nyear+1), ylim = c(0,1), lwd = 2, frame.plot = FALSE, las = 1, main = paste('True occupancy (', nspec, ' species, ', nsite, ' sites)', sep = '') )
  matplot(year, psi.obs, type = "l", lty = 1, xlab = "Year", ylab = "Occupancy prob.", xlim = c(0,nyear+1), ylim = c(0,1), lwd = 2, frame.plot = FALSE, las = 1, main = paste('Observed occupancy (', nspec, ' species, ', nsite, ' sites)', sep = ''))
}

# Return data
return(list(nspec = nspec, nsite = nsite, nrep = nrep, nyear = nyear, mean.psi1 = mean.psi1, sig.lpsi1 = sig.lpsi1, mu.beta.lpsi1 = mu.beta.lpsi1, sig.beta.lpsi1 = sig.beta.lpsi1, range.mean.phi = range.mean.phi, sig.lphi = sig.lphi, mu.beta.lphi = mu.beta.lphi, sig.beta.lphi = sig.beta.lphi, range.mean.gamma = range.mean.gamma, sig.lgamma = sig.lgamma, mu.beta.lgamma = mu.beta.lgamma, sig.beta.lgamma = sig.beta.lgamma, range.mean.p = range.mean.p, sig.lp = sig.lp, mu.beta.lp = mu.beta.lp, sig.beta.lp = sig.beta.lp, range.beta1.season = range.beta1.season, range.beta2.season = range.beta2.season, range.sd.site = range.sd.site, range.sd.survey = range.sd.survey, Xpsi1 = Xpsi1, Xphi = Xphi, Xgamma = Xgamma, Xp = Xp, beta0.lpsi = beta0.lpsi, beta1.lpsi = beta1.lpsi, psi = psi, mean.phi = mean.phi, mean.gamma = mean.gamma, eps.lphi = eps.lphi, eps.lgamma = eps.lgamma, beta0.lphi = beta0.lphi, beta0.lgamma = beta0.lgamma, beta1.lphi = beta1.lphi, beta1.lgamma = beta1.lgamma, phi = phi, gamma = gamma, mean.p = mean.p, eps.lp = eps.lp, beta0.lp = beta0.lp, beta1.lp = beta1.lp, beta1 = beta1, beta2 = beta2, sd.site = sd.site, sd.survey = sd.survey, eps1 = eps1, eps2 = eps2, n.occ = n.occ, psi.fs = psi.fs, mean.psi = mean.psi, z.obs = z.obs, n.occ.obs = n.occ.obs, psi.obs = psi.obs, nyear.pres = nyear.pres, nspec.pres = nspec.pres, nyear.det = nyear.det, nspec.det = nspec.det, z = z, p = p, y = y))
} # ------------------ End function definition ---------------------

