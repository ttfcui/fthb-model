# likelihood and density regeneration functions
# TODO: this likelihood is defined over densities rather than
# probabilities. Fix discrepancy with integration?
tau <- function(mu, sigma) {
  -sum(log(log(lambda)) - log(sigma*sqrt(2*pi*allT^3)) -
         (log(lambda) - mu*allT)^2/(2*allT*sigma^2))
}
tau_mix <- function(mu, sigma1, sigma2, prop) {
    dist1 <- log(lambda)/(sigma1*sqrt(2*pi*allT^3))*
             exp(-(log(lambda) - mu*allT)^2/
                 (2*allT*sigma1^2))
    allT2 <- allT-tLate
    allT2[which(allT2 < 0)] <- 0
    dist2 <- log(lambda)/(sigma2*sqrt(2*pi*allT2^3))*
             exp(-(log(lambda) - mu*allT2)^2/
                 (2*allT2*sigma2^2))
    dist2[which(
          is.na(dist2$randTime))] <- 0
    -sum(log(prop*dist1 + (1-prop)*dist2))
}
tau_dens <- function(x, mu, sigma) {
  data.frame(age_primary=19+x,
             prop=log(lambda)/(sigma*sqrt(2*pi*x^3))*exp(
                 -(log(lambda) - mu*x)^2/(2*x*sigma^2)))
}
tau_mix_dens <- function(x, mu, sigma1, sigma2, prop) {
    x1 <- x[which(x <= tLate)]
    x2 <- x[which(x > tLate)]
    x2a <- x2 - tLate
    dist1 <- prop*log(lambda)/(sigma1*sqrt(2*pi*x1^3))*exp(
             -(log(lambda) - mu*x1)^2/(2*x1*sigma1^2))
    dist2 <- prop*log(lambda)/(sigma1*sqrt(2*pi*x2^3))*exp(
             -(log(lambda) - mu*x2)^2/(2*x2*sigma1^2)) +
             (1-prop)*log(lambda)/(sigma2*sqrt(2*pi*x2a^3))*exp(
               -(log(lambda) - mu*x2a)^2/(2*x2a*sigma2^2))
    data.frame(age_primary=19+x, prop=c(dist1, dist2))
}

library(data.table)
library(ggplot2)

# Load in the original IRS data
ecdfs <- fread('~/dropbox/fthb/data/irs/AGE_COUNTS.csv')
setkey(ecdfs, TAX_YR, age_primary)
totals <- ecdfs[,.(denom=sum(count)),by=TAX_YR]
setkey(totals, TAX_YR)
ecdfs <- merge(ecdfs, totals)
ecdfs[, prop := count/denom]

for (i in 2002:2013) {
  # Prep data for a single year
test <- cumsum(ecdfs[TAX_YR == i, prop])
randNo <- runif(100000)
randTime <- sapply(randNo, function(x) which(1/(test - x) == max(1/(test - x))))
allT <- data.table(randTime) - 0.1 # As a concession to continuous time
ecdf_rep <- allT[,.N/100000,by=randTime]
ecdf_rep[, randTime := randTime + 0.1]
setkey(ecdf_rep, randTime)

# Lambda is p_h*theta
lambda <- if (i != 2009) 15*0.2*exp(1) else 0.9*15*0.2*exp(1)
#lik_params <- list(mu=1/20, sigma=1/10)
lik_params <- list(mu=1/20, sigma1=1/10, sigma2=1/8, prop=0.8)
for (j in 0) {
tLate <- j
# Run MLE estimation
out <- stats4::mle(tau_mix, start = lik_params, method="L-BFGS-B",
                   lower=c(1e-2, 1e-1, 1e-1, 0.1),
                   upper=c(Inf, Inf, Inf, 1))
coefs <- attr(out, "coef")
cat(i, ": ", coefs, "\n")
cat(AIC(out), "\n")

#sim_call <- tau_dens(1:50, coefs[1], coefs[2])
sim_call <- tau_mix_dens(1:50, coefs[1], coefs[2], coefs[3],
                         coefs[4])

# Plot empirical PMF with estimated continuous function
print(ggplot(ecdfs[(TAX_YR == i) & (age_primary < 60), .(age_primary, prop)],
      aes(x=age_primary, y=prop)) + geom_line() +
      geom_line(data=sim_call, colour="red"))
}
}
# NOTES: Fits post-08 seem plausible with just one distribution.
# Fixing lambda to a certain value (from the data?) seems important
# because mu, lambda not easy to identify separately.
# Changing the value of lambda to discount it in only 2009 makes
# endogenous wealth parameters in line with all other years.
# Fits pre-08 are quite bad. Is this where you estimate a mixture of two
# stopping times instead?
