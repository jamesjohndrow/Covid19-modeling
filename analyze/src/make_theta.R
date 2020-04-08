#
# Authors:     KL
# Maintainers: JJ, PB, MG
# Copyright:  
# =========================================
# Covid19-modeling/analyze/src/make_theta.R

if (!require('pacman')) {install.packages('pacman')}
pacman::p_load(here)

args <- list(theta_plot=here::here('analyze/output/theta.png'),
             theta=here::here('analyze/output/theta.Rds'))

#The median time from illness onset ... to death was 18·5 days (15·0–22·0; table 2). 15-22 is IQR
#onset.to.death <- rbinom(100000, 122, 0.157)
scale <- rgamma(100000, 27.75, 1.5)
onset.to.death <- rpois(100000, scale)
hist(onset.to.death)
median(onset.to.death)
quantile(onset.to.death, prob=c(0.25, 0.75))

#The median incubation period was estimated to be 5.1 days (95% CI, 4.5 to 5.8 days),
#and 97.5% of those who develop symptoms will do so within 11.5 days (CI, 8.2 to 15.6 days) of infection.
#mean 5.1, 97.5%ile 11.5
#incubation <- rnbinom(100000, 8, 0.63)
scale <-rgamma(100000, 5.5, 1.1)
incubation <- rpois(100000, scale)
median(incubation)
quantile(incubation, prob=0.975)
quantile(incubation, prob=0.99)

samps <- onset.to.death + incubation
max.time <- quantile(samps, prob=0.99) - 1
theta <- table(c(samps, 0:max.time))[0:max.time + 1]
theta <- theta / sum(theta)

png(args$theta_plot, width=500, height=400)
plot(theta, xlab="T", ylab="probability", cex.axis=1.5, cex.lab=1.5)
dev.off()

saveRDS(theta, file=args$theta)

# done.
