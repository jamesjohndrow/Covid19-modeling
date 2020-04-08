#The median time from illness onset ... to death was 18·5 days (15·0–22·0; table 2). 15-22 is IQR
#onset.to.death <- rbinom(100000,122, .157)
scale <- rgamma(100000, 27.75, 1.5)
onset.to.death <- rpois(100000,scale)
hist(onset.to.death)
median(onset.to.death)
quantile(onset.to.death, prob=c(.25, .75))

#The median incubation period was estimated to be 5.1 days (95% CI, 4.5 to 5.8 days),
#and 97.5% of those who develop symptoms will do so within 11.5 days (CI, 8.2 to 15.6 days) of infection.
#mean 5.1, 97.5%ile 11.5
#incubation <- rnbinom(100000, 8, .63)
scale <-rgamma(100000, 5.5, 1.1)
incubation <- rpois(100000, scale)
median(incubation)
quantile(incubation, prob=.975)
quantile(incubation, prob=.99)

samps <- onset.to.death + incubation
max.time <- quantile(samps, prob=.99)-1
theta <- table(c(samps, 0:max.time))[0:max.time + 1]
theta <- theta/sum(theta)

png(paste('./output/theta.png',sep=''),width=500,height=400)
plot(theta, xlab="T", ylab = "probability", cex.axis=1.5, cex.lab=1.5)
dev.off()

saveRDS(theta, file = "./output/theta.Rds")
