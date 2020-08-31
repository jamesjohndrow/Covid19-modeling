#!/usr/bin/env Rscript --vanilla
# set expandtab ts=4 sw=4 ai fileencoding=utf-8
#
# Author: KL, JJ
# Maintainer(s): PB, MG
# License: (c) HRDAG 2020, GPL v2 or newer
#
# Covid19-undercount/analyze/src/1-simulate_data.R
# -----------------------------------------------------------
#
# simulate data for testing code

if (!require('pacman')) {install.packages('pacman')}
pacman::p_load(ggplot2, dplyr, tidyr, forcats, stringr)
stopifnot(str_ends(getwd(), "Covid19-modeling/analyze"))

set.seed(19481210)

getargs <- function() {
	# sets file paths in one place, could be used in command line call
	# works in RStudio
    parser <- argparse::ArgumentParser()
    parser$add_argument("--countrypops", default="output/country_pops.RData")
    parser$add_argument("--thetagraph", default="output/theta.png")
 	parser$add_argument("--theta", default="output/theta.rds")
	parser$add_argument("--example_w_nu_pdf", default="output/example_realization_w_nu.pdf")
	parser$add_argument("--example_w_nu_png", default="output/example_realization_w_nu.png")
	parser$add_argument("--example1", default="output/example_realization1.png")
	parser$add_argument("--funs", default="src/functions.R")
    parser$parse_args()
}
args <- getargs()
source(args$funs)

country <- 'us'
R0.case <- seq(from=.75, to=2.25, by=0.25)
n.cases <- length(R0.case)

# country pops
pops <- matrix(0, 6, 1)
pops[1] <- 3.3e8
pops[2] <- 6.05e7 #it
pops[3] <- 6.7e7 #fr
pops[4] <- 8.29e7 #de
pops[5] <- 1.72e7 #nl
pops[6] <- 6.64e7 #uk
rownames(pops) <- c('us', 'it', 'fr', 'de', 'nl', 'uk')
saveRDS(pops, file=args$countrypops)

beta <- .2 # infection parameter
gammar <- .05 # recovery parameter
N <- 1000000 # pop size
I <- 10 # initial number of infected
max.time <- 30 # max number of days to die given death

# The median time from illness onset ... to death was 18·5 days (15·0–22·0;
# table 2). 15-22 is IQR
#---onset.to.death <- rbinom(100000,122, .157)
scale <- rgamma(100000, 27.75, 1.5)
onset.to.death <- rpois(100000, scale)
hist(onset.to.death)
median(onset.to.death)
quantile(onset.to.death, prob=c(.25, .75))

# The median incubation period was estimated to be 5.1 days (95% CI, 4.5 to 5.8
# days), and 97.5% of those who develop symptoms will do so within 11.5 days
# (CI, 8.2 to 15.6 days) of infection.  mean 5.1, 97.5%ile 11.5
#---incubation <- rnbinom(100000, 8, .63)
scale <-rgamma(100000, 5.5, 1.1)
incubation <- rpois(100000, scale)
median(incubation)
quantile(incubation, prob=.975)
quantile(incubation, prob=.99)

samps <- onset.to.death + incubation
max.time <- quantile(samps, prob=.99)-1
theta <- table(c(samps, 0:max.time))[0:max.time + 1]
theta <- theta/sum(theta)

png(args$thetagraph, width=500, height=400)
plot(theta, xlab="T", ylab="probability", cex.axis=1.5, cex.lab=1.5)
dev.off()
saveRDS(theta, file=args$theta)

# case fatality rate
p <- .01
p_str <- sub('\\.','_',as.character(p))
T <- 300 # number days in the simulation

gammar <- .1
sim0 <- simulate.sir.synthetic.data(N, beta, gammar, I, T, max.time, theta, p)
Is <- sim0$Is
Ss <- sim0$Ss
Rs <- sim0$Rs
Nus <- sim0$Nus
Ds <- sim0$Ds[1:T]
save.death.times <- sim0$save.death.times


plot(Ss, type = 'l', col='black', ylim=c(0, N))
lines(Is, col = 'blue')
lines(Rs, col = 'green')

# FIXME: remove this deprecated code
# df <- data.frame(S=Ss, I=Is, R=Rs, T=1:T, Nu = Nus, D = Ds[1:T])
# df <- mutate(df, I.scale = I/max(I)*max(Ds), Nu.scale = Nu/max(Nu)*max(Ds))
# #df.long <- melt(df, id="T")
# df.long <- df %>% gather(variable,value,-T)

df.long <- data.frame(S=Ss, I=Is, R=Rs, T=1:T, Nu = Nus, D = Ds[1:T]) %>%
	mutate(I.scale=(I / max(I)) * max(Ds), Nu.scale=(Nu / max(Nu)) * max(Ds)) %>%
	gather(variable, value, -T)

sir_levels <- c("S", "I", "R", "Nu", "Nu.scale", "I.scale", "D")
df.long$variable <- factor(df.long$variable,
						   ordered=TRUE,
						   levels=sir_levels)


df.long.sir <- filter(df.long, variable %in% c("S", "I", "R"))
example.plot <- ggplot(df.long.sir, aes(x=T, y=value, color=variable)) +
  geom_line() +
  theme_bw() +
  ylab("Number") +
  ggtitle("One realization of the SIR model") +
  scale_color_manual(values=c("blue", 'orange', 'green'))
ggsave(args$example1, example.plot, width=400, height=250, units="mm")
# png(paste('./output/example_realization1.png',sep=''),width=800,height=500)
rm(df.long.sir)

df.long2 <- df.long %>%
  filter(variable %in% c("S", "I", "R", "Nu", "D")) %>%
  mutate(variable = fct_recode(variable,
                               nu = "Nu",
                               Deaths = "D",
                               Susceptible = "S",
                               Infected = "I",
                               Recovered = "R"),
         panel = fct_collapse(variable,
                              SIR = c("Susceptible",
                                      "Infected",
                                      "Recovered"),
                              group_other = FALSE),
         panel = fct_rev(panel))

p <- ggplot(df.long2, aes(x = T, y = value, color = variable)) +
  geom_line() +
  scale_color_manual(values=c(Susceptible = 'green',
                              Infected = 'orange',
                              Recovered = 'blue',
                              nu = 'magenta',
                              Deaths = "red"),
                     labels = c("nu" = expression(nu))) +
  facet_grid(panel ~ ., scales = "free_y", labeller = label_parsed) +
  theme_bw() + theme(legend.title = element_blank()) +
  ylab("Number") + xlab("Days during epidemic") +
  theme(text=element_text(size=24))

ggsave(args$example_w_nu_pdf, p, width=400, height=250, units="mm")
ggsave(args$example_w_nu_png, p, width=400, height=250, units="mm")

plot(Ds, col = 'red', type = 'l', ylim=range(c(Ds, Nus/100)))
lines(Nus/100, col = 'black')
legend("topright", fill = c("red", "black"),
	   legend=c("Deaths", "New Infections/100"))
hist(save.death.times)


