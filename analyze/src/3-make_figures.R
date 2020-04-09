#
# Authors:     KL
# Maintainers: JJ, PB, MG
# Copyright:
# =========================================
# Covid19-modeling/analyze/src/2-fit_model_to_real_data.R

#simulate data for testing code
if (!require('pacman')) {install.packages('pacman')}
pacman::p_load(ggplot2, dplyr, tidyr, readr, lubridate,
               glue, future, future.apply, reshape2, argparse, here)

parser <- ArgumentParser()
parser$add_argument('--functions', default=here::here('analyze/src/functions.R'))
parser$add_argument('--theta', default=here::here('analyze/output/theta.Rds'))
parser$add_argument('--pops', default=here::here('analyze/input/country_pops.RData'))
parser$add_argument('--worldometers', default=here::here('analyze/input/worldcounts_deaths_us.csv'))
parser$add_argument('--p', default=0.01)
parser$add_argument('--country', default='us')

args <- parser$parse_args()

args[['mcmc']] <- here::here(glue('analyze/output/mcmc_run_{args$p}.Rdata'))
args[['trace_plots']] <- here::here(glue('analyze/output/trace_plots_fractional_{args$country}_{args$p}.png'))
args[['post_hists']] <- here::here(glue('analyze/output/post_hists_fractional_{args$country}_{args$p}.png'))
args[['bigDs']] <- here::here(glue('analyze/output/bigDs_{args$p}.Rdata'))
args[['bigNus']] <- here::here(glue('analyze/output/bigNus_{args$p}.Rdata'))
args[['num_deaths_cases']] <- here::here(glue('analyze/output/num_deaths_cases_fractional_{args$country}_{args$p}.png'))
args[['num_newly_infected']] <- here::here(glue('analyze/output/num_newly_infected_cases_fractional_{args$country}_{args$p}.png'))
args[['predicted_vs_actual']] <- here::here(glue('analyze/output/Predicted_vs_actual_fractional_{args$country}_{args$p}.png'))
args[['R0']] <- here::here(glue('analyze/output/R0_fractional_{args$country}_{args$p}.png'))
args[['infected_10']] <- here::here(glue('analyze/output/infected_10_days_ago_fractional_{args$country}_{args$p}.png'))
args[['newly_infected_10']] <- here::here(glue('analyze/output/newly_infected_10_days_ago_fractional_{args$country}_{args$p}.png'))
args[['time_to_recovery']] <- here::here(glue('analyze/output/time_to_recovery_fractional_{args$country}_{args$p}.png'))


source(args$functions)

##----- set up
R0.case <- seq(from=0.75, to=2.25, by=0.25)
n.cases <- length(R0.case)

theta <- readRDS(args$theta)
#read in population size data
pops <- readRDS(args$pops)

max.time <- length(theta) - 1
N <- pops[rownames(pops) == args$country]

days.since.mitigation <- 10

burn <- 10000
n.thin <- 1000
disp.int <- 1000

# worldometers data entered manually from website
dat <- read_csv(args$worldometers)
dat <- dat %>% mutate(date=mdy(date))
Ds <- dat$deaths_us
Ds <- c(rep(0, 14), Ds)

dt.all <- c(min(dat$date) + seq(from=-14, to=-1, by=1),
            dat$date,
            max(dat$date) + seq(21))

# In a "normal" SIR model, the R0 is beta/gamma
# The average days to recovery is 1/gamma
# Let's choose
#gammar <- 1/6
#zeta <- c(1/25,1/2) # limits of uniform prior on gammar
#beta <- gammar*2.5
#xi <- c(1.2,5.0) # prior on beta is uniform on this interval

T <- length(Ds)

#load the results of the mcmc
load(args$mcmc)

nmc <- nrow(mcmc.out$BETA)
thinned <- round(seq(burn, nmc, length.out=n.thin))

#make trace plots
df.bg <- data.frame(beta=mcmc.out$BETA,
                    gammar=mcmc.out$GAMMAR,
                    T0=mcmc.out$T0S,
                    iter=seq(nmc))
df.bg <- df.bg %>%
  mutate(R0=beta / gammar) %>%
  gather(variable, value, -iter)

png(args$trace_plots, width=600, height=400)
ggplot(df.bg, aes(x=iter, y=value)) +
  geom_point(size=0.4) +
  facet_wrap(~variable, scales='free') +
  theme(text=element_text(size=24), axis.text.x=element_text(angle=90)) +
  xlab("iteration") +
  theme_bw()
dev.off()
print("traceplot done")

#make histograms of posterior distribution of parameters
png(args$post_hists, width=600, height=400)
ggplot(df.bg %>% filter(iter>=burn), aes(x=value)) +
  geom_histogram() +
  facet_wrap(~variable, scales='free') +
  theme(text=element_text(size=24), axis.text.x=element_text(angle=90)) +
  theme_bw()
dev.off()
print("posterior dists done")

#load post-processed MCMC output
load(file=args$bigDs)
load(file=args$bigNus)

n.days.future <- dim(big.Ds)[1] - T

#get quantiles of these
mean.Ds <- apply(big.Ds, c(1, 3), mean, na.rm=TRUE)
qlo.Ds <- apply(big.Ds, c(1,3), quantile, probs=0.025)
qhi.Ds <- apply(big.Ds, c(1,3), quantile, probs=0.975)
df.Ds <- data.frame(mean.Ds)
names(df.Ds) <- c(R0.case,
                  round(mean(mcmc.out$BETA[burn:nmc] / mcmc.out$GAMMAR[burn:nmc]), 2))
#df.Ds$day <- seq(T+n.days.future)
df.Ds$day <- dt.all
df.Ds <- df.Ds %>% gather(R0, value, -day)

df.qlo <- data.frame(qlo.Ds)
names(df.qlo) <- c(R0.case,
                   round(mean(mcmc.out$BETA[burn:nmc] / mcmc.out$GAMMAR[burn:nmc]), 2))
df.qlo$day <- dt.all
#df.qlo$day <- seq(T+n.days.future)
df.qlo <- df.qlo %>% gather(R0, q025, -day)

df.qhi <- data.frame(qhi.Ds)
names(df.qhi) <- c(R0.case,
                   round(mean(mcmc.out$BETA[burn:nmc]/mcmc.out$GAMMAR[burn:nmc]), 2))
#df.qhi$day <- seq(T+n.days.future)
df.qhi$day <- dt.all
df.qhi <- df.qhi %>% gather(R0, q975, -day)

df.Ds <- df.Ds %>%
  inner_join(df.qlo, by=c('R0', 'day')) %>%
  inner_join(df.qhi, by=c('R0','day'))
df.Ds$R0 <- factor(df.Ds$R0)

#plot future trajectories under different scenarios for R0
png(args$num_deaths_cases, width=600, height=400)

ggplot(df.Ds %>%
         filter(day > date("2020-03-28") & day <= date("2020-03-28") + n.days.future),
       aes(x=day, y=log10(value), col=R0)) +
  geom_line() +
  theme_bw() +
  geom_ribbon(data=df.Ds %>%
                filter(day > date("2020-03-28") & day <= date("2020-03-28") + n.days.future),
              aes(ymin=log10(q025),
                  ymax=log10(q975),
                  fill=R0),
              alpha=0.2) +
  theme(text=element_text(size=24)) +
  ylab(expression(paste(log[10], (deaths)))) +
  ylim(c(0, 5)) +
  xlab("date")
dev.off()
print("future trajectories done")

#get more quantiles
mean.Nus <- apply(big.Nus, c(1, 3), mean)
qlo.Nus <- apply(big.Nus, c(1, 3), quantile, probs=0.025)
qhi.Nus <- apply(big.Nus, c(1, 3), quantile, probs=0.975)
df.Nus <- data.frame(mean.Nus)
names(df.Nus) <- c(R0.case,
                   round(mean(mcmc.out$BETA[burn:nmc] / mcmc.out$GAMMAR[burn:nmc]), 2))
df.Nus$day <- dt.all
df.Nus <- df.Nus %>% gather(R0, value, -day)

df.qlo <- data.frame(qlo.Nus)
names(df.qlo) <- c(R0.case,
                   round(mean(mcmc.out$BETA[burn:nmc] / mcmc.out$GAMMAR[burn:nmc]), 2))
df.qlo$day <- dt.all
#df.qlo$day <- seq(T+n.days.future)
df.qlo <- df.qlo %>% gather(R0, q025, -day)

df.qhi <- data.frame(qhi.Nus)
names(df.qhi) <- c(R0.case,
                   round(mean(mcmc.out$BETA[burn:nmc] / mcmc.out$GAMMAR[burn:nmc]), 2))
df.qhi$day <- dt.all
#df.qhi$day <- seq(T+n.days.future)
df.qhi <- df.qhi %>% gather(R0, q975, -day)

df.Nus <- df.Nus %>%
  inner_join(df.qlo, by=c('R0', 'day')) %>%
  inner_join(df.qhi, by=c('R0', 'day'))
df.Nus <- df.Nus %>%
  mutate(ribbon = day > date("2020-03-20") | (R0 > 2.5 &  day > date("2020-03-15")))
df.Nus$R0 <- factor(df.Nus$R0)

#plot future new infections under different scenarios for reductions in R0
png(args$num_newly_infected, width=600, height=400)
ggplot(df.Nus %>% filter(
						 (day > date("2020-03-20") - days.since.mitigation &
						  day <= date("2020-03-29") + n.days.future)),
       aes(x=day, y=log10(value+1), col=R0)) +
  geom_line() +
  theme_bw() +
  theme(text=element_text(size=24)) +
  geom_ribbon(df.Nus %>%
                filter((day > date("2020-03-20") - days.since.mitigation &
                          day <= date("2020-03-29") + n.days.future)),
              mapping = aes(ymin=log10(q025), ymax=log10(q975), fill=R0), alpha=0.2) +
  xlab("date") +
  ylab(expression(paste(log[10], nu)))
dev.off()
print("future new infections done.")

big.Ds.all <- big.Ds
big.Ds <- big.Ds[1:T,,n.cases+1]
colnames(big.Ds) <- 1:length(thinned)
big.Ds <- data.frame(big.Ds)
#big.Ds$time <- 1:T
big.Ds$date <- dt.all[1:T]
#big.Ds.long <- big.Ds %>% gather(variable,value,-time)  #melt(big.Ds, id.vars=c("time"))
big.Ds.long <- melt(big.Ds, id.vars=c("date"))
big.Ds.long$alpha <- 0.05
big.Ds.long$type <- "simulations"
means <- data.frame(date=dt.all[1:T],
                    variable="mean",
                    value=apply(big.Ds[,1:length(thinned)], 1, mean),
                    alpha=1,
                    type = "mean")
big.Ds.long <- rbind(big.Ds.long, means)
data <- data.frame(date=dt.all[1:T],
                   variable="data",
                   value=Ds,
                   alpha=5,
                   type="actual deaths")
big.Ds.long <- rbind(big.Ds.long, data)

#big.Ds.long$type <- factor(big.Ds.long$type, ordered=TRUE, levels=c("simulations", "mean", "actual deaths"))

png(args$predicted_vs_actual, width=800, height=500)
ggplot(filter(big.Ds.long, date > date("2020-02-15")),
       aes(x=date, y=value, group=variable, alpha=alpha, color=type)) +
  geom_line() +
  scale_color_manual(values=c("red", "blue", "gray")) +
  guides(alpha=FALSE) +
  theme_bw() +
  theme(text=element_text(size=24))
#geom_line(aes(x=time,y=variable), data=filter(big.Ds.long, type=="actual deaths"))
dev.off()
print("pred vs actual done")

R0 <- data.frame(r=mcmc.out$BETA / mcmc.out$GAMMAR)
R0 <- R0 %>% mutate(iter=seq(nmc))

#mean(1/mcmc.out$GAMMAR)

png(args$R0, width=400, height=300)
ggplot(R0 %>% filter(iter > burn), aes(x=r)) +
  geom_histogram() +
  theme(text=element_text(size=24)) +
  theme_bw()
dev.off()
print("R0 histogram done")


df1 <- data.frame(total_infections=apply(mcmc.out$NUS[(burn + 1):nmc, 1:(T - 10)],
                                         1,
                                         sum))
png(args$infected_10, width=600, height=400)
ggplot(df1, aes(x=total_infections)) +
  geom_histogram() +
  theme_bw() +
  xlab("Infected + Recovered") +
  theme(plot.title=element_text(hjust=0.5)) +
  theme(text=element_text(size=24))
dev.off()
print("inf + recovered done")

df1 <- data.frame(new_infections=mcmc.out$NUS[(burn+1):nmc,T-10])

png(args$newly_infected_10, width=600, height=400)
ggplot(df1, aes(x=new_infections)) +
  geom_histogram() +
  theme_bw() +
  xlab("New Infected") +
  theme(plot.title=element_text(hjust=0.5)) +
  theme(text=element_text(size=24))
dev.off()
print("newly infected done")

png(args$time_to_recovery, width=600, height=400)
hist(1 / mcmc.out$GAMMAR[(burn+1):nmc],
     main='posterior distribution of \n average length of the infectious period',
     xlab=expression(gamma[r]^(-1)))
dev.off()
print("time to recovery")

quantile(R0$r[R0$iter > burn], c(0.025, 0.975))

df.Ds %>%
  filter(day %in% c(date("2020-03-29"), date("2020-03-30"), date("2020-03-31"))) %>%
  group_by(day) %>%
  summarise(mn.deaths=mean(value))

# done.
