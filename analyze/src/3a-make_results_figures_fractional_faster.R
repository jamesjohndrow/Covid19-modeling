#
# Authors:     KL, JJ
# Maintainers: PB, MG
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
p <- .01
per_decrease <- seq(from=0,to=0.8,by=0.1)
n.cases <- length(per_decrease)

theta <- readRDS(args$theta)
#read in population size data
pops <- readRDS(args$pops)

max.time <- length(theta) - 1
N <- pops[rownames(pops) == args$country]

days.since.mitigation <- 10
n.days.future <- 100


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

dt.all <- c(min(dat$date)+seq(from=-14,to=-1,by=1),dat$date,max(dat$date)+seq(n.days.future))

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
df.bg <- data.frame(beta=mcmc.out$BETA,gammar=mcmc.out$GAMMAR,T0=mcmc.out$T0S,iter=seq(nmc))
df.bg <- df.bg %>% mutate(R0=beta/gammar) %>% gather(variable,value,-iter)

png(paste('./output/trace_plots_fractional_',p,'.png',sep=''),width=600,height=400)
ggplot(df.bg,aes(x=iter,y=value)) + geom_point(size=0.4) + facet_wrap(~variable,scales='free') +
  theme(text=element_text(size=24),axis.text.x=element_text(angle=90)) + xlab("iteration") + theme_bw()
dev.off()
#note: this one should look bad in the beginning and then look good after adaptation begins at 5k iterations. 

#make histograms of posterior distribution of parameters
png(paste('./output/post_hists_fractional_',p,'.png',sep=''),width=600,height=400)
ggplot(df.bg %>% filter(iter>=burn),aes(x=value)) + geom_histogram() + facet_wrap(~variable,scales='free') +
  theme(text=element_text(size=24),axis.text.x=element_text(angle=90)) + theme_bw()
dev.off()

#load post-processed MCMC output
load(file=paste('output/bigDs_',p,'.RData',sep=''))
load(file=paste('output/bigNus_',p,'.RData',sep=''))

n.days.future <- dim(big.Ds)[1] - T

#get quantiles of these
mean.Ds <- apply(big.Ds,c(1,3),mean,na.rm=T)
qlo.Ds <- apply(big.Ds,c(1,3),quantile,probs=0.025)
qhi.Ds <- apply(big.Ds,c(1,3),quantile,probs=0.975)
df.Ds <- data.frame(mean.Ds)
names(df.Ds) <- per_decrease
#df.Ds$day <- seq(T+n.days.future)
df.Ds$day <- dt.all
df.Ds <- df.Ds %>% gather(per_decrease,value,-day)

df.qlo <- data.frame(qlo.Ds)
names(df.qlo) <- per_decrease
df.qlo$day <- dt.all
#df.qlo$day <- seq(T+n.days.future)
df.qlo <- df.qlo %>% gather(per_decrease,q025,-day)

df.qhi <- data.frame(qhi.Ds)
names(df.qhi) <- per_decrease
#df.qhi$day <- seq(T+n.days.future)
df.qhi$day <- dt.all
df.qhi <- df.qhi %>% gather(per_decrease,q975,-day)

df.Ds <- df.Ds %>% inner_join(df.qlo,by=c('per_decrease','day')) %>% inner_join(df.qhi,by=c('per_decrease','day'))
df.Ds$per_decrease <- factor(df.Ds$per_decrease)

#plot future trajectories under different scenarios for R0
png(paste('./output/num_deaths_cases_fractional_fast_',p,'.png',sep=''),width=600,height=400)
ggplot(df.Ds %>% filter(day>date("2020-03-28") & day <= date("2020-03-28") + n.days.future),aes(x=day,y=log10(value),col=per_decrease)) +
  geom_line() +
  theme_bw() +
  geom_ribbon(data=df.Ds %>% filter(day>date("2020-03-28") & day <= date("2020-03-28") + n.days.future),
              aes(ymin=log10(q025),
                  ymax=log10(q975),
                  fill=per_decrease),
              alpha=0.2) +
  theme(text=element_text(size=24)) +
  ylab(expression(paste(log[10], (deaths)))) +
  ylim(c(0,5)) +
  xlab("date")
dev.off()


#get more quantiles
mean.Nus <- apply(big.Nus,c(1,3),mean)
qlo.Nus <- apply(big.Nus,c(1,3),quantile,probs=0.025)
qhi.Nus <- apply(big.Nus,c(1,3),quantile,probs=0.975)
df.Nus <- data.frame(mean.Nus)
names(df.Nus) <- per_decrease
df.Nus$day <- dt.all
df.Nus <- df.Nus %>% gather(per_decrease,value,-day)

df.qlo <- data.frame(qlo.Nus)
names(df.qlo) <- per_decrease
df.qlo$day <- dt.all
#df.qlo$day <- seq(T+n.days.future)
df.qlo <- df.qlo %>% gather(per_decrease,q025,-day)

df.qhi <- data.frame(qhi.Nus)
names(df.qhi) <- per_decrease
df.qhi$day <- dt.all
#df.qhi$day <- seq(T+n.days.future)
df.qhi <- df.qhi %>% gather(per_decrease,q975,-day)

df.Nus <- df.Nus %>% inner_join(df.qlo,by=c('per_decrease','day')) %>% inner_join(df.qhi,by=c('per_decrease','day'))
df.Nus <- df.Nus %>% mutate(ribbon = day > date("2020-03-20") | (per_decrease == 0 &  day > date("2020-03-15")))
df.Nus$per_decrease <- factor(df.Nus$per_decrease)


#plot future new infections under different scenarios for reductions in R0
png(paste('./output/num_newly_infected_cases_fractional_fast_',p,'.png',sep=''),width=600,height=400)
ggplot(df.Nus %>% filter(day>date("2020-03-20")-days.since.mitigation & day <= date("2020-03-29") + n.days.future),aes(x=day,y=log10(value+1),col=per_decrease)) +
  geom_line() +
  theme_bw() +
  theme(text=element_text(size=24)) +
  geom_ribbon(data=df.Nus %>% filter(day>date("2020-03-20")-days.since.mitigation & day <= date("2020-03-29") + n.days.future),aes(ymin=log10(q025),ymax=log10(q975),fill=per_decrease),alpha=0.2) +
  xlab("date") +
  ylab(expression(paste(log[10], nu)))
dev.off()



big.Ds.all <- big.Ds
big.Ds <- big.Ds[1:T,,n.cases]
colnames(big.Ds) <- 1:length(thinned)
big.Ds <- data.frame(big.Ds)
#big.Ds$time <- 1:T
big.Ds$date <- dt.all[1:T]
#big.Ds.long <- big.Ds %>% gather(variable,value,-time)  #melt(big.Ds, id.vars=c("time"))
big.Ds.long <- melt(big.Ds, id.vars=c("date")) #for some unknown reason, gather fucks this up. just leave it as melt
big.Ds.long$alpha <- .05
big.Ds.long$type <- "simulations"
means <- data.frame(date= dt.all[1:T], variable="mean", value=apply(big.Ds[,1:length(thinned)], 1, mean), alpha = 1, type = "mean")
big.Ds.long <- rbind(big.Ds.long, means)
data <- data.frame(date=dt.all[1:T], variable = "data", value = Ds, alpha = 5, type = "actual deaths")
big.Ds.long <- rbind(big.Ds.long, data)

#big.Ds.long$type <- factor(big.Ds.long$type, ordered=TRUE, levels=c("simulations", "mean", "actual deaths"))

png(paste('./output/Predicted_vs_actual_fractional_fast_', p,'.png',sep=''),width=800,height=500)
ggplot(filter(big.Ds.long, date > date("2020-02-15")), aes(x=date, y=value, group=variable, alpha = alpha, color=type)) +
  geom_line() +
  scale_color_manual(values=c("red", "blue", "gray")) +
  guides(alpha=FALSE) + theme_bw() +
  theme(text=element_text(size=24))
#geom_line(aes(x=time,y=variable), data=filter(big.Ds.long, type=="actual deaths"))
dev.off()


R0 <- data.frame(r=mcmc.out$BETA/mcmc.out$GAMMAR)
R0 <- R0 %>% mutate(iter = seq(nmc))

#mean(1/mcmc.out$GAMMAR)

png(paste('./output/R0_fractional_fast_',p,'.png', sep=''),width=400,height=300)
ggplot(R0 %>% filter(iter>burn), aes(x=r)) +
  geom_histogram() +
  theme(text=element_text(size=24)) +
  theme_bw()
dev.off()


df1 <- data.frame(total_infections=apply(mcmc.out$NUS[(burn+1):nmc,1:(T-10)],1,sum))

png(paste('output/infected_10_days_ago_fractional_fast_',p,'.png',sep=''),width=600,height=400)
ggplot(df1, aes(x=total_infections)) + geom_histogram() + theme_bw() +
  #  ggtitle("posterior distribution of Infected + Recovered \n as of 10 days ago") +
  xlab("Infected + Recovered") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=24))
dev.off()


df1 <- data.frame(new_infections=mcmc.out$NUS[(burn+1):nmc,T-10])

png(paste('output/newly_infected_10_days_ago_fractional_fast_',p,'.png',sep=''),width=600,height=400)
ggplot(df1, aes(x=new_infections)) + geom_histogram() + theme_bw() +
  #  ggtitle("posterior distribution of Infected + Recovered \n as of 10 days ago") +
  xlab("New Infected") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text=element_text(size=24))
dev.off()

png(paste('output/time_to_recovery_fractional_fast_',p,'.png',sep=''),width=600,height=400)
hist(1/mcmc.out$GAMMAR[(burn+1):nmc],main='posterior distribution of \n average length of the infectious period',xlab=expression(gamma[r]^(-1)))
dev.off()



quantile(R0$r[R0$iter>burn],c(0.025,0.975))

df.Ds %>% filter(day %in% c(date("2020-03-29"),date("2020-03-30"),date("2020-03-31"))) %>% group_by(day) %>% summarise(mn.deaths = mean(value))






