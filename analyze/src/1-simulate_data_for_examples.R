#simulate data for testing code
rm(list=ls(all=T))
if (!require('pacman')) {install.packages('pacman')}
pacman::p_load(ggplot2,dplyr,tidyr)
#setwd('/data1/GitHub/Covid19-undercount/analyze/')
#setwd('~/Documents/GitHub/Covid19-undercount/analyze/')
setwd('~/git/Covid19-modeling/analyze/')
source('./src/functions.R')


#set some parameters for simulated example figures
beta <- .2 # infection parameter
gammar <- .1 # recovery parameter
N <- 1000000 # pop size
I <- 1 # initial number of infected

theta <- readRDS('./output/theta.Rds')
max.time <- length(theta)-1 # max number of days to die given death


p <- .01 # set case fatality rate
T <- 300 # number days in the simulation
T0 <- 1 #set day with first infection

#simulate SIR
sim0 <- simulate.sir(N,beta,gammar,I,T,max.time,theta,p, T0=T0)
Is <- sim0$Is; Ss <- sim0$Ss; Rs <- sim0$Rs; Nus <- sim0$Nus; Ds <- sim0$Ds[1:T]; save.death.times <- sim0$save.death.times


#plot the epidemic curves
df <- data.frame(S=Ss, I=Is, R=Rs, T=1:T, Nu = Nus, D = Ds[1:T])
df <- mutate(df, I.scale = I/max(I)*max(Ds), Nu.scale = Nu/max(Nu)*max(Ds))
#df.long <- melt(df, id="T")
df.long <- df %>% gather(variable,value,-T)
df.long$variable <- factor(df.long$variable, ordered=TRUE, levels=c("S", "I", "R", "Nu", "Nu.scale", "I.scale", "D"))


df.long <- mutate(df.long, panel = as.character(variable))%>%
  mutate(panel=ifelse(variable %in% c("S", "I", "R"), "SIR", as.character(variable)))
tmp.gg <- ggplot(filter(df.long, variable %in% c("S", "I", "R", "Nu", "D")), aes(x=T, y=value, color =variable)) + geom_line()

png(paste('./output/example_realization.png',sep=''),width=800,height=500)
tmp.gg + facet_grid(panel~., scales = "free_y") +
  scale_color_manual(values=c( 'green', 'orange', 'blue','magenta', "red")) +
  theme_bw() +
  ylab("Number") +
  theme(text=element_text(size=24)) +
  ggtitle("Modified SIR Model")
dev.off()

