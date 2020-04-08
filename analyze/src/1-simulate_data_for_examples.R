#!/usr/bin/env Rscript --vanilla
# set expandtab ts=4 sw=4 ai fileencoding=utf-8
#
# Author: KL
# Maintainer(s): JJ, MG, PB
# License: (c) HRDAG 2020, GPL v2 or newer
#
# -----------------------------------------------------------
# Covid19-modeling/analyze/src/1-simulate_data_for_examples.R
#
# simulate data for testing code


if (!require('pacman')) {install.packages('pacman')}
p_load(ggplot2, tidyr, dplyr)
stopifnot(str_ends(getwd(), "Covid19-modeling/analyze"))

files <- list(funcs = "src/functions.R",
			  theta = "output/theta.Rds",
			  output = "output/example_realization.png")
source(files$funcs)


#set some parameters for simulated example figures
beta <- .2 # infection parameter
gammar <- .1 # recovery parameter
N <- 1000000 # pop size
I <- 1 # initial number of infected

theta <- readRDS(files$theta)
max.time <- length(theta)-1 # max number of days to die given death

p <- .01 # set case fatality rate
T <- 300 # number days in the simulation
T0 <- 1 #set day with first infection

#simulate SIR
sim0 <- simulate.sir(N, beta, gammar, I, T, max.time, theta, p, T0=T0)
Is <- sim0$Is;
Ss <- sim0$Ss;
Rs <- sim0$Rs;
Nus <- sim0$Nus;
Ds <- sim0$Ds[1:T];
save.death.times <- sim0$save.death.times


# FIXME tidyr
# plot the epidemic curves
df <- data.frame(S=Ss, I=Is, R=Rs, T=1:T, Nu=Nus, D=Ds[1:T])
df <- mutate(df, I.scale=(I/max(I))*max(Ds), Nu.scale=(Nu/max(Nu))*max(Ds))
df.long <- df %>% gather(variable,value,-T)
df.long$variable <- factor(df.long$variable, ordered=TRUE, levels=c("S", "I", "R", "Nu", "Nu.scale", "I.scale", "D"))


# FIXME tidyr
df.long <- mutate(df.long, panel = as.character(variable))%>%
  mutate(panel=ifelse(variable %in% c("S", "I", "R"), "SIR", as.character(variable)))
tmp.gg <- ggplot(filter(df.long, variable %in% c("S", "I", "R", "Nu", "D")),
				 aes(x=T, y=value, color=variable)) + geom_line()

# FIXME: pdf? ggsave()? Tarak's graph?
png(files$output, width=800, height=500)
tmp.gg + facet_grid(panel~., scales = "free_y") +
  scale_color_manual(values=c( 'green', 'orange', 'blue','magenta', "red")) +
  theme_bw() +
  ylab("Number") +
  theme(text=element_text(size=24)) +
  ggtitle("Modified SIR Model")
dev.off()

# done.
