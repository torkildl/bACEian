#
# Some simple estimation experiments with twin models in Stan
#
# The script fits the same Stan rnadom effects mixed model parameterization of the classic twin design to
# simulated data on 2000 twin pairs. (A lot, I know but I am working with relatively large samples.)
#
# It fits the model using Hamiltonian MC and Variational Bayes. I am not an expert in any of these, but 
# have gotten som estimates out at least. It compares the estimates to twin estimates obtained from 
# the Mz/Dz correlations.
#
# I have borrowed the Stan code from Espen M. Eilertsen at FHI, but probably broken it, so don't blame him. 
#
#
library(rstan)
library(here)
library(tidyverse)
library(data.table)

### Read in simulated data on 2000 twinpairs
simtwins <- fread(input = here("simtwins.csv")) %>% sample_n(size=300)


# Split by zygosity
mzs <- filter(simtwins, mz==1)
dzs <- filter(simtwins, dz==1)

# Get the variance components the "manual way"
mz_cor <- cor(mzs$tw1,mzs$tw2)
dz_cor <- cor(dzs$tw1,dzs$tw2)
VC_A <- 2*(mz_cor-dz_cor)
VC_C <- (mz_cor-dz_cor)
VC_E <- 1-VC_A-VC_C
manual_results <- data.frame(vc = c("Asd","Csd","Esd"), prop = c(VC_A,VC_C,VC_E))
sum(VC_A,VC_C,VC_E)


# Stan estimation
#
# put data in list
d <- simtwins %>% 
  gather(key=twinvar, value=value, tw1, tw2) %>% 
  arrange(fam) %>% 
  group_by(fam) %>% 
  mutate(twin = row_number()) %>% 
  select(-twinvar) %>% 
  ungroup %>% 
  mutate(famtw = row_number())

dat = list(n_obs = nrow(d), n_fam = max(d$fam), n_famtw = max(d$famtw),
           fam = d$fam, famtw = d$famtw,
           y = d$value, mz = d$mz, dz = d$dz)

# fit the model using regular Stan MC
hmcfit <- stan(file = "model.stan", 
           data = dat,
           pars = c("mu", "a", "c", "e", "A", "C", "E", "Asd", "Csd", "Esd"),
           iter = 2000)
hmcresults <- rstan::extract(hmcfit) %>% map(as.vector) %>% bind_cols

# fit the model using Stans VB algorithm
themodel <- stan_model(file = "model.stan",verbose = T)
vbfit <- rstan::vb(object = themodel,
           data = dat, 
           pars = c("mu", "a", "c", "e", "A", "C", "E", "Asd", "Csd", "Esd"),
           output_samples = 1000)
vbresults <- rstan::extract(vbfit) %>% map(as.vector) %>% bind_cols


df <- list("Stan-HMC" =hmcresults, "Stan-VB"=vbresults) %>% bind_rows(.id="method") %>% 
  select(method, Asd, Csd, Esd) %>% 
  gather(key=vc, value=prop, -method)


# Look at the distributions for A, C and E
hist <- ggplot(df) + 
  geom_histogram( aes(x=prop,fill=method),bins=100) + 
  facet_wrap(~vc) + 
  geom_segment(data=manual_results, mapping=aes(x=prop,xend=prop,y=0,yend=400)) + 
  scale_fill_brewer(type="qual") + 
  coord_cartesian(ylim = c(0,400)) + 
  theme_bw() + 
  theme(legend.position = "bottom") 
hist

ggsave(filename = "histograms.png", plot = hist, device = "png")


# Checking out the relationship between A and C components
AandC <- list("Stan-HMC" =hmcresults, "Stan-VB"=vbresults) %>% bind_rows(.id="method") %>% 
  select(method, Asd, Csd) %>% 
  ggplot(aes(x=Asd, y=Csd, color=method)) + geom_point(size=.5) + geom_smooth(metho="loess",size=0.7) +
  scale_color_brewer(type="qual")
ggsave(filename = "./acscatter.png",plot = AandC, device = "png")


hist

AandC

