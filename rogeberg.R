#
# Some simple estimation experiments with twin models in Stan
#
# The script fits the Ole Røgeberg's random effects mixed model parameterization of the classic twin design in 
# Stan tp simulated data on 2000 twin pairs. (A lot, I know but I am working with relatively large samples.)
#
# It fits the model using Hamiltonian MC.
#
options(mc.cores = parallel::detectCores())

#
library(rstan)
rstan_options(auto_write = TRUE)
library(here)
library(tidyverse)
library(data.table)
library(tidybayes)
library(mets)
library(lemon)

#legendmovefunction
shift_legend2 <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"
  
  # now we just need a simple call to reposition the legend
  reposition_legend(p, 'center', panel=names)
}
#shift_legend2(p)



### Read in simulated data on 2000 twinpairs
simtwins <- fread(input = here("simtwins.csv"))
head(simtwins)

# Get the variance components the "manual way"
mzs <- filter(simtwins, mz==1)
dzs <- filter(simtwins, dz==1)
mz_cor <- cor(mzs$tw1,mzs$tw2)
dz_cor <- cor(dzs$tw1,dzs$tw2)
VC_A <- 2*(mz_cor-dz_cor)
VC_C <- (mz_cor-dz_cor)
VC_E <- 1-VC_A-VC_C
manual_results <- data.frame(vc = c("Asd","Csd","Esd"), prop = c(VC_A,VC_C,VC_E))
sum(VC_A,VC_C,VC_E)


### Stan estimation
###
### Ole Røgeberg's Dirichlet priors
###
data_stan = list(n_fam = max(simtwins$fam),
                 n_famtw_mz = sum(simtwins$mz == 1),
                 n_famtw_dz = sum(simtwins$mz == 0),
                 y_mz = c(simtwins[mz == 1,]$tw1,
                          simtwins[mz == 1,]$tw2),
                 y_dz = c(simtwins[mz == 0,]$tw1,
                          simtwins[mz == 0,]$tw2),
                 fam_mz = c(simtwins[mz == 1,]$fam,
                            simtwins[mz == 1,]$fam),
                 fam_dz = c(simtwins[mz == 0,]$fam,
                            simtwins[mz == 0,]$fam),
                 outcome_sd = sd(c(simtwins$tw1, simtwins$tw2)),
                 outcome_mean = mean(c(simtwins$tw1, simtwins$tw2)))

dirichlet_fit <- stan(file = here::here("./mixed-ace-dirichlet.stan"), 
           data = data_stan, iter=1000, control = list(adapt_delta = 0.9), pars = c("mu", "a", "c", "e_sigma", "A", "C", "E", "Asd", "Csd", "Esd"))


## Get posterior
posterior <- tidy_draws(dirichlet_fit)


## Histogram of ACE
histoace <- posterior %>% 
  select(Asd, Csd, Esd, .chain) %>% 
  gather(key=parm, value=est, -.chain) %>% 
  ggplot() + geom_density(aes(fill=as.factor(.chain), x=est), alpha=.3) + facet_wrap(~parm, ncol=2)
shift_legend2(histoace)



## Ask some Qs:
p_C_larger_A <- nrow(filter(posterior,Csd>Asd)) / nrow(posterior)
p_A_larger_E <- nrow(filter(posterior,Asd>Esd)) / nrow(posterior)

## Scatter of Asd and Csd
ggplot(posterior) + geom_point(aes(x=Asd, y=Csd), alpha=0.3, color="seagreen") +theme_minimal()





###
### Test with mets pkg's BMI data for twins.
###
data("twinbmi")
bmidf <- as.data.frame(twinbmi) %>% 
  group_by(tvparnr) %>% 
  filter(n()>1) %>% 
  ungroup %>% 
  mutate(var = str_c("tw",num)) %>% 
  select(bmi,age,var,zyg,tvparnr) %>% 
  spread(key=var,value=bmi) %>% 
  mutate(mz = if_else(zyg=="MZ",1,0)) %>% 
  mutate(dz = 1-mz) %>% 
  mutate(age = round(age,0)) %>% 
  mutate(fam = row_number())


data_stan = list(n_fam = max(bmidf$fam),
                 n_famtw_mz = sum(bmidf$mz == 1),
                 n_famtw_dz = sum(bmidf$mz == 0),
                 y_mz = c(pull(filter(bmidf,mz == 1),tw1),
                          pull(filter(bmidf,mz == 1),tw2)),
                 y_dz = c(pull(filter(bmidf,mz == 0),tw1),
                          pull(filter(bmidf,mz == 0),tw2)),
                 fam_mz = rep(pull(filter(bmidf, mz == 1),fam),2),
                 fam_dz = rep(pull(filter(bmidf, mz == 0),fam),2),
                 outcome_mean = mean(c(bmidf$tw1, bmidf$tw2)),
                 outcome_sd = sd(c(bmidf$tw1, bmidf$tw2)))

dirichlet_fit <- stan(file = here::here("./mixed-ace-dirichlet.stan"), 
                      data = data_stan, iter=2000, control = list(adapt_delta = 0.9), 
                      pars = c("mu", "a", "c", "e_sigma", "A", "C", "E", "Asd", "Csd", "Esd"))


posterior2 <- tidy_draws(dirichlet_fit)
posterior2

## Histogram of ACE
posterior2 %>% 
  select(Asd, Csd, Esd, .chain) %>% 
  gather(key=parm, value=est, -.chain) %>% 
  ggplot() + geom_density(aes(fill=as.factor(.chain), x=est), alpha=.3) + facet_wrap(~parm, ncol=2)


## Ask some Qs:
p_C_larger_A <- nrow(filter(posterior2,Csd>Asd)) / nrow(posterior)
p_A_larger_E <- nrow(filter(posterior2,Asd>Esd)) / nrow(posterior)

## Scatter of Asd and Csd
ggplot(posterior2) + geom_point(aes(x=Asd, y=Csd), alpha=0.3, color="seagreen") +theme_minimal()
