stan_dat <- readr::read_rds(here::here("output/rds/stan_dat.rds"))

library(cmdstanr)
# model_NC <- cmdstanr::cmdstan_model("stan/testPoisson_AgeYear_NC.stan")
model_GP <- cmdstanr::cmdstan_model("stan/HierarchicalGP.stan")
fit <- model_GP$sample(data = stan_dat,
                       parallel_chains = 1,
                       chains = 1,
                adapt_delta = 0.95,seed = 56465,
                max_treedepth = 15,
                iter_sampling = 1000,
                iter_warmup = 1000)
f <- fit$summary()
d <- fit$draws()

library(bayesplot)
mcmc_trace(d, pars = glue::glue("vars[{1:3}]"))
mcmc_trace(d, pars = glue::glue("tot_var"))
mcmc_acf(d, regex_pars =  "groups_re")#pars = glue::glue("tot_var"))
mcmc_trace(d, regex_pars =  "groups_re")#pars = glue::glue("tot_var"))
# h <- rstan::stan("stan/HierarchicalGP.stan", data = stan_dat)
# fit <- model_NC$sample(data = stan_dat, parallel_chains = 8,chains = 8,
#                        adapt_delta = 0.95,seed = 56465,max_treedepth = 15,
#                  iter_sampling = 2000,
#                  iter_warmup = 1000)
fit$save_object("/mnt/Storage/LargeRFiles/Dropping-Counts/Heirarchical_AgeYear_NC.rds")
fit <- read_rds("/mnt/Storage/LargeRFiles/Dropping-Counts/testPoisson_AgeYear.rds")
s3 <- fit$summary()
# readr::write_rds(s3, here::here("output/rds/testmodel_varsNC.rds") )
fit$cmdstan_diagnose()
library(dplyr)
library(ggplot2)
library(tidyr)
f %>% filter(grepl("GP_d2s_group\\[", variable)) %>%
  separate(variable, into = c("Var", "D2S", "gr", "ex"),
           sep = "\\[|\\]|\\,", remove = F, convert = T) %>%
  mutate(Year = ifelse(gr<3, 2007, 2008),
         Age = ifelse(gr%%2==0, "Juvenile", "Adults")) %>%
  ggplot(aes(D2S, mean, colour =Age)) +
  geom_pointrange(aes(ymin = q5, ymax = q95, colour = Age))+
  geom_line() + facet_wrap(~Year)


s3 %>% filter(grepl("X1", variable)) %>%
  separate(variable, into = c("Var", "DOY", "ex"),
           sep = "\\[|\\]|\\,", remove = F, convert = T) %>%
  mutate(doy = lubridate::mdy("01-01-2007") +as.numeric(colnames(stan_dat$D_doy)) )%>%
  ggplot(aes(doy, mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2)+
  geom_line()
s3 %>% filter(grepl("alpha", variable)) %>%
  ggplot(aes(variable, mean)) +
  geom_pointrange(aes(ymin = q5, ymax = q95), alpha = 0.2)


model_GP2 <- cmdstanr::cmdstan_model("stan/HierarchicalGP2.stan")
fit2 <- model_GP2$sample(data = stan_dat, parallel_chains = 1,chains = 1,
                       # adapt_delta = 0.95,seed = 56465,
                       # max_treedepth = 15,
                       iter_sampling = 1000,
                       iter_warmup = 1000 )

f2 <- fit2$summary()
d2 <- fit2$draws()#f2
f2 %>% tidyr::separate(variable, sep = "\\[", into = c("a", "b")) %>% pull(a) %>% unique()
library(bayesplot)
mcmc_trace(d2, regex_pars = "d2s_re")
mcmc_trace(d2, regex_pars = "tot_var")
mcmc_trace(d2, regex_pars = "vars")
mcmc_trace(d2, regex_pars = "alpha")
mcmc_trace(d2, regex_pars = "eta")
mcmc_trace(d2, regex_pars = "GP_d2s_group\\[1")

f2 %>% filter(grepl("GP_doy\\[", variable)) %>%
  separate(variable, into = c("Var", "doy", "gr", "ex"),
           sep = "\\[|\\]|\\,", remove = F, convert = T) %>%
  # mutate(Year = ifelse(gr<3, 2007, 2008),
         # Age = ifelse(gr%%2==0, "Juvenile", "Adults")) %>%
  ggplot(aes(doy, mean)) +
  geom_pointrange(aes(ymin = q5, ymax = q95))+
  geom_line()

f2 %>%
  filter(grepl("GP_d2s_group\\[", variable)) %>%
  separate(variable, into = c("Var", "D2S", "gr", "ex"),
           sep = "\\[|\\]|\\,", remove = F, convert = T) %>%
  mutate(Year = ifelse(gr<3, 2007, 2008),
         Age = ifelse(gr%%2==0, "Juvenile", "Adults")) %>%
  ggplot(aes(D2S, mean, colour =Age)) +
  geom_pointrange(aes(ymin = q5, ymax = q95, colour = Age))+
  geom_line() + facet_wrap(~Year)
