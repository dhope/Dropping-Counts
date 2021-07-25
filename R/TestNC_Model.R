stan_dat <- readr::read_rds(here::here("output/rds/stan_dat.rds"))

library(cmdstanr)
model_NC <- cmdstanr::cmdstan_model("stan/testPoisson_AgeYear_NC.stan")
fit <- model_NC$sample(data = stan_dat, parallel_chains = 8,chains = 8,
                       adapt_delta = 0.95,seed = 56465,max_treedepth = 15,
                 iter_sampling = 2000,
                 iter_warmup = 1000)
fit$save_object("/mnt/Storage/LargeRFiles/Dropping-Counts/testPoisson_AgeYear_NC.rds")
fit <- read_rds("/mnt/Storage/LargeRFiles/Dropping-Counts/testPoisson_AgeYear.rds")
s3 <- fit$summary()
readr::write_rds(s3, here::here("output/rds/testmodel_varsNC.rds") )
fit$cmdstan_diagnose()
library(dplyr)
library(ggplot2)
library(tidyr)
s3 %>% filter(grepl("X2", variable)) %>%
  separate(variable, into = c("Var", "group", "D2S", "ex"),
           sep = "\\[|\\]|\\,", remove = F, convert = T) %>%
  mutate(Year = ifelse(group<3, 2007, 2008),
         Age = ifelse(group%%2==0, "Juvenile", "Adults")) %>%
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

