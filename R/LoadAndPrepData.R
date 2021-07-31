library(tidyverse)
l <- list.files("Data", full.names = T)
d <- map_df(l,read_csv) %>%
  transmute(#Age = Age,
    D2S = `Ditance From Shore`,
            DOY = `Julian Date`,
            Year = ifelse(is.na(Year), 2007, Year),
            ExposureTime_HR = `Exposure Time (Hr)`,
            ExposureTime_MM = `Exposure Time (min)`,
            DC = `DC of (15 ) quadrats`) %>%
  mutate(gr = case_when(Year==2007&DOY<212~1,
                        Year==2007~2,
                        Year==2008&DOY<212~3,
                        Year==2008~4)) %>%
  replace_na(list(DC=0)) %>%
  filter(D2S<=1000)
glimpse(d)


# int n_doy;
# int n_year;
# int n_d2s;
# int<lower=0> n; // Number of rows
# int d2s[n]; // Distance to shore
# int doy[n]; // Day of Year
# int year[n]; // Year
# vector[n] exposure; // Minutes of exposure
# int DC[n]; // Counts summed by 15 quadrats
# matrix[n_doy, n_doy] D_doy; // Distance between points in day of year
# matrix[n_d2s, n_d2s] D_d2s;


D_doy <-
expand_grid(X = sort(unique(d$DOY)), Y = sort(unique(d$DOY))) %>%
  mutate(d = abs(X-Y)) %>%
  pivot_wider(id_cols = X, names_from = Y, values_from = d) %>%
  column_to_rownames("X") %>%
  as.matrix()
D_d2s <-
  expand_grid(X = sort(unique(d$D2S)), Y = sort(unique(d$D2S))) %>%
  mutate(d = abs(X-Y)) %>%
  pivot_wider(id_cols = X, names_from = Y, values_from = d) %>%
  column_to_rownames("X") %>%
  as.matrix()



stan_dat <-
  d %>%
  mutate(doy = as.numeric(factor(DOY)),
         d2s = as.numeric(factor(D2S)),
         year = as.numeric(factor(Year)),
         exposure = log(ExposureTime_MM)) %>%
  tidybayes::compose_data(
    n_doy = n_distinct(DOY),
    n_d2s = n_distinct(D2S),
    n_year = n_distinct(Year),
    D_doy = D_doy,
    D_d2s = D_d2s,
    n_groups = 4
  )
write_rds(stan_dat, here::here("output/rds/stan_dat.rds"))
library(cmdstanr)
m <- cmdstanr::cmdstan_model("stan/testMod.stan")
fit <- m$sample(data = stan_dat, chains = 1, iter_sampling = 1000, iter_warmup = 1000)
m2 <- cmdstanr::cmdstan_model("stan/testPoisson.stan")
fit <- m2$sample(data = stan_dat, chains = 1, iter_sampling = 1000, iter_warmup = 1000)
sum_ary <- fit$summary()
sum_ary$variable
sum_ary %>% filter(grepl("X1", variable)) %>%
  separate(variable, into = c("Var", "DOY", "ex"),
           sep = "\\[|\\]", remove = F, convert = T) %>%
  ggplot(aes(DOY, mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2)+
  geom_line()

sum_ary %>% filter(grepl("X2", variable)) %>%
  separate(variable, into = c("Var", "D2S", "ex"),
           sep = "\\[|\\]", remove = F, convert = T) %>%
  ggplot(aes(D2S, mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2)+
  geom_line()


m3 <- cmdstanr::cmdstan_model("stan/testPoisson_AgeYear.stan")
fit <- m3$sample(data = stan_dat, parallel_chains = 8,chains = 8,
                 iter_sampling = 2000,
                 iter_warmup = 1000)
fit$save_object("/mnt/Storage/LargeRFiles/Dropping-Counts/testPoisson_AgeYear_4gp.rds")
fit <- read_rds("/mnt/Storage/LargeRFiles/Dropping-Counts/testPoisson_AgeYear.rds")
 s3 <- fit$summary()
write_rds(s3, here::here("output/rds/testmodel_vars2.rds"))
s3 %>% filter(grepl("X2", variable)) %>%
  separate(variable, into = c("Var", "group", "D2S", "ex"),
           sep = "\\[|\\]|\\,", remove = F, convert = T) %>%
  mutate(Year = ifelse(group<3, 2007, 2008),
         Age = ifelse(group%%2==0, "Juvenile", "Adults")) %>%
  ggplot(aes(D2S, mean, colour =Age)) +
  geom_pointrange(aes(ymin = q5, ymax = q95, colour = Age))+
  geom_line() + facet_wrap(~Year)


s3 %>% filter(grepl("X1", variable)) %>%
  separate(variable, into = c("Var", "year","DOY", "ex"),
           sep = "\\[|\\]|\\,", remove = F, convert = T) %>%
  ggplot(aes(DOY, mean, colour = factor(year))) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.2)+
  geom_line()
s3 %>% filter(grepl("alpha", variable)) %>%
  ggplot(aes(variable, mean)) +
  geom_pointrange(aes(ymin = q5, ymax = q95), alpha = 0.2)+
  geom_line()

