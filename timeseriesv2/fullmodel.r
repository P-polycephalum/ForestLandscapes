library(pacman)
p_load(dclone, MASS, ggplot2, snow, tidyverse, parallel, gridExtra)

day_to_DOY <- function(day, pheno_year, start_month = pheno_start_month) {
  if (is.na(start_month)) start_month <- 1
  format(
    as.Date(paste0(as.numeric(as.character(pheno_year)), "-", sprintf("%02d", start_month), "-01")) + day,
    "%B %d"
  )
}
DOY_to_day <- function(month, day, pheno_year, start_month = pheno_start_month) {
  if (is.na(start_month)) start_month <- 1
  cal_year <- if (month >= start_month) pheno_year else pheno_year + 1
  as.numeric(
    as.Date(paste0(cal_year, "-", sprintf("%02d", month), "-", sprintf("%02d", day))) -
      as.Date(paste0(pheno_year, "-", sprintf("%02d", start_month), "-01"))
  )
}

# Read and prepare data
pheno_start_month <- 8
data <- read.csv("timeseriesv2//data//cavanillesia_leafing_valid.csv")
data <- data %>%
  mutate(
    y_norm = pmin(pmax(leafing / 100, 1e-4), 1 - 1e-4),
    date = as.Date(time),
    DOY = yday(date),
    year = year(date),
    month = month(date),
    pheno_year = if (is.na(pheno_start_month)) year else if_else(month >= pheno_start_month, year, year - 1),
    day = if (is.na(pheno_start_month)) DOY else as.numeric(difftime(
      date,
      as.Date(paste0(pheno_year, "-", sprintf("%02d", pheno_start_month), "-01")),
      units = "days"
    )) + 1L,
    tree = as.factor(tag),
    pheno_year = as.factor(pheno_year),
    tree_year = as.factor(paste0(tree, "_", pheno_year)))

trees1 <- unique(data$tree)
years1 <- unique(data$pheno_year)
all_before_threshold <- data.frame()

for (i in 1:length(trees1)) {
  for (j in 1:length(years1)) {
    subset_data <- data %>%
      filter(tree == trees1[i], pheno_year == years1[j]) %>%
      arrange(day)

    if (nrow(subset_data) == 0 || all(is.na(subset_data$y_norm))) {
      next
    }

    min_indices <- which(subset_data$y_norm == min(subset_data$y_norm, na.rm = TRUE))
    cut_idx <- min_indices[length(min_indices)]
    seg_before <- subset_data[1:cut_idx, ]

    n_obs   <- nrow(seg_before)
    y_range <- max(seg_before$y_norm, na.rm = TRUE) - min(seg_before$y_norm, na.rm = TRUE)

    if (n_obs < 2 || y_range < 0.5) {
      print(paste0("SKIPPED tree=", trees1[i], " year=", years1[j],
                   " n=", n_obs, " range=", round(y_range, 3)))
      print(min(seg_before$y_norm, na.rm = TRUE))
      print(max(seg_before$y_norm, na.rm = TRUE))
      next
    }

    print(paste0("tree=", trees1[i], " year=", years1[j],
                 " n=", n_obs, " range=", round(y_range, 3)))
    all_before_threshold <- bind_rows(all_before_threshold, seg_before)
  }
}

drop_seg <- all_before_threshold  %>% filter(pheno_year != "2017") %>%
mutate(
  pheno_year = droplevels(as.factor(pheno_year)),
  tree = droplevels(as.factor(tree))
)
eps_logit <- 1e-6
drop_prob <- pmin(pmax(drop_seg$y_norm, eps_logit), 1 - eps_logit)
drop_logit <- log((1 - drop_prob) / drop_prob)
drop_year_id <- as.numeric(as.factor(drop_seg$pheno_year))
drop_tree_id <- as.numeric(as.factor(drop_seg$tree))

windows()
ggplot(drop_seg, aes(x = day, y = y_norm, color = tree)) +
  geom_line() +
  geom_point(size = 0.8, alpha = 0.6) +
  facet_wrap(~ pheno_year) + theme_bw() 

# i want to know how many years of data does each tree have in the drop segment
tree_year_counts <- drop_seg %>%
  group_by(tree) %>%
  summarise(n_years = n_distinct(pheno_year)) %>%
  arrange(desc(n_years))
print(tree_year_counts)

for (i in 1:nrow(tree_year_counts)) {
  print(paste0("Tree ", tree_year_counts$tree[i], " has data for ", tree_year_counts$n_years[i], " phenological years."))
}

year_tree_counts <- drop_seg %>%
  group_by(pheno_year) %>%
  summarise(n_trees = n_distinct(tree)) %>%
  arrange(desc(n_trees))

for (i in 1:nrow(year_tree_counts)) {
  print(paste0("Phenological Year ", year_tree_counts$pheno_year[i], " has data for ", year_tree_counts$n_trees[i], " trees."))
}


#first model, drop population level for cava
start_time <- Sys.time()
leaves.drop <- function(){
  kd    ~ dunif(0.01, 1)
  Td    ~ dunif(60, 275)
  base  ~ dunif(0, 1)
  amp_raw ~ dunif(0, 1)
  sigsq ~ dunif(1, 25)
  amp <- amp_raw * (1 - base)
  for (j in 1:n) {
    eta[j] <- kd * (days[j] - Td)
    mu_y[j] <- base + amp / (1 + exp(eta[j]))
    muf[j] <- log((1 - mu_y[j]) / mu_y[j])
  }
  for (k in 1:K) {
    for (i in 1:n) {
      X[i, k] ~ dnorm(muf[i], 1 / sigsq)
    }
  }
}
data_drop <- list(
  K = 1,
  X = dcdim(data.matrix(drop_logit)),
  days = drop_seg$day,
  n = nrow(drop_seg)
)
inits.drop<- list(
  kd = runif(1, 0.01, 1),
  Td = runif(1, 50, 275),
  base = runif(1, 0, 1),
  amp_raw = runif(1, 0, 1),
  sigsq = runif(1, 1, 25)
)
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(5)
model.drop<- dc.parfit(
  cl,
  data_drop,
  params = c("kd", "Td", "base", "amp", "sigsq"),
  model = leaves.drop,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.drop
)
summary(model.drop)
saveRDS(model.drop, "timeseriesv2//data//dipt_drop.rds")
dcdiag(model.drop)
end_time <- Sys.time()
print(paste("Time taken for population-level model:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))

windows()
ggplot(drop_seg, aes(x = day, y = y_norm,color= tree)) +
  geom_point(size = 0.8, alpha = 0.6) +
  facet_wrap(~ pheno_year) +
  theme_bw() +
  labs(x = "Day of Phenological Year", y = "Normalized Leafing (%)") +
  theme(legend.position = "None") +
  scale_color_brewer(palette = "Set1", name = "Phenological Year") +
#---------------------------------------------------------------#################
# Drop: year fixed effects (amp + base)

start_time <- Sys.time()
leaves.drop_year_fe <- function(){
  kd    ~ dunif(0.01, 1)
  base  ~ dunif(0, 1)
  amp_raw ~ dunif(0, 1)
  sigsq ~ dunif(1, 25)

  for (y in 1:nyear) {
    yTd[y] ~ dunif(50, 300)
  }

  amp <- amp_raw * (1 - base)

  for (j in 1:n) {
    eta[j] <- kd * (days[j] - yTd[year[j]])
    mu_y[j] <- (amp / (1 + exp(eta[j])))+base
    muf[j] <- log((1 - mu_y[j]) / mu_y[j])
  }

  for (k in 1:K) {
    for (i in 1:n) {
      X[i, k] ~ dnorm(muf[i], 1 / sigsq)
    }
  }
}
data_drop_year_fe <- list(
  K = 1,
  X = dcdim(data.matrix(drop_logit)),
  days = drop_seg$day,
  year = drop_year_id,
  n = nrow(drop_seg),
  nyear = length(unique(drop_year_id))
)
inits.drop_year_fe <- list(
  kd = runif(1, 0.01, 1),
  sigsq = runif(1, 1, 25),
  base = runif(1, 0, 0.9),
  amp_raw = runif(1, 0, 1),
  yTd_raw = runif(length(unique(drop_year_id)), 50, 300)
)
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(5)
leaves.drop_year_fe <- dc.parfit(
  cl,
  data_drop_year_fe,
  params = c("kd","sigsq","yTd"),
  model = leaves.drop_year_fe,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nyear"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.drop_year_fe
)
summary(leaves.drop_year_fe)
dcdiag(leaves.drop_year_fe)
saveRDS(leaves.drop_year_fe, "timeseriesv2//data//dipt_drop_year_fe.rds")
end_time <- Sys.time()
print(paste("Time taken for year fixed effects model:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))


table<-dctable(leaves.drop_year_fe)
windows()
plot(table, 1:8, type="log.var")
plot(table, 9:12, type="log.var")

df_drop_year_fe<- data.frame(
  day = drop_seg$day,
  year = drop_year_id,
  tree = drop_seg$tree
) %>%
  filter(tree == unique(drop_seg$tree)[1]) %>%
  group_by(year) %>%
  mutate(
    interval = day - lag(day, default = first(day)),
    interval_for_stats = if_else(row_number() == 1L, NA_real_, interval),
    interval_abs = abs(day - 143.5409)
  ) %>%
  ungroup()

View(df_drop_year_fe)
View(drop_seg)
df_drop_year_fe_metrics <- df_drop_year_fe %>%
  group_by(year) %>%
  summarise(
    mean_interval = mean(interval_for_stats, na.rm = TRUE),
    median_interval = median(interval_for_stats, na.rm = TRUE),
    max_interval = max(interval_for_stats, na.rm = TRUE),
    obs_within_991_of_mean = sum(interval_abs <= 32.42129, na.rm = TRUE),
    .groups = "drop"
  )
View(df_drop_year_fe_metrics)

#---------------------------------------------------------------#################
start_time <- Sys.time()
leaves.drop_indv_fe<- function(){
  kd      ~ dunif(0.01, 1)
  base    ~ dunif(0, 1)
  amp_raw ~ dunif(0, 1)
  sigsq   ~ dunif(1, 25)

  for (y in 1:nindv) {
    iTd[y] ~ dunif(1, 365)
  }
  amp <- amp_raw * (1 - base)

  for (j in 1:n) {
    eta[j] <- kd * (days[j] - (iTd[indv[j]])) 
    mu_y[j] <- base + amp / (1 + exp(eta[j]))
    muf[j] <- log((1 - mu_y[j]) / mu_y[j])
  }

  for (k in 1:K) {
    for (i in 1:n) {
      X[i, k] ~ dnorm(muf[i], 1 / sigsq)
    }
  }
}
data4dclone_drop_indv_fe <- list(
  K = 1,
  X = dcdim(data.matrix(drop_logit)),
  days = drop_seg$day,
  indv = as.numeric(as.factor(drop_seg$tree)),
  n = nrow(drop_seg),
  nindv = length(unique(drop_seg$tree))
)
inits.drop_indv_fe <- list(
  kd = runif(1, 0.01, 1),
  sigsq = runif(1, 1, 25),
  base = runif(1, 0, 0.9),
  amp_raw = runif(1, 0, 1),
  iTd = runif(length(unique(drop_seg$tree)), 1, 365)
)
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(3)
leaves.drop_indv_fe <- dc.parfit(
  cl,
  data4dclone_drop_indv_fe,
  params = c("kd", "sigsq", "base", "amp", "iTd"),
  model = leaves.drop_indv_fe,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nindv"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.drop_indv_fe
)
saveRDS(leaves.drop_indv_fe, "timeseriesv2//data//dipt_drop_indv_fe.rds")
dcdiag(leaves.drop_indv_fe)
summary(leaves.drop_indv_fe)
end_time <- Sys.time()
print(paste("Time taken for individual fixed effects model:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))

### Drop: year fixed effects (amp + base) + year random effect on kd###3
leaves.drop_tdkd_fe<- function(){ #year and kd fixed effect
  base  ~ dunif(0, 1)
  amp_raw ~ dunif(0, 1)
  sigsq ~ dunif(1, 30)

  for (y in 1:nyear) {
    yTd[y] ~ dunif(60, 215)
    ykd[y] ~ dunif(0.01, 2)
  }

  amp <- amp_raw * (1 - base)

  for (j in 1:n) {
    eta[j] <- ykd[year[j]] * (days[j] - yTd[year[j]])
    mu_y[j] <- base + amp / (1 + exp(eta[j]))
    muf[j] <- log((1 - mu_y[j]) / mu_y[j])
  }

  for (k in 1:K) {
    for (i in 1:n) {
      X[i, k] ~ dnorm(muf[i], 1 / sigsq)
    }
  }
}
data_drop_tdkd_fe <- list(
  K = 1,
  X = dcdim(data.matrix(drop_logit)),
  days = drop_seg$day,
  year = drop_year_id,
  n = nrow(drop_seg),
  nyear = length(unique(drop_year_id))
)
inits.drop_tdkd_fe <- list(
  base = runif(1, 0, 0.9),
  amp_raw = runif(1, 0, 1),
  sigsq = runif(1, 1, 30),
  yTd = runif(length(unique(drop_year_id)), 60, 215),
  ykd = runif(length(unique(drop_year_id)), 0.01, 2)
)
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(3)
leaves.drop_tdkd_fe <- dc.parfit(
  cl,
  data_drop_tdkd_fe,
  params = c("base", "amp_raw", "sigsq", "yTd", "ykd"),
  model = leaves.drop_tdkd_fe,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nyear"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.drop_tdkd_fe
)
saveRDS(leaves.drop_tdkd_fe, "timeseriesv2//data//leaves_drop_tdkd_fe.rds")
###--------------------------------------------#################



leaves.drop_tdkd_indv_fe<- function(){ #year and kd fixed effect
  base  ~ dunif(0, 1)
  amp_raw ~ dunif(0, 1)
  sigsq ~ dunif(1, 30)

  for (y in 1:nindv) {
    iTd[y] ~ dunif(60, 215)
    ikd[y] ~ dunif(0.01, 2)
  }

  amp <- amp_raw * (1 - base)

  for (j in 1:n) {
    eta[j] <- ikd[indv[j]] * (days[j] - iTd[indv[j]])
    mu_y[j] <- base + amp / (1 + exp(eta[j]))
    muf[j] <- log((1 - mu_y[j]) / mu_y[j])
  }

  for (k in 1:K) {
    for (i in 1:n) {
      X[i, k] ~ dnorm(muf[i], 1 / sigsq)
    }
  }
}
data_drop_tdkd_indv_fe <- list(
  K = 1,
  X = dcdim(data.matrix(drop_logit)),
  days = drop_seg$day,
  indv = as.numeric(as.factor(drop_seg$tree)),
  n = nrow(drop_seg),
  nindv = length(unique(drop_seg$tree))
)
inits.drop_tdkd_indv_fe <- list(
  base = runif(1, 0, 0.9),
  amp_raw = runif(1, 0, 1),
  sigsq = runif(1, 1, 30),
  iTd = runif(length(unique(drop_seg$tree)), 60, 215),
  ikd = runif(length(unique(drop_seg$tree)), 0.01, 2)
)
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(3)
model.drop_tdkd_indv_fe <- dc.parfit(
  cl,
  data_drop_tdkd_indv_fe,
  params = c("base", "amp_raw", "sigsq", "iTd", "ikd"),
  model = leaves.drop_tdkd_indv_fe,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nindv"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.drop_tdkd_indv_fe
)
saveRDS(model.drop_tdkd_indv_fe, "timeseriesv2//data//leaves_drop_tdkd_indv_fe.rds")


leaves.drop_tdkd_year_fe<- function(){
  base  ~ dunif(0, 1)
  amp_raw ~ dunif(0, 1)
  sigsq ~ dunif(1, 30)

  for (y in 1:nyear) {
    yTd[y] ~ dunif(60, 215)
    ykd[y] ~ dunif(0.01, 2)
  }

  amp <- amp_raw * (1 - base)

  for (j in 1:n) {
    eta[j] <- ykd[year[j]] * (days[j] - yTd[year[j]])
    mu_y[j] <- base + amp / (1 + exp(eta[j]))
    muf[j] <- log((1 - mu_y[j]) / mu_y[j])
  }

  for (k in 1:K) {
    for (i in 1:n) {
      X[i, k] ~ dnorm(muf[i], 1 / sigsq)
    }
  }
}
data4dclone_drop_tdkd_year_fe <- list(
  K = 1,
  X = dcdim(data.matrix(drop_logit)),
  days = drop_seg$day,
  year = drop_year_id,
  n = nrow(drop_seg),
  nyear = length(unique(drop_year_id))
)
inits.drop_tdkd_year_fe <- list(
  base = runif(1, 0, 0.9),
  amp_raw = runif(1, 0, 1),
  sigsq = runif(1, 1, 30),
  yTd = runif(length(unique(drop_year_id)), 60, 215),
  ykd = runif(length(unique(drop_year_id)), 0.01, 2)
)
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(3)
leaves.drop_tdkd_year_fe <- dc.parfit(
  cl,
  data4dclone_drop_tdkd_year_fe,
  params = c("base", "amp_raw", "sigsq", "yTd", "ykd"),
  model = leaves.drop_tdkd_year_fe,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nyear"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.drop_tdkd_year_fe
)
saveRDS(leaves.drop_tdkd_year_fe, "timeseriesv2//data//leaves_drop_tdkd_year_fe.rds")

leaves.drop_amp_base_indv_fe <- function(){
  kd      ~ dunif(0.01, 2)
  base    ~ dunif(0, 1)
  amp_raw ~ dunif(0, 1)
  sigsq   ~ dunif(1, 30)

  for (y in 1:nindv) {
    iTd[y] ~ dunif(60, 215)
  }
  amp <- amp_raw * (1 - base)

  for (j in 1:n) {
    eta[j]  <- kd * (days[j] - (iTd[indv[j]])) 
    mu_y[j] <- base + amp / (1 + exp(eta[j]))
    muf[j]  <- log((1 - mu_y[j]) / mu_y[j])
  }

  for (k in 1:K) {
    for (i in 1:n) {
      X[i, k] ~ dnorm(muf[i], 1 / sigsq)
    }
  }
}

#model for dclone year random effect only and also fixed effect. it works for both

cl.seq <- c(1, 4, 8, 16, 32)
n.iter <- 12000
n.update <- 5000
n.adapt <- 2000
n.chains <- 3
thin <- 5
start_time <- Sys.time()
leaves.drop_yearfe_indvre <- function(){
  kd      ~ dunif(0.01, 1)
  base    ~ dunif(0, 1)
  amp_raw ~ dunif(0, 1)
  sigsq   ~ dunif(1, 25)
  sigma_i ~ dunif(1, 40)
  tau_i   <- 1/pow(sigma_i, 2)

  for (y in 1:nyear) {
    yTd[y] ~ dunif(50, 210)
  }

  # Individual random effects
  for (a in 1:nindv) {
    iTd_raw[a] ~ dnorm(0, tau_i)
  }
  iTd_bar <- mean(iTd_raw[])
  for (a in 1:nindv) {
    iTd[a] <- iTd_raw[a] - iTd_bar
  }
  amp <- amp_raw * (1 - base)

  for (j in 1:n) {
    timing[j] <- yTd[year[j]] + iTd[indv[j]]
    eta[j] <- kd * (days[j] - timing[j])
    mu_y[j] <- base + amp / (1 + exp(eta[j]))
    muf[j] <- log((1 - mu_y[j]) / mu_y[j])
  }
  for (k in 1:K) {
    for (j in 1:n) {
      X[j, k] ~ dnorm(muf[j], 1 / sigsq)
    }
  }
}
data_drop_yearfe_indvre <- list(
  K = 1,
  X = dcdim(data.matrix(drop_logit)),
  days = drop_seg$day,
  year = drop_year_id,
  indv = drop_tree_id,
  n = nrow(drop_seg),
  nyear = length(unique(drop_year_id)),
  nindv = length(unique(drop_tree_id))
)
inits.drop_yearfe_indvre <- list(
  kd = runif(1, 0.01, 1),
  base = runif(1, 0, 0.9),
  amp_raw = runif(1, 0, 1),
  sigsq = runif(1, 1, 25),
  yTd = runif(length(unique(drop_year_id)), 50, 210),
  sigma_i = runif(1, 1, 40),
  iTd_raw = rnorm(length(unique(drop_tree_id)), 0, sd=runif(1, 1, 40))
)

try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(5)
model.drop_yearfe_indvre <- dc.parfit(
  cl,
  data_drop_yearfe_indvre,
  params = c("amp","base", "kd", "sigsq", "iTd" , "yTd"),
  model = leaves.drop_yearfe_indvre,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nyear", "nindv"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.drop_yearfe_indvre
)
dcdiag(model.drop_yearfe_indvre)
diag_test(dcdiag(model.drop_yearfe_indvre))
summary(model.drop_yearfe_indvre)
saveRDS(model.drop_yearfe_indvre, "timeseriesv2//data//cava_drop_yearfe_indvre.rds")

end_time <- Sys.time()
print(paste("Time taken for year fixed effects + individual random effects model:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))

diag_test(dcdiag(model.drop_yearfe_indvre))
table <- dctable(model.drop_yearfe_indvre)
summary(model.drop_yearfe_indvre)
windows()
plot(table, 1:6,type="log.var")

# windows()
# plot(table, 1:8,type="log.var")
# plot(table, 9:16,type="log.var")
# plot(table, 17:24,type="log.var")
# plot(table, 25:32,type="log.var")
# plot(table, 33:34,type="log.var")
end_time <- Sys.time()
print(paste("Time taken for year fixed effects + individual random effects model:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))

results_yearfe_indvre <- summary(model.drop_yearfe_indvre)
print(results_yearfe_indvre$statistics)

amp_mean <- results_yearfe_indvre$statistics["amp", "Mean"]
base_mean <- results_yearfe_indvre$statistics["base", "Mean"]
kd_mean <- results_yearfe_indvre$statistics["kd", "Mean"]

iTd_mean <- sapply(1:length(unique(drop_seg$tree)), function(i) {
  results_yearfe_indvre$statistics[paste0("iTd[", i, "]"), "Mean"]
})
yTd_mean <- sapply(1:length(unique(drop_year_id)), function(y) {
  results_yearfe_indvre$statistics[paste0("yTd[", y, "]"), "Mean"]
})

x_range <- seq(min(drop_seg$day), max(drop_seg$day), length.out = 100)
year_lookup <- tibble(
  year    = seq_along(levels(drop_seg$pheno_year)),
  pheno_year = levels(drop_seg$pheno_year)
)
indv_lookup <- tibble(
  indv = seq_along(levels(drop_seg$tree)),
  tree = levels(drop_seg$tree)
)
df_pred <- tidyr::expand_grid(day = x_range, year = year_lookup$year, indv = indv_lookup$indv) %>%
  left_join(year_lookup, by = "year") %>%
  left_join(indv_lookup, by = "indv") %>%
  mutate(
    timing    = yTd_mean[year] + iTd_mean[indv],
    predicted_y = base_mean + amp_mean / (1 + exp(kd_mean * (day - timing))),
    tree_year = paste0(tree, "_", pheno_year)
  )

windows()
ggplot(drop_seg, aes(x = day, y = y_norm, color = tree)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_line(data = df_pred, aes(x = day, y = predicted_y, group = tree_year, color = tree), linewidth = 0.4) +
  facet_wrap(~ pheno_year)

ggplot(drop_seg, aes(x=day, y=y_norm, color=as.factor(year))) +
  geom_point(size=2, alpha=0.6) +
  geom_line(data=df_pred, aes(x=day, y=predicted_y, group=pheno_year, color=pheno_year), linewidth=0.4) +
  facet_wrap(~ tree)
data4dclone_drop_ampbase_indv_fe <- list(
  K = 1,
  X = dcdim(data.matrix(drop_logit)),
  days = drop_seg$day,
  indv = as.numeric(as.factor(drop_seg$tree)),
  n = nrow(drop_seg),
  nindv = length(unique(drop_seg$tree))
)

data4dclone_drop_ampbase_yearkd_fe <- list(
  K = 1,
  X = dcdim(data.matrix(drop_logit)),
  days = drop_seg$day,
  year = drop_year_id,
  n = nrow(drop_seg),
  nyear = length(unique(drop_year_id))
)

# Data lists for dclone
data4dclone_drop_ampbase <- list(
  K = 1,
  X = dcdim(data.matrix(drop_logit)),
  n = nrow(drop_seg),
  days = drop_seg$day
)

#data for dclone year random effect only and also fixed effect. it works for both 
data4dclone_drop_ampbase_year_fe <- list(
  K = 1,
  X = dcdim(data.matrix(drop_logit)),
  days = drop_seg$day,
  year = drop_year_id,
  n = nrow(drop_seg),
  nyear = length(unique(drop_year_id))
)

# Initial values
inits.dropampbase <- list(
  kd = runif(1, 0.01, 2),
  sigsq = runif(1, 1, 30),
  Td = runif(1, 60, 215),
  base = runif(1, 0, 1),
  amp_raw = runif(1, 0, 1)
)

inits.drop_ampbase_year_fe <- list(
  kd = runif(1, 0.01, 2),
  sigsq = runif(1, 1, 30),
  base = runif(1, 0, 0.9),
  amp_raw = runif(1, 0, 1),
  yTd = runif(length(unique(drop_year_id)), 60, 215)
)

inits.drop_ambase_yearkd_fe <- list(
  sigsq = runif(1, 1, 30),
  base = runif(1, 0, 0.9),
  amp_raw = runif(1, 0, 1),
  yTd = runif(length(unique(drop_year_id)), 60, 215),
  ykd = runif(length(unique(drop_year_id)), 0.01, 2)
)
inits.drop_ambase_year_re<- list(
  Td = runif(1, 60, 215),
  kd = runif(1, 0.01, 2),
  base = runif(1, 0, 0.9),
  amp_raw = runif(1, 0, 1),
  sigsq = runif(1, 1, 30),
  sigma_year = runif(1, 1, 60)
)

inits.drop_ampbase_indv_fe <- list(
  kd = runif(1, 0.01, 2),
  sigsq = runif(1, 1, 30),
  base = runif(1, 0, 0.9),
  amp_raw = runif(1, 0, 1),
  iTd = runif(length(unique(drop_seg$tree)), 60, 215)
)

inits.drop_ampbase_yearfe_indvre <- list(
  kd = runif(1, 0.01, 2),
  base = runif(1, 0, 0.9),
  amp_raw = runif(1, 0, 1),
  sigsq = runif(1, 1, 30),
  sigma_indv = runif(1, 1, 60),
  yTd = runif(length(unique(drop_year_id)), 60, 215),
  iTd_raw = rnorm(length(unique(drop_seg$tree)), mean = 0, sd = runif(1, 1, 60))
)

# Parallel MCMC settings
cl.seq <- c(1, 4, 8)
n.iter <- 10000
n.adapt <- 2000
n.update <- 5000
thin <- 5
n.chains <- 3

try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(3)
leaves.drop_ampbase_yearfe_indvre <- dc.parfit(
  cl,
  data4dclone_drop_ampbase_yearfe_indvre,
  params = c("kd", "sigsq", "base", "amp_raw", "yTd", "iTd"),
  model = leaves.drop_amp_base_yearfe_indvre,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nyear", "nindv"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.drop_ampbase_yearfe_indvre
)
saveRDS(leaves.drop_ampbase_yearfe_indvre, "timeseriesv2//data//leaves_drop_ampbase_yearfe_indvre.rds")
summary(leaves.drop_ampbase_yearfe_indvre)
dcdiag(leaves.drop_ampbase_yearfe_indvre)
table<- dctable(leaves.drop_ampbase_yearfe_indvre)

windows()
plot(table,1:6,type='log.var')

results<-summary(leaves.drop_ampbase_yearfe_indvre)
print(results$statistics["kd", "Mean"])
print(results$statistics["kd", "DC SD"])
print(results$statistics["sigsq", "Mean"])
print(results$statistics["sigsq", "DC SD"])
iTd_mean <- sapply(1:length(unique(drop_seg$tree)), function(i) {
  results$statistics[paste0("iTd[", i, "]"), "Mean"]
})
yTd_mean <- sapply(1:length(unique(drop_year_id)), function(y) {
  results$statistics[paste0("yTd[", y, "]"), "Mean"]
})

base_mean <- results$statistics["base", "Mean"]
amp_raw_mean <- results$statistics["amp_raw", "Mean"]
amp_mean <- amp_raw_mean * (1 - base_mean)
kd_mean <- results$statistics["kd", "Mean"]

year_lookup <- tibble(
  year = sort(unique(drop_year_id)),
  pheno_year = levels(as.factor(drop_seg$pheno_year))
)

indv_lookup <- tibble(
  indv = sort(unique(as.numeric(as.factor(drop_seg$tree)))),
  tree = levels(as.factor(drop_seg$tree))
)

x_range <- seq(min(drop_seg$day), max(drop_seg$day))

df_pred <- tidyr::expand_grid(
  day = x_range,
  year = year_lookup$year,
  indv = indv_lookup$indv
) %>%
  left_join(year_lookup, by = "year") %>%
  left_join(indv_lookup, by = "indv") %>%
  mutate(
    timing = yTd_mean[match(year, year_lookup$year)] +
      iTd_mean[match(indv, indv_lookup$indv)],
    y_pred = base_mean + amp_mean / (1 + exp(kd_mean * (day - timing))),
    yTd = yTd_mean[match(year, year_lookup$year)],
    iTd = iTd_mean[match(indv, indv_lookup$indv)]
  )

windows()
ggplot(df_pred, aes(x = day, y = y_pred, group = tree, color = tree)) +
  geom_line(alpha = 0.85) +
  geom_point(
    data = drop_seg,
    aes(x = day, y = y_norm, group = tree, color = tree),
    alpha = 0.35,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ pheno_year, ncol = 1) +
  labs(
    title = "Drop model: year fixed effect + tree random effect",
    x = "Day of Year",
    y = "Leafing"
  ) +
  theme_minimal() +
  theme(legend.position = "right")



df_pred_filtered<- df_pred %>%
  filter(tree %in% c("134166","4468")) %>%
  filter(pheno_year %in% c("2024","2018"))

windows()
ggplot(df_pred_filtered, aes(x = day, y = y_pred, group = tree, color = tree)) +
  geom_line(alpha = 0.85) +
  geom_point(
    data = drop_seg %>% filter(tree %in% c("134166","4468")) %>% filter(pheno_year %in% c("2024","2018")),
    aes(x = day, y = y_norm, group = tree, color = tree),
    alpha = 0.35,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ pheno_year, ncol = 1) +
  labs(
    title = "Drop model: year fixed effect + tree random effect",
    x = "Day of Year",
    y = "Leafing"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

dcdiag(leaves.drop_ampbase_yearfe_indvre)
diag_test(dcdiag(leaves.drop_ampbase_yearfe_indvre))

try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(3)
leaves.drop_ampbase_indv_fe <- dc.parfit(
  cl,
  data4dclone_drop_ampbase_indv_fe,
  params = c("kd", "sigsq", "base", "amp_raw", "iTd"),
  model = leaves.drop_amp_base_indv_fe,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nindv"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.drop_ampbase_indv_fe
)
summary(leaves.drop_ampbase_indv_fe)
dcdiag(leaves.drop_ampbase_indv_fe)

#fit the year random effects model
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(3)
leaves.drop_ampbase_year_re <- dc.parfit(
  cl,
  data4dclone_drop_ampbase_year_fe,
  params = c("kd", "sigsq", "base", "amp_raw", "yTd"),
  model = leaves.drop_amp_base_year_re,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nyear"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.drop_ambase_year_re
)

# Fit drop amp+base global model
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(3)
leaves.drop_ampbase <- dc.parfit(
  cl,
  data4dclone_drop_ampbase,
  params = c("kd", "sigsq", "Td", "base", "amp_raw"),
  model = leaves.drop_amp_base,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.dropampbase
)
saveRDS(leaves.drop_ampbase, "timeseriesv2//results//leaves_drop_ampbase.rds")

# Fit drop amp+base year fixed-effects model
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(3)
leaves.drop_ampbase_year_fe <- dc.parfit(
  cl,
  data4dclone_drop_ampbase_year_fe,
  params = c("kd", "sigsq", "base", "amp_raw", "yTd"),
  model = leaves.drop_amp_base_year_fe,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nyear"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.drop_ampbase_year_fe
)
saveRDS(leaves.drop_ampbase_year_fe, "timeseriesv2//results//leaves_drop_ampbase_year_fe.rds")
# read it back in 
leaves.drop_ampbase_year_fe <- readRDS("timeseriesv2//data//leaves_drop_ampbase_year_fe.rds")
#fit drop amp+base year kd random effects model
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(3)
leaves.drop_ampbase_yearkd_fe <- dc.parfit(
  cl,
  data4dclone_drop_ampbase_yearkd_fe,
  params = c("sigsq", "base", "amp_raw", "yTd", "ykd"),
  model = leaves.drop_amp_base_yearkd_fe,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nyear"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.drop_ambase_yearkd_fe
)



summary(leaves.drop_ampbase_year_re)
dcdiag(leaves.drop_ampbase_year_re)
dctable(leaves.drop_ampbase_year_re)

summary(leaves.drop_ampbase_yearkd_fe)

day_to_DOY(round(mean(summary(leaves.drop_ampbase_yearkd_fe)$statistics[paste0("yTd[", 1:length(unique(drop_year_id)), "]"), "Mean"]),0), pheno_year = 2020, start_month = pheno_start_month)

mean(summary(leaves.drop_ampbase_yearkd_fe)$statistics[paste0("ykd[", 1:length(unique(drop_year_id)), "]"), "Mean"])




# Summaries and diagnostics
result_drop_ampbase <- summary(leaves.drop_ampbase)
result_ampbase <- summary(leaves.drop_ampbase_year_fe)

print(result_drop_ampbase)
print(result_ampbase)
diag_test(dcdiag(leaves.drop_ampbase))
diag_test(dcdiag(leaves.drop_ampbase_year_fe))

# Global model predictions
df_pred_drop_ampbase <- data.frame(
  day = drop_seg$day,
  y_norm = drop_seg$y_norm,
  y_pred = result_drop_ampbase$statistics["base", "Mean"] +
    (result_drop_ampbase$statistics["amp_raw", "Mean"] * (1 - result_drop_ampbase$statistics["base", "Mean"])) /
    (1 + exp(result_drop_ampbase$statistics["kd", "Mean"] * (drop_seg$day - result_drop_ampbase$statistics["Td", "Mean"])))
)

windows()
ggplot(df_pred_drop_ampbase, aes(x = day)) +
  geom_point(aes(y = y_norm), color = "blue", alpha = 0.5) +
  geom_line(aes(y = y_pred), color = "red") +
  labs(title = "Drop Amp/Base Global Model", x = "Day of Year", y = "Leafing") +
  theme_minimal()

# Year FE model predictions by year
results <- summary(leaves.drop_ampbase_indv_fe)
x_range <- seq(min(drop_seg$day), max(drop_seg$day))
indv_lookup <- tibble(
  indv = sort(unique(as.numeric(as.factor(drop_seg$tree)))),
  tree = levels(as.factor(drop_seg$tree))
)

iTd_mean <- sapply(indv_lookup$indv, function(i) {
  results$statistics[paste0("iTd[", i, "]"), "Mean"]
})

mean(iTd_mean)
day_to_DOY(round(mean(iTd_mean),0), pheno_year = 2020, start_month = pheno_start_month)

kd_mean <- results$statistics["kd", "Mean"]


print(results$statistics["kd", "Mean"])
print(results$statistics["kd", "DC SD"])

print(results$statistics["sigsq", "Mean"])
print(results$statistics["sigsq", "DC SD"])

x_range <- seq(min(drop_seg$day), max(drop_seg$day))
year_lookup <- tibble(
  year = sort(unique(drop_year_id)),
  pheno_year = levels(as.factor(drop_seg$pheno_year))
)
base_mean <- results$statistics["base", "Mean"]
amp_raw_mean <- results$statistics["amp_raw", "Mean"]
kd_mean <- results$statistics["kd", "Mean"]

df_pred_year_indv_fe <- tidyr::expand_grid(
  day = x_range,
  indv = indv_lookup$indv
) %>%
  left_join(indv_lookup, by = "indv") %>%
  mutate(
    y_pred = base_mean +
      (amp_raw_mean * (1 - base_mean)) /
      (1 + exp(kd_mean * (day - iTd_mean[match(indv, indv_lookup$indv)])))
  )

windows()
ggplot(df_pred_year_indv_fe, aes(x = day, y = y_pred, color = tree)) +
  geom_line() +
  geom_point(data = drop_seg, aes(x = day, y = y_norm, color = pheno_year), alpha = 0.4, inherit.aes = FALSE) +
  labs(title = "Drop Amp/Base Individual FE Model", x = "Day of Year", y = "Leafing") +
  theme_minimal() +
  theme(legend.position = "right")


yTd_mean <- sapply(year_lookup$year, function(y) {
  result_ampbase$statistics[paste0("yTd[", y, "]"), "Mean"]
})

base_mean <- result_ampbase$statistics["base", "Mean"]
amp_raw_mean <- result_ampbase$statistics["amp_raw", "Mean"]
kd_mean <- result_ampbase$statistics["kd", "Mean"]

df_pred_year_fe <- tidyr::expand_grid(
  day = x_range,
  year = year_lookup$year
) %>%
  left_join(year_lookup, by = "year") %>%
  mutate(
    yTd_mean = yTd_mean[match(year, year_lookup$year)],
    y_pred = base_mean +
      (amp_raw_mean * (1 - base_mean)) /
      (1 + exp(kd_mean * (day - yTd_mean)))
  )

windows()
ggplot(df_pred_year_fe, aes(x = day, y = y_pred, color = pheno_year)) +
  geom_line() +
  geom_point(data = drop_seg, aes(x = day, y = y_norm, color = pheno_year), alpha = 0.4, inherit.aes = FALSE) +
  facet_wrap(~ pheno_year, ncol = 1) +
  labs(title = "Drop Amp/Base Year FE Model", x = "Day of Year", y = "Leafing") +
  theme_minimal() +
  theme(legend.position = "right")