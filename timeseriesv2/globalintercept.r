library(pacman)
p_load(dclone, MASS, ggplot2, snow, tidyverse, parallel, gridExtra)

# Read and prepare data
# 8 for drop in cava, 11 for drop in dipteryx
pheno_start_month <- 10
data <- read.csv("timeseriesv2//data//dipteryx_oleifera_leafing_valid.csv")
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

all_after_threshold <- data.frame()
for (i in 1:length(trees1)) {
  for (j in 1:length(years1)) {
    subset_data <- data %>%
      filter(tree == trees1[i], pheno_year == years1[j]) %>%
      arrange(day)

    if (nrow(subset_data) == 0 || all(is.na(subset_data$y_norm))) {
      next
    }

    min_indices <- which(subset_data$y_norm == min(subset_data$y_norm, na.rm = TRUE))
    cut_idx <- min_indices[1]
    seg_after <- subset_data[cut_idx:nrow(subset_data), ]

    n_obs   <- nrow(seg_after)
    y_range <- max(seg_after$y_norm, na.rm = TRUE) - min(seg_after$y_norm, na.rm = TRUE)

    if (n_obs < 2 || y_range < 0.5) {
      print(paste0("SKIPPED tree=", trees1[i], " year=", years1[j],
                   " n=", n_obs, " range=", round(y_range, 3)))
      print(min(seg_after$y_norm, na.rm = TRUE))
      print(max(seg_after$y_norm, na.rm = TRUE))
      next
    }

    print(paste0("tree=", trees1[i], " year=", years1[j],
                 " n=", n_obs, " range=", round(y_range, 3)))
    all_after_threshold <- bind_rows(all_after_threshold, seg_after)
  }
}

flush_seg <- all_after_threshold %>%
  mutate(
    pheno_year = droplevels(as.factor(pheno_year)),
    tree = droplevels(as.factor(tree))
  )
eps_logit <- 1e-6
flush_prob <- pmin(pmax(flush_seg$y_norm, eps_logit), 1 - eps_logit)
flush_logit <- log((1 - flush_prob) / flush_prob)
flush_year_id <- as.numeric(as.factor(flush_seg$pheno_year))
flush_tree_id <- as.numeric(as.factor(flush_seg$tree))


#sanity check
windows()
ggplot(flush_seg, aes(x = day, y = y_norm, color=tree_year)) +
  geom_line()+
  theme(legend.position = "none")

cl.seq <- c(1, 4, 8, 16, 32)
n.iter <- 10000
n.update <- 5000
n.adapt <- 2000
n.chains <- 3
thin <- 5

#first model, drop population level for cava
start_time <- Sys.time()
leaves.flush <- function(){
  kd    ~ dunif(0.01, 1)
  Td    ~ dunif(100, 365)
  base  ~ dunif(0, 1)
  amp_raw ~ dunif(0, 1)
  sigsq ~ dunif(1, 25)
  amp <- amp_raw * (1 - base)
  for (j in 1:n) {
    eta[j] <- (-1 * kd) * (days[j] - Td)
    mu_y[j] <- base + amp / (1 + exp(eta[j]))
    muf[j] <- log((1 - mu_y[j]) / mu_y[j])
  }
  for (k in 1:K) {
    for (i in 1:n) {
      X[i, k] ~ dnorm(muf[i], 1 / sigsq)
    }
  }
}
data_flush <- list(
  K = 1,
  X = dcdim(data.matrix(flush_logit)),
  days = flush_seg$day,
  n = nrow(flush_seg)
)
inits.drop<- list(
  kd = runif(1, 0.01, 1),
  Td = runif(1, 100, 365),
  base = runif(1, 0, 1),
  amp_raw = runif(1, 0, 1),
  sigsq = runif(1, 1, 25)
)
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(5)
model.flush<- dc.parfit(
  cl,
  data_flush,
  params = c("kd", "Td", "base", "amp", "sigsq"),
  model = leaves.flush,
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
summary(model.flush)
saveRDS(model.flush, "timeseriesv2//data//dipt_flush.rds")
dcdiag(model.flush)
end_time <- Sys.time()
print(paste("Time taken for population-level model:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))


x_range <- seq(min(flush_seg$day), max(flush_seg$day), length.out = 100)
y_pred <- sapply(x_range, function(x) {
  kd <- summary(model.flush)$statistics["kd", "Mean"]
  Td <- summary(model.flush)$statistics["Td", "Mean"]
  base <- summary(model.flush)$statistics["base", "Mean"]
  amp <- summary(model.flush)$statistics["amp", "Mean"]
  eta <- (-1 * kd) * (x - Td)
  mu_y <- base + amp / (1 + exp(eta))
  return(mu_y)
})
df_pred <- data.frame(x_range = x_range, y_pred = y_pred)
windows()
ggplot(flush_seg, aes(x = day, y = y_norm, color=tree_year)) +
  geom_line() +
  geom_line(data = df_pred, aes(x = x_range, y = y_pred), color = "black", linewidth = 1) +
  theme(legend.position = "none")


