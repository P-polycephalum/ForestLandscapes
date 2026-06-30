library(pacman)
p_load(dclone, MASS, ggplot2, snow, tidyverse, parallel, gridExtra)

# Read and prepare data
# 8 for drop in cava, 10 for flush in cava
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
    cut_idx<- min_indices[1]
    seg_after <- subset_data[cut_idx:nrow(subset_data), ]
    #loop for drop
    n_obs_before   <- nrow(seg_after)
    y_range_after <- max(seg_after$y_norm, na.rm = TRUE) - min(seg_after$y_norm, na.rm = TRUE)

    
    y_range_after <- max(seg_after$y_norm, na.rm = TRUE) - min(seg_after$y_norm, na.rm = TRUE)
    if (n_obs_before < 2 || y_range_after < 0.5) {
      print(paste0("SKIPPED tree=", trees1[i], " year=", years1[j],
                   " n=", n_obs_before, " range=", round(y_range_after, 3)))
      next
    }
    print(paste0("tree=", trees1[i], " year=", years1[j],
                 " n=", n_obs_before, " range=", round(y_range_after, 3)))
    all_after_threshold <- bind_rows(all_after_threshold, seg_after)
  }
}
flush_seg <- all_after_threshold %>%
  filter(!tree %in% c("3811", "4250","2619")) %>%
  mutate(
    pheno_year = droplevels(as.factor(pheno_year)),
    tree = droplevels(as.factor(tree))
  )

flush_logit <- log((1 - flush_seg$y_norm) / flush_seg$y_norm)
flush_year_id <- as.numeric(as.factor(flush_seg$pheno_year))
flush_tree_id <- as.numeric(as.factor(flush_seg$tree))

max(flush_seg$y_norm)
min(flush_seg$y_norm)
windows()
ggplot(flush_seg, aes(x = day, y = y_norm, color=tree_year)) +
  geom_point()+
  facet_wrap(~tree)+
  theme(legend.position = "none")

cl.seq <- c(1, 4, 8, 16, 32)
n.iter <- 10000
n.update <- 5000
n.adapt <- 2000
n.chains <- 3
thin <- 5
#first model, drop population level for cava
start_time <- Sys.time()
leaves.flush_indv_fe <- function(){
  kd    ~ dunif(0.01, 1)
  base  ~ dunif(0, 1)
  amp_raw ~ dunif(0, 1)
  sigsq ~ dunif(1, 25)

  for (i in 1:nindv) {
    iTd[i] ~ dunif(1, 365)
  }
  amp <- amp_raw * (1 - base)

  for (j in 1:n) {
    eta[j] <-  (-1*kd) * (days[j] - iTd[indv[j]])
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
  X = dcdim(data.matrix(flush_logit)),
  days = flush_seg$day,
  n = nrow(flush_seg),
  indv = flush_tree_id,
  nindv = length(unique(flush_tree_id))
)
inits.flush<- list(
  kd = runif(1, 0.01, 1),
  iTd = runif(length(unique(flush_tree_id)), 1, 365),
  base = runif(1, 0, 1),
  amp_raw = runif(1, 0, 1),
  sigsq = runif(1, 1, 25)
)
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(5)
model.flush_indv_fe<- dc.parfit(
  cl,
  data_drop,
  params = c("kd", "iTd", "base", "amp", "sigsq"),
  model = leaves.flush_indv_fe,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nindv"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.flush
)
summary(model.flush_indv_fe)
saveRDS(model.flush_indv_fe, "timeseriesv2//data//dipt_flush_indv_fe.rds")
end_time <- Sys.time()
print(paste("Time taken for individual-level model:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))



model.flush_indv_fe <- readRDS("timeseriesv2//data//dipt_flush_indv_fe.rds") # we skip trees 3811 and 4250
diag_test(dcdiag(model.flush_indv_fe))
summary(model.flush_indv_fe)
base_hat <- summary(model.flush_indv_fe)$statistics["base", "Mean"]
amp_hat  <- summary(model.flush_indv_fe)$statistics["amp", "Mean"]
kd_hat   <- summary(model.flush_indv_fe)$statistics["kd", "Mean"]

tree_levels <- levels(flush_seg$tree)

iTd_hat <- summary(model.flush_indv_fe)$statistics[
  paste0("iTd[", seq_along(tree_levels), "]"),
  "Mean"
]

mean_iTd <- mean(iTd_hat)
sd_iTd <- sd(iTd_hat)
print("The mean timing of the drop transition for each tree is:")
print(day_to_DOY(mean(iTd_hat), 2019, start_month = pheno_start_month))
print("The max difference in timing of the drop transition for each tree is:")
print(round(max(iTd_hat) - min(iTd_hat), 4))


print("The rate of the flush transition for each tree is:")
print(round(summary(model.flush_indv_fe)$statistics["kd", "Mean"], 4))
print("The DC SD of the rate is:")
print(round(summary(model.flush_indv_fe)$statistics["kd", "DC SD"], 4))
print("The residual variance for each tree is:")
print(round(summary(model.flush_indv_fe)$statistics["sigsq", "Mean"], 4))
print("The DC SD of the residual variance is:")
print(round(summary(model.flush_indv_fe)$statistics["sigsq", "DC SD"], 4))
print("The duration of the leaf flush transition for each tree is:")
print(round(dur_trans(model.flush_indv_fe), 4))
print("The sd of the timing parameters for individual trees is:")
print(round(sd_iTd, 4))

names(iTd_hat) <- tree_levels


days <- seq(1, 365, by = 1)

df_drop_indv <- bind_rows(lapply(tree_levels, function(tree_name) {

  data.frame(
    day = days,
    tree = tree_name,
    y_norm = base_hat + amp_hat / (
      1 + exp((-1*kd_hat) * (days - iTd_hat[tree_name]))
    ),
    transition = "drop"
  )

})) %>%
  mutate(tree = factor(tree, levels = tree_levels))
windows()
ggplot(df_drop_indv, aes(x = day, y = y_norm, color=tree)) +
  geom_line()+
  geom_point(data = flush_seg, aes(x = day, y = y_norm, color=pheno_year), alpha = 0.5)+
  # geom_segment(x=min(iTd_hat), xend=max(iTd_hat), y=0.5, yend=0.5, color="black", size=1)+
  # geom_segment(x=min(iTd_hat), xend=min(iTd_hat), y=0.45, yend=0.55, color="black", size=1)+
  # geom_segment(x=max(iTd_hat), xend=max(iTd_hat), y=0.45, yend=0.55, color="black", size=1)+
  # geom_text(x=mean(iTd_hat)+60, y=0.50, label=paste(round(max(iTd_hat) - min(iTd_hat), 0)," days" ), color="black", size=5)+
  facet_wrap(~tree)+
  theme(legend.position = "none")

table <- dctable(model.flush_indv_fe)

windows()
plot(table, 1:8, type="log.var")
plot(table, 9:16, type="log.var")
plot(table, 17:24, type="log.var")
plot(table, 25:28, type="log.var")
