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
  filter(pheno_year != "2019") %>%
  mutate(
    pheno_year = droplevels(as.factor(pheno_year)),
    tree = droplevels(as.factor(tree))
)
eps_logit <- 1e-6
flush_prob <- pmin(pmax(flush_seg$y_norm, eps_logit), 1 - eps_logit)
flush_logit <- log((1 - flush_prob) / flush_prob)
flush_year_id <- as.numeric(as.factor(flush_seg$pheno_year))
flush_tree_id <- as.numeric(as.factor(flush_seg$tree))

windows()
ggplot(flush_seg, aes(x = day, y = y_norm, color=tree_year)) +
  geom_line()+
  facet_wrap(~pheno_year)+
  theme(legend.position = "none")

tree_year_counts <- flush_seg %>%
  group_by(tree) %>%
  summarise(n_years = n_distinct(pheno_year)) %>%
  arrange(desc(n_years))
print(tree_year_counts)

for (i in 1:nrow(tree_year_counts)) {
  print(paste0("Tree ", tree_year_counts$tree[i], " has data for ", tree_year_counts$n_years[i], " phenological years."))
}

year_tree_counts <- flush_seg %>%
  group_by(pheno_year) %>%
  summarise(n_trees = n_distinct(tree)) %>%
  arrange(desc(n_trees))

for (i in 1:nrow(year_tree_counts)) {
  print(paste0("Phenological Year ", year_tree_counts$pheno_year[i], " has data for ", year_tree_counts$n_trees[i], " trees."))
}

cl.seq <- c(1, 2)
n.iter <- 10000/2
n.update <- 5000/2
n.adapt <- 2000/2
n.chains <- 3
thin <- 2
#first model, drop population level for cava
start_time <- Sys.time()
leaves.flush_year_fe <- function(){
  kd    ~ dunif(0.01, 1)
  base  ~ dunif(0, 1)
  amp_raw ~ dunif(0, 1)
  sigsq ~ dunif(1, 25)

  for (i in 1:nyear) {
    yTd[i] ~ dunif(120, 360)
  }

  amp <- amp_raw * (1 - base)

  for (j in 1:n) {
    eta[j] <- (-1 * kd) * (days[j] - yTd[year[j]])
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
  n = nrow(flush_seg),
  year = flush_year_id,
  nyear = length(unique(flush_year_id))
)
inits.flush <- list(
  kd = runif(1, 0.01, 1),
  yTd = runif(length(unique(flush_year_id)), 120, 360),
  base = runif(1, 0, 1),
  amp_raw = runif(1, 0, 1),
  sigsq = runif(1, 1, 25)
)
try(stopCluster(cl), silent = TRUE)
cl <- makePSOCKcluster(2)
model.flush_yearfe_dipt<- dc.parfit(
  cl,
  data_flush,
  params = c("amp", "base", "kd", "yTd", "sigsq"),
  model = leaves.flush_year_fe,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = c("n", "nyear"),
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits.flush
)
summary(model.flush_yearfe_dipt)
dcdiag(model.flush_yearfe_dipt)
saveRDS(model.flush_yearfe_dipt, "timeseriesv2//data//dipt_flush_year_fe.rds")
end_time <- Sys.time()
print(paste("Time taken for year-level model:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))
try(stopCluster(cl), silent = TRUE)

table<-dctable(model.flush_yearfe_dipt)

windows()
plot(table, 1:8, type="log.var")
diag_test(dcdiag(model.flush_yearfe_dipt))
summary(model.flush_yearfe_dipt)
