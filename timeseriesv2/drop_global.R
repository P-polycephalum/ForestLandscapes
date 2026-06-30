library(pacman)
p_load(dclone, MASS, ggplot2, snow, tidyverse, parallel, gridExtra)

# Read and prepare data
# 8 for drop in cava, 10 for flush in cava
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
    # cut_idx <- min_indices[1]
    cut_idx_end <- min_indices[length(min_indices)]
    seg <- subset_data[1:cut_idx_end, ]
    #seg_after <- subset_data[cut_idx:nrow(subset_data), ]
    
    #loop for flush
    n_obs   <- nrow(seg)
    y_range <- max(seg$y_norm, na.rm = TRUE) - min(seg$y_norm, na.rm = TRUE)

    if (n_obs < 2 || y_range < 0.5) {
      print(paste0("SKIPPED tree=", trees1[i], " year=", years1[j],
                   " n=", n_obs, " range=", round(y_range, 3)))
      print(min(seg$y_norm, na.rm = TRUE))
      print(max(seg$y_norm, na.rm = TRUE))
      next
    }
    print(paste0("tree=", trees1[i], " year=", years1[j],
                 " n=", n_obs, " range=", round(y_range, 3)))
    all_before_threshold <- bind_rows(all_before_threshold, seg)
  }
}

windows()
ggplot(all_before_threshold, aes(x = day, y = y_norm, color=tree_year)) +
  geom_jitter()+
  theme(legend.position = "none")
drop_seg <- all_before_threshold %>% 
  filter(!tree %in% c("3811","4250")) %>%
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
ggplot(drop_seg, aes(x = day, y = y_norm, color=tree_year)) +
  geom_jitter()+
  facet_wrap(~pheno_year)+
  theme(legend.position = "none")

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

cl.seq <- c(1, 4, 8, 16, 32)
n.iter <- 10000
n.update <- 5000
n.adapt <- 2000
n.chains <- 3
thin <- 5
#first model, drop population level for cava
start_time <- Sys.time()
leaves.drop_indv_fe <- function(){
  kd    ~ dunif(0.01, 1)
  base  ~ dunif(0, 1)
  amp_raw ~ dunif(0, 1)
  sigsq ~ dunif(1, 25)
  Td ~ dunif(50, 250)

  amp <- amp_raw * (1 - base)

  for (j in 1:n) {
    eta[j] <-  kd * (days[j] - Td)
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
  Td = runif(1, 50, 250),
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
  model = leaves.drop_indv_fe,
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
summary(model.dropt)
diag_test(dcdiag(model.drop))
saveRDS(model.drop, "timeseriesv2//data//cava_drop.rds")
end_time <- Sys.time()
print(paste("Time taken for individual-level model:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))


