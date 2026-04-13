library(pacman)
p_load(dclone, MASS, ggplot2, snow, tidyverse, parallel, lubridate)

# linear response prediction function
logit.pf <- function(kd,Td,x){
  out <- kd*(x-Td)
  return(out)
}
# JAGS model for intercepts
leaves <- function(){
  kd~ dunif(0,2)
  Td~ dunif(61,210)
  sigsq~ dunif(10,25)
  for(j in 1:n){
    muf[j] <-  kd*(days[j]-Td)
  }
  for(k in 1:K){
    for(i in 1:n){
      X[i,k] ~ dnorm(muf[i],1/sigsq)
    } 
  }
}
# Function to generate sampling dates with variability
generate_sampling_dates <- function(min_date, n_years, interval_days, interval_var) {
  max_date <- min_date + years(n_years) - days(1)
  n_approx <- ceiling(n_years * 365 / interval_days * 1.5)
  intervals <- pmax(1, rnorm(n_approx, interval_days, interval_var))
  dates <- c(min_date, min_date + days(cumsum(round(intervals))))
  dates[dates <= max_date]
}
# Function to simulate phenology data with interannual and intraspecific variability
simulate_phenology_data <- function(sampling_dates, n_individuals,
                                    Td, kd,
                                    interannual_var_td, interannual_var_kd,
                                    intraspecific_var_td, intraspecific_var_kd,
                                    noise) {
  n_years  <- length(unique(year(sampling_dates)))
  all.days <- rep(sampling_dates, n_individuals)
  cal_year <- year(all.days)
  year_idx <- as.numeric(as.factor(cal_year))
  indv.id  <- rep(1:n_individuals, each = length(sampling_dates))

  uTd <- rnorm(n_years, 0, interannual_var_td)[year_idx]
  iTd <- rnorm(n_individuals, 0, intraspecific_var_td)[indv.id]
  uKd <- rnorm(n_years, 0, interannual_var_kd)[year_idx]
  iKd <- rnorm(n_individuals, 0, intraspecific_var_kd)[indv.id]
  yiTd <- Td + uTd + iTd
  yiKd <- kd + uKd + iKd

  doy        <- yday(all.days)
  norm.samps <- rnorm(length(all.days), logit.pf(yiKd, yiTd, doy), sqrt(noise))

  data.frame(
    days = all.days, doy, indv = indv.id, year = cal_year,
    indv_year = as.factor(paste0(indv.id, "_", cal_year)),
    y = 1 / (1 + exp(norm.samps)), samps = norm.samps,
    yiTd, uTd, iTd, yiKd, uKd, iKd
  )
}

# -------------------------
# SETTINGS
# -------------------------
n.years <- 10
n.individuals <- 5
sigsq <- 0
kd <- 0.1
Td <- 120

inteval_days <- 30
interval_var <- 2
min_date <- as.Date("2018-01-01")

vector_inter_var <- c(0, 5, 10, 20, 40)
vector_intra_var <- c(0, 5, 10, 20, 40)
n.sims <- 3

all.days <- generate_sampling_dates(
  min_date,
  n_years = n.years,
  interval_days = inteval_days,
  interval_var = interval_var
)

matrix_sigsq <- matrix(NA,
  nrow = length(vector_inter_var),
  ncol = length(vector_intra_var)
)
rownames(matrix_sigsq) <- paste0("InterVar_", vector_inter_var)
colnames(matrix_sigsq) <- paste0("IntraVar_", vector_intra_var)

cl.seq <- c(1,5,10)
n.iter <- 3000
n.adapt <- 500
n.update <- 1000
thin <- 1
n.chains <- 3
cl <- makePSOCKcluster(3)

combo_grid <- expand.grid(inter_var = vector_inter_var, intra_var = vector_intra_var)
combo_grid <- combo_grid[sample(nrow(combo_grid)), ]
for (i in seq_len(nrow(combo_grid))) {
  inter_var <- combo_grid$inter_var[i]
  intra_var <- combo_grid$intra_var[i]
  cat("Inter:", inter_var, "Intra:", intra_var, "\n")
  tile_sigsq <- numeric(n.sims)
  for (sim in 1:n.sims) {
      iter_start <- proc.time()
      sim_data <- simulate_phenology_data(
        sampling_dates = all.days,
        n_individuals = n.individuals,
        Td = Td,
        kd = kd,
        interannual_var_td = inter_var,
        interannual_var_kd = 0,
        intraspecific_var_td = intra_var,
        intraspecific_var_kd = 0,
        noise = sigsq
      )

      data4dclone <- list(
        K = 1,
        X = dcdim(data.matrix(sim_data$samps)),
        n = nrow(sim_data),
        days = sim_data$doy
      )

      inits <- list(
        list(kd=runif(1,0,15), Td=runif(1,1,365), sigsq=runif(1,0.01,30)),
        list(kd=runif(1,0,15), Td=runif(1,1,365), sigsq=runif(1,0.01,30)),
        list(kd=runif(1,0,15), Td=runif(1,1,365), sigsq=runif(1,0.01,30))
      )

      leaves.dclone <- dc.parfit(
        cl,
        data4dclone,
        params = c("kd","Td","sigsq"),
        model = leaves,
        n.clones = cl.seq,
        multiply = "K",
        unchanged = "n",
        n.chains = n.chains,
        n.adapt = n.adapt,
        n.update = n.update,
        n.iter = n.iter,
        thin = thin,
        inits = inits
      )

      print(summary(leaves.dclone))
      print(dcdiag(leaves.dclone))
      summ <- summary(leaves.dclone)$statistics['sigsq','Mean']
      tile_sigsq[sim] <- summ
      iter_elapsed <- proc.time() - iter_start
      cat(sprintf("  [sim %d/%d] Inter: %g  Intra: %g  -> %.1f sec elapsed\n",
                  sim, n.sims, inter_var, intra_var, iter_elapsed["elapsed"]))
    }

    matrix_sigsq[
      paste0("InterVar_", inter_var),
      paste0("IntraVar_", intra_var)
    ] <- mean(tile_sigsq)

    print(matrix_sigsq)
}
stopCluster(cl)

df <- as.data.frame(matrix_sigsq) %>%
  rownames_to_column("InterVar") %>%
  pivot_longer(-InterVar, names_to = "IntraVar", values_to = "sigsq")%>%
  mutate(InterVar = as.numeric(gsub("InterVar_", "", InterVar)),
         IntraVar = as.numeric(gsub("IntraVar_", "", IntraVar)))

windows()
ggplot(df,
       aes(x=factor(IntraVar), y=factor(InterVar), fill=sigsq)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=round(sigsq, 2)), color="white", size=5, fontface="bold") +
  scale_fill_gradient(low="lightblue", high="blue") +
  labs(title="Estimated Observational Variance across Intraspecific and Interannual Variability",
       x="Intraspecific variability in Td (days)",
       y="Interannual variability in Td (days)",
       fill="Observational\nVariance") +
  theme_minimal()
  
#########################################################################
# Try with real data
########################################################################
data<- read.csv("cavallinesia_leafing_timeseries.csv")
head(data)

data <- data %>%
  mutate(
    y_norm= pmin(pmax(leafing / 100, 1e-4), 1 - 1e-4),
    date = as.Date(date),
    date_num = as.numeric(difftime(date, as.Date("2018-04-04"), units = "days")),
    DOY= yday(date),
    year= year(date),
    month= month(date),
    pheno_year = if_else(month >= 9, year, year - 1),
    day = as.numeric(difftime(date, as.Date(paste0(pheno_year, "-09-01")), units = "days")),
    tree= as.factor(tag),
    pheno_year= as.factor(pheno_year),
    tree_year= as.factor(paste0(tree, "_", pheno_year))
  )

trees1<- unique(data$tree)
years1<- unique(data$pheno_year)
all_before_threshold <- data.frame()

for (i in 1:length(trees1)) {
  for (j in 1:length(years1)) {
    subset_data <- data %>% filter(tree == trees1[i], pheno_year == years1[j])
    subset_data<- subset_data %>% arrange(day)
    cat("Processing Tree:", trees1[i], "Year:", print(years1[j]), "N rows:", nrow(subset_data), "\n")
    if (nrow(subset_data) == 0) {
      cat("  → No data for this tree-year combination\n\n")
      next
    }

    threshold <- 0.1
    below_threshold <- which(subset_data$y_norm < threshold)
    
    cat("  → Values below threshold:", length(below_threshold), "\n")
    cat("  → Range of y_norm:", min(subset_data$y_norm, na.rm=T), "to", max(subset_data$y_norm, na.rm=T), "\n")
    
    if (length(below_threshold) == 0) {
      cat("  → No values below threshold; skipping\n\n")
      print(below_threshold)
      next
    }
    
    last_item <- max(below_threshold, na.rm=TRUE)
    cat("  → Last below-threshold index:", last_item, "of", nrow(subset_data), "\n")
    
    #subset all values before last_item
    if (last_item > 0 && last_item <= nrow(subset_data)) {
      subset_before_threshold <- subset_data[1:last_item, ]
      all_before_threshold <- bind_rows(all_before_threshold, subset_before_threshold)
      cat("  → Added", nrow(subset_before_threshold), "rows\n\n")
    } else {
      cat("  → Last below-threshold item at end or at start; skipping\n\n")
      print(last_item)
      print(nrow(subset_data))
      print(last_item > 0 && last_item < nrow(subset_data))
    }
  }
}

windows()
ggplot(all_before_threshold, aes(x=day, y=y_norm, group=tree_year, color=as.factor(pheno_year))) +
  geom_line() +
  labs(title="Cavallinesia phenology data",
       y="Predicted leafing",
       x="Day of year",
       color="Phenological Year") +
  theme_minimal()


cl.seq <- c(1,2)
n.iter <- 10000
n.adapt <- 2000
n.update <- 5000
thin <- 5
n.chains <- 3
cl <- makePSOCKcluster(3)

test.data <- log((1 - all_before_threshold$y_norm) / all_before_threshold$y_norm)
data4dclone <- list(
  K = 1,
  X = dcdim(data.matrix(test.data)),
  n = nrow(all_before_threshold),
  days = all_before_threshold$day
)
inits <- list(
  list(kd=runif(1,0,2), Td=runif(1,61,210), sigsq=runif(1,10,25)),
  list(kd=runif(1,0,2), Td=runif(1,61,210), sigsq=runif(1,10,25)),
  list(kd=runif(1,0,2), Td=runif(1,61,210), sigsq=runif(1,10,25))
)

leaves.dclone <- dc.parfit(
  cl,
  data4dclone,
  params = c("kd","Td","sigsq"),
  model = leaves,
  n.clones = cl.seq,
  multiply = "K",
  unchanged = "n",
  n.chains = n.chains,
  n.adapt = n.adapt,
  n.update = n.update,
  n.iter = n.iter,
  thin = thin,
  inits = inits
)
dctable(leaves.dclone)

results <- summary(leaves.dclone)
table_summary <- results$statistics
windows()
grid.table(round(as.data.frame(table_summary), 4))
dcdiag(leaves.dclone)

table <- dcdiag(leaves.dclone)
windows()
grid.table(round(as.data.frame(table), 4))
# -------------------------
#plot


results$statistics[, "Mean"]
x_range<- seq(1,365, by=1)
y_predicted <- 1 / (1 + exp(results$statistics['kd','Mean'] * (x_range - results$statistics['Td','Mean'])))
df_pred <- data.frame(day = x_range, predicted_leafing = y_predicted)

day_to_DOY <- function(day, pheno_year) {
  format(as.Date(paste0(as.numeric(as.character(pheno_year)), "-09-01")) + day, "%B %d")
}
DOY_to_day <- function(month, day, pheno_year) {
  cal_year <- if (month >= 9) pheno_year else pheno_year + 1
  as.numeric(as.Date(paste0(cal_year, "-", month, "-", day)) -
               as.Date(paste0(pheno_year, "-09-01")))
}

SOF <- min(df_pred$day[df_pred$predicted_leafing < 0.9])  #day in which predicted leafing falls below 0.9
MOF <- min(df_pred$day[df_pred$predicted_leafing < 0.5])  #day in which predicted leafing falls below 0.5
EOF <- min(df_pred$day[df_pred$predicted_leafing < 0.1])  #day in which predicted leafing falls below 0.11

day_to_DOY(SOF, 2020)
day_to_DOY(MOF, 2020)
day_to_DOY(EOF, 2020)

windows()
ggplot(all_before_threshold, aes(x=day, y=y_norm, color=as.factor(tree_year))) +
  geom_point() +
  geom_line(data=df_pred, aes(x=day, y=predicted_leafing), color="black", size=1) +
  geom_vline(xintercept = SOF, linetype="dashed", color="black", size=0.5) +
  geom_vline(xintercept = MOF, linetype="dashed", color="black", size=0.5) +
  geom_vline(xintercept = EOF, linetype="dashed", color="black", size=0.5) +
  scale_x_continuous(
    breaks = {
      ref <- as.Date("2020-09-01")
      as.numeric(seq(ref, ref + 365, by = "month") - ref)
    },
    labels = format(seq(as.Date("2020-09-01"),
                        as.Date("2020-09-01") + 365,
                        by = "month"), "%b")
  ) +
  labs(title="Cavallinesia platanifolia: Global intercept model fit",
       y="Predicted leafing",
       x="Day of year") +
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("plots/cavallinesia_phenology_fit.png", width=10, height=6)


## the real thing happens between november and may as described by condit. we buffer by one month in each side. so my new prior should be oct 1 to june 30. 
DOY_to_day(11, 1, 2020) # -30 november 1 
DOY_to_day(3, 30, 2020) # 303  # we will run with an appropiate priorm  march 30
