library(pacman)
p_load(dclone, MASS, ggplot2, snow, tidyverse, parallel)
logit.pf <- function(kd,Td,x){
  out <- kd*(x-Td)
  return(out)
}

prepare_real_phenology_data <- function(file_path) {
  raw_data <- read.csv(file_path)

  raw_data %>%
    mutate(
      y_norm = pmin(pmax(leafing / 100, 1e-4), 1 - 1e-4),
      date = as.Date(date),
      date_num = as.numeric(difftime(date, as.Date("2018-04-04"), units = "days")),
      DOY = yday(date),
      year = year(date),
      month = month(date),
      pheno_year = if_else(month >= 9, year, year - 1),
      day = as.numeric(difftime(date, as.Date(paste0(pheno_year, "-09-01")), units = "days")),
      tree = as.factor(tag),
      pheno_year = as.factor(pheno_year),
      tree_year = as.factor(paste0(tree, "_", pheno_year))
    )
}

keep_before_last_threshold <- function(subset_data, threshold = 0.1, verbose = TRUE) {
  subset_data <- subset_data %>% arrange(day)

  if (nrow(subset_data) == 0) {
    return(subset_data)
  }

  below_threshold <- which(subset_data$y_norm < threshold)

  if (verbose) {
    cat(
      "Processing Tree:", as.character(subset_data$tree[1]),
      "Year:", as.character(subset_data$pheno_year[1]),
      "N rows:", nrow(subset_data), "\n"
    )
    cat("  → Values below threshold:", length(below_threshold), "\n")
    cat("  → Range of y_norm:", min(subset_data$y_norm, na.rm = TRUE), "to", max(subset_data$y_norm, na.rm = TRUE), "\n")
  }

  if (length(below_threshold) == 0) {
    if (verbose) cat("  → No values below threshold; skipping\n\n")
    return(subset_data[0, ])
  }

  last_item <- max(below_threshold, na.rm = TRUE)

  if (verbose) {
    cat("  → Last below-threshold index:", last_item, "of", nrow(subset_data), "\n")
  }

  if (last_item > 0 && last_item <= nrow(subset_data)) {
    subset_before_threshold <- subset_data[1:last_item, ]
    if (verbose) cat("  → Added", nrow(subset_before_threshold), "rows\n\n")
    return(subset_before_threshold)
  }

  if (verbose) cat("  → Invalid last index; skipping\n\n")
  subset_data[0, ]
}

extract_pre_threshold_data <- function(data, threshold = 0.1, verbose = TRUE) {
  groups <- split(data, list(data$tree, data$pheno_year), drop = TRUE)
  purrr::map_dfr(groups, keep_before_last_threshold, threshold = threshold, verbose = verbose)
}

plot_phenology_series <- function(data) {
  windows()
  ggplot(data, aes(x = day, y = y_norm, color = as.factor(tree_year))) +
    geom_line() +
    labs(
      title = "Cavallinesia phenology data",
      y = "Normalized leafing",
      x = "Date"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

plot_threshold_subset <- function(all_before_threshold) {
  windows()
  ggplot(all_before_threshold, aes(x = day, y = y_norm, group = tree_year, color = as.factor(pheno_year))) +
    geom_line() +
    labs(
      title = "Cavallinesia phenology data",
      y = "Predicted leafing",
      x = "Day of year",
      color = "Phenological Year"
    ) +
    theme_minimal()
}


generate_sampling_dates <- function(min_date, n_years, interval_days, interval_var) {
  # Generate sampling dates with specified sampling interval and variability.
  # min_date: start date.
  # n_years: number of years to simulate.
  # interval_days: average days between samples.
  # interval_var: standard deviation of sampling interval.
  max_date <- min_date + years(n_years) - days(1)
  print(paste("Sampling dates from", min_date, "to", max_date))

  sampling_dates <- min_date
  current_date <- min_date

  while (current_date < max_date) {
    interval <- rnorm(1, mean = interval_days, sd = interval_var)
    interval <- max(1, interval)
    current_date <- current_date + days(round(interval))

    if (current_date <= max_date) {
      sampling_dates <- c(sampling_dates, current_date)
    }
  }

  print(paste("Number of sampling dates:", length(sampling_dates)))
  sampling_dates
}

simulate_phenology_data <- function(sampling_dates,
                                    n_years,
                                    n_individuals,
                                    Td,
                                    kd,
                                    interannual_var_td,
                                    interannual_var_kd,
                                    intraspecific_var_td,
                                    intraspecific_var_kd,
                                    noise){
  # Simulate phenology data with interannual and intraspecific variability.
  all.days <- rep(sampling_dates, n_individuals)
  calendar_years <- year(all.days)
  year_indices <- as.numeric(as.factor(calendar_years))
  indv.id <- rep(1:n_individuals, each = length(sampling_dates))

  uTd <- rnorm(n = n_years, mean = 0, sd = interannual_var_td)
  iTd <- rnorm(n = n_individuals, mean = 0, sd = intraspecific_var_td)
  yiTd <- Td + uTd[year_indices] + iTd[indv.id]

  uKd <- rnorm(n = n_years, mean = 0, sd = interannual_var_kd)
  iKd <- rnorm(n = n_individuals, mean = 0, sd = intraspecific_var_kd)
  yiKd <- kd + uKd[year_indices] + iKd[indv.id]

  sampling_dates_doy <- yday(all.days)
  mu.true <- logit.pf(kd = yiKd, Td = yiTd, x = sampling_dates_doy)
  norm.samps <- rnorm(n = length(all.days), mean = mu.true, sd = sqrt(noise))
  y.sims <- 1 / (1 + exp(norm.samps))

  sim_df <- data.frame(
      days = all.days,
      indv = indv.id,
      year = calendar_years,
      y = y.sims,
      samps = norm.samps,
      yiTd = yiTd,
      uTd = uTd[year_indices],
      iTd = iTd[indv.id],
      yiKd = yiKd,
      uKd = uKd[year_indices],
      iKd = iKd[indv.id]
    ) %>%
      mutate(
        indv_year = as.factor(paste0(indv, "_", year)),
        doy = yday(days)
      )

  return(sim_df)
}

plot_simulated_phenology <- function(df, n_years) {
  windows()
  ggplot(df, aes(x = doy, y = y, group = as.factor(indv_year))) +
    geom_line(aes(color = as.factor(year))) +
    labs(
      title = paste("Simulated phenology data with", n_years, "years"),
      y = "Simulated y",
      x = "Day of year"
    ) +
    theme_minimal()
}

summarize_effects <- function(model_fit, param_name, param_levels = NULL, param_label = NULL) {
  # Extract and format fixed or random effects from JAGS model fit
  # param_name: parameter pattern to search for (e.g., "yTd", "iTd")
  # param_levels: optional labels for each effect level
  # param_label: descriptive label for output (e.g., "Phenological Year", "Individual")
  
  if (is.null(param_label)) {
    param_label <- param_name
  }
  
  fit_summary <- summary(model_fit)
  
  param_rows <- grep(param_name, rownames(fit_summary$statistics))
  param_means <- fit_summary$statistics[param_rows, "Mean"]
  param_sds <- fit_summary$statistics[param_rows, "SD"]
  
  cat("\n========================================\n")
  cat("Empirical", param_label, "Effects (", param_name, ")\n")
  cat("========================================\n\n")
  
  effect_table <- data.frame(
    index = seq_along(param_means),
    mean = round(param_means, 2),
    sd = round(param_sds, 2)
  )
  colnames(effect_table) <- c("index", paste0(param_name, "_mean"), paste0(param_name, "_sd"))
  
  if (!is.null(param_levels)) {
    effect_table <- cbind(
      level = param_levels[seq_along(param_means)],
      effect_table
    )
  }
  
  print(effect_table)
  
  cat("\n--- Summary Statistics ---\n")
  cat("Mean of", param_name, "estimates:", round(mean(param_means), 2), "\n")
  cat("SD of", param_name, "estimates:  ", round(sd(param_means), 2), "\n")
  cat("Min", param_name, ":              ", round(min(param_means), 2), "\n")
  cat("Max", param_name, ":              ", round(max(param_means), 2), "\n\n")
  
  invisible(list(means = param_means, sds = param_sds, table = effect_table))
}



#one simple model like this one with fixed effects is misspecified
#because it doesnt accounts for variability in the Td given the individual
leaves <- function(){
  kd~ dunif(0,2)
  Td~ dunif(1,365)
  sigsq~ dunif(0.01,100)
  for(j in 1:n){
    muf[j] <-  kd*(days[j]-Td)
  }
  for(k in 1:K){
    for(i in 1:n){
      X[i,k] ~ dnorm(muf[i],1/sigsq)
    } 
  }
}

leaves_year_fe_td <- function(){
  kd    ~ dunif(0, 2)
  sigsq ~ dunif(0.01, 100)

  for (y in 1:nyear) {
    yTd[y] ~ dunif(1, 365)
  }

  for (j in 1:n) {
    muf[j] <- kd * (days[j] - yTd[year[j]])
  }

  for (k in 1:K) {
    for (i in 1:n) {
      X[i,k] ~ dnorm(muf[i], 1/sigsq)
    }
  }
}

leaves_fe_kdtd <- function(){
  sigsq ~ dunif(0.01, 100)

  for (y in 1:nyear) {
    yKd[y] ~ dunif(0.01, 2)
  }
  for (y in 1:nyear) {
    yTd[y] ~ dunif(1, 365)
  }

  for (j in 1:n) {
    muf[j] <- yKd[year[j]] * (days[j] - yTd[year[j]])
  }

  for (k in 1:K) {
    for (i in 1:n) {
      X[i,k] ~ dnorm(muf[i], 1/sigsq)
    }
  }
}

leaves_year_tdkd <- function(){
  sigsq ~ dunif(0.01, 100)

  for (y in 1:nyear) {
    yTd[y] ~ dunif(1, 365)
  }

  for (y in 1:nyear) {
    yKd[y] ~ dunif(0, 2)
  }

  for (j in 1:n) {
    muf[j] <- yKd[year[j]] * (days[j] - yTd[year[j]])
  }

  for (k in 1:K) {
    for (i in 1:n) {
      X[i,k] ~ dnorm(muf[i], 1/sigsq)
    }
  }
}

# what if i did the opposite and monitor the individuals instead of the years? 
leaves_indv_fe <- function(){
  kd    ~ dunif(0, 2)
  sigsq ~ dunif(0.01, 100)

  for (i in 1:nindv) {
    iTd[i] ~ dunif(1, 365)
  }

  for (j in 1:n) {
    muf[j] <- kd * (days[j] - iTd[indv[j]])
  }

  for (k in 1:K) {
    for (i in 1:n) {
      X[i,k] ~ dnorm(muf[i], 1/sigsq)
    }
  }
}


leaves_re_both <- function(){

  # Fixed parameters
  kd    ~ dunif(0, 2)
  sigsq ~ dunif(0.01, 100)

  Td0 ~ dunif(1, 365)

  # Variance components
  sigma_year ~ dunif(0, 50)
  sigma_indv ~ dunif(0, 50)

  tau_year <- 1 / pow(sigma_year, 2)
  tau_indv <- 1 / pow(sigma_indv, 2)

  # Year random effects
  for (y in 1:nyear){
    u_raw[y] ~ dnorm(0, tau_year)
  }

  u_bar <- mean(u_raw[])

  for (y in 1:nyear){
    u[y] <- u_raw[y] - u_bar
  }

  # Individual random effects
  for (i in 1:nindv){
    v_raw[i] ~ dnorm(0, tau_indv)
  }

  v_bar <- mean(v_raw[])
  for (i in 1:nindv){
    v[i] <- v_raw[i] - v_bar
  }

  # Phenology model
  for (j in 1:n){
    muf[j] <- kd * (days[j] - (Td0 + u[year[j]] + v[indv[j]]))
  }


  for (k in 1:K){
    for (j in 1:n){
      X[j,k] ~ dnorm(muf[j], 1/sigsq)
    }
  }
}

leaves_fe <- function() {
  kd    ~ dunif(0, 2)
  sigsq ~ dunif(0.01, 100)

  for (yy in 1:nyear) {
    for (ii in 1:nindv) {
      yiTd[yy, ii] ~ dunif(1, 365)
    }
  }

  for (j in 1:n) {
    muf[j] <- kd * (days[j] - yiTd[year[j], indv[j]])
  }

  for (k in 1:K) {
    for (j in 1:n) {
      X[j, k] ~ dnorm(muf[j], 1 / sigsq)
    }
  }
}


############################################################################
# Simulate data for intraspecific variability model
##################################################################
min_date <- as.Date("2018-01-01")
interval_days <- 30
interval_var <- 2
n.years <- 4
interannual_var_td <- 5
interannual_var_kd<- 0.01
intraspecific_var_td <- 2
intraspecific_var_kd<- 0.01
noise <- 2
Td <- 120
kd <- 0.1
n.individuals <-4

sampling_dates <- generate_sampling_dates(
  min_date = min_date,
  n_years = n.years,
  interval_days = interval_days,
  interval_var = interval_var
)

df <- simulate_phenology_data(
  sampling_dates = sampling_dates,
  n_years = n.years,
  n_individuals = n.individuals,
  Td = Td,
  kd = kd,
  interannual_var_td = interannual_var_td,
  intraspecific_var_kd = intraspecific_var_kd,
  intraspecific_var_td = intraspecific_var_td,
  interannual_var_kd = interannual_var_kd,
  noise = noise
)
View(df)
plot_simulated_phenology(df, n.years)

#input data
data <- prepare_real_phenology_data("cavallinesia_leafing_timeseries.csv")
plot_phenology_series(data)
all_before_threshold <- extract_pre_threshold_data(
  data = data,
  threshold = 0.1,
  verbose = TRUE
)
all_before_threshold <- all_before_threshold %>%
filter(!pheno_year %in% c("2017"))
plot_threshold_subset(all_before_threshold)
View(all_before_threshold)
#extract all the vectors necessary
# Re-factor to ensure consecutive indices 1:n after filtering
all_before_threshold <- all_before_threshold %>%
  mutate(
    pheno_year = droplevels(as.factor(pheno_year)),
    tree = droplevels(as.factor(tree))
  )

test.data<- log(1-all_before_threshold$y_norm) - log(all_before_threshold$y_norm)
n<- nrow(all_before_threshold)
all.days<- all_before_threshold$day
year.id<- as.numeric(all_before_threshold$pheno_year)
indv.id <- as.numeric(all_before_threshold$tree)
n.indv<- length(unique(indv.id))
n.year<- length(unique(year.id))

#################################################################################################
#Jags PARAMETERS
#################################################################################################
cl.seq <- c(1,4,8,16);
n.iter<-5000;n.adapt<-500;n.update<-200;n.thin<-5;n.chains<-3;


###################################################################################################
#Run the global intercept model with the real data
###################################################################################################

data4dclone_intercept <- list(K=1, X=dcdim(data.matrix(test.data)), n=n, days=all.days)
data4dclone_sim <- list(K=1, X=dcdim(data.matrix(df$samps)), n=nrow(df), days=yday(df$days))
cl<- makePSOCKcluster(3)
out.parms <- c("kd", "Td", "sigsq")
inits.list <- list(
    list(kd=runif(1, min=0, max=2), Td=runif(1, min=0, max=365), sigsq=runif(1, min=0.01, max=100)),
    list(kd=runif(1, min=0, max=2), Td=runif(1, min=0, max=365), sigsq=runif(1, min=0.01, max=100)),
    list(kd=runif(1, min=0, max=2), Td=runif(1, min=0, max=365), sigsq=runif(1, min=0.01, max=100))
)
leaves.intercept <- dc.parfit(cl,data4dclone_sim, params=out.parms, model=leaves, n.clones=cl.seq,
                        multiply="K",unchanged="n",
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=n.thin,
                        inits=inits.list
                        )

summary(leaves.intercept)
dcdiag(leaves.intercept)

table<-dctable(leaves.intercept)
windows()
plot(table, 1:3, type='log.var')

library(coda)

# Combine all chains
post <- as.matrix(leaves.intercept)

# Extract parameters
Td <- post[,"Td"]
kd <- post[,"kd"]
sigsq <- post[,"sigsq"]

# Plot
windows()
hist(Td,
     breaks=40,
     col="gray",
     border="white",
     main="Posterior distribution of Td",
     xlab="Td")
     xlab="sigsq")
windows()
plot(Td, kd,
     pch=16,
     col=rgb(0,0,0,0.2),
     xlab="Td",
     ylab="kd",
     main="Posterior samples: Td vs kd")
#######################################################################################
#run the leaves year fe model with the real data
#######################################################################################

data4dclone_year_fe <- list(K=1, X=dcdim(data.matrix(test.data)), n=n, days=all.days, year=year.id, nyear=n.year)
data4dclone_year_sim<- list(K=1, X=dcdim(data.matrix(df$samps)), n=nrow(df), days=yday(df$days), year=as.numeric(as.factor(df$year)), nyear=length(unique(df$year)))
cl<- makePSOCKcluster(3)
out.parms_year_fe <- c("kd", "sigsq", "yTd")
inits.list_year_fe <- list(
    list(kd=0.1, sigsq=1, yTd=runif(n.year, min=1, max=365)),
    list(kd=0.3, sigsq=10, yTd=runif(n.year, min=1, max=365)),
    list(kd=0.5, sigsq=30, yTd=runif(n.year, min=1, max=365))
)
leaves.year_fe <- dc.parfit(cl,data4dclone_year_fe, params=out.parms_year_fe, model=leaves_year_fe, n.clones=cl.seq,
                        multiply="K",unchanged=c("n","nyear"),
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=n.thin,
                        inits=inits.list_year_fe
                        )

summary(leaves.year_fe)
year_levels <- levels(as.factor(all_before_threshold$pheno_year))
year_effects <- summarize_effects(leaves.year_fe, "yTd", param_levels = year_levels, param_label = "Phenological Year")

#######################################################################################
# MODEl: Kd by year, Td by year fixed effects
#######################################################################################

data4dclone_year_fe <- list(K=1, X=dcdim(data.matrix(test.data)), n=n, days=all.days, year=year.id, nyear=n.year)
data4dclone_year_sim<- list(K=1, X=dcdim(data.matrix(df$samps)), n=nrow(df), days=yday(df$days), year=as.numeric(as.factor(df$year)), nyear=length(unique(df$year)))
cl<- makePSOCKcluster(3)
out.parms_year_fe <- c("kd", "sigsq", "yTd", "yKd")
inits.list_year_fe <- list(
    list(kd=0.1, sigsq=1, yTd=runif(n.year, min=1, max=365), yKd=runif(n.year, min=0.01, max=2)),
    list(kd=0.3, sigsq=10, yTd=runif(n.year, min=1, max=365), yKd=runif(n.year, min=0.01, max=2)),
    list(kd=0.5, sigsq=30, yTd=runif(n.year, min=1, max=365), yKd=runif(n.year, min=0.01, max=2))
)
leaves.year_kdtd_fe <- dc.parfit(cl,data4dclone_year_fe, params=out.parms_year_fe, model=leaves_fe_kdtd, n.clones=cl.seq,
                        multiply="K",unchanged=c("n","nyear"),
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=n.thin,
                        inits=inits.list_year_fe
                        )

summary(leaves.year_kdtd_fe)
dcdiag(leaves.year_kdtd_fe)


#######################################################################################
#MODEL: Kd by year, Td Fixed
#######################################################################################

data4dclone_year_fe <- list(K=1, X=dcdim(data.matrix(test.data)), n=n, days=all.days, year=year.id, nyear=n.year)
data4dclone_year_sim<- list(K=1, X=dcdim(data.matrix(df$samps)), n=nrow(df), days=yday(df$days), year=as.numeric(as.factor(df$year)), nyear=length(unique(df$year)))
cl<- makePSOCKcluster(3)
out.parms_year_fe <- c("Td", "sigsq", "yKd")
inits.list_year_fe <- list(
    list(Td=50, sigsq=1, yKd=runif(length(unique(df$year)), min=0.01, max=2)),
    list(Td=150, sigsq=10, yKd=runif(length(unique(df$year)), min=0.01, max=2)),
    list(Td=250, sigsq=30, yKd=runif(length(unique(df$year)), min=0.01, max=2))
)
model_kdyear_tdfixed <- dc.parfit(cl,data4dclone_year_sim, params=out.parms_year_fe, model=leaves_year_fe_kd, n.clones=cl.seq,
                        multiply="K",unchanged=c("n","nyear"),
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=n.thin,
                        inits=inits.list_year_fe
                        )

summary(model_kdyear_tdfixed)
dcdiag(model_kdyear_tdfixed)
year_levels <- levels(as.factor(all_before_threshold$pheno_year))
year_effects <- summarize_effects(model_kdyear_tdfixed, "yKd", param_levels = year_levels, param_label = "Phenological Year")
#######################################################################################
#run the Td-Kd year fe model with the real data
#######################################################################################
data4dclone_year_tdkd <- list(K=1, X=dcdim(data.matrix(test.data)), n=n, days=all.days, year=year.id, nyear=n.year)
cl<- makePSOCKcluster(3)
out.parms_year_tdkd <- c("sigsq", "yTd", "yKd")
inits.list_year_tdkd <- list(
    list(sigsq=1, yTd=runif(n.year, min=1, max=365), yKd=runif(n.year, min=0, max=2)),
    list(sigsq=10, yTd=runif(n.year, min=1, max=365), yKd=runif(n.year, min=0, max=2)),
    list(sigsq=30, yTd=runif(n.year, min=1, max=365), yKd=runif(n.year, min=0, max=2))
)
leaves.year_tdkd <- dc.parfit(cl,data4dclone_year_tdkd, params=out.parms_year_tdkd, model=leaves_year_tdkd, n.clones=cl.seq,
                        multiply="K",unchanged=c("n","nyear"),
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=n.thin,
                        inits=inits.list_year_tdkd
                        )
summary(leaves.year_tdkd)
dcdiag(leaves.year_tdkd)
################################################################################
#fit the leaves_indv_fe model with the real data
################################################################################

data4dclone_indv_fe <- list(K=1, X=dcdim(data.matrix(test.data)), n=n, days=all.days, indv=indv.id, nindv=n.indv)
cl<- makePSOCKcluster(3)
out.parms_indv_fe <- c("kd", "sigsq", "iTd")
inits.list_indv_fe <- list(
    list(kd=0.1, sigsq=1, iTd=runif(n.indv, min=1, max=365)),
    list(kd=0.3, sigsq=10, iTd=runif(n.indv, min=1, max=365)),
    list(kd=0.5, sigsq=30, iTd=runif(n.indv, min=1, max=365))
)

leaves.indv_fe <- dc.parfit(cl,data4dclone_indv_fe, params=out.parms_indv_fe, model=leaves_indv_fe, n.clones=cl.seq,
                        multiply="K",unchanged=c("n","nindv"),
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=n.thin,
                        inits=inits.list_indv_fe
                        )
summary(leaves.indv_fe)

indv_levels <- levels(as.factor(all_before_threshold$tree))
indv_effects <- summarize_effects(leaves.indv_fe, "iTd", param_levels = indv_levels, param_label = "Individual")


##########################################################################3
#run the leaves_re_both model with the simulated data
###########################################################################
head(df)
data4dclone_re_both <- list(K=1, X=dcdim(data.matrix(df$samps)), n=nrow(df), days=yday(df$days), indv=as.factor(df$indv), year=as.factor(df$year), nindv= length(unique(df$indv)), nyear=length(unique(df$year)))

cl<- makePSOCKcluster(3)
out.parms_re_both <- c("kd", "sigsq", "Td0", "sigma_year", "sigma_indv", "u", "v")

inits.list_re_both <- list(
  list(
    kd=0.10, sigsq=18, Td0=120,
    sigma_year=15, sigma_indv=8,
    u_raw=rnorm(n.year, mean=0, sd=15),
    v_raw=rnorm(n.indv, mean=0, sd=8)
  ),
  list(
    kd=0.12, sigsq=22, Td0=110,
    sigma_year=20, sigma_indv=10,
    u_raw=rnorm(n.year, mean=0, sd=20),
    v_raw=rnorm(n.indv, mean=0, sd=10)
  ),
  list(
    kd=0.08, sigsq=15, Td0=130,
    sigma_year=10, sigma_indv=6,
    u_raw=rnorm(n.year, mean=0, sd=10),
    v_raw=rnorm(n.indv, mean=0, sd=6)
  )
)

leaves.all <- dc.parfit(cl,data4dclone_re_both, params=out.parms_re_both, model=leaves_re_both, n.clones=cl.seq,
                        multiply="K",unchanged=c("n","nindv","nyear"),
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=n.thin,
                        inits=inits.list_re_both
                        )

summary(leaves.all)
dcdiag(leaves.all)

dctable_all <- dctable(leaves.all)
windows()
plot(dctable_all, 1:4)

#######################################################################################
#run fixed effects model with the simulated data
#######################################################################################
data4dclone_real<- list(K=1, X=dcdim(data.matrix(test.data)), n=n, days=all.days, indv=indv.id, year=year.id, nindv=n.indv, nyear=n.year)
data4dclone_fe <- list(K=1, X=dcdim(data.matrix(df$samps)), n=nrow(df), days=yday(df$days), indv=as.numeric(as.factor(df$indv)), year=as.numeric(as.factor(df$year)), nindv= length(unique(df$indv)), nyear=length(unique(df$year)))

cl<- makePSOCKcluster(3)
out.parms_fe <- c("kd", "sigsq", "yiTd")

inits <- list(
  list(kd=0.1, sigsq=20, yiTd=matrix(runif(n.year*n.indv, 1, 365), n.year, n.indv)),
  list(kd=0.2, sigsq=15, yiTd=matrix(runif(n.year*n.indv, 1, 365), n.year, n.indv)),
  list(kd=0.05, sigsq=30, yiTd=matrix(runif(n.year*n.indv, 1, 365), n.year, n.indv))
)

leaves.fe <- dc.parfit(cl,data4dclone_real, params=out.parms_fe, model=leaves_fe, n.clones=cl.seq,
                        multiply="K",unchanged=c("n","nindv","nyear"),
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=n.thin,
                        inits=inits
                        )

summary(leaves.fe)

dcdiag(leaves.fe)

# Extract year and individual effects from simulated data
yiTd_effects_sim <- df %>%
  distinct(year, indv, yiTd) %>%
  arrange(year, indv)

# Create year index mapping (calendar year to index)
year_mapping <- tibble(
  calendar_year = sort(unique(df$year)),
  year_index = seq_along(unique(df$year))
)

# Add year_index to simulated effects
yiTd_effects_sim <- yiTd_effects_sim %>%
  left_join(year_mapping, by = c("year" = "calendar_year"))

# Extract fitted yiTd effects from JAGS output
post <- summary(leaves.fe)$statistics
yiTd_rows <- grep("^yiTd\\[", rownames(post))
yiTd_fitted <- post[yiTd_rows, c("Mean", "SD")]
yiTd_fitted <- as.data.frame(yiTd_fitted)
yiTd_fitted$param <- rownames(yiTd_fitted)

# Parse year and indv indices from parameter names (e.g., "yiTd[1,3]" -> year=1, indv=3)
yiTd_fitted <- yiTd_fitted %>%
  mutate(
    year_index = as.numeric(sub("yiTd\\[([0-9]+),.*", "\\1", param)),
    indv = as.numeric(sub(".*,([0-9]+)\\]", "\\1", param)),
    fitted_mean = Mean,
    fitted_sd = SD
  ) %>%
  select(year_index, indv, fitted_mean, fitted_sd)

# Merge simulated with fitted effects for comparison
comparison <- yiTd_effects_sim %>%
  left_join(yiTd_fitted, by = c("year_index", "indv")) %>%
  mutate(
    bias = fitted_mean - yiTd,
    rel_error = abs(bias) / yiTd * 100
  )

cat("\n--- Bias Summary ---\n")
cat("Mean bias:", round(mean(comparison$bias, na.rm = TRUE), 2), "\n")
cat("SD bias:  ", round(sd(comparison$bias, na.rm = TRUE), 2), "\n")
cat("Mean relative error:", round(mean(comparison$rel_error, na.rm = TRUE), 2), "%\n\n")

head(comparison)
comparison %>%
  group_by(year_index) %>%
  summarise(
    mean_bias = mean(bias, na.rm = TRUE),
    mean_rel_error = mean(rel_error, na.rm = TRUE)
  ) %>%
  print()

df %>%
    group_by(year) %>%
    summarise(mean_yTd = mean(yiTd, na.rm = TRUE)) %>%
    print()  