library(pacman)
p_load(dclone, MASS, ggplot2, snow, tidyverse, parallel,coda)
logit.pf <- function(kd,Td,x){
  out <- kd*(x-Td)
  return(out)
}

leaves <- function(){
  kd~ dunif(0,15)
  Td~ dunif(1,365)
  sigsq~ dunif(0.01,30)
  for(j in 1:n){
    muf[j] <-  kd*(days[j]-Td)
  }
  for(k in 1:K){
    for(i in 1:n){
      X[i,k] ~ dnorm(muf[i],1/sigsq)
    } 
  }
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


############################################################################
# Simulate data for intraspecific variability model
##################################################################
min_date <- as.Date("2018-01-01")
interval_days <- 30
interval_var <- 2
interannual_var_td <- 20
interannual_var_kd<- 0
intraspecific_var_td <- 20
intraspecific_var_kd<- 0
noise <- 2
Td <- 120
kd <- 0.1

#vectors
n.years_vect <- c(2,6,10,14)
n.individuals_vect <- c(1,5,15,30)


n.chains <- 3;n.adapt <- 50;n.iter <- 1000;n.thin <- 1;n.update <- 100
#for identifiability we will evaluate
cov_matrix<- matrix(nrow=length(n.years_vect), ncol=length(n.individuals_vect))
colnames(cov_matrix) <- paste0("n.ind", n.individuals_vect)
rownames(cov_matrix) <- paste0("n.y_", n.years_vect)
lambdamax_matrix <- matrix(nrow=length(n.years_vect), ncol=length(n.individuals_vect))
colnames(lambdamax_matrix) <- paste0("n.ind", n.individuals_vect)
rownames(lambdamax_matrix) <- paste0("n.y_", n.years_vect)
td_spread_matrix <- matrix(nrow=length(n.years_vect), ncol=length(n.individuals_vect))
colnames(td_spread_matrix) <- paste0("n.ind", n.individuals_vect)
rownames(td_spread_matrix) <- paste0("n.y_", n.years_vect)
td_bias_matrix <- matrix(nrow=length(n.years_vect), ncol=length(n.individuals_vect))
colnames(td_bias_matrix) <- paste0("n.ind", n.individuals_vect)
rownames(td_bias_matrix) <- paste0("n.y_", n.years_vect)


for (n.years in n.years_vect) {
  for (n.individuals in n.individuals_vect) {
    print(paste("Simulating data for", n.years, "years and", n.individuals, "individuals"))
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

    # for estimability we estimate the expected value of yiTd
    n_tree_year_combos <- dplyr::n_distinct(df$indv, df$year)
    n_unique_yiTd <- dplyr::n_distinct(df$yiTd)
    cat("Tree-year combos:", n_tree_year_combos,
      "| Unique yiTd values:", n_unique_yiTd, "\n")

    emp_td <- mean(unique(df$yiTd))
    
    data4dclone_sim <- list(K=1, X=dcdim(data.matrix(df$samps)), n=nrow(df), days=yday(df$days))
    out.parms <- c("kd", "Td", "sigsq")
    cl.seq <- c(1,4,8,16);
    cl<- makePSOCKcluster(3)
    listinnits<-list(list(kd=runif(1, min=0, max=15), Td=runif(1, min=0, max=365), sigsq=runif(1, min=0.01, max=30)),
                     list(kd=runif(1, min=0, max=15), Td=runif(1, min=0, max=365), sigsq=runif(1, min=0.01, max=30)),
                     list(kd=runif(1, min=0, max=15), Td=runif(1, min=0, max=365), sigsq=runif(1, min=0.01, max=30))
    )
    leaves.sim <- dc.parfit(cl,data4dclone_sim, params=out.parms, model=leaves, n.clones=cl.seq,
                            multiply="K",unchanged="n",
                            n.chains = n.chains, 
                            n.adapt=n.adapt, 
                            n.iter=n.iter, 
                            n.thin=n.thin,
                            n.update=n.update,
                            inits = listinnits
    )
    stopCluster(cl)
    print(summary(leaves.sim))
    print(dcdiag(leaves.sim))
    
    
    diag_table<- data.frame(dcdiag(leaves.sim)) %>% mutate(klambdamax = n.clones * lambda.max)
    cov<- sd(diag_table$klambdamax)/mean(diag_table$klambdamax)
    cov_matrix[paste0("n.y_", n.years), paste0("n.ind", n.individuals)] <- cov
    print(paste("Coefficient of variation of klambdamax:", cov))
    print(cov_matrix)

    lambdamax_matrix[paste0("n.y_", n.years), paste0("n.ind", n.individuals)] <- diag_table$lambda.max[1]
    print(paste("Maximum eigenvalue of the variance-covariance matrix:", diag_table$lambda.max[1]))
    print(lambdamax_matrix)

    td_spread_in_days <- confint(leaves.sim)["Td","97.5 %"] - confint(leaves.sim)["Td","2.5 %"] 
    td_spread_matrix[paste0("n.y_", n.years), paste0("n.ind", n.individuals)] <- td_spread_in_days
    print(paste("Spread in Td (days):", td_spread_in_days))
    print(td_spread_matrix)

    td_bias_in_days <- abs(emp_td - mean(summary(leaves.sim)$statistics["Td","Mean"]))
    td_bias_matrix[paste0("n.y_", n.years), paste0("n.ind", n.individuals)] <- td_bias_in_days
    print(paste("Empirical Td:", emp_td))
    print(paste("Estimated Td:", mean(summary(leaves.sim)$statistics["Td","Mean"])))
    print(td_bias_matrix)
    
  }}

############################################################################
# demostrate sigsq inflation with increasing variability in the shifts
##################################################################
min_date <- as.Date("2018-01-01")
interval_days <- 30
interval_var <- 2

interannual_var_td_vector <- c(0,5,10,20)
interannual_var_kd<- 0
intraspecific_var_td_vector <- c(0,5,10,20)
intraspecific_var_kd<- 0
noise <- 0
Td <- 120
kd <- 0.1

n.years <- 5
n.individuals <- 5


#for identifiability we will evaluate
cov_matrix<- matrix(nrow=length(interannual_var_td_vector), ncol=length(intraspecific_var_td_vector))
colnames(cov_matrix) <- paste0("Y.var", intraspecific_var_td_vector)
rownames(cov_matrix) <- paste0("I.var", interannual_var_td_vector)

lambdamax_matrix <- matrix(nrow=length(interannual_var_td_vector), ncol=length(intraspecific_var_td_vector))
colnames(lambdamax_matrix) <- paste0("Y.var", intraspecific_var_td_vector)
rownames(lambdamax_matrix) <- paste0("I.var", interannual_var_td_vector)

sigsq_matrix <- matrix(nrow=length(interannual_var_td_vector), ncol=length(intraspecific_var_td_vector))
colnames(sigsq_matrix) <- paste0("Y.var", intraspecific_var_td_vector)
rownames(sigsq_matrix) <- paste0("I.var", interannual_var_td_vector)

for (indv_var in intraspecific_var_td_vector) {
  for (year_var in interannual_var_td_vector) {
    print(paste("Simulating data for", n.years, "years and", n.individuals, "individuals"))
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
        interannual_var_td = year_var,
        intraspecific_var_kd = intraspecific_var_kd,
        intraspecific_var_td = indv_var,
        interannual_var_kd = interannual_var_kd,
        noise = noise
        )

    # for estimability we estimate the expected value of
    emp_td <- mean(unique(df$uTd+120))
    data4dclone_sim <- list(K=1, X=dcdim(data.matrix(df$samps)), n=nrow(df), days=yday(df$days))
    out.parms <- c("kd", "Td", "sigsq")
    cl.seq <- c(1,4,8,16);
    cl<- makePSOCKcluster(3)
    listinnits<-list(list(kd=runif(1, min=0, max=15), Td=runif(1, min=0, max=365), sigsq=runif(1, min=0.01, max=30)),
                     list(kd=runif(1, min=0, max=15), Td=runif(1, min=0, max=365), sigsq=runif(1, min=0.01, max=30)),
                     list(kd=runif(1, min=0, max=15), Td=runif(1, min=0, max=365), sigsq=runif(1, min=0.01, max=30))
    )
    leaves.sim <- dc.parfit(cl,data4dclone_sim, params=out.parms, model=leaves, n.clones=cl.seq,
                            multiply="K",unchanged="n",
                            n.chains = n.chains, 
                            n.adapt=n.adapt, 
                            n.iter=n.iter, 
                            n.thin=n.thin,
                            n.update=n.update,
                            inits = listinnits
    )
    stopCluster(cl)
    print(summary(leaves.sim))
    print(dcdiag(leaves.sim))
    
    
    diag_table<- data.frame(dcdiag(leaves.sim)) %>% mutate(klambdamax = n.clones * lambda.max)
    cov<- sd(diag_table$klambdamax)/mean(diag_table$klambdamax)
    cov_matrix[paste0("I.var", year_var), paste0("Y.var", indv_var)] <- cov
    print(paste("Coefficient of variation of klambdamax:", cov))
    print(cov_matrix)

    lambdamax_matrix[paste0("I.var", year_var), paste0("Y.var", indv_var)] <- diag_table$lambda.max[1]
    print(paste("Maximum eigenvalue of the variance-covariance matrix:", diag_table$lambda.max[1]))
    print(lambdamax_matrix)

    sigsq_estimate <- summary(leaves.sim)$statistics["sigsq","Mean"]
    sigsq_matrix[paste0("I.var", year_var), paste0("Y.var", indv_var)] <- sigsq_estimate
    print(paste("Estimated sigsq:", sigsq_estimate))
    print(sigsq_matrix)

  }}


sigsq_matrix

sigsq_df <- as.data.frame(sigsq_matrix) %>%
  rownames_to_column("interannual_var") %>%
  pivot_longer(-interannual_var, names_to="intraspecific_var", values_to="sigsq") %>%
  mutate(
    interannual_var_num = as.numeric(gsub("I.var", "", interannual_var)),
    intraspecific_var_num = as.numeric(gsub("Y.var", "", intraspecific_var))
  )
  
windows()
ggplot(sigsq_df,
       aes(x=factor(intraspecific_var_num), y=factor(interannual_var_num), fill=sigsq)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=round(sigsq, 2)), color="white", size=4) +
  scale_fill_gradient(low="lightblue", high="blue") +
  labs(title="Estimated Observational Variance across Intraspecific and Interannual Variability",
       x="Intraspecific variability in Td (days)",
       y="Interannual variability in Td (days)",
       fill="Observational\nVariance") +
  theme_minimal()
ggsave("plots/sigsq_heatmap.png", width=8, height=6)