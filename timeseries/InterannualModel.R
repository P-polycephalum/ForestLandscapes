library(pacman)
p_load(dclone, MASS, ggplot2, snow, tidyverse, parallel)
logit.pf <- function(kd,Td,x){
  out <- kd*(x-Td)
  return(out)
}

leaves_year_fe <- function(){
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

windows()
ggplot(data, aes(x=day, y=y_norm, color=as.factor(tree_year))) +
  geom_line() +
  labs(title="Cavallinesia phenology data",
       y="Normalized leafing",
       x="Date") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")



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
    
    # windows()
    # plot(subset_data$day, subset_data$y_norm, type='l', main=paste("Tree:", trees1[i], "Year:", years1[j]),
    #      xlab="Day", ylab="Normalized Leafing")
    
    #find values below threshold
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
      #points(subset_before_threshold$day, subset_before_threshold$y_norm, col='red', pch=19)
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

View(all_before_threshold)

test.data<- log(1-all_before_threshold$y_norm) - log(all_before_threshold$y_norm)
n<- nrow(all_before_threshold)
all.days<- all_before_threshold$day
year.id<- as.numeric(as.factor(all_before_threshold$pheno_year))
n.years<- length(unique(year.id))

data4dclone <- list(K=1, X=dcdim(data.matrix(test.data)), n=n, days=all.days, year=year.id, nyear=n.years)

cl.seq <- c(1,4,8,16,32);
n.iter<-10000;n.adapt<-5000;n.update<-100;thin<-5;n.chains<-3;

cl <- makePSOCKcluster(3) 
inits.list <- list(
  list(kd=0.1, sigsq=1, yTd= runif(n=length(unique(year.id)), 1, 365)),
  list(kd=0.2, sigsq=20, yTd= runif(n=length(unique(year.id)), 1, 365)),
  list(kd=0.9, sigsq=40, yTd= runif(n=length(unique(year.id)), 1, 365))
)
annual_model<- dc.parfit(cl, data4dclone, params=c("kd","sigsq","yTd"), model=leaves_year_fe, n.clones=cl.seq,
                        multiply="K",unchanged=c("n","nyear"),
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=thin
)

summary(annual_model)
dcdiag(annual_model)

dctable <- dctable(annual_model)

windows()
plot(dctable, 1:4)

results <- summary(annual_model)
results$statistics[, "Mean"]

x_range<- seq(1,365, by=1)

# Get mapping from year index to actual phenological year
pheno_year_levels <- levels(as.factor(all_before_threshold$pheno_year))

# Generate predictions for all years
all_predictions <- data.frame()
for (y in 1:n.years) {
  kd <- results$statistics["kd", "Mean"]
  Td <- results$statistics[paste0("yTd[", y, "]"), "Mean"]
  predicted <- logit.pf(kd, Td, x_range)
  predicted <- 1 / (1 + exp(predicted))
  
  this_year <- data.frame(
    day = x_range,
    y_norm = predicted,
    pheno_year = pheno_year_levels[y]
  )
  all_predictions <- bind_rows(all_predictions, this_year)
}

windows()
ggplot(all_before_threshold, aes(x=day, y=y_norm, group=tree_year)) +
  geom_point(aes(color=as.factor(tag))) +
  geom_line(data=all_predictions, aes(x=day, y=y_norm, group=pheno_year), 
            color="red", size=0.5) +
  facet_wrap(~pheno_year,ncol=1) +
  labs(title="Cavallinesia phenology data",
       y="Predicted leafing",
       x="Day of year") +
  theme_minimal()

#what is the sd of the year effects?
sd(results$statistics[grep("yTd", rownames(results$statistics)), "Mean"])
mean(results$statistics[grep("yTd", rownames(results$statistics)), "Mean"])
summary(annual_model)


#when does Av becomes identifiable? how many years we need
#n.years<- c(7,14,21,28)
min_date<- as.Date("2018-01-01")
interval_days<- 30
interval_var<- 2
n.years<-7


interannual_var_td<- 15
interannual_vect<- c(5,10,15,20)
intraspecific_var_td<- 5
intraspecific_var_vect<- c(3,5,7,10)
noise<- 5
Td<- 100
kd<- 0.1
n.individuals<- 10


for (year_var in interannual_vect) {
  for (indv_var in intraspecific_var_vect) {

    print(paste("Running model with", year_var, "interannual variance and", indv_var, "intraspecific variance"))
    max_date <- min_date + years(n.years) - days(1)
    print(paste("Sampling dates from", min_date, "to", max_date))
    
    # Generate sampling dates with variable intervals
    sampling_dates <- min_date
    current_date <- min_date
    
    while (current_date < max_date) {
      interval <- rnorm(1, mean=interval_days, sd=interval_var)
      interval <- max(1, interval)  # Ensure positive interval
      current_date <- current_date + days(round(interval))
      if (current_date <= max_date) {
        sampling_dates <- c(sampling_dates, current_date)
      }
    }
    
    print(paste("Number of sampling dates:", length(sampling_dates)))

    all.days <- rep(sampling_dates, n.individuals)
    calendar_years <- year(all.days)
    year_indices <- as.numeric(as.factor(calendar_years))  # Convert to sequential indices 1, 2, 3, ...
    indv.id <- rep(1:n.individuals, each=length(sampling_dates))
    
    if (length(all.days) != length(year_indices) || length(all.days) != length(indv.id)) {
      stop("Length mismatch: all.days, year_indices, and indv.id must have equal length")
    }
    
    # Simulate data
    uTd <- rnorm(n=n.years, mean=0, sd=year_var)
    print(paste("Simulated interannual effects (uTd):", round(Td+uTd, 2)))
    iTd<- rnorm(n=n.individuals, mean=0, sd=indv_var)
    yTd <- Td + uTd[year_indices] + iTd[indv.id]

    # for (uTd_val in uTd) {
    #   for (iTd_val in iTd) {
    #     print(paste("  → Simulating with uTd:", round(uTd_val, 2), "and iTd:", round(iTd_val, 2)))
    #     print(paste("resulting yTd :", round(Td + uTd_val + iTd_val, 2)))
    #   }
    # }

    sampling_dates_doy <- yday(all.days)  # Get day of year for all observations
    mu.true <- logit.pf(kd=kd, Td=yTd[year_indices], x=sampling_dates_doy)
    norm.samps <- rnorm(n=length(all.days), mean=mu.true, sd=sqrt(noise))
    y.sims <- 1/(1+exp(norm.samps))

    df <- data.frame(
      days=all.days,
      indv=indv.id,
      year=calendar_years,
      y=y.sims,
      samps=norm.samps
    ) %>% mutate(indv_year= as.factor(paste0(indv,"_",year)),
                doy= yday(days))

    data4dclone <- list(K=1,
                        X=dcdim(data.matrix(df$samps)),
                        n=nrow(df),
                        days=df$doy,
                        year=as.numeric(as.factor(df$year)),
                        nyear=length(unique(df$year))
    )

    cl <- makePSOCKcluster(3)
    inits.list <- list(
      list(kd=0.1, sigsq=1, yTd= runif(n=length(unique(df$year)), 1, 365)),
      list(kd=0.2, sigsq=6, yTd= runif(n=length(unique(df$year)), 1, 365)),
      list(kd=0.05, sigsq=10, yTd= runif(n=length(unique(df$year)), 1, 365))
    )
    annual_iter<- dc.parfit(cl, data4dclone, params=c("kd","sigsq","yTd"), model=leaves_year_fe, n.clones=cl.seq,
                          multiply="K",unchanged=c("n","nyear"),
                          n.chains = n.chains, 
                          n.adapt=n.adapt, 
                          n.update=n.update,
                          n.iter = n.iter, 
                          thin=thin,
                          inits=inits.list
    )
    print(summary(annual_iter))
    print(mean(summary(annual_iter)$statistics[grep("yTd", rownames(summary(annual_iter)$statistics)), "Mean"]))
    print(sd(summary(annual_iter)$statistics[grep("yTd", rownames(summary(annual_iter)$statistics)), "Mean"]))
    print(dcdiag(annual_iter))
    windows()
    plot(dctable(annual_iter), 1:4)
    windows()
    plot(dctable(annual_iter), 5:8)
    }
    
}

windows()
ggplot(df, aes(x=doy, y=y, group=as.factor(indv_year))) +
      geom_line(aes(color=as.factor(year))) +
      facet_wrap(~year, ncol=1) +
      labs(title=paste("Simulated phenology data with", year.n, "years"),
          y="Simulated y",
          x="Day of year") +
      theme_minimal()
View(df)
dcdiag(annual_model) 
yTd <- rnorm(n=n.years, mean=Td, sd=aV)
mu.true <- logit.pf(kd=kd,Td=yTd[year.id],x=all.days)
norm.samps <- rnorm(n=n, mean=mu.true, sd=sqrt(sigsq))
y.sims <- 1/(1+exp(norm.samps))

df<- data.frame(
  days=all.days,
  indv=indv.id,
  year=year.id,
  y=y.sims
) %>% mutate(indv_year= as.factor(paste0(indv,"_",year)))

windows()
ggplot(df, aes(x=days, y=y, group=as.factor(indv_year))) +
  geom_line(aes(color=as.factor(year))) +
  labs(title="Simulated phenology data",
       y="Simulated y",
       x="Day of year") +
  theme_minimal()
