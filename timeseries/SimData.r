library(pacman)
p_load(dclone, MASS, ggplot2, snow, tidyverse, parallel, gridExtra, lubridate, reshape2)

logit.pf <- function(kd,Td,x){
  out <- kd*(x-Td)
  return(out)
}
logit.pf_inv <- function(kd,Td,x){
  out <- 1/(1+exp(-kd*(x-Td)))
  return(out)
}

#lets try real data with the same model
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



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
    threshold <- 0.2
    below_threshold <- which(subset_data$y_norm < threshold)
    
    
    cat("  → Values below threshold:", length(below_threshold), "\n")
    cat("  → Range of y_norm:", min(subset_data$y_norm, na.rm=T), "to", max(subset_data$y_norm, na.rm=T), "\n")
    
    if (length(below_threshold) == 0) {
      cat("  → No values below threshold; skipping\n\n")
      print(below_threshold)
      next
    }
    
    first_item <- min(below_threshold, na.rm=TRUE)
    cat("  → First below-threshold index:", first_item, "of", nrow(subset_data), "\n")
    
    #subset all values from first_item onwards
    if (first_item > 1 && first_item <= nrow(subset_data)) {
      subset_after_threshold <- subset_data[first_item:nrow(subset_data), ]
      #points(subset_after_threshold$day, subset_after_threshold$y_norm, col='red', pch=19)
      all_before_threshold <- bind_rows(all_before_threshold, subset_after_threshold)
      cat("  → Added", nrow(subset_after_threshold), "rows\n\n")
    } else {
      
      cat("  → First below-threshold item at start or beyond range; skipping\n\n")
      print(first_item)
      print(nrow(subset_data))
      print(first_item > 1 && first_item <= nrow(subset_data))
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

nrow(all_before_threshold)

leaves <- function(){
  kd ~ dunif(0, 2)
  Td ~ dunif(90, 360)
  sigsq ~ dunif(1.0E-6, 40)

  for(j in 1:n){
    muf[j] <- -kd * (days[j] - Td)
  }

  for(k in 1:K){
    for(i in 1:n){
      X[i,k] ~ dnorm(muf[i], 1/sigsq)
    }
  }
}

#logit transform 
test.data <- log(1-all_before_threshold$y_norm) - log(all_before_threshold$y_norm)
data4dclone <- list(K=1, X=dcdim(data.matrix(test.data)), n=nrow(all_before_threshold), days=all_before_threshold$day)
cl.seq <- c(1,4,8,16);
n.iter<-10000;n.adapt<-5000;n.update<-100;thin<-5;n.chains<-3;
out.parms <- c("kd", "Td", "sigsq")
cl<- makePSOCKcluster(3)
cava.intercept <- dc.parfit(cl,data4dclone, params=out.parms, model=leaves, n.clones=cl.seq,
                        multiply="K",unchanged="n",
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=thin,
                        inits = list(
                          list(kd = 0.1, Td = 90, sigsq = 10),
                          list(kd = 0.2, Td = 200,  sigsq = 15),
                          list(kd = 0.3, Td = 300, sigsq = 12)
                        )
)


table<-summary(cava.intercept)
windows()
grid.table(round(table$statistics[, c("Mean", "SD",'DC SD', "R hat")], 2))

windows()
plot(cava.intercept)

mcmc_recortado <- cava.intercept[1:100, ]

# Create one plot with all parameters and chains in different colors
windows()

plot(1:100, mcmc_recortado[[1]][, 'Td'], type='l', col='red', lwd=2,
     main="Trace plot (first 100 iterations) - All parameters and chains",
     xlab="Iteration", ylab="Parameter value", na.rm=T)
lines(1:100, mcmc_recortado[[2]][, 'Td'], col='blue', lwd=2, lty=1)
lines(1:100, mcmc_recortado[[3]][, 'Td'], col='green', lwd=2, lty=1)

dcdiagtable <- dcdiag(cava.intercept)
windows()
grid.table(round(dcdiagtable, 2))
cavatable <- dctable(cava.intercept)
windows()
plot(cavatable)
windows()
plot(cavatable, type="log.var")

#generate predicted values
logit.pf(kd=-coef(cava.intercept)["kd"],Td=coef(cava.intercept)["Td"],x=seq(1,365,1))->LC_values
LC_values<- 1/(1+exp(LC_values))

x <- seq(1, 365)

calendar_doy <- ((x + 243 - 1) %% 365) + 1

results_df <- data.frame(
  x = x,
  y = LC_values,
  doy = calendar_doy,
  date= as.Date(calendar_doy - 1, origin = "2018-01-01")
)
head(results_df)

#what date is 90%
date10 <- results_df %>% filter(y <= 0.1) %>% slice(n())
date50 <- results_df %>% filter(y <= 0.5) %>% slice(n()) 
date90 <- results_df %>% filter(y <= 0.9) %>% slice(n()) 

day_date<- seq(from=as.Date("2018-09-01"), to=as.Date("2019-08-31"), by="month")
day_date_lookup <- data.frame(
  day = as.numeric(difftime(day_date, as.Date("2018-09-01"), units = "days")),
  date = day_date
)


windows()
ggplot(all_before_threshold,
       aes(x = day, y = y_norm, color = as.factor(tree_year))) +
  geom_point() +
  geom_line(data = results_df, aes(x = x, y = y),
            color = "black", linewidth = 1) +
  scale_x_continuous(
    breaks = day_date_lookup$day,
    labels = format(day_date_lookup$date, "%b %d")
  ) +
  labs(
    title = "Cavallinesia phenology data",
    y = "Predicted leafing",
    x = "Date"
  ) +theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")+
  geom_vline(xintercept = c(date10 %>% pull(x),
                             date50 %>% pull(x),
                             date90 %>% pull(x)),
             linetype = "dashed", color = "grey")

print(date50)
print(date90)
print(date10)

######################
#interannual simulation
######################
##Simulate year effects

n.years  <- 7
one.year <- seq(from=1,to=365,by=30)
samp.days <- rep(one.year,n.years)
n.inds <- 10
all.days <- rep(samp.days,n.inds)
n <- length(all.days)
year.id  <- rep(rep(1:n.years, each = length(one.year)), n.inds)
indv.id  <- rep(1:n.inds, each = length(samp.days))

sigsq <- 0.45  #noise levels
kd <- 0.1      #leaf drop rate
Td <- 100      #day of year when 50% of leaves are dropped
aV<- 15        #anual variability in Td
                    
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


leaves_year <- function(){
  lkd ~ dnorm(0,0.4)
  kd <- exp(lkd)

  ls ~ dnorm(0,1)
  sigsq <- pow(exp(ls),2)

  for (y in min(year):max(year)) {
  yTd[y] ~ dnorm(90, 0.5)
  }
  
  for(j in 1:n){
    muf[j] <-  kd*(days[j]-yTd[year[j]])
  }

  for(k in 1:K){
    for(i in 1:n){
      X[i,k] ~ dnorm(muf[i],1/sigsq)
    } 
  }
}
unique(all_before_threshold$year)
test.data <- log(1-y.sims) - log(y.sims)
data4dclone <- list(K=1, X=dcdim(data.matrix(test.data)), n=n, days=all.days, year=year.id, nyear=n.years)

cl.seq <- c(1,4,8,16);
n.iter<-10000;n.adapt<-5000;n.update<-100;thin<-10;n.chains<-3;

cl <- makePSOCKcluster(3) 
annual_model<- dc.parfit(cl, data4dclone, params=c("kd","sigsq","yTd"), model=leaves_year, n.clones=cl.seq,
                        multiply="K",unchanged=c("n","nyear"),
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=thin,
                        inits=list(lkd=log(0.2))
)

summary(annual_model)
yTd

mean(yTd)
sd(yTd)
 
mean(coef(annual_model)[grep("yTd", names(coef(annual_model)))])
sd(coef(annual_model)[grep("yTd", names(coef(annual_model)))])



dcdiag(annual_model)
annual_table <- dctable(annual_model)
windows()
plot(annual_table,1:6)
windows()
plot(annual_table, type="log.var",6:8)


######################
#interannual simulation
######################
##Simulate year effects, one global intercept and a random effect for the shift

n.years  <- 7
one.year <- seq(from=1,to=365,by=30)
samp.days <- rep(one.year,n.years)
n.inds <- 16

all.days <- rep(samp.days,n.inds)

all.days<- all_before_threshold$day
n <- length(all.days)
year.id  <- rep(rep(1:n.years, each = length(one.year)), n.inds)
indv.id  <- rep(1:n.inds, each = length(samp.days))

sigsq <- 11  #noise levels
kd <- 0.05      #leaf drop rate
Td <- 120      #day of year when 50% of leaves are dropped
aV<- 21        #anual variability in Td
                    
uTd <- rnorm(n=n.years, mean=0, sd=aV)
yTd <- Td + uTd[year.id]
mu.true <- logit.pf(kd=kd,Td=yTd,x=all.days)
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


leaves_year_global <- function(){
  lkd ~ dnorm(log(0.1),0.4)
  kd <- exp(lkd)

  ls ~ dnorm(log(15),1)
  sigsq <- pow(exp(ls),2)

  ltd ~ dnorm(log(100),1)
  Td <- exp(ltd)

  log.aV ~ dnorm(log(15),4)
  aV <- exp(log.aV)
  tau <- 1/pow(aV,2) 

  for (y in min(year):max(year)) {
    uRaw[y] ~ dnorm(0, tau)
  }
  
  u_bar <- mean(uRaw[min(year):max(year)])

  for (y in min(year):max(year)) {
    uY[y] <- uRaw[y] - u_bar
    yTd[y] <- Td + uY[y]
  }
  
  for(j in 1:n){
    muf[j] <-  kd*(days[j]-yTd[year[j]])
  }

  for(k in 1:K){
    for(i in 1:n){
      X[i,k] ~ dnorm(muf[i],1/sigsq)
    } 
  }
}

test.data <- log(1-y.sims) - log(y.sims)
data4dclone <- list(K=1, X=dcdim(data.matrix(test.data)), n=n, days=all.days, year=year.id, nyear=n.years)

cl.seq <- c(1,4,8,16);
n.iter<-1000;n.adapt<-500;n.update<-100;thin<-10;n.chains<-3;

cl <- makePSOCKcluster(3) 
annual_model<- dc.parfit(cl,data4dclone, params=c("Td","kd","sigsq","uY"), model=leaves_year_global, n.clones=cl.seq,
                        multiply="K",unchanged=c("n","nyear"),
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=thin,
                        partype= "balancing")

sum_model <- summary(annual_model)
coef(annual_model)
dcdiag(annual_model)
table_summary <- round(cbind(True_Value = c(Td, kd, sigsq, uTd), 
                             sum_model$statistics[, c("Mean", "SD", "R hat")], 
                             sum_model$quantiles[, c("2.5%", "97.5%")]), 2)

windows()
grid.table(table_summary)


dcdiag(annual_model)
annual_table <- dctable(annual_model)
windows()
plot(annual_table,1:6)
windows()
plot(annual_table, type="log.var",7:10)
windows()
plot(annual_table, type="log.var",1:6)

#######################################
## real data with year effects model
#####################################
test.data <- log(1-all_before_threshold$y_norm) - log(all_before_threshold$y_norm)
all_before_threshold$pheno_year <- as.numeric(as.factor(all_before_threshold$pheno_year))

data4dclone <- list(K=1,
                    X=dcdim(data.matrix(test.data)),
                    n=nrow(all_before_threshold),
                    days=all_before_threshold$day,
                    year=all_before_threshold$pheno_year,
                    nyear=length(unique(all_before_threshold$pheno_year)))

cl.seq <- c(1,4);
n.iter<-10000;n.adapt<-5000;n.update<-100;thin<-10;n.chains<-3;
cl <- makePSOCKcluster(3)
cava.year <- dc.parfit(cl, data4dclone, params=c("Td","kd","sigsq","uY","tau"), model=leaves_year_global, n.clones=cl.seq,
                        multiply="K",unchanged=c("n","nyear"),
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=thin)


#summmary of the model
summary(cava.year)
coef(cava.year)[["tau"]]^(-1/2)
dcdiag(cava.year)
quantile(cava.year, probs = c(0.025))[["Td"]]
sum_model <- summary(cava.year)
table_summary <- round(cbind(sum_model$statistics[, c("Mean", "SD", "R hat")], 
                             sum_model$quantiles[, c("2.5%", "97.5%")]), 2)

windows()
grid.table(table_summary)

results_df <- data.frame(
  x=numeric(),
  y=numeric(),
  doy=numeric(),
  date=as.Date(character())
)
for (year in 1:length(unique(all_before_threshold$pheno_year))) {

  cat("Estimated uY for year", year, ":", 
      coef(cava.year)[paste0("uY[", year, "]")], "\n")

  x <- seq(1, 365)

  mu <- coef(cava.year)[["Td"]] + coef(cava.year)[[paste0("uY[", year, "]")]]
  mu_low <- coef(cava.year)[["Td"]] + 
            quantile(cava.year, probs = c(0.025))[[paste0("uY[", year, "]")]]
  mu_high <- coef(cava.year)[["Td"]] + 
             quantile(cava.year, probs = c(0.975))[[paste0("uY[", year, "]")]]

  # Mean
  LC_mean <- logit.pf(kd = coef(cava.year)["kd"], Td = mu, x = x)
  y_mean <- 1 / (1 + exp(LC_mean))

  # CI
  LC_low <- logit.pf(kd = coef(cava.year)["kd"], Td = mu_low, x = x)
  y_low <- 1 / (1 + exp(LC_low))

  LC_high <- logit.pf(kd = coef(cava.year)["kd"], Td = mu_high, x = x)
  y_high <- 1 / (1 + exp(LC_high))

  calendar_doy <- ((x + 243 - 1) %% 365) + 1

  temp_df <- data.frame(
    x = x,
    y = y_mean,
    y_low = y_low,
    y_high = y_high,
    doy = calendar_doy,
    date = as.Date(calendar_doy - 1, origin = "2018-01-01"),
    year = as.factor(year)
  )

  results_df <- bind_rows(results_df, temp_df)
}
head(results_df)


day_date_lookup <- data.frame(
  day = as.numeric(difftime(seq(from=as.Date("2018-09-01"), to=as.Date("2019-08-31"), by="month"),
                            as.Date("2018-09-01"), units = "days")),
  date = seq(from=as.Date("2018-09-01"), to=as.Date("2019-08-31"), by="month")
)
year_labels <- c(
  "1" = "2017–2018",
  "2" = "2018–2019",
  "3" = "2019–2020",
  "4" = "2020–2021",
  "5" = "2021–2022",
  "6" = "2022–2023",
  "7" = "2023–2024"
)

windows()
ggplot(results_df, aes(x = x)) +
  geom_ribbon(
    aes(ymin = y_low, ymax = y_high, fill = year),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(
    aes(y = y, color = year),
    linewidth = 1
  ) +
  geom_point(
    data = all_before_threshold %>% 
      filter(pheno_year %in% unique(results_df$year)),
    aes(x = day, y = y_norm, color = as.factor(pheno_year)),
    size = 1.5,
    alpha = 0.7
  ) +
  scale_x_continuous(
    breaks = day_date_lookup$day,
    labels = format(day_date_lookup$date, "%b %d")
  ) +
  scale_color_discrete(labels = year_labels) +
  scale_fill_discrete(labels = year_labels) +
  labs(
    title = "Cavallinesia planatifolia - Year Effects Model",
    y = "Predicted Leaf Coverage",
    x = "Date",
    color = "Year",
    fill = "Year"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/cava_year_effects_model.png", width = 10, height = 6)





cavayeartable <- dctable(cava.year)
windows()
plot(cavayeartable,1:6)
windows()
plot(cavayeartable,6:9)

windows()
plot(cavayeartable,1:4, type="log.var")
windows()
plot(cavayeartable,5:9, type="log.var")


