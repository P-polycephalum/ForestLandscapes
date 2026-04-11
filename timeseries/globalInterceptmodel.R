library(pacman)
p_load(dclone, MASS, ggplot2, snow, tidyverse, parallel)
logit.pf <- function(kd,Td,x){
  out <- kd*(x-Td)
  return(out)
}

n.years  <- 7
# Sample 12 dates per year with 25-35 day gaps.
make_year_days <- function(n_dates = 12, gap_min = 25, gap_max = 35, max_day = 365) {
  repeat {
    gaps <- sample(gap_min:gap_max, size = n_dates - 1, replace = TRUE)
    start_day <- sample(1:5, size = 1)
    days <- start_day + c(0, cumsum(gaps))
    if (max(days) <= max_day) {
      return(days)
    }
  }
}

year.days <- lapply(seq_len(n.years), function(x) make_year_days())
samp.days <- unlist(year.days, use.names = FALSE)
n.inds <- 10
all.days <- rep(samp.days, n.inds)
n <- length(all.days)
year.id <- rep(seq_len(n.years), times = vapply(year.days, length, 1L))
year.id <- rep(year.id, n.inds)
indv.id <- rep(seq_len(n.inds), each = length(samp.days))

sigsq <- 0.45  #noise levels 
kd <- 0.1
Td <- 100
mu.true <- logit.pf(kd=kd,Td=Td,x=all.days)

norm.samps <- rnorm(n=n, mean=mu.true, sd=sqrt(sigsq))
y.sims <- 1/(1+exp(norm.samps))

df<- data.frame(
  days=all.days,
  indv=indv.id,
  year=year.id,
  y=y.sims
)

windows()
ggplot(df, aes(x=days, y=y, color=as.factor(year))) +
  geom_point() +
  labs(title="Simulated phenology data",
       y="Simulated y",
       x="Day of year") +
  theme_minimal()


##JAGS model for intercepts
leaves <- function(){
  lkd ~ dnorm(0,0.4)
  kd <- exp(lkd)
  ltd ~ dnorm(0,4)
  Td <- exp(ltd)
  ls ~ dnorm(0,1)
  sigsq <- pow(exp(ls),2)
  for(j in 1:n){
    muf[j] <-  kd*(days[j]-Td)
  }
  for(k in 1:K){
    for(i in 1:n){
      X[i,k] ~ dnorm(muf[i],1/sigsq)
    } 
  }
}



test.data <- log(1-y.sims) - log(y.sims)
data4dclone <- list(K=1, X=dcdim(data.matrix(test.data)), n=n, days=all.days)

cl.seq <- c(1,4,8,16);
n.iter<-10000;n.adapt<-5000;n.update<-100;thin<-10;n.chains<-3;

cl<- makePSOCKcluster(3)
out.parms <- c("kd", "Td", "sigsq")
leaves.dclone <- dc.parfit(cl,data4dclone, params=out.parms, model=leaves, n.clones=cl.seq,
                        multiply="K",unchanged="n",
                        n.chains = n.chains, 
                        n.adapt=n.adapt, 
                        n.update=n.update,
                        n.iter = n.iter, 
                        thin=thin
                        #inits=list(lkd=log(0.2), ltd=log(40))
                        )

dcdiag(leaves.dclone)
summary(leaves.dclone)

dctable <- dctable(leaves.dclone)
windows()
plot(dctable)
windows()
plot(dctable, type="log.var")

