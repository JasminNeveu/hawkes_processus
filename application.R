library(ETAS)
library(chron)
library(stats)

data(japan.quakes)

temp<-chron(as.character(japan.quakes[,1]),
            as.character(japan.quakes[,2]),
            format=c(date="y-m-d",times="h:m:s"))
temp<-as.numeric(temp)/365
temp<-temp-min(temp)

# estimations des parametres lamba alpha beta

k = length(temp)
A = function(temp, beta){
  A_vals = numeric(k)
  A_vals[1] = 0
  for(i in 2:k){
    A_vals[i] = exp(-beta * (temp[i] - temp[i-1])) * (1 + A_vals[i-1])
  }
  return(A_vals)
}

log_vrais = function(theta, temp){
  lamb  = theta[1]
  alpha = theta[2]
  beta  = theta[3]
  A_vals = A(temp, beta)
  s1 = sum(log(lamb + alpha * A_vals))
  s2 = (alpha / beta) * sum(exp(-beta * (temp[k] - temp)) - 1)
  return(-(s1 - lamb * temp[k] + s2))
}
startVals <- c(1,1,2)

res = optim(
  par = startVals,
  fn  = log_vrais,
  temp = temp,
  method = "L-BFGS-B",
  control = list(maxit=1000))

lamb_hat <- res$par[1]
alpha_hat <- res$par[2]
beta_hat <- res$par[3]


# 

lamb_star_i = function(t){
  lamb[i] + sum()
}
gen_hawkes_process = function(lamb_star,temp){
  epsilon <- 10**(-10)
  T <- k
  t <- 0
  P <- numeric(k)
  while(t <= T){
    M <- lamb_star(t+epsilon)
    E <- exp_distrib(M)
    t <- t + E
    U <- inifor(0,M)
    if(t <= T & U <= lamb_star(t)){
      P <- [P,t]
    }
  }
  return(P)
}
