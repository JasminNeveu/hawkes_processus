# Notes:
# 2) La fonction de génération d'un processus de Hawkes est calquée sur l'algo
#    écrit dans le papier étudié
# Les codes fonctionnent mais il peut y avoir des erreurs quand même ....

# Introduction =================================================================

# install.packages("ETAS")
# install.packages("crhon")
library(ETAS)
library(chron)
library(ggplot2)

gc()
# cat("\014")
# dev.off()
rm(list=ls())

data(japan.quakes)

temp<-chron::chron(as.character(japan.quakes[,1]),
                   as.character(japan.quakes[,2]),
                   format=c(date="y-m-d",times="h:m:s"))

temp<-as.numeric(temp)/365
temp<-temp-min(temp)

plot(temp)


# I - Likelihood function ======================================================
l <- function(teta, ti=temp){
  lambda = teta[1]
  alpha = teta[2]
  beta = teta[3]
  
  k = length(ti)
  Ai <- c(0)
  for(i in 2:k){
    t11 <- exp( -beta*(ti[i]-ti[i-1]) )
    t12 <- 1 + Ai[i-1]
    Ai <- c(Ai, t11*t12)
  }
  t21 <- sum(log(lambda + alpha*Ai))
  t22 <- -lambda*ti[k]
  t23 <- (alpha/beta) * sum( exp(-beta*(ti[k]-ti)) -1 )
  
  return(-(t21+t22+t23))
}
# II - hawkes process gen ======================================================
gen_process <- function(eps=10e-10,
                        t0=0, b_T=80, lbda, alpha, beta){
  processus <- c(0)
  t <- t0
  
  f_mu <- mu_hawkes <- function(a=alpha, b=beta, t){
    return(a*exp(-b*t))
  }
  st_l <- function(lambda=lbda, t, ti, mu){
    return(lambda + sum( mu(t=t-ti)*(ti<t) ) )
  }
  
  while(t < b_T){
    
    # "Find a new upper bound"
    M <- st_l(t=t+eps, ti=processus, mu=f_mu)
    
    # "Generate next candidate point"
    E <- rexp(n=1, rate=M)
    t <- t+E
    
    # "Keep it with some probs"
    U <- runif(n=1, min=0, max=M)
    
    if( (t < b_T) && (U <= st_l(t=t, ti=processus, mu = f_mu)) ){
      processus <- c(processus, t)  
    }
  }
  return(processus)
}  

# III -  data likelihhod computation ===========================================

opti <- optim(par=c(1, 1, 2), fn=l, control=list(maxit=1000))

res <- opti$par
st_lambda <- res[1]
st_alpha <- res[2]
st_beta <- res[3]


set.seed(1149)
p <- c()

for(i in 1:4){
  st_hawkes_proc <- gen_process(lbda = st_lambda, alpha=st_alpha,
                                beta = st_beta)
  pi <- ggplot()+
    geom_step(aes(x=temp, y=1:length(temp), col="séismes"))+
    geom_step(aes(x=st_hawkes_proc, y=1:length(st_hawkes_proc),
                  col="Hawkes model realisation"))+
    xlim(c(1,2))+
    ylim(100,200)+
    theme_bw()
  p <- c(p, pi)
}

p1 <- p[1]
p2 <- p[2]
p3 <- p[3]
p4 <- p[4]

p1
p2
p3
p4
