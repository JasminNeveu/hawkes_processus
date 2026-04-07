# Introduction =================================================================

# install.packages("ETAS")
# install.packages("chron")
library(ETAS)
library(chron)
library(ggplot2)

setwd("~/work/hawkes_processus/simulations")
gc()
# cat("\014")
# dev.off()
rm(list=ls())

data(japan.quakes)

# donné par M.Lavancier => variable de temps facilement exploitable
temp<-chron::chron(as.character(japan.quakes[,1]),
                   as.character(japan.quakes[,2]),
                   format=c(date="y-m-d",times="h:m:s"))

temp<-as.numeric(temp)/365
temp<-temp-min(temp)

plot(temp)

# 0 - production de graphes des données étuduiées ==============================
### les données sous forme de step function =====
n <- length(temp)-1
plot1 <- ggplot()+
  geom_step(aes(x=temp, y=0:n))+
  labs(x = "temps (t)", y="N(t)", title="")+
  theme_bw(base_size = 10)

plot1

ggsave(filename = "plot_1_japan_earthq.pdf", plot=plot1, width = 6,         
       height = 5,        
       units = "cm",       
       dpi = 300,          
       device = "pdf")
### les données sous forme des m 1ers instants de sauts =====
m <-100


plot2 <- ggplot()+
  geom_point(aes(x=temp[1:m], y=rep(1, m)))+
  labs(x = "temps (t)", y="", title="")+
  theme_bw()

plot2

ggsave(filename = "plot_2_japan_earthq.pdf", plot=plot2, width = 6,         
       height = 5,        
       units = "cm",       
       dpi = 300,          
       device = "pdf")

# I - Likelihood function ======================================================
l <- function(teta, ti=temp){
  lambda = teta[1]
  alpha = teta[2]
  beta = teta[3]
  
  k = length(ti)
  Ai <- numeric(k)
  Ai[1] <- 0
  
  for(i in 2:k){
    t11 <- exp( -beta*(ti[i]-ti[i-1]) )
    t12 <- 1 + Ai[i-1]
    A[i] <-  t11*t12
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
### calcul du maximum likelihood =====
set.seed(1149)
st_opti <- optim(par=c(1, 1, 2), fn=l, control=list(maxit=10000))

res <- st_opti$par
st_lambda <- res[1]
st_alpha <- res[2]
st_beta <- res[3]


### évolution de la log_vraissamblance avec l'algo =====
# N <- 1000
# opti <- matrix(ncol=3, nrow=N)
# for(i in 1:N){
#   optim_i <- optim(par=c(1, 1, 2), fn=l, control=list(maxit=i),
#                    method="L-BFGS-B",
#                    lower = c(1e-5, 1e-5, 1e-5),
#                    upper = c(2000, 2000, 2000))
#   opti[i, ] <- optim_i$par
#}

# v2:
N <- 1000
opti <- matrix(ncol=3, nrow=N+1)

opti[1, ] <- c(1, 1, 2)
for(i in 1:N+1){
  optim_i <- optim(par=opti[i-1, ], fn=l, control=list(maxit=1),
                   method="L-BFGS-B", 
                   lower = c(1e-5, 1e-5, 1e-5),
                   upper = c(2000, 2000, 2000))
  opti[i, ] <- optim_i$par
}
# 

ggplot()+
  geom_point(aes(x=1:N, y=opti[, 1]))+
  labs(main="évolution du maximum de vraisemblance de lambda", 
       x="nombre d'itérations", y="maximum de vraisemblance en lambda")+
  theme_bw()

ggplot()+
  geom_point(aes(x=1:N, y=opti[, 2]))+
  labs(main="évolution du maximum de vraisemblance d'alpha", 
       x="nombre d'itérations", y="maximum de vraisemblance en alpha")+
  theme_bw()

ggplot()+
  geom_point(aes(x=1:N, y=opti[, 3]))+
  labs(main="évolution du maximum de vraisemblance de beta", 
       x="nombre d'itérations", y="maximum de vraisemblance en beta")+
  theme_bw()

### exemple avec 4 graines différentes de l'estimateur =====
M <- 100
n <- length(temp)-1

p <- c()
set.seed(1149)
for(i in 1:M){
  st_hawkes_proc <- gen_process(lbda = st_lambda, alpha=st_alpha,
                                beta = st_beta)
  ni <- length(st_hawkes_proc)-1
  pi <- ggplot()+
    geom_step(aes(x=st_hawkes_proc, y=0:ni, color=paste0("simulation", i)))+
    geom_step(aes(x=temp, y=0:n, color="data"))+
    labs(x = "temps (t)", y="N(t)", title="")+
    theme_bw()
  
  p <- c(p, pi)
}

p1 <- p[1]
p2 <- p[2]
p3 <- p[3]
p4 <- p[4]

ggsave(filename = "plot_data_model_4.pdf", plot=p4, width = 6,         
       height = 5,        
       units = "cm",       
       dpi = 300,          
       device = "pdf")



