# But: étudier les variations du paramètre alpha en jouant sur les paramètres
# Remarque: même code que dans "simulations_beta.R"

# Hawkes processus generation ==================================================

# install.packages("patchwork")
library(ggplot2)
library(patchwork)

gc()
# cat("\014")
# dev.off()
rm(list=ls())


# Parameters ===================================================================
lambda <- 1
beta_0 <- 2
epsilon <- 10e-10

# paramètres alpha à fair varier

a0 <- 1
a1 <- 0.5
a2 <- 1.5
a3 <- 2


# nombre d'itérations
big_T <- 10


# mu function that define the Hawkes process ===================================
mu_hawkes <- function(a, b, t){
  return(a*exp(-b*t))
}

# main algorithm ===============================================================

gen_process <- function(alpha, beta, mu_h,
                        seed=1910, eps=epsilon,
                        t0=0, b_T=big_T, lbda=lambda){
  set.seed(seed)
  processus <- c(0)
  t <- t0
  
  st_l <- function(lambda=lbda, t, ti, mu){
    return(lambda + sum( mu(a=alpha, b=beta, t=t-ti)*(ti<t) ) )
  }
  
  while(t < b_T){
    
    # "Find a new upper bound"
    M <- st_l(t=t+eps, ti=processus, mu=mu_h)
    
    # "Generate next candidate point"
    E <- rexp(n=1, rate=M)
    t <- t+E
    
    # "Keep it with some probs"
    U <- runif(n=1, min=0, max=M)
    
    if( (t < b_T) && (U <= st_l(t=t, ti=processus, mu = mu_h)) ){
      processus <- c(processus, t)  
    }
  }
  return(processus)
}  

# Results ======================================================================

hawkes_0 <- gen_process(alpha=a0, beta=beta_0, mu_h=mu_hawkes)
hawkes_1 <- gen_process(alpha=a1, beta=beta_0, mu_h=mu_hawkes)
hawkes_2 <- gen_process(alpha=a2, beta=beta_0, mu_h=mu_hawkes)
hawkes_3 <- gen_process(alpha=a3, beta=beta_0, mu_h=mu_hawkes)

max1 = length(hawkes_0)
max2 = length(hawkes_3)

# Création des 4 graphiques
p1 <- ggplot() + 
  geom_step(aes(x=hawkes_0, y=0:(length(hawkes_0)-1)))+
  labs(title=paste0("Processus de Hawkes, alpha = ", a0), x="t", y="N(t)") + 
  ylim(c(0, max1)) +
  theme_bw()

p2 <- ggplot() + 
  geom_step(aes(x=hawkes_1, y=0:(length(hawkes_1)-1))) +
  labs(title = paste0("Processus de Hawkes, alpha = ", a1), x="t", y="N(t)") + 
  ylim(c(0, max1)) +
  theme_bw()

p3 <- ggplot() + 
  geom_step(aes(x=hawkes_2, y=0:(length(hawkes_2)-1))) +
  labs(title = paste0("Processus de Hawkes, alpha = ", a2), x="t", y="N(t)") + 
  ylim(c(0, max2)) +
  theme_bw()

p4 <- ggplot() + 
  geom_step(aes(x=hawkes_3, y=0:(length(hawkes_3)-1))) +
  labs(title = paste0("Processus de Hawkes, alpha = ", a3), x="t", y="N(t)") +
  ylim(c(0, max2)) +
  theme_bw()

# Assemblage 2x2
(p1 + p2) / (p3 + p4)