Graf <- function(u,d,n=25){
  A<-matrix(0,n,n)
  for(i in 1:n){
    for(j in i:n){
      A[i,j]<-u^(j-i)*d^(i-1)
    }
  }
  return(A)
}
#Graf(1.1,0.9,5)

Vertex_Value <- function(p,r,dt,Cu,Cd){
  return(exp(-r*dt)*(p*Cu+(1-p)*Cd))
}



EOP_tree <- function(S0=50,r=0.05,K=52,Time=2,dt=1,sigma=0.3, call=1, EU=1){
  
  u <- exp(sigma*sqrt(dt))
  d <- exp(-sigma*sqrt(dt))
  n <- Time/dt+1
  p <- (exp(r*dt)-d)/(u-d)
  
  u <- 1.2
  d <- 0.8
  p <- (exp(r*dt)-d)/(u-d)
  
  #tree <- matrix(rep(0,625), nrow = 25)
  #first_line <- c(1, cumprod(c(rep(u,24))))
  #work_tree <- 
  graph <- (Graf(u,d,n) * S0 - K)*(-1)^(call + 1)
  graph[graph < 0] <- 0
  
  Cofanie <- matrix(0,Time/dt+1,Time/dt+1)
  Cofanie[,Time/dt+1] <- graph[,Time/dt+1]
  
  if (EU == 1) {
    for(j in (Time/dt+1):2){
      for(k in 1:j-1){
        Cofanie[k,j-1] <- Vertex_Value(p,r,dt,Cofanie[k,j],Cofanie[k+1,j])
      }
    }
  }
  else {
    for(j in (Time/dt+1):2){
      for(k in 1:j-1){
        Cofanie[k,j-1] <- max(Vertex_Value(p,r,dt,Cofanie[k,j],Cofanie[k+1,j]), graph[k,j-1])
        #print(Vertex_Value(p,r,dt,Cofanie[k,j],Cofanie[k+1,j]))
        #print(graph[k,j-1])
      }
    }
  }
  
  
  return(Cofanie[1,1])
}

EOP_tree(call=0, EU=0)

###################### Testy wrażliwości:

#test dla K

vector_K <- seq(40, 60, 1)

values_K_eu_call <- sapply(vector_K, EOP_tree, S0=50,r=0.05,dt=1/12,Time=2,sigma=0.3, call=1, EU=1)
values_K_eu_put <- sapply(vector_K, EOP_tree, S0=50,r=0.05,dt=1/12,Time=2,sigma=0.3, call=0, EU=1)
values_K_usa_call <- sapply(vector_K, EOP_tree, S0=50,r=0.05,dt=1/12,Time=2,sigma=0.3, call=1, EU=0)
values_K_usa_put <- sapply(vector_K, EOP_tree, S0=50,r=0.05,dt=1/12,Time=2,sigma=0.3, call=0, EU=0)

plot(vector_K, values_K_eu_call, type='l', xlab='K', ylab='Cena', main='Cena opcji europejskiej call')
plot(vector_K, values_K_eu_put, type='l', xlab='K', ylab='Cena', main='Cena opcji europejskiej put')
plot(vector_K, values_K_usa_call, type='l', xlab='K', ylab='Cena', main='Cena opcji amerykańskiej call')
plot(vector_K, values_K_usa_put, type='l', xlab='K', ylab='Cena', main='Cena opcji amerykańskiej put')

#test dla T

vector_T <- seq(0.5, 5, 0.5)

values_T_eu_call <- sapply(vector_T, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,sigma=0.3, call=1, EU=1)
values_T_eu_put <- sapply(vector_T, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,sigma=0.3, call=0, EU=1)
values_T_usa_call <- sapply(vector_T, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,sigma=0.3, call=1, EU=0)
values_T_usa_put <- sapply(vector_T, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,sigma=0.3, call=0, EU=0)

plot(vector_T, values_T_eu_call, type='l', xlab='T', ylab='Cena', main='Cena opcji europejskiej call')
plot(vector_T, values_T_eu_put, type='l', xlab='T', ylab='Cena', main='Cena opcji europejskiej put')
plot(vector_T, values_T_usa_call, type='l', xlab='T', ylab='Cena', main='Cena opcji amerykańskiej call')
plot(vector_T, values_T_usa_put, type='l', xlab='T', ylab='Cena', main='Cena opcji amerykańskiej put')

#test dla S0

vector_S0 <- seq(40, 60, 1)

values_S0_eu_call <- sapply(vector_S0, EOP_tree, Time=2,r=0.05,dt=1/12,K=48,sigma=0.3, call=1, EU=1)
values_S0_eu_put <- sapply(vector_S0, EOP_tree, Time=2,r=0.05,dt=1/12,K=48,sigma=0.3, call=0, EU=1)
values_S0_usa_call <- sapply(vector_S0, EOP_tree, Time=2,r=0.05,dt=1/12,K=48,sigma=0.3, call=1, EU=0)
values_S0_usa_put <- sapply(vector_S0, EOP_tree, Time=2,r=0.05,dt=1/12,K=48,sigma=0.3, call=0, EU=0)

plot(vector_S0, values_S0_eu_call, type='l', xlab=expression('S'['0']), ylab='Cena', main='Cena opcji europejskiej call')
plot(vector_S0, values_S0_eu_put, type='l', xlab=expression('S'['0']), ylab='Cena', main='Cena opcji europejskiej put')
plot(vector_S0, values_S0_usa_call, type='l', xlab=expression('S'['0']), ylab='Cena', main='Cena opcji amerykańskiej call')
plot(vector_S0, values_S0_usa_put, type='l', xlab=expression('S'['0']), ylab='Cena', main='Cena opcji amerykańskiej put')

#test dla sigma

vector_sigma <- seq(0.1, 2, 0.1)

values_sigma_eu_call <- sapply(vector_sigma, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,Time=2, call=1, EU=1)
values_sigma_eu_put <- sapply(vector_sigma, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,Time=2, call=0, EU=1)
values_sigma_usa_call <- sapply(vector_sigma, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,Time=2, call=1, EU=0)
values_sigma_usa_put <- sapply(vector_sigma, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,Time=2, call=0, EU=0)

plot(vector_sigma, values_sigma_eu_call, type='l', xlab=expression(sigma), ylab='Cena', main='Cena opcji europejskiej call')
plot(vector_sigma, values_sigma_eu_put, type='l', xlab=expression(sigma), ylab='Cena', main='Cena opcji europejskiej put')
plot(vector_sigma, values_sigma_usa_call, type='l', xlab=expression(sigma), ylab='Cena', main='Cena opcji amerykańskiej call')
plot(vector_sigma, values_sigma_usa_put, type='l', xlab=expression(sigma), ylab='Cena', main='Cena opcji amerykańskiej put')

#test dla r

vector_r <- seq(0.01, 0.2, 0.01)

values_r_eu_call <- sapply(vector_r, EOP_tree, S0=50,Time=2,dt=1/12,K=48,sigma=0.3, call=1, EU=1)
values_r_eu_put <- sapply(vector_r, EOP_tree, S0=50,Time=2,dt=1/12,K=48,sigma=0.3, call=0, EU=1)
values_r_usa_call <- sapply(vector_r, EOP_tree, S0=50,Time=2.05,dt=1/12,K=48,sigma=0.3, call=1, EU=0)
values_r_usa_put <- sapply(vector_r, EOP_tree, S0=50,Time=2,dt=1/12,K=48,sigma=0.3, call=0, EU=0)

plot(vector_r, values_r_eu_call, type='l', xlab='r', ylab='Cena', main='Cena opcji europejskiej call')
plot(vector_r, values_r_eu_put, type='l', xlab='r', ylab='Cena', main='Cena opcji europejskiej put')
plot(vector_r, values_r_usa_call, type='l', xlab='r', ylab='Cena', main='Cena opcji amerykańskiej call')
plot(vector_r, values_r_usa_put, type='l', xlab='r', ylab='Cena', main='Cena opcji amerykańskiej put')

#test dla dt

vector_dt <- seq(0.05, 1, 0.05)

values_dt_eu_call <- sapply(vector_dt, EOP_tree, S0=50,r=0.05,K=52,Time=2,sigma=0.3, call=1, EU=1)
values_dt_eu_put <- sapply(vector_dt, EOP_tree, S0=50,r=0.05,K=52,Time=2,sigma=0.3, call=0, EU=1)
values_dt_usa_call <- sapply(vector_dt, EOP_tree, S0=50,r=0.05,K=52,Time=2,sigma=0.3, call=1, EU=0)
values_dt_usa_put <- sapply(vector_dt, EOP_tree, S0=50,r=0.05,K=52,Time=2,sigma=0.3, call=0, EU=0)

plot(vector_dt, values_dt_eu_call, type='l', xlab='dt', ylab='Cena', main='Cena opcji europejskiej call')
plot(vector_dt, values_dt_eu_put, type='l', xlab='dt', ylab='Cena', main='Cena opcji europejskiej put')
plot(vector_dt, values_dt_usa_call, type='l', xlab='dt', ylab='Cena', main='Cena opcji amerykańskiej call')
plot(vector_dt, values_dt_usa_put, type='l', xlab='dt', ylab='Cena', main='Cena opcji amerykańskiej put')
