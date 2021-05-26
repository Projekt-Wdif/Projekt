Graf<- function(u,d,n=25){
  A<-matrix(0,n,n)
  for(i in 1:n){
    A[1:i,i]<-u^(seq(i,-i+2,-2)-1)
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
  
  #tree <- matrix(rep(0,625), nrow = 25)
  #first_line <- c(1, cumprod(c(rep(u,24))))
  #work_tree <- 
  
  graph <- (Graf(u,d,n) * S0 - K)*(-1)^(call + 1)
  graph[graph < 0] <- 0
  
  Cofanie <- matrix(0,n,n)
  Cofanie[,n] <- graph[,n]
  
  if (EU == 1) {
    for(j in n:2){
      for(k in 1:(j-1)){
        Cofanie[k,j-1] <- Vertex_Value(p,r,dt,Cofanie[k,j],Cofanie[k+1,j])
      }
    }
  }
  else {
    for(j in n:2){
      for(k in 1:(j-1)){
        Cofanie[k,j-1] <- max(Vertex_Value(p,r,dt,Cofanie[k,j],Cofanie[k+1,j]), graph[k,j-1])
      }
    }
  }
  
  
  return(Cofanie[1,1])
}

EOP_tree(call=0, EU=0)

#Opłacalność wykonania opcji amerykańskiej
Worth_American <- function(S0=50,r=0.05,K=52,Time=2,dt=1,sigma=0.3, call=1){

  u <- exp(sigma*sqrt(dt))
  d <- exp(-sigma*sqrt(dt))
  n <- Time/dt+1
  p <- (exp(r*dt)-d)/(u-d)
  
  graph <- (Graf(u,d,n) * S0 - K)*(-1)^(call + 1)
  graph[graph < 0] <- 0
  
  Cofanie <- matrix(0,n,n)
  Cofanie[,n] <- graph[,n]
  
  Worth <- matrix(0,n,n)
  Worth[,n] <- (graph[,n]>0)
  
  #wiemy że opcję call opłaca sie wykonać tylko na końcu
  #Dlatego wystarczy spojrzeć na kolumne ostatnich momentów
  if(call == 1){
    return(Worth)
  }
  
  for(j in n:2){
    for(k in 1:(j-1)){
      v1 <- Vertex_Value(p,r,dt,Cofanie[k,j],Cofanie[k+1,j])
      v2 <- graph[k,j-1]
      Cofanie[k,j-1] <- max(v1,v2)
      if(v2>v1){
        Worth[k,j-1] <- 1
      }
    }
  }
  
  return(Worth)
}

###### Delta i alpha

delta_akcji <- function(Vb, Vc, Stb, Stc) {
  return((Vb - Vc)/(Stb - Stc))
}

alpha_gotowki <- function(delta_akcji, Stc, r, dt) {
  return(-delta_akcji * Stc * exp(-r * dt))
}


portfel <- function(S0=20, dt=3/12, Time=1/2, K=21, r=0.12, sigma=0.3, EU=1, call=1, akcje=1) {
  
  u <- exp(sigma*sqrt(dt))
  d <- exp(-sigma*sqrt(dt))
  #u = 1.1
  #d = 0.9
  n <- Time/dt+1
  p <- (exp(r*dt)-d)/(u-d)
  
  macierz <- Graf(u, d, n) * S0
  
  graph <- (Graf(u, d, n) * S0 - K)*(-1)^(call + 1)
  graph[graph < 0] <- 0
  
  Cofanie <- matrix(0,Time/dt+1,Time/dt+1)
  Cofanie[,Time/dt+1] <- graph[,Time/dt+1]
  
  M_akcje <- matrix(0, n, n)
  M_gotowki <- matrix(0, n, n)
  if (EU == 1) {
    for(j in n:2){
      for(k in 1:(j-1)){
        Cofanie[k,j-1] <- Vertex_Value(p,r,dt,Cofanie[k,j],Cofanie[k+1,j])
        M_akcje[k,j-1] <- delta_akcji(Cofanie[k,j], Cofanie[k+1,j], macierz[k,j], macierz[k+1,j]) 
        M_gotowki[k,j-1] <- alpha_gotowki(delta_akcji(Cofanie[k,j], Cofanie[k+1,j], macierz[k,j], macierz[k+1,j]), macierz[k+1,j], r, dt)
      }
    }
  }
  else {
    for(j in n:2){
      for(k in 1:(j-1)){
        Cofanie[k,j-1] <- max(Vertex_Value(p,r,dt,Cofanie[k,j],Cofanie[k+1,j]), graph[k,j-1])
        M_akcje[k,j-1] <- delta_akcji(Cofanie[k,j], Cofanie[k+1,j], macierz[k,j], macierz[k+1,j])
        M_gotowki[k,j-1] <- alpha_gotowki(delta_akcji(Cofanie[k,j], Cofanie[k+1,j], macierz[k,j], macierz[k+1,j]), macierz[k+1,j], r, dt)
      }
    }
  }
  if(akcje == 1){
    return(M_akcje)
  }
  else {
    return(M_gotowki)
  }

}
akcje <- portfel(50, 1/12, 2, 48, 0.02, 0.3, 1, 1, 1)
gotowka <- portfel(50, 1/12, 2, 48, 0.02, 0.3, 1, 1, 2)
#akcje


for (i in 1:sqrt(length(akcje))) {
  for (j in 1:sqrt(length(akcje))) {
    if(akcje[i,j] == 0) {
      akcje[i,j] <- NA
    }
  }
}
#akcje
for (i in 1:sqrt(length(gotowka))) {
  for (j in 1:sqrt(length(gotowka))) {
    if(gotowka[i,j] == 0) {
      gotowka[i,j] <- NA
    }
  }
}
#gotowka

X <- melt(akcje, na.rm = TRUE)
ggplot(X, aes(x = Var2, y = Var1)) + 
  geom_point(aes(col=value), size=6) +
  xlab("X") + ylab("Y")
#X

Y <- melt(gotowka, na.rm = TRUE)
ggplot(Y, aes(x = Var2, y = Var1)) + 
  geom_point(aes(col=value), size=6) +
  xlab("X") + ylab("Y")

akcje <- portfel(50, 1/12, 2, 48, 0.02, 0.3, 1, 0, 1)
gotowka <- portfel(50, 1/12, 2, 48, 0.02, 0.3, 1, 0, 2)
#akcje


for (i in 1:sqrt(length(akcje))) {
  for (j in 1:sqrt(length(akcje))) {
    if(akcje[i,j] == 0) {
      akcje[i,j] <- NA
    }
  }
}
#akcje
for (i in 1:sqrt(length(gotowka))) {
  for (j in 1:sqrt(length(gotowka))) {
    if(gotowka[i,j] == 0) {
      gotowka[i,j] <- NA
    }
  }
}
#gotowka

X <- melt(akcje, na.rm = TRUE)
ggplot(X, aes(x = Var2, y = Var1)) + 
  geom_point(aes(col=value), size=6) +
  xlab("X") + ylab("Y")
#X

Y <- melt(gotowka, na.rm = TRUE)
ggplot(Y, aes(x = Var2, y = Var1)) + 
  geom_point(aes(col=value), size=6) +
  xlab("X") + ylab("Y")

###################### Testy wrażliwości:

#test dla K

vector_K <- seq(40, 60, 1)

values_K_eu_call <- sapply(vector_K, EOP_tree, S0=50,r=0.05,dt=1/12,Time=2,sigma=0.3, call=1, EU=1)
values_K_eu_put <- sapply(vector_K, EOP_tree, S0=50,r=0.05,dt=1/12,Time=2,sigma=0.3, call=0, EU=1)
values_K_usa_call <- sapply(vector_K, EOP_tree, S0=50,r=0.05,dt=1/12,Time=2,sigma=0.3, call=1, EU=0)
values_K_usa_put <- sapply(vector_K, EOP_tree, S0=50,r=0.05,dt=1/12,Time=2,sigma=0.3, call=0, EU=0)

# grid <- matrix(c(1,2,3,4) , ncol=2, nrow=2, byrow = T)
# layout(grid)
# 
# plot(vector_K, values_K_eu_call, type='l', xlab='K', ylab='Cena', main='Cena opcji europejskiej call')
# plot(vector_K, values_K_eu_put, type='l', xlab='K', ylab='Cena', main='Cena opcji europejskiej put')
# plot(vector_K, values_K_usa_call, type='l', xlab='K', ylab='Cena', main='Cena opcji amerykanskiej call')
# plot(vector_K, values_K_usa_put, type='l', xlab='K', ylab='Cena', main='Cena opcji amerykanskiej put')
# 
# grid <- matrix(c(1,2,3,4) , ncol=2, nrow=2, byrow = T)
# layout(grid)

plot(vector_K, values_K_eu_call, type = "o", col = 3, pch = 1, lty = 1, lwd = 2, cex = 1, xlab='K', ylab='Cena', main='Cena opcji w zaleznosci od ceny wykonania K', ylim = c(0,17))
lines(vector_K, values_K_eu_put, type = "o", col = 2, pch = 2, lty = 1, lwd = 2, cex = 1)
lines(vector_K, values_K_usa_put, type = "o", col = 4, pch = 4, lty = 1, lwd = 2, cex = 1)

par(xpd = TRUE)

legend(39.8,11, c("EU/AM call", "EU put", "AM put"), pch = c(1,2,4), col = c(3,2,4), lty = 1, bg = "white", pt.bg = "white", cex = 1, lwd = 2)




#test dla T

vector_T <- seq(0.5, 5, 0.5)

values_T_eu_call <- sapply(vector_T, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,sigma=0.3, call=1, EU=1)
values_T_eu_put <- sapply(vector_T, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,sigma=0.3, call=0, EU=1)
values_T_usa_call <- sapply(vector_T, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,sigma=0.3, call=1, EU=0)
values_T_usa_put <- sapply(vector_T, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,sigma=0.3, call=0, EU=0)

# plot(vector_T, values_T_eu_call, type='l', xlab='T', ylab='Cena', main='Cena opcji europejskiej call')
# plot(vector_T, values_T_eu_put, type='l', xlab='T', ylab='Cena', main='Cena opcji europejskiej put')
# plot(vector_T, values_T_usa_call, type='l', xlab='T', ylab='Cena', main='Cena opcji amerykanskiej call')
# plot(vector_T, values_T_usa_put, type='l', xlab='T', ylab='Cena', main='Cena opcji amerykanskiej put')
# 
# grid <- matrix(c(1,2,3,4) , ncol=2, nrow=2, byrow = T)
# layout(grid)

plot(vector_T, values_T_eu_call, type = "o", col = 3, pch = 1, lty = 1, lwd = 2, cex = 1, xlab='T', ylab='Cena', main='Cena opcji w zaleznosci od zapadalnosci T', ylim = c(0,20))
lines(vector_T, values_T_eu_put, type = "o", col = 2, pch = 2, lty = 1, lwd = 2, cex = 1)
lines(vector_T, values_T_usa_put, type = "o", col = 4, pch = 4, lty = 1, lwd = 2, cex = 1)

par(xpd = TRUE)

legend(0.4,20, c("EU/AM call", "EU put", "AM put"), pch = c(1,2,4), col = c(3,2,4), lty = 1, bg = "white", pt.bg = "white", cex = 1, lwd = 2)



#test dla S0

vector_S0 <- seq(40, 60, 1)

values_S0_eu_call <- sapply(vector_S0, EOP_tree, Time=2,r=0.05,dt=1/12,K=48,sigma=0.3, call=1, EU=1)
values_S0_eu_put <- sapply(vector_S0, EOP_tree, Time=2,r=0.05,dt=1/12,K=48,sigma=0.3, call=0, EU=1)
values_S0_usa_call <- sapply(vector_S0, EOP_tree, Time=2,r=0.05,dt=1/12,K=48,sigma=0.3, call=1, EU=0)
values_S0_usa_put <- sapply(vector_S0, EOP_tree, Time=2,r=0.05,dt=1/12,K=48,sigma=0.3, call=0, EU=0)

# plot(vector_S0, values_S0_eu_call, type='l', xlab=expression('S'['0']), ylab='Cena', main='Cena opcji europejskiej call')
# plot(vector_S0, values_S0_eu_put, type='l', xlab=expression('S'['0']), ylab='Cena', main='Cena opcji europejskiej put')
# plot(vector_S0, values_S0_usa_call, type='l', xlab=expression('S'['0']), ylab='Cena', main='Cena opcji amerykanskiej call')
# plot(vector_S0, values_S0_usa_put, type='l', xlab=expression('S'['0']), ylab='Cena', main='Cena opcji amerykanskiej put')
# 
# grid <- matrix(c(1,2,3,4) , ncol=2, nrow=2, byrow = T)
# layout(grid)

plot(vector_S0, values_S0_eu_call, type = "o", col = 3, pch = 1, lty = 1, lwd = 2, cex = 1, xlab='spot', ylab='Cena', main='Cena opcji w zaleznosci od poczatkowej ceny spot', ylim = c(0,20))
lines(vector_S0, values_S0_eu_put, type = "o", col = 2, pch = 2, lty = 1, lwd = 2, cex = 1)
lines(vector_S0, values_S0_usa_put, type = "o", col = 4, pch = 4, lty = 1, lwd = 2, cex = 1)

par(xpd = TRUE)

legend(39.5,20, c("EU/AM call", "EU put", "AM put"), pch = c(1,2,4), col = c(3,2,4), lty = 1, bg = "white", pt.bg = "white", cex = 1, lwd = 2)



#test dla sigma

vector_sigma <- seq(0.1, 2, 0.1)

values_sigma_eu_call <- sapply(vector_sigma, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,Time=2, call=1, EU=1)
values_sigma_eu_put <- sapply(vector_sigma, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,Time=2, call=0, EU=1)
values_sigma_usa_call <- sapply(vector_sigma, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,Time=2, call=1, EU=0)
values_sigma_usa_put <- sapply(vector_sigma, EOP_tree, S0=50,r=0.05,dt=1/12,K=48,Time=2, call=0, EU=0)

# plot(vector_sigma, values_sigma_eu_call, type='l', xlab=expression(sigma), ylab='Cena', main='Cena opcji europejskiej call')
# plot(vector_sigma, values_sigma_eu_put, type='l', xlab=expression(sigma), ylab='Cena', main='Cena opcji europejskiej put')
# plot(vector_sigma, values_sigma_usa_call, type='l', xlab=expression(sigma), ylab='Cena', main='Cena opcji amerykanskiej call')
# plot(vector_sigma, values_sigma_usa_put, type='l', xlab=expression(sigma), ylab='Cena', main='Cena opcji amerykanskiej put')
# 
# grid <- matrix(c(1,2,3,4) , ncol=2, nrow=2, byrow = T)
# layout(grid)

plot(vector_sigma, values_sigma_eu_call, type = "o", col = 3, pch = 1, lty = 1, lwd = 2, cex = 1, xlab=expression(paste(sigma)), ylab='Cena', main=expression(paste('Cena opcji w zaleznosci od parametru ', sigma)), ylim = c(0,50))
lines(vector_sigma, values_sigma_eu_put, type = "o", col = 2, pch = 2, lty = 1, lwd = 2, cex = 1)
lines(vector_sigma, values_sigma_usa_put, type = "o", col = 4, pch = 4, lty = 1, lwd = 2, cex = 1)

par(xpd = TRUE)

legend(0.1,50, c("EU/AM call", "EU put", "AM put"), pch = c(1,2,4), col = c(3,2,4), lty = 1, bg = "white", pt.bg = "white", cex = 1, lwd = 2)




#test dla r

vector_r <- seq(0.01, 0.2, 0.01)

values_r_eu_call <- sapply(vector_r, EOP_tree, S0=50,Time=2,dt=1/12,K=48,sigma=0.3, call=1, EU=1)
values_r_eu_put <- sapply(vector_r, EOP_tree, S0=50,Time=2,dt=1/12,K=48,sigma=0.3, call=0, EU=1)
values_r_usa_call <- sapply(vector_r, EOP_tree, S0=50,Time=2.05,dt=1/12,K=48,sigma=0.3, call=1, EU=0)
values_r_usa_put <- sapply(vector_r, EOP_tree, S0=50,Time=2,dt=1/12,K=48,sigma=0.3, call=0, EU=0)

# plot(vector_r, values_r_eu_call, type='l', xlab='r', ylab='Cena', main='Cena opcji europejskiej call')
# plot(vector_r, values_r_eu_put, type='l', xlab='r', ylab='Cena', main='Cena opcji europejskiej put')
# plot(vector_r, values_r_usa_call, type='l', xlab='r', ylab='Cena', main='Cena opcji amerykanskiej call')
# plot(vector_r, values_r_usa_put, type='l', xlab='r', ylab='Cena', main='Cena opcji amerykanskiej put')
# 
# grid <- matrix(c(1,2,3,4) , ncol=2, nrow=2, byrow = T)
# layout(grid)

plot(vector_r, values_r_eu_call, type = "o", col = 3, pch = 1, lty = 1, lwd = 2, cex = 1, xlab='r', ylab='Cena', main='Cena opcji w zaleznosci od stopy wolnej od ryzyka r', ylim = c(0,22))
lines(vector_r, values_r_eu_put, type = "o", col = 2, pch = 2, lty = 1, lwd = 2, cex = 1)
lines(vector_r, values_r_usa_put, type = "o", col = 4, pch = 4, lty = 1, lwd = 2, cex = 1)

par(xpd = TRUE)

legend(0.005,22, c("EU/AM call", "EU put", "AM put"), pch = c(1,2,4), col = c(3,2,4), lty = 1, bg = "white", pt.bg = "white", cex = 1, lwd = 2)




#test dla dt

vector_dt <- c(1/1460,1/730,1/365,1/200,1/168,1/144,1/120,1/96,1/72,1/60,1/48,1/32,1/24,1/16,1/12,1/8,1/6,1/4,1/2,1)

values_dt_eu_call <- sapply(vector_dt, EOP_tree, S0=50,r=0.05,K=52,Time=2,sigma=0.3, call=1, EU=1)
values_dt_eu_put <- sapply(vector_dt, EOP_tree, S0=50,r=0.05,K=52,Time=2,sigma=0.3, call=0, EU=1)
values_dt_usa_call <- sapply(vector_dt, EOP_tree, S0=50,r=0.05,K=52,Time=2,sigma=0.3, call=1, EU=0)
values_dt_usa_put <- sapply(vector_dt, EOP_tree, S0=50,r=0.05,K=52,Time=2,sigma=0.3, call=0, EU=0)



plot(vector_dt, values_dt_eu_call, type = "o", col = 3, pch = 1, lty = 1, lwd = 2, cex = 1, xlab=expression(paste(Delta, 't (skala logarytmiczna)')), ylab='Cena', main=expression(paste('Cena opcji w zaleznosci od ', Delta, 't')), ylim = c(6,10.1), log = "x")
lines(vector_dt, values_dt_eu_put, type = "o", col = 2, pch = 2, lty = 1, lwd = 2, cex = 1)
lines(vector_dt, values_dt_usa_put, type = "o", col = 4, pch = 4, lty = 1, lwd = 2, cex = 1)

par(xpd = TRUE)

legend(0.0006,9.1, c("EU/AM call", "EU put", "AM put"), pch = c(1,2,4), col = c(3,2,4), lty = 1, bg = "white", pt.bg = "white", cex = 1, lwd = 2)

grid <- matrix(c(1,1,2,3) , ncol=2, nrow=2, byrow = T)
layout(grid)

plot(vector_dt, values_dt_eu_call, xlab=expression(paste(Delta, 't (skala logarytmiczna)')), ylab='Cena',type = "o", col = 3, pch = 1, lty = 1, lwd = 2, cex = 1, main=expression(paste('Cena EU/AM call w zaleznosci od ', Delta, 't')), log = "x")
plot(vector_dt, values_dt_eu_put, type = "o", col = 2, pch = 2, lty = 1, lwd = 2, cex = 1, xlab=expression(paste(Delta, 't (skala logarytmiczna)')), ylab='Cena', main=expression(paste('Cena EU put w zaleznosci od ', Delta, 't')), log = "x")
plot(vector_dt, values_dt_usa_put, type = "o", col = 4, pch = 4, lty = 1, lwd = 2, cex = 1, xlab=expression(paste(Delta, 't (skala logarytmiczna)')), ylab='Cena', main=expression(paste('Cena AM put w zaleznosci od ', Delta, 't')), log = "x")

