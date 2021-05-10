#Ktoś może od nas kupić
EOP_call <- function(S0,r,K,Time,dt,sigma){
  
  u <- exp(sigma*sqrt(dt))
  d <- exp(-sigma*sqrt(dt))
  n <- Time/dt
  p <- (exp(r*dt)-d)/(u-d)
  
  Suma <- 0
  
  for(j in 0:n){
    Suma <- Suma + choose(n,j)*p^j*(1-p)^(n-j)*max(S0*u^j*d^(n-j)-K,0)      
  }
  
  return(exp(-r*Time)*Suma)
}

EOP_call(20,0.12,21,24,1,0.3)

#Ktoś może nam sprzedać
EOP_put <- function(S0,r,K,Time,dt,sigma){
  
  u <- exp(sigma*sqrt(dt))
  d <- exp(-sigma*sqrt(dt))
  n <- Time/dt
  p <- (exp(r*dt)-d)/(u-d)
  
  Suma <- 0
  
  for(j in 0:n){
    Suma <- Suma + choose(n,j)*p^j*(1-p)^(n-j)*max(K-S0*u^j*d^(n-j),0)      
  }
  
  return(exp(-r*Time)*Suma)
}

EOP_put(20,0.12,21,24,1,0.3)


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
  
  
  #Cu <-
  #Cd <- 
  #vertex <- Vertex_Value(p,r,dt,Cu,Cd)
}
EOP_tree(call=0, EU=0)

###test dla dt


vector_dt <- seq(0.05, 1, 0.05)
sapply(vector_dt, EOP_tree, S0=50,r=0.05,K=52,Time=2,sigma=0.3, call=0, EU=0)



