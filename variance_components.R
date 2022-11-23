library(dplyr)
library(BalancedSampling)

load("results_43.Rdata")
start_date <- 43
mse_fpps <- (0.5*y$y)^2
mse_lp <- (0.47*y$y)^2
mse_lc_v <- (0.42*y$y)^2
mse_lc <- (0.41*y$y)^2


map <- map05
N <- 400
Ni <- map$map$p
n <- 80
pi <- map$map$p/20000*n
X <- cbind(map$map$x, map$map$y)


# pi_il
u <- dataset_1s$pred - dataset_1s$inf_rate
mean(u)
var(u)
sigma_u <- var(u)

gen_sample <- function(x){
  k <- sample(N, n, prob = pi)
  return(as.numeric((1:N) %in% k))
}

simul <- 100000
simulation <- lapply(1:simul, gen_sample)
m <- matrix(data=unlist(simulation), ncol=N, nrow=simul, byrow = T)
pi_il <- (t(m) %*% m)/simul
rm(m, simulation)
dist <- as.matrix(dist(cbind(map$map$x[1:400], map$map$y[1:400]), ))

# rho

dist2 <- map05$data %>%
  filter(t == start_date)

u2 <- dataset$pred - dataset$inf_rate
sigma_u2 <- var(u2)
round(sigma_u2, 6)
dist2 <- as.matrix(dist(cbind(dist2$x, dist2$y), ))
u2_mat <- u2 %*% t(u2)
dist2 <- as.numeric(dist2)
u2_mat <- log(as.numeric(u2_mat)^2)
y <- u2_mat - 2*log(sigma_u)
x <- dist2*2

data2 <- data.frame(x=x, y = y)
data2_sub <- data2[sample(1:nrow(data2), 100000), ]
coef <- coef(lm(data = data2_sub, y~x-1))
rho <- as.numeric(exp(coef))
rho
rho <- 1
rho
# componenti

A <- sigma_u2*(20000/80)^2*sum(pi_il*rho^dist*((diag(1, 400)*-1)+1))
B <- mse_lp - mse_lc
E <- sum(1/pi*Ni^2/3*sigma_u2)
D <- sum(1/pi*Ni*((Ni-3)/3)*pred_sigma)
K <- mse_fpps-(A + B + D + E)
rho2 <- 1 
C <- sigma_u2*sum(1/pi*Ni^2/3*(3-1)*rho2)
VM <- C - K


(A + B + C + D + E)
(A + B + C + D + E - VM)/(mse_fpps)
round(c(A, B, C, D, E)/(A+B+C+D+E), 2)
c(A, B, C, D, E, VM, rho, rho2)



rho
rho2

VM1 <- sum(Ni*sigma_u2)
VM2 <- sum(Ni*(Ni-1))*sigma_u2*rho2
VM3 <- (t(Ni) %*% rho^(dist*((diag(1, 400)*-1)+1)) %*% Ni)*sigma_u2
VM1  
VM2
VM3
VM1 + VM2 + VM3

c(VM1, E)
c(VM2, C)
c(VM3, A)
VM


A + B + C + D + E - VM
A + B + C + D + E - VM1 -VM2 - VM3
mse_fpps




# Ra
gen_sample_lpm2 <- function(x){
  k <- lpm2(pi,X)
  return(as.numeric((1:N) %in% k))
}

simulation <- lapply(1:simul, gen_sample_lpm2)
m_lpm2 <- matrix(data=unlist(simulation), ncol=N, nrow=simul, byrow = T)
pi_il_lpm2 <- (t(m_lpm2) %*% m_lpm2)/simul
rm(m_lpm2, simulation)

A_star <- sigma_u2*(20000/80)^2*sum(pi_il_lpm2*rho^dist*((diag(1, 400)*-1)+1))
R_a <- sigma_u2*(20000/80)^2*sum((pi_il-pi_il_lpm2)*rho^dist*((diag(1, 400)*-1)+1))

# Rb
B_star <- mse_lp - mse_lc_v
R_b <- B - B_star
c(R_a, R_b)

