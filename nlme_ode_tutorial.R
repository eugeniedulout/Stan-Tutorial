install.packages(c('dplyr', 'tidyverse','cmdstan', 'ggplot2','gridExtra','bayesplot','posterior','deSolve','msm','magrittr','dplyr')) #Install all packages

library(dplyr)    #Load all packages
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(bayesplot)
library(cmdstanr)
library(posterior)
library(deSolve)
library(msm)
library(magrittr)
library(dplyr)
options(mc.cores = parallel::detectCores())

set.seed(123)

# I - Data Simulation

## Simulate the ODE System

ode.model <- function (t, x, params) {
  # State variables
  S <- x[1]
  I <- x[2]
  V <- x[3]
  # Parameters
  alpha_S <- params["alpha_S"]
  delta_S <- params["delta_S"]
  beta <- params["beta"]
  K <- params["K"]
  h <- params["h"]
  pi_ <- params["pi_"]
  gamma <- 23
  # Equations 
  dSdt <- alpha_S-delta_S*S-beta*S*V
  dIdt <- beta*S*V- K*I*I^(h)
  dVdt <- pi_*I-gamma*V-beta*S*V
  # Return
  dxdt <- c(dSdt,dIdt,dVdt)
  list(dxdt)
}

## Sample the parameters for each subject

N = 30 # number of subjects
Ttot = 90 # number of days
all_params = matrix(0, nrow = N, ncol = 7) # matrix containing all final parameters



# Set the hierarchical means 
mu_alpha_S <- 63.7
mu_delta_S <- 0.000751
mu_lbeta <-  -4.18
mu_k <- 0.356
mu_h <- 0.148
mu_lpi_ <- 1.09

# Set the hierarchical variances 
tau_alpha_S <- 0.1
tau_delta_S <- 0.0023
tau_lbeta <- 0.1
tau_k <- 0.1
tau_h <- 0.01
tau_lpi_ <- 0.1

for (i in 1:N){
  
  # Sample each parameter
  params <- list()
  params[1] <- rnorm(1,mu_alpha_S,tau_alpha_S)
  params[2] <- rtnorm(1,mu_delta_S,tau_delta_S, lower=0)
  params[3] <- 10^(rnorm(1,mu_lbeta,tau_lbeta))
  params[4] <- rtnorm(1,mu_k,tau_k, lower=0)
  params[5] <- rtnorm(1,mu_h,tau_h, lower=0)
  params[6] <- 10^(rnorm(1,mu_lpi_,tau_lpi_))
  
  
  
  # Add to matrix containing all parameters
  
  for (j in 1:6){
    all_params[i,j] <- all_params[i,j]+ as.numeric(params[j][1])
  }
  all_params[i,7] = 23 #gamma is not variable between individuals
}

## Simulate the data for each subject and time point

ytot = matrix(0, nrow = Ttot, ncol = N) # matrix containing all the simulated data for the viral load

for (j in 1:N){ 
  individual_params <- c(alpha_S=all_params[j,1], # get the individual parameters for each subject
                         delta_S=all_params[j,2], 
                         beta=all_params[j,3],
                         K=all_params[j,4],
                         h=all_params[j,5],
                         pi_=all_params[j,6], 
                         gamma=all_params[j,7])
  
  times <- seq(from=1,to=Ttot,by=1) # time range for the simulation
  xstart <- c(S=8.388815e+04,I=1.869510e-02,V=1.0e-02) # initial state of the system
  
  
  ode( # using the R general solver for ODEs
    func=ode.model,
    y=xstart,
    times=times,
    parms=individual_params
  ) %>%
    as.data.frame() -> out 
  
  for (i in 1:Ttot){
    ytot[i,j] <- ytot[i,j]+ as.numeric(out$V[i]) # add the viral load for each subject and each time point
  }
}

## Add noise

sigma_noise = 0.1
ytot_test = matrix(0, nrow = Ttot, ncol = N)
for (j in 1:N){
  for (i in 2:Ttot){
    pos = 0
    while (pos==0){
      ytot_test[i,j] <- 10^(log10(ytot[i,j]) + rnorm(1,0,sigma_noise))
      if(ytot_test[i,j]>0.01){
        pos=1
      }
    }
    ytot[i,j] <- ytot_test[i,j]
  }
}

## Clean simulated data

time <- list(1:Ttot)
names(time)[1] <- "days"
RV217_sim <- as.list(as.data.frame(ytot))
RV217_sim <- append(RV217_sim, time, 0)
RV217_simdf <- as.data.frame(RV217_sim)

## Save simulated data in an external file

write.table(RV217_simdf, file = "RV217_simulated_data", append = FALSE, sep = ",", dec = ".",
            row.names = TRUE, col.names = TRUE)

## Plot the simulated data

days = list()
VL = list()
ID = list()

for (i in 1:N){
  for (j in 1:Ttot){
    days <- append(days, RV217_sim[[1]][j])
    VL <- append(VL,RV217_sim[[i+1]][j])
    ID <- append(ID, i)
  }
}

days <- as.numeric(unlist(days))
VL <- as.numeric(unlist(VL))
ID <- as.numeric(unlist(ID))


RV217_sim_unpivot = list(days, VL, ID)
names(RV217_sim_unpivot)[1] <- "days" 
names(RV217_sim_unpivot)[2] <- "VL" 
names(RV217_sim_unpivot)[3] <- "ID" 
RV217_sim_unpivot <- data.frame(RV217_sim_unpivot)

plt <- ggplot(RV217_sim_unpivot,aes(y = log10(1000*VL),x = days,group = ID, col = VL)) + geom_line()
plt +scale_color_gradient(low="blue", high="red")


# II - Stan Inference

## Read the data


data <- read.table("RV217_simulated_data",sep=',', header=TRUE)


## Time optimization

data_red=data[c(1,2,7,14,21,28,35,42,49,56,63,70,77,84),]


## Format our data for Stan use

n_sub <- 30 #number of subjects
T <- length(data_red$days)-1 #number of time points
y0_V <- as.numeric(data_red[data_red$days == 1, -1]) #initial value
t0 <- 1 #first time point
times <- as.numeric(data_red$days)[-1] #vector of all time points without first one
y <- data_red[data_red$days > 1, -1] #vector of all values without first one
model_data_red <- list(n_sub = n_sub, T = T, y = y, y0_V = y0_V, times = times, t0 = t0)
print(model_data_red)


## Inference and optimization

mod_CH_fixed <- cmdstan_model("Reeves_model_fixed.stan")
mod_CH_fixed$exe_file()
fit_CH_red <- mod_CH_fixed$sample(data = model_data_red, chains = 2, num_warmup =250, num_samples  =200, cores=2)

# III - Fit Analysis

fit_CH_red$time()
print(fit_CH_red, max_rows = 200, digits = 5)
fit_CH_red$draws()

stanfit <- rstan::read_stan_csv(fit_CH_red$output_files())
stan_trace(stanfit)

mcmc_trace(fit_CH_red$draws("k"), pars = c("k[27,1]")) #trace plots
mcmc_dens_overlay(fit_CH_red$draws("alpha_S"), pars = c("alpha_S[1,1]","alpha_S[2,1]","alpha_S[3,1]","alpha_S[4,1]")) #separates the Markov chains
mcmc_areas(fit_CH_red$draws("logbeta"), pars = c("logbeta[27,1]")) #uncertainty intervals as shaded areas under the estimated posterior density curves




