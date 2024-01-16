// Stan code for a mixed effect model with ODEs

functions{
  
    // The function block is for any user-defined function. It must be the first block that appears.
    // The ODE systems must be decsribed here with the following signature:
    
    vector ode_system(real t, vector y, real alpha_S, real delta_S, real logbeta, real K, real h, real logpi_){

        // Parameters
        
        // Write here the fixed parameters. For instance, in that model gamma = 23 for all subject.
        // Write here all the transformed parameters. For instance, the parameters beta and pi_ follow a lognormal distribution,
        // so we need to transform them before using them in the set of equations.

        int gamma = 23;
        real beta = 10^(logbeta);
        real pi_ = 10^(logpi_);
        
        // Equations
        
        vector[3] dydt;

        dydt[1] = alpha_S - delta_S*y[1] - beta*y[1]*y[3]; //S
        dydt[2] = beta*y[1]*y[3] - K*y[2]*y[2]^h;          //I
        dydt[3] = pi_*y[2]-gamma*y[3]-beta*y[1]*y[3];      //V
        
        return dydt;

    }
} 


data{
  
  //The data block is where all the required data for the model is described. It must be the same as the list used for the stan input.
  
    int<lower=0> n_sub;                //number of subjects
    int<lower=0> T;                    //number of time points
    array[T, n_sub] real<lower=-10> y; //vector of all values without first one
    array[n_sub] real<lower=0> y0_V;   //initial value for V
    array[T] real<lower=0> times;      //vector of all time points without first one
    real<lower=0> t0;                  //first time point
}
transformed data{
  
    array[n_sub] vector[3] y0; //vector of initial values for the ODE system
    
    for(n in 1:n_sub){
      
      y0[n,1]=8.388815e+04;  //S
      y0[n,2]=1.869510e-02;; //I
      y0[n,3]=y0_V[n];       //V
    }
}

parameters{
  
  //The parameters block declares the model’s parameters and their distribution.
  //As this is a mixed effect model, each parameter for each subject is sampled from a group sampled mean and a group sampled variance.
  //Here, all the parameters are separated in different arrays so that we can use lower and upper bounds if necessary,
  //to avoid divergence of the model after sampling.
  
  //alpha_S
    array[1] real<lower=50, upper = 70> mu_alpha_S;        //Mean for alpha_S
    array[1] real<lower=0> tau_alpha_S;     //Variance for alpha_S
    array[n_sub,1] real<lower=60, upper = 70> alpha_S;     //alpha_S, sampled below using mu_alpha_S and tau_alpha_S
    
  //delta_S 
    array[1] real<lower=0> mu_delta_S;                     //Mean for delta_S
    array[1] real<lower=0> tau_delta_S;     //Variance for delta_S
    array[n_sub,1] real<lower=0> delta_S;                  //delta_S
  
  //logbeta
    array[1] real<lower=-5, upper = -3> mu_logbeta;        //Mean for logbeta
    array[1] real<lower=0> tau_logbeta;     //Variance for logbeta
    array[n_sub,1] real<lower=-4.3, upper = -4.1> logbeta; //Logbeta
    
  //k
    array[1] real<lower=-0> mu_k;                          //Mean for k
    array[1] real<lower=0> tau_k;                          //Variance for k
    array[n_sub,1] real<lower=0> k;                        //k
    
  //h
    array[1] real<lower=-0> mu_h;                          //Mean for h
    array[1] real<lower=0> tau_h;                          //Variance for h
    array[n_sub,1] real<lower=0> h;                        //h
    
  //logpi_
    array[1] real<lower=-0> mu_logpi_;                     //Mean for logpi_
    array[1] real<lower=0> tau_logpi_;                     //Variance for logpi_
    array[n_sub,1] real<lower=0> logpi_;                   //logpi_
    
  //variance
    real<lower=0> mu_sigma;                               //group variance's mean
    real<lower=0> tau_sigma;                              //group variance's variance
    array[n_sub] real<lower=0> sigma;                     //group variance
}

model{
  
  //The model block is where the log-probability function is defined.
  //Here, at each step the ODE solver and the parameter estimation are applied simultaneously.
  
  
  // ODE solver
  array[T] vector[n_sub] yhat;
  
  for (n in 1:n_sub){
        array[T] vector[3] yhat_temp3;
        
        //Solve ode using one of stan's 2 solvers: ode_bdf or ode_rk45. They must respect the following signature.
        yhat_temp3 = ode_bdf(ode_system, y0[n,:], t0, times, alpha_S[n,1], delta_S[n,1], logbeta[n,1],k[n,1],h[n,1],logpi_[n,1]);
        for(t in 1:T){
            yhat[t,n]=yhat_temp3[t,3]; 
            //the ODE solver estimate all the variables of the model,
            //but we only fit on the variable for which we have data, here it's V
        }
    }
  
  //Priors
  
  //As we are using simulated data, our priors are well defined. In some cases, the priors might be unknown,
  //and the use of weakly informative priors is recommended (such as normal, student, or cauchy distributions).
  //The use of uniform distributions is possible but not recommended, unless the dataset is very large, as they are non-informative.
  
  //Hierarchical parameters sampled
    mu_alpha_S[1] ~ normal(63.7,1);
    mu_delta_S[1] ~ normal(0.000751,0.001);
    mu_logbeta[1] ~ normal(-4.18,1); 
    mu_k[1] ~ normal(0.356,0.5); 
    mu_h[1] ~ normal(0.148,0.05);
    mu_logpi_[1] ~ normal(1.09,0.5);
    

    tau_alpha_S[1] ~ normal(0,1);
    tau_delta_S[1] ~ normal(0.0023,0.1);
    tau_logbeta[1] ~ normal(0,1);
    tau_k[1] ~ normal(0,1);
    tau_h[1] ~ normal(0,1);
    tau_logpi_[1] ~ normal(0,1);
    

    mu_sigma ~ normal(0,1);  //group variance's mean
    tau_sigma ~ normal(0,1); //group variance's variance

    
    for (n in 1:n_sub){
      //Sample each parameter for each subject
          alpha_S[n,1] ~ normal(mu_alpha_S[1],tau_alpha_S[1]);
          delta_S[n,1] ~ normal(mu_delta_S[1],tau_delta_S[1]);
          logbeta[n,1] ~ normal(mu_logbeta[1],tau_logbeta[1]);
          k[n,1] ~ normal(mu_k[1],tau_k[1]);
          h[n,1] ~ normal(mu_h[1],tau_h[1]);
          logpi_[n,1] ~ normal(mu_logpi_[1],tau_logpi_[1]);
          
      //Sample the group variance for each subject  
          sigma[n] ~ normal(mu_sigma,tau_sigma); 
    }
    
    for(n in 1:n_sub){
        for(t in 1:T){
            y[t,n] ~ normal(yhat[t,n], sigma[n]);
        }
    }
}
