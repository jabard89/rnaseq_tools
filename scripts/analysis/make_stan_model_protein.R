make_stan_model <- function(column_list) {
  model <- paste0(
    '// -*- mode: C -*-
  data {
    // Number of proteins
    int<lower=1> NPROT;     
    
    // Note: These are all real numbers
    // columns tot, sup, pellet
    real<lower=1> tot_obs[NPROT];\n',
    paste(paste0("\treal<lower=1> ",column_list,"_obs[NPROT];\n"),collapse=""),
    '}
  parameters {\n',
    paste(paste0("\treal<lower=0> mixing_",column_list,";\n"),collapse=""),
    '  
    // dispersion parameter for counts
    real<lower=0> sigma;
  }
  model{
    // mixing ratios\n',
    paste(paste0("\tmixing_",column_list," ~gamma(1,1);\n"),collapse=""),
    '
    // Cauchy prior for negbin dispersion parameter
    sigma ~ cauchy(0,3);
    
    for(idx in 1:NPROT){ 
      // count dist log-normal with specified means
      // Total
      log(tot_obs[idx]) ~ normal(log(\n',
    paste("\t",paste0("mixing_",column_list[1:length(column_list)-1]," * ",
                      column_list[1:length(column_list)-1],"_obs[idx] + "),collapse=""),
    paste0("mixing_",column_list[length(column_list)]," * ",
           column_list[length(column_list)],"_obs[idx]),sigma);\n"),
    ' }
  
  }')
  return(model)
}