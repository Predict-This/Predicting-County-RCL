# Predicting-County-RCL
We merged data from the 2010 Census, 2012 County Presidential Data from the MIT Elections Lab, 
and Small Area Estimates from the National Surveys on Drug Use and Health (NSDUH) from 2010 to 2012 at the county level. 
This data is cleaned and avaialble in the final_data.rds file

The rcl_glm_lag_2014.R file uses this data and stacks an ensemble of logistic regressions fit to random realizations of 
the observed data using a minority oversampling technique. One county where recreational cannabis could be legally sold in 2014 is
sampled for every two counties where it could not. County-level probability estimates are weighted averages of predictions from each logistic model. 
The optimal cut-off can be varied to either minimize the false positive (Type I error) or false negative (Type II error) rate. By default,
accuracy weighted probabilities are used with a cut-off of .333 (equal to the prevalence of legalzied counties within each resample).
