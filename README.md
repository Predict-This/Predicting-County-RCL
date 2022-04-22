# Predicting-County-RCL
This project uses publicly available data from 2010, 2011, and 2012 to predict which counties would make recreational cannabis sales legal in 2014. 

We merged data from the 2010 Census, 2012 County Presidential Data from the MIT Elections Lab, 
and Small Area Estimates from the National Surveys on Drug Use and Health (NSDUH, 2010-2012) at the county level. 
This data is cleaned and avaialble in the final_data.rds file.

The rcl_glm_lag_2014.R file uses this data and stacks an ensemble of logistic regressions fit to random realizations of 
the observed data using a minority oversampling technique. One county where recreational cannabis could be legally sold in 2014 is
sampled for every two counties where it could not. County-level probability estimates are weighted averages of predictions from each logistic model. 
The optimal cut-off and weighting mechnaism can be varied to either minimize the false positive (Type I error) or false negative (Type II error) rate. By default,
probabilities are weighted by the overall accuracy of the model from which they were dervied with a cut-off of .333 (equal to the prevalence of legalzied counties within each resample).
