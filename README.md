# Simulation_SAS
The purpose of this project is to investigate what happens to the estimation of fixed effect parameters, the standard errors estimates of the fixed effect parameters and the Type I and Type II error for the fixed effects 
hypothesis testing when the correlation structure is over- and under-specified.

Fitting data with over and under specified variance covariance stuctures tends to be a common problem for researchers. Hence, it is insightful to see what errors we have to deal with, if and when that happens. 

A simulation study is the best possible way to look into this. Here, I will start by simulating data with a compound symmetry (random intercept) correlation structure, an auto-regressive process of order one (AR(1)) 
correlation structure and an unstructured (UN) correlation structure. Then fit models to this simulated data with different correlation structures. This helps to compare the results of different scenarios. 

There have been many similar simulation studies done using R. But as of my knowledge, there have been not been many simulation studies using SAS. 

In this project, I will mainly use SAS to achieve this. 


