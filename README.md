# Simulation_SAS
The purpose of this project is to investigate what happens to the estimation of fixed effect parameters, the standard error estimates of the fixed effect parameters and the Type I and
Type II error for the fixed effects hypothesis testing when the correlation structure is over- and under-specified. 

Fitting data with over and under specified variance covariance stuctures tend to be a common mistake that researchers can make. Hence, it is insightful to see what errors researchers have to deal with, if and when that happens. 

A simulation exercise is the best possible way to look into this. Here, I have started by simulating data with a compound symmetry (random intercept) correlation structure, an auto-regressive process of order one (AR(1)) correlation structure and an unstructured (UN) correlation structure. Then, I have fitted models to this simulated data with different correlation structures. This helps to compare the results of different case scenarios. For now I have only fitted models with an Uncorrelated (UC), Random Intercept (CS) and Unstructured (UN) Correlation structure. 

There have been many similar simulation exercises done using R. But as of my knowledge, I have not seen much similar simulation exercises being conducted using SAS. Therefore, in this project, I am mainly using SAS to achieve this. This started as a class project for one of my Bio-Statistics class but I have and will continue to update this and add more case scenarios as I work on this further. Actually, no one in my class used SAS for this other than myself and that is one of key reasons why, I felt like sharing this with anyone who is interested to do similar exercises in SAS.

Also, I did find a few important articles with regards to this most notably written by Dr. Rick Wicklin, who is a distinguished researcher in computational statistics at SAS and is a principal developer of PROC IML and SAS/IML Studio. I have referenced all of these articles in my report.  

This will be an ongoing exercise as I simulate data with other correlation structures and then fit models to these simulated data with different correlation structures to analyze the subsequent results. The code, results and corresponding report for this will also be updated accordingly.

 


