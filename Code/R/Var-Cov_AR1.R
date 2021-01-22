# Computing the AR(1) variance-covariance structure based on the assumed values of the parameters 

#set seed for reproducibility
set.seed(4747)
#number datasets
iter<-10000
#list for simulated dataframe
sim.df.ar<- list()
#number of participants
n<-20
#participants per group
n1<-5
#number of times
p<-5
# generate AR(1) covariance matrix
#matrix of distance between times
dist<- 1:p 
H<- abs(outer(dist, dist, "-"))
#set AR(1) parameters
#sigma^2e = sigma^2z/ (1- phi^2)
#we know sigma^2z =1 and phi^2=0.6
#thus, sigma^2e = 1/(1-0.6) = 2.5
#this is the same as 1 + 1.5 from CS, where 1 is from sigma^2e and 1.5 is from sigma^2b
sig2<- 2.5
phi<-sqrt(0.6)
#variance covariance matrix
V<- sig2*phi^H
V
