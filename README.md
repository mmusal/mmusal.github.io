#On this page we will discuss various steps in developing a Hidden Markov Model Spatial Model. 
First we need to discuss each model separately. 
Hidden Markov Models:
Assume there is a random variable Y, indexed by time t, $Y_{t}$. We would like to learn about the factors that effect $Y_{t}$. To do this we will assume a distribution over Y.  
To decide on the distribution we need to know about what it is Y actually represents. 
For the sake of this discussion Y will represent the biweekly mortality in California. 
Since we have information on each of California's counties, we will model the biweekly mortality $Y_{it}$ where $t=1 \cdots T=77$ and $i=1 \cdots N=58$.  
A discrete probability distribution which can be used to model mortality is Poisson. 

 $\Pr(Y{=}y)={\frac {\lambda ^{y}e^{-\lambda }}{y!}}$ where $\lambda$ is the mean and variance of the distribution. 
Once we have decided on the distribution of mortality, what separates out one distribution from the other is (are) the parameter(s) of the distribution. There are countless variables we can use instead of mortality but what changes uncertainty evaluation (probability) for Y is not what Y stands for but the parameter(s) used in the equation. 
For instance, in the equation for poisson distribution, the equation itself is fixed and you are trying to evaluate the probability that Y (what the random variable stands for) is y (the number itself). What changes the result (P(Y=y)) is really $\lambda$.   
To demonstrate the distribution of Y with different $\lambda$ values, we will use R:

library(ggplot2)
library(dplyr)

#values that lambda is going to take

lambda=c(1,3,5,10,20)

#how many variates we will generate

n=10000

#the two dimensional array that is going to be populated via the rpois 

#(randomly generate poisson distributions) function.

poissonvariates=as.data.frame(matrix(nrow=n,ncol=length(lambda)))

#Below is the for loop that is going to generate the variates. The number of times 
#that the statements in the curly bracket is going to run is the length of the lambda vector (5).
#the rpois has 2 required parameters that we need to provide values for; n and lambda. 
#n is the number of values we will generate, lambda is the mean/variance of the distribution. 

for(i in 1:length(lambda)){
poissonvariates[,i]=rpois(n=n,lambda[i])}

####give names to the dataframe columns
names(poissonvariates)=c("lambda=1","lambda=3","lambda=5","lambda=10","lambda=20")

#we need to stack the columns for ggplot function

stackedpoisson=stack(poissonvariates)

#once the columns are stacked they have the names: values ind

#We need to summarize the data to create the barchart.

stackedpoisson2 <- stackedpoisson %>% 
  group_by(values, ind) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/n) 
  
#and here is the plot that represents poisson variates simulated with 

#different lambda values

ggplot(data=stackedpoisson2, aes(x=values,y=perc, fill=ind))+
  geom_bar(stat="identity")+
  facet_grid(. ~ ind)+
  labs(x = "Sim. Values", y = "Percent", fill = "Lambda")
  
  
![Rplot](https://github.com/mmusal/mmusal.github.io/assets/11746560/32b2a741-5c22-4a1c-8898-d7d3b8b6e442)
