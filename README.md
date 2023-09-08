#On this page we will discuss various steps in developing a Hidden Markov Model Spatial Model. 
First we need to discuss each model separately. These descriptions are going to be very shallow, which we will elaborate as we progress. 
Spatial Models:
A spatial model is based on the premise that things that are closer to each other are more similar than things that are further away. There are different frameworks on describing distances between things. Here we will assume that distances are described by a neighbourhood structure between counties of California. A county will effect the counties that it shares a boundry with and vice versa.   

Hidden Markov Models:
A hidden markov model assumes there is an unobserved variable that has an effect on the distribution of the observed variable. For instance in finance we can conceptualize an up/down trend for the stock market based on whether there is a bull or a bear market. In this illustration we will be utilizing up/down stages of the disease to inform us of the distribution of biweekly mortality. The unobsorved trend will have a probability structure that ideally we should further model. 

General Introduction:
Assume there is a random variable Y, indexed by time t, $Y_{t}$. We would like to learn about the factors that effect $Y_{t}$. To do this we will assume a distribution over Y.  
To decide on the distribution we need to know about what it is Y actually represents. 
For the sake of this discussion Y will represent the biweekly mortality in California. 
Since we have information on each of California's counties, we will model the biweekly mortality $Y_{it}$ where $t=1 \cdots T=77$ and $i=1 \cdots N=58$.  
A discrete probability distribution which can be used to model mortality is Poisson. 

 $\Pr(Y{=}y)={\frac {\lambda ^{y}e^{-\lambda }}{y!}}$ where $\lambda$ is the mean and variance of the distribution. 
Once we have decided on the distribution of mortality, what separates out one distribution from the other is (are) the parameter(s) of the distribution. There are countless variables we can use instead of mortality but what changes uncertainty evaluation (probability) for Y is not what Y stands for but the parameter(s) used in the equation. 
For instance, in the equation for poisson distribution, the equation itself is fixed and you are trying to evaluate the probability that Y (what the random variable stands for) is y (the number itself). What changes the result (P(Y=y)) is really $\lambda$.   
To demonstrate the distribution of Y with different $\lambda$ values, we will use R:
#Libraries to be used:
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

#different lambda values. 

ggplot(data=stackedpoisson2, aes(x=values,y=perc, fill=ind))+
  geom_bar(stat="identity")+
  facet_grid(. ~ ind)+
  labs(x = "Sim. Values", y = "Percent", fill = "Lambda")
  
  
![Rplot](https://github.com/mmusal/mmusal.github.io/assets/11746560/32b2a741-5c22-4a1c-8898-d7d3b8b6e442)

\\

The point of the plot is to show you how changing the value lambda, changes the uncertainty evaluations regarding the values of Y.

In this presentation we assume that the biweekly mortality in the counties of California is Poisson distributed. 
$Y_{it}\sim Poisson(\lambda_{it})$
To eloborate on what I mentioned earlier, this by itself is not really interesting. We can assign probabilities to values other than the ones we have seen at a particular time $t$ in a particular county $i$. However this would rely very heavily on the assumption of biweekly mortality being Poisson distributed and not be of particular use by itself.  

The goal of the project is to model the biweekly mortality in the counties of California in order to answer the following questions:
1) How did socio-demographic factors effect biweekly mortality
2) How did these change through the progression of the pandemic
3) Can we detect regions (counties or collection of counties) that performed better than others.
4) Did specific policies at county or state level effect the unobserved up/down trends.

In order to answer these questions we will need some background so that we can develop models to help answer those questions:
A) Data/tools needed for spatial models in R.
B) Some neccessary Probability concepts.
C) Spatial Models (We will use Besag York Mollie).
D) Hidden Markov Models.
E) Bayesian Modeling Tool (STAN).
6) Combine Bayesian Spatial and Hidden Markov Models  
