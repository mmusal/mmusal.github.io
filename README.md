# General Project Description:
#In this project we will discuss various steps in developing a Hidden Markov Model Spatial Model. 
First we need to discuss each model separately. These descriptions are going to be very shallow, which we will elaborate as we progress. 
## Spatial Models:
A spatial model is based on the premise that things that are closer to each other are more similar than things that are further away. There are different frameworks on describing distances between things. Here we will assume that distances are described by a neighbourhood structure between counties of California. A county will effect the counties that it shares a boundry with and vice versa.   

## Hidden Markov Models:
A hidden markov model assumes there is an unobserved variable that has an effect on the distribution of the observed variable. For instance in finance we can conceptualize an up/down trend for the stock market based on whether there is a bull or a bear market. In this illustration we will be utilizing up/down stages of the disease to inform us of the distribution of biweekly mortality. The unobsorved trend will have a probability structure that ideally we should further model. 

## Detailed Introduction:
Assume there is a random variable Y, indexed by time t, $Y_{t}$. We would like to learn about the factors that effect $Y_{t}$. To do this we will assume a distribution over Y.  
To decide on the distribution we need to know about what it is Y actually represents. 
For the sake of this discussion Y will represent the biweekly mortality in California. 
Since we have information on each of California's counties, we will model the biweekly mortality $Y_{it}$ where $t=1 \cdots T=77$ and $i=1 \cdots N=58$.  
A discrete probability distribution which can be used to model mortality is Poisson. 

 $\Pr(Y{=}y)={\frac {\lambda ^{y}e^{-\lambda }}{y!}}$ where $\lambda$ is the mean and variance of the distribution. 
Once we have decided on the distribution of mortality, what separates out one distribution from the other is (are) the parameter(s) of the distribution. There are countless variables we can use instead of mortality but what changes uncertainty evaluation (probability) for Y is not what Y stands for but the parameter(s) used in the equation. 
For instance, in the equation for poisson distribution, the equation itself is fixed and you are trying to evaluate the probability that Y (what the random variable stands for) is y (the number itself). What changes the result (P(Y=y)) is really $\lambda$.   


In this presentation we assume that the biweekly mortality in the counties of California is Poisson distributed. 
$Y_{it}\sim Poisson(\lambda_{it})$
To eloborate on what I mentioned earlier, this by itself is not really interesting. We can assign probabilities to values other than the ones we have seen at a particular time $t$ in a particular county $i$. However this would rely very heavily on the assumption of biweekly mortality being Poisson distributed and not be of particular use by itself.  

The goal of the project is to model the biweekly mortality in the counties of California in order to answer the following questions:
## 1) How did socio-demographic factors effect biweekly mortality
## 2) How did these change through the progression of the pandemic
## 3) Can we detect regions (counties or collection of counties) that performed better than others.
## 4) Did specific policies at county or state level effect the unobserved up/down trends.

In order to answer these questions we will need some background so that we can develop models to help answer those questions:

## A) Data/tools for spatial models in R.

## B) [Some neccessary Probability concepts:](./SomeProbabilityConcepts.md) 

## C) Spatial Models (We will use Besag York Mollie).

## D) Hidden Markov Models.

## E) Bayesian Modeling Tool (STAN).

## F) Combine Bayesian Spatial and Hidden Markov Models
   
