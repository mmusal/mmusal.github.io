# Independence: 
If there are two random variables Y and $\theta$ (in Bayesian statistics parameters are random variables) 
and if we declare Y and $\theta$ to be probabilistically independent (there can be different interpretations of independence, in this illustration we specifically mean probabilistic independence) what we mean by that notationally is:

$P(Y|\theta)=P(Y)$ and $P(\theta|Y)=P(\theta)$

and in plain terms we are saying that the knowledge of $theta$ does not inform us as to the uncertainty evaluation of Y.
Of course this would lead to 

$P(Y)*P(\theta)=P(Y,\theta)$ if Y and $\theta$ are independent.

# Conditional Independence: 
Assume there are two random variables $Y_{1}$ and $Y_{2}$. If $Y_{1}$ and $Y_{2}$ are conditionally independent given the random variable $\theta$:

$P(Y_{1}|\theta)*P(Y_{2}|\theta)=P(Y_{1},Y_{2}|\theta)$

# Likelihood: 
The likelihood is the joint probability (function) of the data given the parameters. In many straightforward models with N data points, the likelihood can be calculated as

$P(Y_{1}|\theta) * \cdots* P(Y_{N}|\theta)=$

$P(Y_{1},\cdots,Y_{N}|\theta)$

This set up of multiplying conditonal probabilities requires $P(Y_{i}|\theta), i = 1\cdots N$ to be independent (conditonally on $\theta$) and identically distributed.

# Inverting Probability
In statistical modeling we use probability distributions which have parameters that have explicit meanings. For instance if Y is assumed to follow a Poisson distribution, the distribution has a single parameter $\lambda$. This is the mean and variance of the distribution. If Y is assumed to follow a Normal distribution it has two parameters $\mu$ and $\sigma^{2}$ that stands for the mean and variance of Y. In modeling we would like to learn about the parameters of the random variable after observing data. In other words we would like to learn

$P(\theta|Y_{1}\cdots Y_{N})$
The formulation to calculate what we refer to as posterior distribution of the parameter(S) is (are) straightforward.

$P(\theta|Y_{1},\cdots,Y_{N})=
\frac{P(Y_{1},\cdots,Y_{N}|\theta)*P(\theta)}{\int_{\theta} P(Y_{1},\cdots,Y_{N}|\theta)*P(\theta)}$

Note that the denominator in the ratio above is simply the joint probability of the data.
$P(Y_{1},\cdots,Y_{N})$ 

when we integrate out the parameter. It should be remembered that choosing the right distribution to represent the uncertainty about Y have important implications since it determines the likelihood and all the calculations that follow. 
As simple as it is to write the formulation it is not usually straightforward to actually carry out the integration in the denominator. This is the reason why we use Monte Carlo Markov Chains.  

# Probability distributions

To demonstrate the distribution of Y with different $\lambda$ values, we will use R:
#Libraries to be used:

```
library(ggplot2)
library(dplyr)
```

#values that lambda is going to take

```
lambda=c(1,3,5,10,20)
```

#how many variates we will generate

``` 
n=10000
```

#the two dimensional array that is going to be populated via the rpois 

#(randomly generate poisson distributions) function.

```
poissonvariates=as.data.frame(matrix(nrow=n,ncol=length(lambda)))
```

#Below is the for loop that is going to generate the variates. The number of times 
#that the statements in the curly bracket is going to run is the length of the lambda vector (5).
#the rpois has 2 required parameters that we need to provide values for; n and lambda. 
#n is the number of values we will generate, lambda is the mean/variance of the distribution. 

```
for(i in 1:length(lambda)){
poissonvariates[,i]=rpois(n=n,lambda[i])}
```

####give names to the dataframe columns
```
names(poissonvariates)=c("lambda=1","lambda=3","lambda=5","lambda=10","lambda=20")
```

#we need to stack the columns for ggplot function

```
stackedpoisson=stack(poissonvariates)
```

#once the columns are stacked they have the names: values ind

#We need to summarize the data to create the barchart.

```
stackedpoisson2 <- stackedpoisson %>% 
  group_by(values, ind) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/n) 
```

#and here is the plot that represents poisson variates simulated with 

#different lambda values. 

```
ggplot(data=stackedpoisson2, aes(x=values,y=perc, fill=ind))+
  geom_bar(stat="identity")+
  facet_grid(. ~ ind)+
  labs(x = "Sim. Values", y = "Percent", fill = "Lambda")
```

  
![Rplot](https://github.com/mmusal/mmusal.github.io/assets/11746560/32b2a741-5c22-4a1c-8898-d7d3b8b6e442)



The point of the plot is to show you how changing the value lambda, changes the uncertainty evaluations regarding the values of Y.
