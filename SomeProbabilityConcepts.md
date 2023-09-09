# Independence: 
If there are two random variables Y and $\theta$ (in Bayesian statistics parameters are random variables) 
and if we declare Y and $\theta$ to be probabilistically independent (there can be different interpretations of independence, in this illustration we specifically mean probabilistic independence) what we mean by that notationally is:

$P(Y|\theta)=P(Y)$ and/or $P(\theta|Y)$

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
The formulation is straightforward

$P(\theta|Y_{1},\cdots,Y_{N})=
\frac{P(Y_{1},\cdots,Y_{N}|\theta)*P(\theta)}{\P(Y_{1},\cdots,Y_{N}|\theta)*P(\theta)}{\int_{\theta} P(Y_{1},\cdots,Y_{N}|\theta)*P(\theta)}{\P(Y_{1},\cdots,Y_{N}|\theta)*P(\theta)}$
