Independence: 
If there are two random variables Y and $\theta$ (in Bayesian statistics parameters are random variables) 
and if we declare Y and $\theta$ to be probabilistically independent (there can be different interpretations of independence, in this illustration we specifically mean probabilistic independence) what we mean by that notationally is:

$P(Y|\theta)=P(Y)$ and/or $P(\theta|Y)$

and in plain terms we are saying that the knowledge of $theta$ does not inform us as to the uncertainty evaluation of Y.
Of course this would lead to 

$P(Y)*P(\theta)=P(Y,\theta)$ if Y and $\theta$ are independent.

Conditional Independence: 
Assume there are two random variables $Y_{1}$ and $Y_{2}$. If $Y_{1}$ and $Y_{2}$ are conditionally independent given the random variable $\theta$:

$P(Y_{1}|\theta)*P(Y_{2}|\theta)=P(Y_{1},Y_{2}|\theta)$

Likelihood: 
The likelihood is the joint probability (function) of the data given the parameters. In many straightforward models with N data points, the likelihood can be calculated as

$P(Y_{1}|\theta) * \cdots* P(Y_{N}|\theta)=$

$P(Y_{1},\cdots,Y_{N}|\theta)$

