#On this page we will discuss various steps in developing a Hidden Markov Model Spatial Model. 
First we need to discuss each model separately. 
Hidden Markov Models:
Assume there is a random variable Y, indexed by time t, $Y_{t}$. We would like to learn about the factors that effect $Y_{t}$. To do this we will assume a distribution over Y.  
To decide on the distribution we need to know about what it is Y actually represents. 
For the sake of this discussion Y will represent the biweekly mortality in California. 
Since we have information on each of California's counties (N=58) we will model the biweekly mortality $Y_{it}$ where $t=1 \cdots T=77$ and $i=1 \cdots N=58$.  

