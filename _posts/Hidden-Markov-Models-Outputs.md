---
layout: post
title: "Hidden Markov Models Output"
author: "Rasim M Musal"
date: "2023-10-18"
output:
  html_document:
   theme: darkly
   highlight: espresso
   toc: true
   keep_md: yes
   toc_float: true
   toc_collapsed: false
   toc_depth: 1
   number_sections: true
   usemathjax: true
tags: [MCMC OUtput, Convergence, Hidden Markov Model]
always_allow_html: true
---
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  TeX: { 
      equationNumbers: {
 
            autoNumber: "all",
            formatNumber: function (n) {return +n}
      } 
  }
});
</script>



# Description of MCMC Output for Analysis
Reading in the three chains, each one containing 10,000 samples for each parameters. Each row represents a joint density of 18,294 parameters. Some of these parameters are calculated based on the others. [The two STAN objects Lambda and Logalpha](https://mmusal.github.io/blog/2023/Hidden-Markov-Models-with-STAN/#trpr) will have the majority of these parameters.  For instance the T (77) length array of K (2) rows and N (58) columns matrix of Lambda has a total of 77*2*58=8,932 parameters recorded. Logalpha has the same number of parameters. Since both Lambda and Logalpha are calculated based on other parameters in the model itself we do not necessarily need to read them into R but it will allow us to ascertain convergence of all the parameters. 

# Reading in the data and checking for convergence


```r
options(mc.cores=parallel::detectCores(),auto_write = TRUE)
library(rstan)
library(xtable)
library(ggspatial)
#Specify the locations for the 3 separate chains and record them in strings

chain1='C:/Users/rm84/Documents/w2to77mcs_1.csv'
chain2='C:/Users/rm84/Documents/w2to77mcs_2.csv'
chain3='C:/Users/rm84/Documents/w2to77mcs_3.csv'

# create a vector of strings
csvfiles=c(chain1,chain2,chain3)
#create a stan object to query 
fit=read_stan_csv(csvfiles)
#Summary statistics for all the parameters in the object fitsummary
fitsummary=summary(fit)$summary

#Find out whether any of the Rhat convergence measures are over 1.1 which will flag non-convergence 
Rhats=fitsummary[,"Rhat"]
#How many parameters did not converge?
sum(Rhats>1.01,na.rm=TRUE)
```

```
## [1] 0
```

```r
#Which Rhats did not converge
Rhats[Rhats>1.01]
```

```
## named numeric(0)
```

A couple of things to note here regarding the parameters convergence. First, just because there are no Rhats less than 1.01 does not necessarily mean the posterior distribution that MCMC elicited converged to the correct posterior distribution. Second is that even though it is tempting to ignore convergence results of parameters that you are not interested, this is not suggested.
Also note that the new [Rhat convergence limit is now tightened to 1.01](https://projecteuclid.org/journals/bayesian-analysis/volume-16/issue-2/Rank-Normalization-Folding-and-Localization--An-Improved-R%cb%86-for/10.1214/20-BA1221.full) unlike previous suggestions of 1.1.



```r
#Create a dataframe from the stan object by extracting the parameters.
listofdraws=as.data.frame(rstan::extract(fit))
#Remove the stan object
rm(fit)
#Garbage collection.
gc()
```

```
##             used   (Mb) gc trigger    (Mb)   max used    (Mb)
## Ncells   1359418   72.7    2558270   136.7    2558270   136.7
## Vcells 551428238 4207.1 1725881378 13167.5 1725878099 13167.5
```

```r
#attach the dataframe to use the parameter names
attach(listofdraws)

#Number of counties in CA
N=58
#Number of time periods
T=77
t=1
#This is the R object lambda not to be confused with the STAN simulation values
lambda=array(dim=c(T,N,nrow(listofdraws)))
```

In order to assess whether the model can make any claim of model validity we will compare the actual values to the predicted averages.    
For the first time period we only need to use marginal probabilities of regime 1 or 2. The marginal and conditional probabilities are time invariant so that $P(k_{t})$ = $P(k=1)$ and $P(k_{t}|k_{t-1})$ is the same for all t. This property is also referred to as time invariance. 


\[E(Y_{1,n})=\lambda_{1,n}=(P(k=1)*\lambda_{1,1,n}+P(k=2)*\lambda_{1,2,n})\]

The eval, parse and paste "trick" is [not necessarily a good option](https://yihui.org/en/2023/02/eval-parse/). However in this instance it serves our needs, we are able to loop over the parameter names explicitly. 


```r
if(t==1){
  for(n in 1:N){
lambda[1,n,]= eval(parse(text=(paste0("lambda.",t,".",1,".",n))))*eval(parse(text=(paste0("pi1.",1,".1"))))+
eval(parse(text=(paste0("lambda.",t,".",2,".",n))))*eval(parse(text=(paste0("pi1.",1,".2"))))  
}}
```
Note that in this particular model we have added a county index to the conditional probabilities to the likelihood we [introduced in the earlier post](https://mmusal.github.io/blog/2023/Hidden-Markov-Models/#mjx-eqn-eqlikespec77). Therefore the new likelihood will be
\begin{equation}
\Pi_{i=1}^{i=58} P(k_{0}) \times P(k_{1,i}|k_{0,i}) \ldots \times P(k_{76,i}|k_{75,i}) \times P(Y_{0,i}|k_{0,i}) \times P(Y_{1,i}|k_{1,i}) \ldots \times P(Y_{76,i}|k_{76,i}) \label{eq:likespec77wa} 
\end{equation}

Therefore the expected value can be calculated via;

\[
\begin{aligned}
E(Y_{t,n})= & (P(k_{t,n}=1|k_{t-1,n}=1)*P(k_{t-1}=1))*\lambda_{t,1,n}+\\
& (P(k_{t,n}=2|k_{t-1,n}=1)*P(k_{t-1}=1))*\lambda_{t,2,n}+\\
& (P(k_{t,n}=1|k_{t-1,n}=2)*P(k_{t-1}=2))*\lambda_{t,1,n}+\\
&(P(k_{t,n}=2|k_{t-1,n}=2)*P(k_{t-1}=2))*\lambda_{t,2,n}
\end{aligned}
\]


```r
for(t in 1:T){
  for(n in 1:N){
lambda[t,n,]=(
      
      eval(parse(text=(paste0("lambda.",t,".",1,".",n))))*eval(parse(text=(paste0("A.",1,".",n,".",1))))*eval(parse(text=(paste0("pi1.",1,".1"))))+
        eval(parse(text=(paste0("lambda.",t,".",2,".",n))))*eval(parse(text=(paste0("A.",1,".",n,".",2))))*eval(parse(text=(paste0("pi1.",1,".1"))))+ 
        eval(parse(text=(paste0("lambda.",t,".",1,".",n))))*eval(parse(text=(paste0("A.",2,".",n,".",1))))*eval(parse(text=(paste0("pi1.",1,".2"))))+
        eval(parse(text=(paste0("lambda.",t,".",2,".",n))))*eval(parse(text=(paste0("A.",2,".",n,".",2))))*eval(parse(text=(paste0("pi1.",1,".2"))))
    )
  }}
```
where $\lambda_{t,k,n}$ is equal to 
\begin{equation}
E_{t,i}\times exp(\mu_{i,k}+\epsilon_t)
\label{eq:lambda1} 
\end{equation}

Now that we calculated the posterior mean of Covid-19 mortality in each time period and county for each simulated posterior distribution value, we will average them out in order to use in the root mean square error calculation. 


```r
meansoflambda=array(dim=c(T,N))

for(t in 1:T){
  for(n in 1:N){
meansoflambda[t,n]=mean(lambda[t,n,])
}
  }

setwd('C:/Users/rm84/Desktop/research/HMM/data/')

o=source(file='R2023noICUcoviddata1.r')
cadata=as.list(o$value)
names(cadata)=c('K','R','T','N','N_edges','node1',
                'node2','log_E','y','x',
                'VacPop','medage','medsex','race',
                'scaling_factor','START','END')

MSEarray=array(dim=c(T))
for(t in 1:T){
  MSEarray[t]=mean((cadata$y[t,]-meansoflambda[t,])^2)
}
MSE=mean(MSEarray)
rMSE=sqrt(MSE)

rMSE
```

```
## [1] 38.9525
```
As we can see from the Root Mean Square value the fit is not good however this is due to the \eqref{eq:lambda1} having only $\mu_{k}$ HMM component and $\epsilon_{t}$ which can be conceptualized as a time variant intercept. 
We can leave a serious discussion of output interpretation to a later post when we have more variables that explain the variation in Covid-19 mortality. We will still present some results. 


```r
mean(pi1.1.1)
```

```
## [1] 0.8963761
```

```r
mean(pi1.1.2)
```

```
## [1] 0.1036239
```
The pi1.1.1 and pi1.1.2 shows the posterior mean probability of a county being in regime 1 or 2 respectively. 


```r
mean(A.1.1.1)
```

```
## [1] 0.8352265
```

```r
mean(A.1.1.2)
```

```
## [1] 0.1647735
```

```r
mean(A.2.1.1)
```

```
## [1] 0.1972098
```

```r
mean(A.2.1.2)
```

```
## [1] 0.8027902
```

Shows the posterior mean transition probabilities of Alameda (county 1). It would be interesting to model these transition probabilities. Furthermore it will be important to relax the time invariance assumption.  
