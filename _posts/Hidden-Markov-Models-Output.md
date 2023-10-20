---
layout: post
title: "Hidden Markov Models Output"
author: "Rasim M Musal"
date: "2023-10-14"
output:
  html_document:
   theme: darkly
   highlight: espresso
   toc: true
   keep_md: yes
   toc_float: true
   toc_collapsed: false
   toc_depth: 3
   number_sections: true
   usemathjax: true
tags: [rshiny, maps,SMR]
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



## Description of Output
Reading in the three chains, each one containing 10,000 samples for each parameters. Each row represents a joint density of 18,294 parameters. Some of these parameters are calculated based on the others. For instance the T length array of K rows and N columns matrix has a total of 77*2*58=8,932 parameters recorded. Logalpha has the same number of parameters. Since both Lambda and Logalpha are calculated based on other parameters in the model itself we do not necessarily need to read them into R but it will allow us to ascertain convergence of all the parameters.        

## Reading in the data and checking for convergence

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
## Ncells   1359328   72.6    2558345   136.7    2558345   136.7
## Vcells 551427459 4207.1 1725874928 13167.4 1725862315 13167.3
```

```r
#attach the dataframe to use the parameter names
attach(listofdraws)

#Number of counties in CA
N=58
t=1
#This is the R object lambda not to be confused with the STAN simulation values
lambda=array(dim=c(T,N,nrow(listofdraws)))
```

For the first time period we only need to use marginal probabilities.

$E(Y_{1,n})=\lambda_{1,n}=P(k=1)*\lambda_{1,1,n}+P(k=2)*\lambda_{1,2,n}$

The eval, parse and paste "trick" is [not necessarily a good option](https://yihui.org/en/2023/02/eval-parse/). However in this instance it serves our needs, we are able to loop over the parameter names explicitly. 


```r
if(t==1){
  for(n in 1:N){
lambda[1,n,]= eval(parse(text=(paste0("lambda.",t,".",1,".",n))))*eval(parse(text=(paste0("pi1.",1,".1"))))+
eval(parse(text=(paste0("lambda.",t,".",2,".",n))))*eval(parse(text=(paste0("pi1.",1,".2"))))  
}}
```

\[
\begin{aligned}
E(Y_{t,n})= & (P(k_{t}=1|k_{t-1}=1)*P(k_{t-1}=1))*\lambda_{t,1,n}+\\
& (P(k_{t}=2|k_{t-1}=1)*P(k_{t-1}=1))*\lambda_{t,2,n}+\\
& (P(k_{t}=1|k_{t-1}=2)*P(k_{t-1}=2))*\lambda_{t,1,n}+\\
&(P(k_{t}=2|k_{t-1}=2)*P(k_{t-1}=2))*\lambda_{t,2,n}
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

print(meansoflambda[1,])
```

```
##  [1] 3.067203e-02 1.475992e-05 1.150977e-03 6.519284e-03 1.411355e-03
##  [6] 4.372385e-04 2.113913e-02 8.363753e-04 1.200759e-03 4.038168e-02
## [11] 5.987910e-04 2.004447e-03 1.433597e-02 7.876607e-04 3.329927e-02
## [16] 6.536105e-03 2.079092e-03 9.623964e-04 4.015781e-01 4.726630e-03
## [21] 4.536395e-03 3.507588e-04 1.650067e-03 1.129069e-02 1.306685e-04
## [26] 1.426345e-04 8.107381e-03 1.909962e-03 1.269306e-03 8.585135e-02
## [31] 7.157347e-03 1.870238e-04 6.518520e-02 4.010567e-02 1.595545e-03
## [36] 8.939659e-02 6.955137e-02 1.799225e-02 3.017880e-02 4.142631e-03
## [41] 8.848652e-03 7.421953e-03 2.737208e-02 2.507848e-03 8.372591e-03
## [46] 9.242664e-05 9.298740e-04 5.933132e-03 7.579390e-03 2.306701e-02
## [51] 2.686078e-03 3.087017e-03 2.833474e-04 2.194817e-02 2.589200e-03
## [56] 1.673767e-02 6.047600e-03 1.624215e-03
```


