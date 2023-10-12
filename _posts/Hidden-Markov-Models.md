---
layout: post
title: "Hidden Markov Models"
author: "Rasim M Musal"
date: "2023-10-10"
output:
  html_document:
   theme: darkly
   highlight: espresso
   toc: true
   keep_md: yes
   toc_float: true
   toc_collapsed: true
   toc_depth: 3
   number_sections: true
tags: [rshiny, maps,SMR]
always_allow_html: true

---


A good introduction to Hidden Markov Models is [illustrated here](https://nipunbatra.github.io/hmm/). In this post we will explain the application of the model we are developing using notation that will be consistent across the project.
We would like to make inferences about effects on $Y_{ti}$, the number of Covid-19 mortality recorded at time $t$ in county $i$. The time index $t$ is between $(t=0 \ldots (T-1)=76)$. There are 58 counties in California $(i=1, \ldots ,58)$ For the moment we will ignore covariates and focus on "regime" $k$ that define the characteristics of the uncertainty that defines the uncertainty of $Y_{ti}$. In the discussion below we are going to ignore hyper-priors for clarity. These will be made explicit in the post explaining the STAN code. 

Assume $i=1$
In order to make inferences via our HMM model we would need to obtain the joint probability density of mortality and regimes at every $t$.   
\[P(Y_{1:T},k_{1:T})\]
To obtain this joint probability we will need to make assumptions.
1) A regime determines the parameters of the uncertainty distribution.
2) Given the parameters the regime does not provide additional information about $Y_{t}$.  
\[P(Y_{t}|\lambda_{tk})=P(Y_{t}|\lambda_{tk},k_{t})=P(Y_{t}|k_{t})\]
3) For the purposes of this post we will assume that at time $t$ the distribution of Covid-19 mortality follows a Poisson distribution.
\[P(Y_{t}|\lambda_{tk},k_{t})=Poisson(Y_{t}|\lambda_{tk})\]

If we did not have different regimes across a time component we could have used mixture models due to their [close relationship](https://www.utstat.toronto.edu/~rsalakhu/sta4273/notes/Lecture11.pdf) with hidden markov models. Assuming we had a mixture of 2 Poisson distributions (K=2) we could write the mixture distribution of Covid-19 mortality at time period $t$ county $i$ as 

\[P(Y_{ti}) \sim \omega_{k=1} \times Pois.(\lambda_{t1i})+\omega_{k=2} \times Pois.(\lambda_{t2i}) \]

This setup assumes that $Y_{ti}$ is independent and identically distributed both across $t$ and $i$. If we wrote the likelihood for this model it would look like
\[P(Y_{0,1},Y_{0,2},\ldots Y_{0,N},Y_{1,1} \ldots Y_{T-1,N}|\lambda_{0,1}\ldots\lambda_{T-1,N},\omega_{1},\omega_{2})=\]
\[\Pi_{i=1}^{i=N}\Pi_{t=0}^{t=T-1}  (\omega_{k=1} \times Pois.(\lambda_{t1i})+\omega_{k=2} \times Pois.(\lambda_{t2i})),\]

where $\lambda_{tki}$ represents the mean/variance of Covid-19 mortality at time $t$, regime k, county i.

We would like to relax this assumption about $t$ and instead make the claim that the previous observations inform us as to $Y_{t}$ and build Hidden Markov Models (HMM).

In formulating the HMM we again assume K=2, however we will assume first order Markov model (markov property). Below, k_{t} is 1 or 2.
\[P(k_{t}|k_{t-1})=P(k_{t}|k_{t-1},\ldots,k_{0}).\]
In an order one Markov model the only thing that matters to inform the uncertainty of $k_{t}$ is at  $k_{t-1}$.
In addition 

\[P(Y_{tki}|k_{t})=P(Y_{tki}|k_{t},Y_{(t-1),k,i},\ldots Y_{0,k,i} )\]

The matrix created from the probabilities between the states is called the transition matrix. The format of the matrix ensures that it is a stochastic matrix such that the rows sum to 1. 
                                                  

| From/To | $k_{t}=1$ | $k_{t}=2$ | 
|:-------|:------|:-----------|
| $k_{t-1}=1$ | $\pi_{11}$ | 1-$\pi_{11}$ |
| $k_{t-1}=2$ | $\pi_{21}$ | 1-$\pi_{21}$ |

We assume that the transition matrix is time invariant, meaning that regardless of the value $t$, the value of $P(k_{t}|k_{t-1})$ does not change. There are models which relax this assumption.  

In order to make inferences we will need to obtain the likelihood of the model

\[P(Y_{01},\ldots,Y_{(T-1)N})\]

Assume T=2 and ignore for the moment county index i.
\[P(Y_{0}|k_{0}) \times P(k_{0}) \times P({k_{1}|k_{0}}) \times P(Y_{1}|k_{1})\]
This is equivalent to 
\[P(Y_{0}|k_{0}) \times P(Y_{1}|k_{1}) \times P(k_{0},k_{1}) \]

\[P(Y_{0},Y_{1}|k_{0},k_{1}) \times P(k_{0},k_{1}) \]
The first component $P(Y_{0},Y_{1}|k_{0},k_{1})$, is calculated by the first multiplication of conditionals due to the conditional independence of the observations given the states. 


Putting it all together for the project and including the county index, for the 77 observed time periods in the 58 counties of CA we get; 
\[\Pi_{i=1}^{i=58} P(k_{0}) \times P(k_{1}|k_{0}) \ldots \times P(k_{76}|k_{75}) \times P(Y_{0 i}|k_{0}) \times P(Y_{1 i}|k_{1}) \ldots P(Y_{76 i}|k_{76}) \]




                                             







