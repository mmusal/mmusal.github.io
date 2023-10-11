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


A good introduction to Hidden Markov Models is [illustrated](https://nipunbatra.github.io/hmm/). In this post we will explain the application of the model we are developing using notation that will be consistent across the project.
We would like to make inferences about effects on $Y_{ti}$, the number of Covid-19 mortality recorded at time $t$ in county $i$. The time index $t$ is between $(t=0 \ldots (T-1)=76)$. There are 58 counties in California $(i=1, \ldots ,58)$ For the moment we will ignore covariates and focus on "regime" $k$ that define the characteristics of the uncertainty that defines the uncertainty of $Y_{ti}$.

Assume $i=1$
In order to make inferences via our HMM model we would need to obtain the joint probability density of mortality and regimes at every $t$.   
\[P(Y_{1:T},k_{1:T})\]
To obtain this joint probability we will need to make assumptions.
1) A regime determines the parameters of the uncertainty distribution.
2) Given the parameters the regime does not provide additional information about $Y_{t}$.  
\[P(Y_{t}|\lambda_{tk})=P(Y_{t}|\lambda_{tk},k_{t})=P(Y_{t}|k_{t})\]
3) For the purposes of this post we will assume that at time $t$ the distribution of Covid-19 mortality follows a Poisson distribution.
\[P(Y_{t}|\lambda_{tk},k_{t})=Poisson(Y_{t}|\lambda_{tk})\]


