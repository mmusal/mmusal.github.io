---
layout: post
title: "Spatial Relations and Besag-York-Mollie Framework"
author: "Rasim M Musal"
date: "2023-11-14"
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
bibliography: bibliography.bib 
tags: [Besag-York-Mollie, Areal Data]
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




# Models Using Areal Data

When determining the presence of spatial relationships, a first step will be to determine whether there is reason to suspect that data that is in close proximity shows a recognizable pattern such that data which are closer to each other are more similar than data that are further apart. There are different frameworks in constructing spatial models. In this project we are interested in using areal data to create this framework.

Areal data, can be defined as data which is observed within well defined boundaries. For instance, square grids, zip code tabulation areas, or as in our project, counties in a particular state are all examples of well defined boundaries. Once we decide to use Areal data we need to assume a neighborhood structure. Figure 1 below will help illustrate this process. 




<figure><img src="assets/img/countiesofcalifornia.png"><figcaption></figcaption></figure>
<a name="cal_map"></a><span style="color:white; text-align:center; display:inline-block; width:100%;">Figure 1: FIPS numbers are recorded on each county</span>

Turn your attention to San Bernardino. It is surrounded by counties with Federal Information Processing Standard (FIPS) codes 27, 29, 37, 59, 65. These in turn are effected by other counties. Therefore there is a need to define what constitute a neighbor for a county and the weight of that neighborhood. These are all questions that can be quantified via a neighborhood matrix. There are different schemes in literature. However the Bayesian model that will follow is rather rich in terms of parameters we will use a binary neighborhood structure where if a county shares a border with another county its' effect is to be estimated. We will explain how to create and represent this matrix in STAN. Based on the neighborhood structure one of the first things to do would be to determine whether there is any spatial effects to model.

# Moran's I 
Moran's I is an exploratory index for spatial correlation and has a range between -1 and 1.   
Per @banerjee2003hierarchical Moran's I and its variance is;
\[
I=\frac{N}{W}\frac{\sum_{i}\sum_{j}(y_{i}-\bar y)(y_{j}-\bar y)}{(y_{i}-\bar y)^{2}}
\label{eq:I}
\]
\[
Var(I)=\frac{n^{2}(n-1)S_{1}-n(n-1)S_{2}-2S_{0}^{2}}{(n+1)(n-1)^{2}S_{0}^{2}}
\label{eq:VarI}
\]
where
\[
S_{0}=\sum_{i \ne j} w_{ij}, S_{1}=\frac{1}{2}\sum_{i \ne j}(w_{ij}+w_{ji})^{2},S_{2}=\sum_{k} \left( \sum_{j} w_{kj} + \sum_{i} w_{ik} \right)^{2}
\]
where;

N is the number of areal components. (58 in California)

W is the sum of weights in the neighborhood weight matrix (58).

$y_{i}$ is the quantity of interest in areal unit i.

$\bar y$ is the mean quantity.
Under the null hypothesis that there is no spatial correlation the expected Moran's I value would be equal to
\[E(I)=\frac{-1}{n-1}
\label{eq:EI}
\]
Once we have \eqref{eq:I}, \eqref{eq:VarI}, and \eqref{eq:EI} we can calculate a z score via
\[z=\frac{(I-E(I))}{\sqrt{(Var(I))}}\]

and assuming Central Limit Theorem applies a hypothesis test can be performed. On the other hand a set of permutations involving areal units can have more robust results. 

In understanding how this all works we need to look at some r code.


```r
library(spdep)
library(spatialreg)
library(spatial)
setwd("C:/Users/rm84/Desktop/research/HMM/data")
load("workspacewithbasedata.RData")
```


```r
#shapeanddata is a simple features class object that contains 
#spatial information and county level feature information. 
#The poly2nb function creates the neighbor indexes 
nb_CA = poly2nb(shapeanddata)
#The neighbour indexes of county 1 in California (Alameda)
nb_CA[1]
```

```
## [[1]]
## [1]  7 38 39 41 43 50
```

```r
#The nb2listw function will create a the weight list of neighbors for each county 
#Each neighboring county's weight is equally weighted. 
#Therefore the sum of all weights in each row is 1 and the sum of all weights 
#is the number of counties in California (58).     
colw <- nb2listw(nb_CA, style="W")
#Gives the list of weights for Alameda 1/N_{1} where N_{1} is the number of 
#Alameda neighbors. 1/length(nb_CA[1]) == 1/6
colw$weights[1]
```

```
## [[1]]
## [1] 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667 0.1666667
```


# The Spatial Model

Once we determine what counts as a neighbor and that there is a spatial effect to account for, we will need to decide how to account for it. In doing so we will need to define 

\begin{equation}
P(Y_{1},\ldots,Y_{N}) 
\label{eq:joint} 
\end{equation}

which is not trivial. If we want to take into account spatial relations, it is clear that this will not happen via $P(Y_{1}) \times \ldots \times P(Y_{N})$. It is standard in Bayesian applications to obtain this joint form from full conditional probabilities. On the other hand as @banerjee2003hierarchical discusses that, these conditional probabilities can not be arbitrary and demonstrates the form they need to take to create a valid joint probability distribution \eqref{eq:joint}.    

