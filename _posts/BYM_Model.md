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
   toc_depth: 3
   number_sections: true
   usemathjax: true
bibliography: bibliography.bib 
tags: [Besag-York-Mollie, Areal Data, Moran's I, Spatial Effect]
always_allow_html: true
editor_options: 
  markdown: 
    wrap: 72
---

```{=html}
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
```




# Models Using Areal Data

In this section we will give a brief description of spatial models, and
describe processes which will eventually allow us to describe the
Besag-York-Mollie model and its extensions. An in depth introduction to
many of the R related concepts we will discuss here is in [Crime Mapping
Textbook](https://maczokni.github.io/crime_mapping_textbook/)

We will first discuss Moran's I which provides an exploratory test for
spatial relations. This is to be followed by a discussion involving the
modeling of joint probability of spatial effects $\phi$. In addition we
will illustrate how to prepare the data for STAN.

# Moran's I

When determining the presence of spatial relationships, a first step
will be to determine whether there is reason to suspect that data that
is in close proximity shows a recognizable pattern such that data which
are closer to each other are more similar than data that are further
apart. There are different frameworks in constructing spatial models. In
this project we are interested in using areal data to create this
framework.

Areal data, can be defined as data which is observed within well defined
boundaries. For instance, square grids, zip code tabulation areas, or as
in our project, counties in a particular state are all examples of well
defined boundaries. Once we decide to use Areal data we need to assume a
neighborhood structure. Figure 1 below will help
illustrate this process.



<figure><img src="assets/img/countiesofcalifornia.png"><figcaption></figcaption></figure>

<a name="cal_map"></a><span style="color:white; text-align:center; display:inline-block; width:100%;">Figure 1: FIPS numbers are recorded on each county</span>

Turn your attention to San Bernardino. It is surrounded by counties with
Federal Information Processing Standard (FIPS) codes 27, 29, 37, 59, 65.
These in turn are effected by other counties. Therefore there is a need
to define what constitute a neighbor for a county and the weight of that
neighborhood. These are all questions that can be quantified via a
neighborhood matrix. There are different schemes in literature. However
the Bayesian model that will follow is rather rich in terms of
parameters we will use a binary neighborhood structure where if a county
shares a border with another county its' effect is to be estimated. We
will explain how to create and represent this matrix in STAN. Based on
the neighborhood structure one of the first things to do would be to
determine whether there is any spatial effects to model.

Moran's I is an exploratory index for spatial correlation and has a
range between -1 and 1.\
Per @banerjee2003hierarchical Moran's I and its variance is; $$
I=\frac{n}{W}\frac{\sum_{i}\sum_{j}w_{ij}(y_{i}-\bar y)(y_{j}-\bar y)}{(y_{i}-\bar y)^{2}}
\label{eq:I}
$$ $$
Var(I)=\frac{n^{2}(n-1)S_{1}-n(n-1)S_{2}-2S_{0}^{2}}{(n+1)(n-1)^{2}S_{0}^{2}}
\label{eq:VarI}
$$ where $$
S_{0}=\sum_{i \ne j} w_{ij}, S_{1}=\frac{1}{2}\sum_{i \ne j}(w_{ij}+w_{ji})^{2},S_{2}=\sum_{k} \left( \sum_{j} w_{kj} + \sum_{i} w_{ik} \right)^{2}
$$ where;

-   n is the number of areal components. (58 in California).

-   $w_{ij}$ is the weight between areal unit i and j which for our
    purposes is either 0 or 1.

-   W is the sum of weights in the neighborhood weight matrix (58).

$y_{i}$ is the quantity of interest in areal unit i.

$\bar y$ is the mean quantity. Under the null hypothesis that there is
no spatial correlation the expected Moran's I value would be equal to
$$E(I)=\frac{-1}{n-1}
\label{eq:EI}
$$ Once we have \eqref{eq:I}, \eqref{eq:VarI}, and \eqref{eq:EI} we can
calculate a z score via $$z=\frac{(I-E(I))}{\sqrt{(Var(I))}}
\label{eq:z}$$

and assuming Central Limit Theorem applies a hypothesis test can be
performed. Per @banerjee2003hierarchical as convergence to normal
distribution is harder for ratio of quadratic values, a large set of
permutations involving areal units that lead to simulated Moran's I
values can be obtained for more robust results. A new I value would be
calculated for each areal permutation and the original I value can be
compared against the randomly permuted areal units' I values.

In understanding how this all works we need to look at some r code.


```r
library(spdep)
library(spatialreg)
library(spatial)
library(tidyr)
library(ggplot2)
setwd("C:/Users/rm84/Desktop/research/HMM/data")
load("workspacewithbasedata.RData")
cnames=read.table("countynames.txt",sep="\t",header=TRUE)
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

We are going to set up some objects to calculate not just the I value
but p value via

1.  Hypothesis testing based on Central Limit Theorem.
2.  Areal Unit Permutation.

Note that the quantity of interest in our application is SMR of Covid-19
mortality rather than Covid-19 mortality itself. We find that this is
more realistic to determine whether spatial structure exists as there
are large fluctuations in Covid-19 mortality due to differences in
population size.


```r
#Number of permutations for 2. Areal Unit Permutation.
nsim <- 5000
#Random number seed
set.seed(1234)
# Number of Biweeks below N is the total number of counties in CA 58.
NumBiWeeks=77
#SMR is standardized mortality ratio observed/expected
SMR=array(dim=c(NumBiWeeks,N))
#Since E is calculated based on aggregate mortality in California and population size of each county SMR is better for investigating whether there is a spatial structure.
SMR=y/exp(log_E)
```

We already discussed how to calculate I in \eqref{eq:I}. Z score can be
calculated via \eqref{eq:z}. We have 77 biweeks and for each biweek we
are going to calculate 2 sets of p values associated with the Moran's I
scores and visualize them across this time series. The objects,
sim_pvalue and asympt_pvalue are going to contain the p values
associated with permuted spatial neighborhoods and central limit theorem
is assumed to hold, respectively. Since we hypothesize that like values
are clustered together the p values are going to be calculated via
$P(Z>z)$.


```r
#Object where p values are calculated based on spatial neighborhood permutations
sim_pvalue=array(dim=c(NumBiWeeks))
#Object where p values are calculated based on assuming Central Limit Theorem applies
asympt_pvalue=array(dim=c(NumBiWeeks))
#Object that holds Moran's I values.
Ivalue=array(dim=c(NumBiWeeks))
#Hypothesis tests are directional as discussed above. The hardcoded values 58 can be parameterized to generalize. n stands for number of units S0 stands for total weight in the neighbourhood matrix. W and S_0 will have the same values.
for(i in 1:NumBiWeeks){
sim_pvalue[i]<-moran.mc(SMR[i,], listw=colw, nsim=nsim, alternative="greater")$p.value
asympt_pvalue[i]<-moran.test(SMR[i,],colw, alternative="greater", zero.policy=TRUE)$p.value
Ivalue[i]<-moran(SMR[i,],colw,n=58,S0=58)$I
}
#Number of hypothesis tests that have p values less than 0.05 (our pick of alpha)
sum(sim_pvalue<0.05)
```

```
## [1] 28
```

```r
sum(asympt_pvalue<0.05)
```

```
## [1] 31
```

```r
#Create a data frame to contain the created values 
moransdf=as.data.frame(matrix(nrow=NumBiWeeks,ncol=4))
moransdf[,1]=sim_pvalue
moransdf[,2]=asympt_pvalue
moransdf[,3]=Ivalue
moransdf[,4]=c(1:NumBiWeeks)  
#Name each column
names(moransdf)=c("p_value_mc","p_value_asymp","I_Value","Biweek")

#Create a data frame which will have 3 columns
#Biweek, the names which has the labels in quotes and the numbers recorded under #value 
moransdf_long <-moransdf %>% pivot_longer(c("p_value_mc","p_value_asymp","I_Value"),names_to = "names")
#Visualizing the values
ggplot(moransdf_long, aes(x=Biweek, y=value))+
geom_point(stat="identity")+
facet_wrap(~names,  ncol=1)+
scale_x_discrete(limits=c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,77))+
xlab("Biweek t")
```

<figure><img src="BYM_Model_files/figure-html/unnamed-chunk-7-1.png"><figcaption></figcaption></figure>

As we can see from the visualization of p values and the number of them
below the arbitrary but often used alpha value of 0.05, that for
majority of the time periods we can not refuse the null hypothesis that
there is no spatial correlation in the areal units. On the other hand in
nearly 30 of the biweeks the p values are less than 0.05 and after all,
the I value is an exploratory measure. We will describe building the
spatial model in the next section.

# The Joint Spatial Effect $\phi$

Once we determine what counts as a neighbor and that there is a spatial
effect to account for, we will need to decide how to account for it. In
doing so we will need to define

\begin{equation}
P(\phi_{1},\ldots,\phi_{N}), 
\label{eq:joint} 
\end{equation} where N is the number of areal units (i.e.: 58 counties
in CA). Writing this form is not trivial. Furthermore, if we want to
take into account spatial relations, it is clear that this will not
happen via $P(\phi_{1}) \times \ldots \times P(\phi_{N})$. It is
standard in Bayesian applications to obtain this joint form from full
conditional probabilities. However as @banerjee2003hierarchical
discusses, these conditional probabilities can not be arbitrary and
demonstrates the form they need to take to create a valid joint
probability distribution \eqref{eq:joint}.

## Markov Random Fields (MRFs)

Markov Random Fields (MRFs), a special stochastic process, provides an
important method to create this joint distribution from conditional
distributions. A good description of these mechanisms are described in
@besag74 and @cressie2015statistics and an introductory set of
illustrations can be [found
here.](https://ermongroup.github.io/cs228-notes/representation/undirected/)
MRFs, a set of random numbers, can be represented as undirected graphs
with specific conditional independence properties. The main advantage of
using MRFs is being able to use spatial effect probabilities conditional
on neighborhood structure rather than using all areal units in order to
create the joint probability distribution.

We will give a very brief introduction to these mechanism by first
describing the term clique within the context of spatial analysis. We
can define a clique as a set of areal units where each unit is a
neighbor to all in the clique. For instance in figure
Figure 1 Kern, Los Angeles, Ventura (County FIPS == 111)
form a clique since they all share a boundary with each other. This is
also the maximal clique since it is not a subset of another clique.

Using MRFs'
[properties](https://www.cs.cmu.edu/~16831-f14/notes/F11/16831_lecture07_bneuman.pdf)
and their relationship to Gibbs distributions (Fields) we can use
quadratic potentials between areal units to come up with the joint
probability of all areal units.

```{=tex}
\begin{equation}
p(\mathbf{\phi})=p(\phi_{1},\ldots,\phi_{N})=
\frac{1}{Z}exp\bigg[-\sum_{c_{i}\epsilon C}f_{i}(c_{i})\bigg],
\label{eq:cliques}
\end{equation}
```
where $c_{i}$ is 2 member cliques, $i$ is going to go from 1 to the
total number of unique pairwise cliques and Z is a normalizing constant.
We will expand more on this as we detail the code. For instance there is
clique of county FIPS 29 (Kern) and 71 (San Bernardino). The function
$f_{i}$, potentials between spatial areas will be based on squared
differences on $\phi$ terms leading to the result; \begin{equation}
-\frac{1}{2} \sum (\phi_{i}-\phi_{j})^2
\label{eq:sqrddif}
\end{equation}

where the index $i$ and $j$ are neighbors that form the clique. The
detailed assumptions and restrictions given by @besag74 allow the normal
distribution to describe the uncertainty in the joint distribution of
$\phi$.

In the
[link](https://mc-stan.org/users/documentation/case-studies/icar_stan.html)
we can see the linear algebra that shows the derivation of the potential
function via the assumption of multivariate normal distribution with a
mean of 0 for the joint probability of $\phi$. Note that the mean of 0
for $\phi$ is to ensure an identifiable joint distribution. Further
discussion will await the code.

## Preparing the Areal Data for STAN

First we are going to create a neighborhood list using the poly2nb
function from the spdep package.


```r
nb_CA = poly2nb(shapeanddata)
print(nb_CA)
```

```
## Neighbour list object:
## Number of regions: 58 
## Number of nonzero links: 288 
## Percentage nonzero weights: 8.561237 
## Average number of links: 4.965517
```

We obtain some information regarding our spatial object (graph object),
namely that there are 58 areal units (counties) 144
shared borders. This means in a neighborhood 58 by 58
matrix, with (3364) cells, 144 of the cells will
have a non zero value. Furthermore the nb_CA object will contain
neighborhood structure between counties. For instance;


```r
nb_CA[[1]]
```

```
## [1]  7 38 39 41 43 50
```

Shows that the first county 1,(Alameda) shares borders with
counties 7 (Contra Costa), 38 (San Francisco), 39
(San Joaquin), (41 San Mateo), (43
Santa Clara) and (50 Stanislaus).

The pairwise cliques (neighbors) is obtained in r via the nb_data_funs
functions provided by Mitzi Morris via the [github
site](https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/nb_data_funs.R).
Among the many functions the author provides, an important one we will
use is nb2graph.


```r
nbs=nb2graph(nb_CA)

N = nbs$N
node1 = nbs$node1
node2 = nbs$node2
N_edges = nbs$N_edges
```

This function converts the spdep object into a form that will allow STAN
to model the joint distribution of $\phi$. We will pass N as the number
of areal units (counties) which is `r`nbs\$N\`. In addition node1 and
node2 will contain information which will be used in computing
$(\phi_{i} - \phi_{j})^2$.

There are going to be 144 unique cliques (neighbor pairs),
we will look at the first 10 such cliques below.


```r
node1[1:10]
```

```
##  [1] 1 1 1 1 1 1 2 2 2 2
```

```r
node2[1:10]
```

```
##  [1]  7 38 39 41 43 50  3  5  9 26
```

As described above county 1,(Alameda), has
6 listed neighbors in node 2 corresponding to
county 1 in node 1, which can be seen in the R output above.

Once the spatial and covariate data is prepared for STAN by combining
them into a list, it is exported via the stan_rdump function as a
text-file via rstan package.


```r
dump=c('K','R','T','N','N_edges','node1','node2','log_E','y','x',
       'VacPop','medage','medsex','race',
       'scaling_factor','START','END')
#################################################
stan_rdump(dump, file = "C:/Users/rm84/Desktop/research/2023coviddata.r", append = FALSE,
           envir = parent.frame(),
           width = options("width")$width,
           quiet = FALSE)
```

# Besag York Mollie (BYM) Model

BYM model was first implemented in @besag1991bayesian which has a log
linear Poisson likelihood. The original form of the distribution can be
stated as;

```{=tex}
\begin{equation}
Y_{i} \sim Poisson(\lambda_{i}) \\
log(\lambda_{i}) = log(E_{i})  \times (\beta*X_{i}+\theta_{i}+\phi_{i}),  
\end{equation}
```
where $\phi$, areal spatial effect, has the form as described via
\eqref{eq:sqrddif}. It will be a vector with N elements same as
$\theta$, the specific random effects of the areal units. You can conceptualize $\theta_{i}+\phi_{i}$ as the model errors.
In this post
we will ignore the original model and focus on the extended version which circumvents the convergence issues arising from decomposing the model errors into two separate components. Per the modified BYM (BYM2) of @morris2019bayesian, using \eqref{eq:lambda} ; 

```{=tex}
\begin{equation}
Y_{it} \sim Poisson(\lambda_{i}) \\
log(\lambda_{i}) = log(E_{i})  \times \left(\beta*X_{i}+\left(\sqrt{\frac{\rho}{s}}\times \phi_{i}+\sqrt{(1-\rho)}\times \theta_{i}\right)\times \sigma \right),
\label{eq:lambda}
\end{equation}
```
-   $\rho \in [0,1]$ controls the contribution weight of spatial vs random
    effects. If $\rho$ is 0 (1), the contribution of the errors are due solely to $\phi$ ($\theta$). 
-   s is the scaling factor which can be computed using the 
    function provided by once again by Mitzi Morris at the [github
    site](https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/nb_data_funs.R).
    The constant value is evaluated from the matrix $Q$ and for ICAR
    models it is equal to $$Q=D(I-A)$$ where
-   D is the N by N diagonal matrix with the number of neighbors for
    each areal unit.
-   I is the N by N identity matrix.
-   A is the N by N adjacency matrix.\
    Using these matrices as input, the function "scale_nb_components"
    calculates the geometric mean of the variances which are on the
    diagonal of $Q^{-1}$ matrix diagonal. This in turn ensures that the
    prior of $\phi$ will have a standard deviation of 1. We will
    set $\theta$'s prior standard deviation as 1. This is necessary to be able to 
    interpret $\sigma$ as the overall standard deviation of combined
    $\phi$ and $\theta$ terms.

We will include a time component t which will modify the equation. For
now we will focus on STAN program.

## STAN program and BYM Original and Extended Version

In this section we explain how to setup and execute the Besag York
Mollie model. Most of this section is paraphrased from
[here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6830524/), where the
analysis is based on a single slice of time. The authors explain the
original formulation for BYM which have convergence issues and propose a
new solution.

## STAN Function for $\phi$

An elaborate explanation of the function created in STAN for $\phi$
exists
[here](https://mc-stan.org/users/documentation/case-studies/icar_stan#adding-an-icar-component-to-a-stan-model).


```stan
functions {
  real icar_normal_lpdf(row_vector phi, int N, int[] node1, int[] node2) 
  {
    return -0.5 * dot_self(phi[node1] - phi[node2])
    + normal_lpdf(sum(phi) | 0, 0.0001 * N);
  }
          }
```

The dot_self function, multiplied with -0.5 is the application of
equation \eqref{eq:sqrddif}. It should be noted here that the variable
declared as sum of $\phi$ on the last line of the puts a soft constraint
on the mean of $\phi$ vector of 0 with a standard deviation of 0.0001.

## STAN Data/Transformed Section

The data section will include the time index extending the research
explained so far. The STAN program has comments next to these constants
that are self explanatory.

### Mortality Matrix: y[T,N];

The matrix y are the mortality values aggregated in biweeks[^1] from the
daily data obtained from California Health and Human Services is between
2020-02-01 and 2023-01-17. We will remove the last 16 days of
observation from the analysis. The first day of 2023 (Sunday) to is
contained in the last observed biweek of 2022.

[^1]: <https://data.chhs.ca.gov/dataset/covid-19-time-series-metrics-by-county-and-state>

### Gini Inequality: matrix [T,N] gini;

The Gini inequality is obtained from American Community Survey data[^2].
Gini inequality is a number that goes between 0 and 1, 0 indicating
every individual element (household) having the same income and 1
indicates one single household have all the income in the population of
interest. [This number is calculated by calculating the percentage of
area under the Lorenz curve vs the total area under the curve under the
assumption of complete income
equality](https://www.census.gov/topics/income-poverty/income-inequality/about/metrics/gini-index.html).
Gini inequality calculated by [Census](https://www.census.gov/) has 1
and 5 year estimates, we use the 5 year estimates in this research.

[^2]: <https://data.census.gov/table/ACSDT5Y2022.B19083?q=inequality%20in%20counties%20of%20california&g=040XX00US06$0500000>


```r
library(ggplot2)
library(viridis)
library(cowplot)
library(ggspatial)
library(ggpubr) 

wd="C:/Users/rm84/Desktop/research/HMM/data/"
setwd(wd)
GINI=read.table("gini20and21and22.csv",header=TRUE,sep=",")
attach(GINI)
#Once we attach GINI we can refer to the Gini_Index_XXXX variables. The values in GINI dataset are in the order of the county FIPS numbers which is how they are also listed in shapeanddata object. This is why aes parameter can use the variable directly rather than having to be attached to shapeanddata. 
map_GINI20=ggplot() +
  annotation_spatial(shapeanddata) +
  ggtitle("Gini 2020")+
  theme(plot.title = element_text(hjust = 0.5))+
  layer_spatial(shapeanddata, aes(fill = (Gini_Index_2020)))+
  theme(legend.title= element_blank())+
  labs(fill = "GINI. 20")+scale_fill_viridis(limits = c(0.39,0.6),direction=-1)+
  theme(legend.position = "none")

map_GINI21=ggplot() +
  annotation_spatial(shapeanddata) +
  ggtitle("Gini 2021")+
  theme(plot.title = element_text(hjust = 0.5))+
  layer_spatial(shapeanddata, aes(fill = (Gini_Index_2021)))+
  theme(legend.title= element_blank())+
  labs(fill = "GINI. 21")+scale_fill_viridis(limits = c(0.39,0.6),direction=-1)+
  theme(legend.position = "none")

map_GINI22=ggplot() +
  annotation_spatial(shapeanddata) +
  ggtitle("Gini 2022")+
  theme(plot.title = element_text(hjust = 0.5))+
  layer_spatial(shapeanddata, aes(fill = (Gini_Index_2022)))+
  theme(legend.title= element_blank())+
  labs(fill = "GINI. 22")+scale_fill_viridis(limits = c(0.39,0.6),direction=-1)+
  theme(legend.position = "none")
              
ggarrange(map_GINI20,map_GINI21,map_GINI22,nrow=1,ncol=3,common.legend = TRUE, legend="bottom")
```

<figure><img src="BYM_Model_files/figure-html/unnamed-chunk-14-1.png"><figcaption></figcaption></figure>

In this particular plot we see the changes in the GINI inequality
measures collected across the counties of California across 5 years.
This measure has been calculated for the years 2020, 2021 and 2022
therefore they do not provide independent information. [On the other
hand they are more reliable than one year
estimates](https://www.census.gov/programs-surveys/acs/guidance/estimates.html).

### Vaccinations vector [N] VacPop[T] ;

We also obtain daily fully vaccinated individuals' data in the counties
and convert it to biweekly rates of vaccinated population. This assumes
the location people get vaccinated is in the county of their residence.
The Covid-19 vaccine has been available only since the second half of December 2020^[https://www.latimes.com/california/story/2020-12-14/covid-19-vaccines-arrive-in-california] in California. The data from the Californian health services has a minor issue that it recorded a partially vaccinated person in January 8 2020 and some minor number of additional individuals after this date until the date of official start of vaccinations. We ignore this record, and only include this covariate from the biweek 25.

Most of the visualization code we present below is from earlier
[post](https://mmusal.github.io/blog/2023/Intro_to_Spatial_Tools_and_Files/#_32_The_final_multi-line_plot_via_ggplot)


```r
library(tidyr)
library(dplyr)
VacPop=read.table('C:/Users/rm84/Documents/VacPop.csv',header = TRUE,sep=",")
#Sort the name of the counties to make sure it merges to the correct counties.
names=sort(shape$NAMELSAD)
gvacPop1=as.data.frame(cbind(names,VacPop))
names(gvacPop1)[1]="names"
names(gvacPop1)[2]="fips"
T=ncol(VacPop)-1
for(c in 1:(T)){
  names(gvacPop1)[c+2]=paste0("Time_",c)
} 
gvacPop2=gvacPop1 %>% 
  pivot_longer(
    cols = starts_with("Time"), 
    names_to = "Time", 
    values_to = "Vac",
    values_drop_na = TRUE
  )

gvacPop2=select(gvacPop2, -c(fips))
gvacPop2=as.data.frame(gvacPop2)

gvacPop2$Time=rep(c(1:T),58)

maxtime=ncol(VacPop)-1
labels = data.frame(Vac = as.numeric(gvacPop1[which(gvacPop1[,maxtime+1]==max(gvacPop1[,maxtime+1])), maxtime+1])+0.02,Time = maxtime,text = paste0(names[which(gvacPop1[,maxtime]==min(gvacPop1[,maxtime]))]))

gvacPop2$Vac=as.numeric(gvacPop2$Vac)
maxlabel = data.frame(Vac = as.numeric(gvacPop1[which(gvacPop1[,maxtime+1]==max(gvacPop1[,maxtime+1])), maxtime+1])+0.02,
                      Time = maxtime,text = paste0(names[which(gvacPop1[,maxtime+1]==max(gvacPop1[,maxtime+1]))]))
#minlabel and maxlabel are created in order to demonstrate how to partially parameterize the annotations on the ggplot. minlabel will have the information regarding the county with the least amount of vaccination. maxlabel will have the information regarding the county with the most amount of vaccination

minlabel = data.frame(Vac = as.numeric(gvacPop1[which(gvacPop1[,maxtime+1]==min(gvacPop1[,maxtime+1])), maxtime+1])+0.02,Time = maxtime,text = paste0(names[which(gvacPop1[,maxtime]==min(gvacPop1[,maxtime]))]))

i3=ggplot(gvacPop2, aes(x = Time, y = Vac,color=names)) +
  geom_line()+  
  theme(legend.position = "none")+
  theme(axis.title.y =element_blank())+
  xlab("Biweek t")+
  annotate("text", x = maxtime-15, y = maxlabel[,1], label =maxlabel[,3])+
  annotate("segment", color="blue", x=maxtime-2, xend = maxtime-7, y=maxlabel[,1], 
           yend=maxlabel[,1], arrow=arrow(length=unit(0.2,"cm")))+
  annotate("text", x = maxtime-15, y = minlabel[,1], label = minlabel[,3])+
  annotate("segment", color="blue", x=maxtime-2, xend = maxtime-7, y=minlabel[,1], 
           yend=minlabel[,1], arrow=arrow(length=unit(0.2,"cm")))+
  ggtitle("Percent of Population Vaccinated")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_discrete(limits=c(1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,77))

i3
```

<figure><img src="BYM_Model_files/figure-html/unnamed-chunk-15-1.png"><figcaption></figcaption></figure>

As we can see in the labels of this plot there is a wide difference between the counties in terms of percentage of population vaccinated. 

### Poverty and Income

In this section we read in the poverty and income data between the years 2020 and 2022, downloaded from [Small Area Income and Poverty Estimates](https://www.census.gov/data-tools/demo/saipe/#/)(SAIPE)


```r
wd='C:/Users/rm84/Desktop/research/HMM/data'
setwd(wd)
#Read Downloaded Data - poverty 2019-2022
poverty<-read.table('saipepoverty.csv',header=TRUE,sep=',')
#First set of rows
head(poverty)
```

```
##   Year   ID             Name Percent.in.Poverty
## 1 2019 6001   Alameda County                8.9
## 2 2019 6003    Alpine County               17.2
## 3 2019 6005    Amador County                9.8
## 4 2019 6007     Butte County               16.1
## 5 2019 6009 Calaveras County               12.1
## 6 2019 6011    Colusa County               12.0
```

```r
#Read Downloaded Data - Income 2019-2022
income<-read.table('saipeincome.csv',header=TRUE,sep=',')
#First set of rows
head(income)
```

```
##   Year   ID             Name Median.Household.Income
## 1 2019 6001   Alameda County                  107589
## 2 2019 6003    Alpine County                   58112
## 3 2019 6005    Amador County                   62640
## 4 2019 6007     Butte County                   58394
## 5 2019 6009 Calaveras County                   68248
## 6 2019 6011    Colusa County                   59048
```

```r
# Specify number of characters to extract from county field that has leading 0s.
n_last <- 5   
# Number of Counties
N<-length(shape$GEOID)
#remove the year 2019 from the first N rows and combine . 
povertyandincome<-cbind(poverty[-c(1:N),],income[-c(1:N),])
#remove dupilicate and unnecessary columns
povertyandincome<-povertyandincome[-c(5,6,7)]

#since county ID values had leading 0s they were lost in reading
#we attach the 0s to the left of the ID values and use the substring function to obtain 5 columns of the ID values from the right.
povertyandincome$ID=paste('0',povertyandincome$ID,sep="")
povertyandincome$ID=substr(povertyandincome$ID, nchar(povertyandincome$ID) - n_last + 1, nchar(povertyandincome$ID))
#Quality assurance that all county ID s correspond to counties.
povertyandincome=povertyandincome[nchar(povertyandincome$ID)==5,]

#putting them in a form to export for STAN
povertyandincome1=povertyandincome %>%
  pivot_wider(names_from = Year, values_from = c(Percent.in.Poverty,Median.Household.Income))
head(povertyandincome1)
```

```
## # A tibble: 6 × 8
##   ID    Name             Percent.in.Poverty_2020 Percent.in.Poverty_2021
##   <chr> <chr>                              <dbl>                   <dbl>
## 1 06001 Alameda County                       8.6                     9.4
## 2 06003 Alpine County                       14.3                    15.8
## 3 06005 Amador County                       10.3                    11.1
## 4 06007 Butte County                        17.3                    16.6
## 5 06009 Calaveras County                    11.6                    13.5
## 6 06011 Colusa County                       10.3                    11.4
## # ℹ 4 more variables: Percent.in.Poverty_2022 <dbl>,
## #   Median.Household.Income_2020 <int>, Median.Household.Income_2021 <int>,
## #   Median.Household.Income_2022 <int>
```

```r
attach(povertyandincome1)

#Use ggplot as we did in GINI for poverty values facet wrap would be a better option to work with but this works well enough
map_Pov20=ggplot() +
  annotation_spatial(shapeanddata) +
  ggtitle("Poverty 2020")+
  theme(plot.title = element_text(hjust = 0.5))+
  layer_spatial(shapeanddata, aes(fill = (Percent.in.Poverty_2020)))+
  theme(legend.title= element_blank())+
  labs(fill = "Pov. 20")+scale_fill_viridis(limits = c(5,25),direction=-1)+
  theme(legend.position = "none")

map_Pov21=ggplot() +
  annotation_spatial(shapeanddata) +
  ggtitle("Poverty 2021")+
  theme(plot.title = element_text(hjust = 0.5))+
  layer_spatial(shapeanddata, aes(fill = (Percent.in.Poverty_2021)))+
  theme(legend.title= element_blank())+
  labs(fill = "Pov. 21")+scale_fill_viridis(limits = c(5,25),direction=-1)+
  theme(legend.position = "none")

map_Pov22=ggplot() +
  annotation_spatial(shapeanddata) +
  ggtitle("Poverty 2022")+
  theme(plot.title = element_text(hjust = 0.5))+
  layer_spatial(shapeanddata, aes(fill = (Percent.in.Poverty_2022)))+
  theme(legend.title= element_blank())+
  labs(fill = "Pov. 22")+scale_fill_viridis(limits = c(5,25),direction=-1)+
  theme(legend.position = "none")

ggarrange(map_Pov20,map_Pov21,map_Pov22,nrow=1,ncol=3,common.legend = TRUE, legend="bottom")
```

<figure><img src="BYM_Model_files/figure-html/unnamed-chunk-16-1.png"><figcaption></figcaption></figure>

```r
map_Inc20=ggplot() +
  annotation_spatial(shapeanddata) +
  ggtitle("2020")+
  theme(plot.title = element_text(hjust = 0.5))+
  layer_spatial(shapeanddata, aes(fill = (Median.Household.Income_2020/1000)))+
  scale_fill_viridis(limits = c(35,160),direction=-1)+labs(fill="Median Household Income")
  

map_Inc21=ggplot() +
  annotation_spatial(shapeanddata) +
  ggtitle("2021")+
  theme(plot.title = element_text(hjust = 0.5))+
  layer_spatial(shapeanddata, aes(fill = (Median.Household.Income_2021/1000)))+
  theme(legend.title= element_blank())+
  scale_fill_viridis(limits = c(35,160),direction=-1)+
  theme(legend.position = "none")

map_Inc22=ggplot() +
  annotation_spatial(shapeanddata) +
  ggtitle("2022")+
  theme(plot.title = element_text(hjust = 0.5))+
  layer_spatial(shapeanddata, aes(fill = (Median.Household.Income_2022/1000)))+
  theme(legend.title= element_blank())+
  scale_fill_viridis(limits = c(35,160),direction=-1)+
  theme(legend.position = "none")

ggarrange(map_Inc20,map_Inc21,map_Inc22,nrow=1,ncol=3,common.legend = TRUE, legend="bottom")
```

<figure><img src="BYM_Model_files/figure-html/unnamed-chunk-16-2.png"><figcaption></figcaption></figure>

### Median Age,Old-age dependency and Sex Ratio, vector [N] medage ; vector [N] medsex; 


```r
#American Community Survey. Does not have 2020 year therefore we have to use 5 year
#Median age (years)
#Sex ratio (males per 100 females)
#Age dependency ratio
#Old-age dependency ratio
#Child dependency ratio
collabels=c("Cnames","Year","Median_Age","Sex_Ratio","Age_Dependence","Old_Age_Dep","Child_Dep")

library(readxl)
```

```
## Warning: package 'readxl' was built under R version 4.3.2
```

```r
wd="C:/Users/rm84/Desktop/research/HMM/data/"
setwd(wd)
SummaryAgeSex20 <- read_excel("sexagesummary20.xlsx", 
                            sheet = "Data")
SummaryAgeSex21 <- read_excel("sexagesummary21.xlsx", 
                              sheet = "Data")
SummaryAgeSex22 <- read_excel("sexagesummary22.xlsx", 
                              sheet = "Data")

sagesex20=as.data.frame(matrix(nrow=nrow(shape),ncol=nrow(SummaryAgeSex20)))
sagesex21=as.data.frame(matrix(nrow=nrow(shape),ncol=nrow(SummaryAgeSex20)))
sagesex22=as.data.frame(matrix(nrow=nrow(shape),ncol=nrow(SummaryAgeSex20)))

I=nrow(shape)
j=0  
for(i in 1:I){
  i=6*(i-1)+1
  j=1+j
  sagesex20[j,]=as.numeric(t(SummaryAgeSex20[,i]))
  sagesex21[j,]=as.numeric(t(SummaryAgeSex21[,i]))
  sagesex22[j,]=as.numeric(t(SummaryAgeSex22[,i]))
}

sagesex20=cbind("2020",sagesex20)
sagesex21=cbind("2021",sagesex21)
sagesex22=cbind("2022",sagesex22)
cnames=read.table("countynames.txt",header = TRUE,sep="\t")[1:58,1]
sagesex20=cbind(cnames,sagesex20)
sagesex21=cbind(cnames,sagesex21)
sagesex22=cbind(cnames,sagesex22)
names(sagesex20)=collabels
names(sagesex21)=collabels
names(sagesex22)=collabels

sagesex=rbind(sagesex20,sagesex21,sagesex22)
sagesex$Year=as.numeric(sagesex$Year)
sagesex$Median_Age=as.numeric(sagesex$Median_Age)

maxchange_county_name=cnames[which(sagesex20$Median_Age/sagesex21$Median_Age==
max(sagesex20$Median_Age/sagesex21$Median_Age))]

maxchange_county_name
```

```
## [1] "Sierra"
```

```r
ggplot(aes(x = as.character(Year), y = Median_Age),data=sagesex)+geom_boxplot()+xlab("Year")
```

<figure><img src="BYM_Model_files/figure-html/unnamed-chunk-17-1.png"><figcaption></figcaption></figure>

```r
ggplot(aes(x = Year, y = Median_Age,color=Cnames),data=sagesex)+
  geom_line()+  
  theme(legend.position = "none")+
  theme(axis.title.y =element_blank())+
  scale_x_continuous(breaks=c(2020,2021, 2022))+
  xlab("Year")+ylab("Median_Age")+
  annotate("text", x = 2020.5, y = sagesex21[sagesex21[,1]==maxchange_county_name,]$Median_Age-2, label =maxchange_county_name)+
  annotate("segment", color="black", x=2020.5, xend = 2021, y=sagesex21[sagesex21[,1]==maxchange_county_name,]$Median_Age-1, 
           yend=sagesex21[sagesex21[,1]==maxchange_county_name,]$Median_Age, arrow=arrow(length=unit(0.2,"cm")))
```

<figure><img src="BYM_Model_files/figure-html/unnamed-chunk-17-2.png"><figcaption></figcaption></figure>

### Race matrix [R,N] race; row_vector[N] white=race[2];

Below we show the distribution of 6 races in the 58 counties of California. Note that due to sparsity of these races across all of the counties we will use only the White only race as a covariate in the analysis.  


```r
wd="C:/Users/rm84/Desktop/research/HMM/data/"
setwd(wd)
library(reshape)
race=read.table(file="DECENNIALPL2020.P2-2023-02-15T024354race.txt",sep="\t",header=TRUE)
head(race)[,1:2]
```

```
##                                              Label Alameda.County.California
## 1                               Hispanic or Latino                    393749
## 2                                      White alone                    472277
## 3                  Black or African American alone                    159499
## 4          American Indian and Alaska Native alone                      4131
## 5                                      Asian alone                    540511
## 6 Native Hawaiian and Other Pacific Islander alone                     13209
```

```r
#race is going to have 6 main races and 1 remainder
racenames=c(race[1:6,1],"rest")
race1=as.data.frame(matrix(nrow=nrow(race)-2,ncol=ncol(race)-1))

for(j in 2:ncol(race)){
  race1[1:6,j-1]=race[1:6,j]/race[nrow(race),j]
  race1[7,j-1]=1-sum(race1[1:6,j-1])
}

head(race1)[,1:3]
```

```
##            V1          V2          V3
## 1 0.234046600 0.069767442 0.148589218
## 2 0.280724081 0.665282392 0.734422098
## 3 0.094807095 0.008305648 0.030019272
## 4 0.002455489 0.177740864 0.014256066
## 5 0.321282751 0.009966777 0.013687800
## 6 0.007851503 0.000000000 0.001803627
```

```r
#transpose the dataframe
race2=t(race1)
head(race2)
```

```
##          [,1]      [,2]        [,3]        [,4]        [,5]        [,6]
## V1 0.23404660 0.2807241 0.094807095 0.002455489 0.321282751 0.007851503
## V2 0.06976744 0.6652824 0.008305648 0.177740864 0.009966777 0.000000000
## V3 0.14858922 0.7344221 0.030019272 0.014256066 0.013687800 0.001803627
## V4 0.18953655 0.6598766 0.015687609 0.014411809 0.048825319 0.002400393
## V5 0.12949307 0.7654332 0.007374371 0.010973240 0.015587742 0.001655922
## V6 0.61706122 0.3178259 0.008333715 0.012821100 0.011538990 0.003205275
##          [,7]
## V1 0.05883248
## V2 0.06893688
## V3 0.05722192
## V4 0.06926174
## V5 0.06948247
## V6 0.02921379
```

```r
#These labels are not official labels of the Decennial Census but are selected to #fit the visual screen.
racenames1=c("Hisp.","White","Black","Indian","Asian","Haw. or other Pac.","Rest")
colnames(race2)=racenames1
#We will be showing the distribution of 6 races, Hawaiian or other Pacific Islanders are going to be added to the other category.
race2[,7]=(race2[,6]+race2[,7])
race2=race2[,-6]

race3=melt(as.data.frame(race2))


ggplot(data=race3,aes(x=variable,y=value))+
  geom_boxplot()+
  ylab("Percentage")+
  theme(axis.title.x =element_blank())+
  ggtitle("Boxplots of Races 2020")+
  theme(plot.title.position = 'plot', 
      plot.title = element_text(hjust = 0.5))
```

<figure><img src="BYM_Model_files/figure-html/unnamed-chunk-18-1.png"><figcaption></figcaption></figure>



```stan
data { 
  int<lower=1> T;//number of time periods
  int<lower=1> R;//number of races specified in dataset
  int START; //time period analysis will start from 
  int END; // time period analysis will stop at
  int<lower=0> N;//number of spatial areas
  int<lower=0> N_edges; //number of shared borders between counties
  int<lower=1, upper=N> node1[N_edges];  // node1[i], node2[i] neighbors
  int<lower=1, upper=N> node2[N_edges];  // node1[i] < node2[i]
  int<lower=0>  y[T,N];              // count outcomes across T time periods
  real  x[167,N];              // covariates which also include interaction effects
  matrix [T,N] log_E; // logarithm of the expected mortality values if mortality risk was uniformly applied 
  vector [N] VacPop[T] ; // Vaccinated population percentage in the county. Changes with biweeks 
  real<lower=0> scaling_factor;  // scales the variance of the spatial effects
  vector [N] medage; // standardized median age
  vector [N] medsex; // sex (coded as dummy)
  matrix [R,N] race; // percentage of races in counties
}

transformed data {
  //Providing a matrix with a single index returns the specified row. For instance, if m is a matrix, then m[2] is the second row. 
//using the covariates matrix X to create specific matrices   
//t = 1 to 25 is 2020
//t = 26 to 51 is 2021
//t = 52 to 77 is 2022
matrix [T,N] poverty;
matrix [T,N] income;
matrix [T,N] povinc;
matrix [T,N] gini;
int ITER=END-START+1; 
  row_vector[N] white=race[2];
  row_vector[N] nonwhite=1-race[2];
poverty[26:51] =rep_matrix(to_row_vector(x[2]), 26);
poverty[52:77] =rep_matrix(to_row_vector(x[3]), 26);
income[1:25] =rep_matrix(to_row_vector(x[4]), 25);
income[26:51] =rep_matrix(to_row_vector(x[5]), 26);
income[52:77] =rep_matrix(to_row_vector(x[6]), 26);
povinc[1:25] =rep_matrix(to_row_vector(x[8]), 25);
povinc[26:51] =rep_matrix(to_row_vector(x[9]), 26);
povinc[52:77] =rep_matrix(to_row_vector(x[10]), 26);
gini[1:25] =rep_matrix(to_row_vector(x[11]), 25);
gini[26:51] =rep_matrix(to_row_vector(x[12]), 26);
gini[52:77] =rep_matrix(to_row_vector(x[13]), 26);
}

```

## STAN Parameter/Transformed Section

In this section we specify the parameters that we will keep track of.
The parameters that are to be dynamically estimated is in transformed
parameters section. We define $\sigma$ of \eqref{eq:lambda} as the
single dimensional array sigma with ITER elements,
ITER==(END-START)+1==77-1+1, making explicit the choice that we shall estimate
the parameter at every t. The parameter gamma_vac is going to be used to
estimate beta_vac the coefficient effect of vaccinations. The spatial
effect $\phi$ is a row vector with N (58) elements. The parameter which
quantifies the trade-off between spatial effect $\phi$ and $\theta$ is
created as an array of size ITER, rho ($\rho$), which will also be
estimated at every iteration t. The intercept parameter array, beta_int
has the same number of elements as rho array. The rest of the
coefficients are associated the covariates written next to "beta".\
The ITER (77) sized array theta contains a matrix of size 1 and N
(1,58). It represents $\theta$, random effects.

In the transformed parameters section, the first object created is the
object alpha. This is an array of ITER size where each element is a
matrix of "1,N" (1,58) dimensions. This object is where the equation
\eqref{eq:lambda} is going to be evaluated. The beta_vac array is an
array of size (T-V+1) and is calculated as 
$$
\beta_{vac_{1}}=\gamma_{1} \\
\beta_{vac_{t}}=\beta_{vac_{t-1}}+\gamma_{t}
\label{eq:vac}
$$

The for loop involving beta_vac is the one associated with
\eqref{eq:vac}. The next for loop calculates convolved_re

```{=tex}
\begin{equation}
\left(\sqrt{\frac{\rho}{s}}\times \phi_{i}+\sqrt{(1-\rho)}\times \theta_{i}\right)\times \sigma
\end{equation}
```
which is the combination of spatial and random effect in the
\eqref{eq:lambda} and is evaluated in to the matrix object with ITER and
N rows, columns respectively.

In the final stage alpha array of ITER size where every element is a
matrix of 1,N (1,58) dimensions is evaluated.


```stan
parameters {
  real <lower = 0> sigma [ITER];
 real gamma_vac[END-V+1];
   row_vector[N] phi;
  real<lower=0,upper=1> rho[ITER];
  matrix [ITER,N] theta;
real beta_int[ITER];
real beta_pov;
real beta_inc;
real beta_dens;
real beta_gini;
real beta_age;
real beta_sex;
real beta_white;
}
transformed parameters {
  matrix[1,N] alpha [ITER] ;
  real beta_vac[T-V+1];


beta_vac[1]=gamma_vac[1];    
  
    matrix [ITER,N] convolved_re;
   
    for(t in (2):(T-V+1)) {
      beta_vac[t]=beta_vac[t-1]+gamma_vac[t];}
      
for(t in START:V){
      convolved_re[t-START+1] = (sqrt(rho[t-START+1] / scaling_factor) * phi + sqrt(1 - rho[t-START+1]) * theta[t-START+1])*sigma[t-START+1]; 
      alpha[t-START+1][1,1:N]=exp(
                                  log_E[t,1:N]+
                                  to_row_vector(convolved_re[t-START+1])+
                                    beta_int[t-START+1]+
                                    beta_pov*poverty[t]+
                                    beta_inc*income[t]+
                                    beta_dens*to_row_vector(x[7])+
                                    beta_gini*gini[t]+
                                    beta_age*to_row_vector(median_age[t])+
                                    beta_sex*to_row_vector(sex[t])+
                                    beta_white*to_row_vector(white)
                                    ); 
                                   
    }
    for(t in (V+1):END){
      convolved_re[t-START+1] = (sqrt(rho[t-START+1] / scaling_factor) * phi + sqrt(1 - rho[t-START+1]) * theta[t-START+1])*sigma[t-START+1]; 
      alpha[t-START+1][1,1:N]=exp(
                                  log_E[t,1:N]+
                                  to_row_vector(convolved_re[t-START+1])+
                                    beta_int[t-START+1]+
                                    beta_pov*poverty[t]+
                                    beta_inc*income[t]+
                                    beta_dens*to_row_vector(x[7])+
                                    beta_gini*gini[t]+
                                    beta_vac[t-V+1]*to_row_vector(VacPop[t-1])+
                                    beta_age*to_row_vector(median_age[t])+
                                    beta_sex*to_row_vector(sex[t])+
                                    beta_white*to_row_vector(white)
                                    ); 
                                   
    }
}
 
```

## STAN Model Section  

STAN utilizes the model section to write the likelihood and the priors for the parameters. The likelihood itself is based on \eqref{eq:lambda}. Due to the lower limit posted in the [parameters section](https://mmusal.github.io/blog/2023/Joint_Spatial_Effects_BYM/#_44_STAN_ParameterTransformed_Section) the normal prior will evaluate to half normal. The prior for the parameter rho is a relatively informative beta distribution with the hyper-parameters 2 and 2.  The covariate coefficients priors all have a normal distribution with 0 mean and 1 standard deviation. These hyper-priors are not as informative as they might look at a first glance since the model is a log-linear one and the covariates are standardized.   

If we did not have vaccination variable the multiplication of the likelihood and the priors would be written as \eqref{eq:olikelihoodprior}.

```{=tex}
\begin{equation}
\prod_{t=1}^{t=T} \bigg[  \prod_{n=1}^{n=N} \bigg[ Poisson(y_{t,n}|\alpha_{t,n}) \times p(\theta_{t,n}) \bigg] \times p(\rho_{t}) \times p(\sigma_{t}) \times p(\gamma_{t})\bigg] \\ \times p(\phi_{1}, \ldots , \phi_{N}) \times p(\beta_{Pov}) \times \ldots \times p(\beta_{white}) 
\label{eq:olikelihoodprior}
\end{equation}
```

Since we have a vaccination covariate which became a factor at the end of December 2020, we will separate out the likelihood in two adding the effect of vaccines after t=24. Since we want to make our program as flexible as possible  we use the constant V (t=24), START (t=1) and END (t=77) as part of data.   

```{=tex}
\begin{equation}
\prod_{t=START}^{t=V} \bigg[  \prod_{n=1}^{n=N} \bigg[ Poisson(y_{t,n}|\alpha_{t,n}) \times p(\theta_{t,n}) \bigg] \bigg] \\
\prod_{t=V+1}^{t=END} \bigg[  \prod_{n=1}^{n=N} \bigg[ Poisson(y_{t,n}|\alpha_{t,n}) \times p(\theta_{t,n}) \bigg] \times p(\rho_{t}) \times p(\sigma_{t}) \times p(\gamma_{t})\bigg] \\ \times p(\phi_{1}, \ldots , \phi_{N}) \times p(\beta_{Pov}) \times \ldots \times p(\beta_{white}) 
\label{eq:newlikelihoodprior}
\end{equation}
```
Note that sigma and Gamma is an array with (T-V+1) and T elements respectively where V is 24. This allows us to write the priot for these parameters using the $\sim normal(0,1)$ notation rather than writing a loop as we end up doing for theta. The parameter theta is a matrix of ITER rows and N columns as declared in the parameters section.      


```stan
model {

target += poisson_lpmf(y[1,1:N] | alpha[1][1,1:N]);
.
.
.
target += poisson_lpmf(y[77,1:N] | alpha[77][1,1:N]);


phi~icar_normal(N,node1,node2);

beta_int~normal(0,1);
beta_pov~normal(0,1);
beta_inc~normal(0,1);
beta_dens~normal(0,1);
beta_gini~normal(0,1);
beta_age~normal(0,1);
beta_sex~normal(0,1);
beta_white~normal(0,1);
gamma_vac~normal(0,1);
sigma~normal(0,1);
phi~icar_normal(N,node1,node2);
rho ~ beta(2, 2);

  for(t in START:END){
    theta[t-START+1] ~ normal(0,1);
    }

```


# References
