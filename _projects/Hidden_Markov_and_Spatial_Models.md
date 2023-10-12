---
layout: page
title: Hidden Markov and Spatial Models
description: This project explores the multiple neccessary components to create a composite of Bayesian Spatial and Hidden Markov Models. Furthermore the project is useful to see more of STAN programming language applications.
toc: true
keep_md: no
toc_float: true
toc_collapsed: true
toc_depth: 2
number_sections: true
img:
importance: 1
category: work
---

In this project we would like to investigate the affect of socio-demographic and other variables on Covid-19 mortality in the counties of California. As we can see from the [application created by rshiny, authored by Anjali Goel,](https://mmusal.shinyapps.io/timeseriesofSMR/) counties seem to be following a degree of spatial structure where counties with low/high SMR values have counties with similar SMR characteristics. We discuss the [details](https://mmusal.github.io/blog/2023/Explaining_rshinyapp/) of this application. We conduct Moran's I as we describe the data. However before doing so we should learn about [some of the capabilities](https://mmusal.github.io/blog/2023/Intro_to_Spatial_Tools_and_Files/) that R has in visualizing maps using SHP data. The questions we would like to eventually answer are listed below:

## Project Questions
0. Can we represent the spatial dependence in the counties of California.

1. Among the available socio-demographic and other factors recorded for each county in California, which of them affected Covid-19 mortality.

2. To which extent were these socio-demographics and other factors affected the Covid-19 mortality counts.

3. How did these effects change through time (time is measured in biweeks).

4. Can we make models more robust by reducing the number of effects estimated per variable and tie it to policies. 

5. Can we represent different pandemic regimes in California and tie it to particular policies. 


To accomplish this we will need to discuss a set of concepts outlined [below](#set_of_concepts). 
First we need to have an understanding of the scripting language R to represent the problem visually and prepare spatial neighborhood structures. We will discuss Hidden Markov Models and Spatial models separately, and give necessary details and background to understand each necessary to create the combination.   

## <a name="set_of_concepts"></a> Set of concepts 

### [Mapping/Spatial Tools](https://mmusal.github.io/blog/2023/Intro_to_Spatial_Tools_and_Files/)    

### [Hidden Markov Model (HMM)](https://mmusal.github.io/blog/2023/Hidden-Markov-Models/)

#### Introduction to some probabilty concepts

#### Using probability concepts to describe HMM

### Besag York Mollie Model

### STAN Programming Language

#### Besag York Mollie Spatial Models with STAN

#### Hidden Markov Models with STAN

### Combining the Models
