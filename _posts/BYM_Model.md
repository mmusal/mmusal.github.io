---
layout: post
title: "Spatial Relations and Besag-York-Mollie Framework"
author: "Rasim M Musal"
date: "2023-10-07"
output:
  html_document:
   theme: darkly
   highlight: espresso
   toc: true
   keep_md: yes
   toc_float: true
   toc_collapsed: true
   toc_depth: 1
   number_sections: true
   fig_caption: yes
tags: [Besag-York-Mollie, Areal Data]
always_allow_html: true
---




# Models Using Areal Data

When determining the presence of spatial relationships, a first step will be to determine whether there is reason to suspect that data that is in close proximity shows a recognizable pattern such that data which are closer to each other are more similar than data that are further apart. There are different frameworks in constructing spatial models. In this project we are interested in using areal data to create this framework.

Areal data, can be defined as data which is observed within well defined boundaries. For instance, square grids, zip code tabulation areas, or as in our project, counties in a particular state are all examples of well defined boundaries. Once we decide to use Areal data we need to assume a neighborhood structure. Figure 1 below will help illustrate this process. 




<figure><img src="assets/img/countiesofcalifornia.png"><figcaption></figcaption></figure>
<a name="cal_map"></a><span style="color:white; text-align:center; display:inline-block; width:100%;">Figure 1: FIPS numbers are recorded on each county</span>

Turn your attention to San Bernardino. It is surrounded by counties with Federal Information Processing Standard (FIPS) codes 27, 29, 37, 59, 65. These in turn are effected by other counties. Therefore there is a need to define what constitute a neighbor for a county and the weight of that neighborhood. These are all questions that can be quantified via a neighborhood matrix. There are different schemes in literature. However the Bayesian model that will follow is rather rich in terms of parameters we will use a binary neighborhood structure where if a county shares a border with another county its' effect is to be estimated. We will explain how to create and represent this matrix in STAN.

# The Spatial Model

There will be 

