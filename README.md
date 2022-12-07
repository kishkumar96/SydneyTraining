# SydneyTraining

The aim of this document is to share a landcover classification methodology

# Table of Contents
- Introduction
- Creating composite images 
- Composite classification
- Accuracy assessment
- Post classification fitering 
- Deriving area measurements 

# Introduction 
The resilient boundaries for the blue pacific project is a research project that is looking at the implications that sea level rise will have on the maritime boundaries of the pacific. Remote sensing is one of the main methodologies used to extract data for the project. 

The project uses a fusion of optical (Sentinel 2) and synthetic aperture radar (sentinel 1) to conduct landcover classificaiton wherever the data is available. The rationale behind this is that fusion of the two datasets produces higher accuracy in classification for sand and urban classes whereby it also solves the cloud problem of classification using optical data only.

# Creating composite images 
###### [Working example of creating image composites using Sentinel 2 and Sentinel 1] (https://code.earthengine.google.com/?scriptPath=users%2Fkishan2196%2FRBBPScripts%3ALandcover%20Scripts%2FScript1_OpticalAndSAR%20indice%20calculation )

'''javascript
var place = 'Tonga';
var group = '';

//Select roi
var roi = Tongatapu;

//Add roi to the map. Using the inspector tab allows name and id to be verified by clicking.
Map.addLayer(roi,{},'ROI');

//Center the map view on the ROI
Map.centerObject(roi);
'''
