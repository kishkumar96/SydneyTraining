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

```javascript
var place = 'test';
var group = '';

//Select roi
var roi = test;

//Add roi to the map. Using the inspector tab allows name and id to be verified by clicking.
Map.addLayer(roi,{},'ROI');

//Center the map view on the ROI
Map.centerObject(roi);
```


# Spatial and temporal filtering 

```javascript
// Select date 
var year = ['2021','2021'];
print('Year: ',year[0]);
print(roi);
var start = ee.Date(year[0]+'-01-01');
var end = ee.Date(year[1]+'-12-31');

var LS8fmaskcol =  ee.ImageCollection('COPERNICUS/S2_SR')
  .filterDate(start,end)
  .filterBounds(roi);
  //.filter(ee.Filter.neq('LANDSAT_SCENE_ID', 'LC81220482016022LGN00'));
print(LS8fmaskcol);
print('Precollection TOA fmask size: ', LS8fmaskcol.size());


/* T1 to account for newer imagery only being available within the new collection system.
RT ensures that the latest possible imagery is always being used*/
var ls8t1col = ee.ImageCollection('COPERNICUS/S2') 
  .filterDate(start,end)
  .filterBounds(roi);
print('T1 size',ls8t1col.size());

// Load Sentinel-1 C-band SAR Ground Range collection (log scale, VV,
//descending)
var collectionVV = ee.ImageCollection('COPERNICUS/S1_GRD')
.filter(ee.Filter.eq('instrumentMode', 'IW'))
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
.filterMetadata('resolution_meters', 'equals' , 10)
.filterBounds(roi)
.select('VV');
print(collectionVV, 'Collection VV'); 



// Load Sentinel-1 C-band SAR Ground Range collection (log
//scale, VH, descending)
var collectionVH = ee.ImageCollection('COPERNICUS/S1_GRD')
.filter(ee.Filter.eq('instrumentMode', 'IW'))
.filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
.filterMetadata('resolution_meters', 'equals' , 10)
.filterBounds(roi)
.select('VH');
print(collectionVH, 'Collection VH');

//Filter by date
var SARVV = collectionVV.filterDate('2019-01-01', '2019-12-10').mosaic();
var SARVH = collectionVH.filterDate('2019-01-01', '2019-12-10').mosaic();
```
