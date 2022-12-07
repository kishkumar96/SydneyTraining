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


# Spatial and Temporal filtering 

Being a cloud based remote sensing platform google earth engine makes available to the user the historical data for the individual collections. In order to ensure that you are working with the right data it is important to carry out filtering. Earth engine enables spatial and temporal level filtering. The code below queries Sentinel 2 data particularly its S2 and S2_SR collections. The data is filtered at a temporal level and at a spatial level. 
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
print(LS8fmaskcol);
print('Precollection TOA fmask size: ', LS8fmaskcol.size());


/* T1 to account for newer imagery only being available within the new collection system.
RT ensures that the latest possible imagery is always being used*/
var ls8t1col = ee.ImageCollection('COPERNICUS/S2') 
  .filterDate(start,end)
  .filterBounds(roi);
print('T1 size',ls8t1col.size());

```

The available varied collections on earth engine allow filtering unique to there individual metadata. The code below queries Senitnel 1, in terms of polarizations focusing on VV and VH at 10 meter level resolution. 
```javascript
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
# Preprocessing 
```javascript
function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000);
}

//Apply filter to reduce speckle
var SMOOTHING_RADIUS = 50;
var SARVV_filtered = SARVV.focal_mean(SMOOTHING_RADIUS, 'circle', 'meters');
var SARVH_filtered = SARVH.focal_mean(SMOOTHING_RADIUS, 'circle', 'meters');
```

# Combining Image Collections 
```javascript

//Perform required filtering, collection reductions and band selection.
var fmaskFilt = LS8fmaskcol.map(maskS2clouds).select('B2', 'B3', 'B4', 'B8',  'B11', 'B12');
var fmaskUnfilt = LS8fmaskcol.select('B2', 'B3', 'B4', 'B8',  'B11', 'B12');
var medianFmasked = fmaskFilt.median().clip(roi).select('B2', 'B3', 'B4', 'B8',  'B11', 'B12');
var t1Filt = ls8t1col.map(maskS2clouds).select('B2', 'B3', 'B4', 'B8',  'B11', 'B12');
var t1Unfilt = ls8t1col.select('B2', 'B3', 'B4', 'B8',  'B11', 'B12');
var t1median = t1Filt.median().clip(roi).select('B2', 'B3', 'B4', 'B8',  'B11', 'B12');
var minUnmasked = ls8t1col.min().clip(roi).select('B2', 'B3', 'B4', 'B8',  'B11', 'B12');


var mergedFilt = ee.ImageCollection(t1Filt.merge(fmaskFilt));
print('merged, filtered T1 + Fmask TOA size: ',mergedFilt.size());

var mergedUnfilt = ee.ImageCollection(fmaskUnfilt.merge(t1Unfilt));
var mergedMedian = mergedFilt.median().clip(roi);
var mergedUnfilt = mergedUnfilt.median().clip(roi);

//This code fills in gaps with the median (or minimum) unfiltered pixels 
var fillerMask = mergedMedian.unmask().not();
var filler = mergedUnfilt.updateMask(fillerMask);
var medT1 = (mergedMedian.unmask().add(filler.unmask()));


//They are filled in the next step */
medT1 = medT1//.clip(roi
//.filter(ee.Filter.neq('Id',9))
//.filter(ee.Filter.neq('Id',10))
//.filter(ee.Filter.neq('Id',11))
//);

// This code fills in gaps with the t2 median pixels 
var fillerMask = medT1.unmask().not();
//var filler = t2median.updateMask(fillerMask);
var final = (medT1.unmask().add(filler.unmask())).clip(roi);
```

# Visualising results 
```javascript
var ls8viz = {bands: 'B4,B3,B2', gamma: 2};
var visualization = {
  min: 0.0,
  max: 0.3,
  bands: ['B4', 'B3', 'B2'],
};

Map.addLayer(final, visualization, 'As below, GF with T2 median (final)');
var image = final//.clip(roi.filter(ee.Filter.eq('Id',24)));

```

# Normalised Difference Spectral Vector
```javascript
var toNDSVLS8 = function(image){
image = image.select(
  ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'],
  ['B1','B2', 'B3', 'B4', 'B5', 'B6']
  );
var ndsv = image.normalizedDifference(['B1','B2']).rename('R1');
ndsv = ndsv.addBands(image.normalizedDifference(['B1','B3']).rename('R2'));
ndsv = ndsv.addBands(image.normalizedDifference(['B1','B4']).rename('R3'));

ndsv = ndsv.addBands(image.normalizedDifference(['B2','B3']).rename('R4'));
ndsv = ndsv.addBands(image.normalizedDifference(['B2','B4']).rename('R5'));
ndsv = ndsv.addBands(image.normalizedDifference(['B5','B1']).rename('R7'));
ndsv = ndsv.addBands(image.normalizedDifference(['B5','B2']).rename('R8'));
ndsv = ndsv.addBands(image.normalizedDifference(['B5','B3']).rename('R9'));
ndsv = ndsv.addBands(image.normalizedDifference(['B5','B4']).rename('R10'));
ndsv = ndsv.addBands(image.normalizedDifference(['B5','B6']).rename('R11'));
ndsv = ndsv.addBands(image.normalizedDifference(['B6','B1']).rename('R12'));
ndsv = ndsv.addBands(image.normalizedDifference(['B6','B2']).rename('R13'));
ndsv = ndsv.addBands(image.normalizedDifference(['B6','B3']).rename('R14'));
ndsv = ndsv.addBands(image.normalizedDifference(['B6','B4']).rename('R15'));
ndsv = ndsv.addBands(image.normalizedDifference(['B4','B3']).rename('NDVI')); //NDVI
  // normalized difference moisture index
ndsv = ndsv.addBands(image.normalizedDifference(['B6','B2']).rename('NDMI'));
ndsv = ndsv.addBands(image.normalizedDifference(['B2','B5']).rename('MNDWI'));// Modified Normalized Difference Water Index
ndsv = ndsv.addBands(image.select('B4').divide(image.select('B3')).rename('SR')); // Simple Ratio
ndsv = ndsv.addBands(image.select('B5').divide(image.select('B4')).rename('R54'));// Ratio 54
ndsv = ndsv.addBands(image.select('B3').divide(image.select('B5')).rename('R35'));
ndsv = ndsv.addBands(image.expression('(NIR/GREEN)-1',{ 
  'NIR':image.select('B4'), 
  'GREEN':image.select('B2')}));
ndsv = ndsv.addBands(image.normalizedDifference(['B3','B4']).rename('R6'));

return ndsv.clip(roi)
};

var nd = toNDSVLS8(image);
```
# Data Fusion 
```javascript
var opt_sar = ee.Image.cat(nd, SARVV_filtered,SARVH_filtered);
```

# Image Export 
```javascript

Export.image.toAsset({
  image: opt_sar, 
  description: year[0]+'opt_sar'+place+'_clipped',
  assetId: place+'/NDSV/'+year[0]+'opt_sar'+place+'_clipped',
  region: roi,//.geometry().bounds(), 
  scale: 10, 
  maxPixels: 1e13
});

Export.image.toAsset({
  image: image.select('B2', 'B3', 'B4', 'B8'),
  description: year[0]+'_S2_Full_median_'+place+'_clipped',
  assetId: place+'/full_composite/S2_'+year[0]+'_new_median_full_composite_'+place+'_clipped', 
  region: roi,//.geometry().bounds(), 
  scale: 10, 
  maxPixels: 1e13});
  
  Export.image.toAsset({
  image: nd, 
  description: year[0]+'_S2_NDSV_'+place+'_clipped',
  assetId: place+'/NDSV/'+year[0]+'_S2_NDSV_'+place+'_clipped',
  region: roi,//.geometry().bounds(), 
  scale: 10, 
  maxPixels: 1e13
});
```
