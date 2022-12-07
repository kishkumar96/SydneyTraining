# SydneyTraining

This document aims to share a landcover classification methodology

# Table of Contents
- Introduction
- Creating composite images 
- Composite classification
- Accuracy assessment
- Post-classification filtering 
- Deriving area measurements 

# Introduction 
The resilient boundaries for the blue pacific project is a research project that is looking at the implications that sea level rise will have on the maritime boundaries of the pacific. Remote sensing is one of the main methodologies used to extract data for the project. 

The project uses a fusion of optical (Sentinel 2) and synthetic aperture radar (sentinel 1) to conduct landcover classification wherever the data is available. The rationale behind this is that fusion of the two datasets produces higher accuracy in classification for sand and urban classes whereby it also solves the cloud problem of classification using optical data only.

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

Being a cloud-based remote sensing platform google earth engine makes available to the user the historical data for the individual collections. To ensure that you are working with the right data it is important to carry out filtering. Earth engine enables spatial and temporal level filtering. The code below queries Sentinel 2 data particularly its S2 and S2_SR collections. The data is filtered at a temporal level and at a spatial level. 
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

The available varied collections on earth engine allow filtering unique to their metadata. The code below queries Sentinel 1, in terms of polarizations focusing on VV and VH at 10-meter level resolution. 
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
# Composite Classification
###### [Script to classify] (https://code.earthengine.google.com/?scriptPath=users%2Fkishan2196%2FRBBPScripts%3ALandcover%20Scripts%2FScript2_Image%20classification)
```javascript
var year ='2021';
print('Year: '+year);
var roi = test
//Select the image to be classified, and the true-colour image to use as a reference
var place = 'test'
var toClassify = image6
toClassify = toClassify.clip(roi);
```
# Creating Training Samples
Training samples generally take the form of polygonal geometries within which pixels are sampled. As such, they can be created within GEE using the built-in geometry tools, or imported from existing shapefiles etc. The important thing here is that the training samples have a property which describes which class they represent. When using the GEE geometry tools, this can be achieved by clicking the cog icon next to the layer in the geometry list at the top left corner of the map view and selecting to import as a feature or featureCollection (I found that featureCollection works well in this instance). This will allow you to add a property to the training sample, such as 'class' etc. Generally, classes are coded by number (i.e. water = 0, vegetation = 1, urban = 2). This allows easy visualisation and later accuracy assessment. The image below shows a composite image with several training polygons representing different classes (e.g. black for water, green for vegetation).

#  Sampling the Image
```javascript

// Produce median images to sample from
var NDmedian= ee.ImageCollection([
    //SuvaModified, OnoILau1, TuvananIRa1, Qelelevu1, Sigatoka
    OnoSar, QelelevuSar, SigatokaSar, SuvaSar, TuvanaIRaSar, FijiSar2019
  ]).median()


 
 var bands = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'B4', 'R10', 'R11', 'R12', 'R13', 'R14', 'R15', 'NDVI', 'NDMI', 'MNDWI', 'SR', 'R54', 'R35']//, 'VH', 'VV']
 
var classes = Water.merge(CoralReef).merge(Urban).merge(Vegetation).merge(Sand)
var image = NDmedian.select(bands)
// Sample images using training samples
var samples = image.sampleRegions({
	collection: classes,
	properties: ['landcover'],
	scale: 10
}).randomColumn('random');

//var training2 = NDmed
```
# Choosing and training a classifier
```javascript

// Combine required training samples into one combined sample set
//var join = training.merge(training2)

var classifier = ee.Classifier.libsvm();
var split = 0.8;
var training = samples.filter(ee.Filter.lt('random', split)); // subset training data
var testing = samples.filter(ee.Filter.gte('random', split)); // subset testing data 

print('Samples n =', samples.aggregate_count('.all'));
print('Training n =', training.aggregate_count('.all'));
print('Testing n = ', testing.aggregate_count('.all'));

// Train the chosen classifier 
var fullClassifier = classifier.train({
  features:training,//.select(['R1', 'R2', 'R3', 'R4', 'R6', 'R7', 'R10', 'landcover']), 
  classProperty: 'landcover', 
  inputProperties: bands //['R1', 'R2', 'R3', 'R4', 'R6', 'R7', 'R10']//, 'R11', 'R12', 'R13', 'R14', 'R15']
});

var validation = testing.classify(fullClassifier);
var testAccuracy = validation.errorMatrix('landcover', 'classification');
print('Validation error matrix RF:', testAccuracy);
print('Validation overall accuracy RF:', testAccuracy.accuracy());
// Classify the images.
var classified = toClassify.classify(fullClassifier);
```

# Reducing Noise and Masking 
```javascript

//The model results may be "noisy". To reduce noise, create a mask to mask
// unconnected pixels
var pixelcount = classified.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask = pixelcount.select(0).gt(25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask = classified.select('classification').gt(1)
var land= classified.updateMask(countmask).updateMask(classMask)
var CoralclassMask = classified.select('classification').lt(2)
var Coral= classified.updateMask(countmask).updateMask(CoralclassMask)
```
Visualising the classification

# Accuracy assessment
```javascript

var validation = testing.classify(fullClassifier);
var testAccuracy = validation.errorMatrix('landcover', 'classification');
print('Validation error matrix RF:', testAccuracy);
print('Validation overall accuracy RF:', testAccuracy.accuracy());
// Classify the images.
var classified = toClassify.classify(fullClassifier);

//The model results may be "noisy". To reduce noise, create a mask to mask
// unconnected pixels
var pixelcount = classified.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask = pixelcount.select(0).gt(25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask = classified.select('classification').gt(1)
var land= classified.updateMask(countmask).updateMask(classMask)
var CoralclassMask = classified.select('classification').lt(2)
var Coral= classified.updateMask(countmask).updateMask(CoralclassMask)
// Compare the ground truth and the working classification via an error matrix 
var errorMatrix = validation.errorMatrix('landcover','classification');

// Print the accuracy results to the console
print('Confusion table:', errorMatrix);
print('Accuracy: (correct/total)', errorMatrix.accuracy());
print('Consumer\'s accuracy (comission) (across):', errorMatrix.consumersAccuracy());
print('Producer\'s accuracy (omission) (down):', errorMatrix.producersAccuracy());
// Produce and print an error matrix of the binary results
var errorMatrixBinary = binarySample.errorMatrix('landcover', 'remapped')
print('Binary Accuracy:', errorMatrixBinary.accuracy());
print('Binary Confusion table:', errorMatrixBinary);
// Visualise the resulting classification
```
```javascript

Export.image.toAsset({
  image: land, 
  description: year+place+'Land',
  assetId: place+'/landcover/'+year+place+'Land',
  region: roi,//.roi().bounds(), 
  scale: 10, 
  maxPixels: 1e13,
  pyramidingPolicy: {".default": "mode"},
});

Export.image.toAsset({
  image: Coral, 
  description: year+place+'Coral',
  assetId: place+'/landcover/'+year+place+'Coral',
  region: roi,//.roi().bounds(), 
  scale: 10, 
  maxPixels: 1e13,
  pyramidingPolicy: {".default": "mode"},
});
Export.table.toDrive({
	collection: aaPoints, 
	description:'aa_points', 
	folder:'seperate_outputs', 
	fileFormat :'SHP'})
  ```
  
