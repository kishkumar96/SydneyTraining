# SydneyTraining

This document aims to share a landcover classification methodology. 

Holdaway, A., Ford, M., & Owen, S. (2021). Global-scale changes in the area of atoll islands during the 21st century. Anthropocene, 33(2020), 100282. https://doi.org/10.1016/j.ancene.2021.100282
![image](https://user-images.githubusercontent.com/53712176/206316052-0ec32b55-df72-45a8-b33e-e278cf914378.png)


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
###### [Working example of creating image composites using Sentinel 2 and Sentinel 1] (https://code.earthengine.google.com/393ffde48d79e37d6beaf8e5be09dd16)

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
Preprocessing enables the removal of noise from the data and enables it to conform to the preferred state. Hence, since the collections have been defined and filtered we can now preprocess them. For Sentinel 2 data preprocessing is done through cloud masking whereas for Sentinel 1 data speckle filtering is applied. 
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
Now it starts getting a bit more technical. While you may wish to use only one image collection in your workflow for the sake of simplicity, better results may be possible by combining multiple collections (i.e. S2_SR, S2, masked and unmasked) to achieve maximum coverage and image quality. In the first code snippet the gaps in the filtered S2 collection (where there have been clouds detected for a pixel representing the same geographic location in every image in the collection) are filled by using an unfiltered median or min composite. This ensures no gaps persist, but the trade-off is that cloud artefacts may persist in the final composite.
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
The normalised difference spectral vector composite is a computed image where every band of the image is a computed remote sensing indices. Using the normalised difference spectral vector (NDSV) approach improved classification performance in this study. This involves producing a pseudo multispectral image from all possible unique band ratios.
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
To enable classification using Sentinel 2 and Sentinel 1 together the merging of the two into one entity is important. Hence the NDSV composite is combined into one composite using ee.Image.cat() function. 
```javascript
var opt_sar = ee.Image.cat(nd, SARVV_filtered,SARVH_filtered);
```

# Image Export 
Depending on the size of your study area and the number of scenes being reduced, it make take some considerable time for GEE to process the final composite. You may also notice that the composite takes a while to reload when the zoom level is changed - this is because GEE processes at the scale set by the zoom level - essentially a level in the pyramid approach common to many GIS platforms.

Any subsequent calculations that rely on the final composite will also be slow, since it will need to be computed beforehand. Some more complex calculations, such as classification, may not work at all, timing out or running over the GEE user memory limit. This issue becomes magnified when trying to deal with multiple composites covering different date ranges - clearly processing multiple composites within the same script would be difficult and highly inefficient. To address this issue, images and features that are created within GEE scripts may be exported as a GEE asset. After export, assets may be imported into a script from the assets tab (on the left of the window by default). Imported assets perform much better in complex calculations, as less processing needs to be done 'on the fly. Both images and features may also be exported to the drive of the google account associated with GEE, allowing data produced in GEE to be used outside of it in the traditional matter. The following code snippet shows how the final composite (finalComp) is exported as a GEE asset. There are a few things to note here: the Export.image.toAsset() function takes several arguments, but not all are required, and some (such as scale) make others redundant. In cases like these, it can be easier to use curly brackets { } within the function brackets ( ). This allows the arguments to be specifically called by name and followed by: before answering the argument as normal. This can also make it clearer what each argument in a complex function is doing, improving readability.
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
Now that the composite imagery has been generated and saved as assets, they can be classified. Classification involves using a special algorithm (a classifier) to determine which of a user defined group of classes each pixel is most likely to represent. In this case, decisions are based upon the spectral values of each pixel (per band) after the classifier has been trained using a labelled dataset of representative pixels. For more information on classification within GEE, see this GEE Classification video tutorial. For the purposes of this tutorial, a single date classification (training data sampled from one image) will be prepared initially, then multi-date classification will be discussed.

For ease of use, I created a new script for classification, keeping it separate from the code which produces the composites detailed above. To begin, in a new script the previously generated composite NDSV image was defined as the variable toClassify and clipped to the ROI (the same as the previous script). Note that clipping does not carry over in exported assets: the area clipped out when generating the composite will be all black (i.e. null), but for visualistion purposes it is best to clip this off by clipping the image again. This code also adds the image to be classified to the map view (Map.addLayer(...)).
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
```javascript

var Training = 
    /* color: #d63000 */
    /* shown: false */
    ee.FeatureCollection(
        [ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45345806569367, -18.15103689128319],
                  [178.46174026506827, -18.141126823720917],
                  [178.46505706043231, -18.13372691935029],
                  [178.47016991101458, -18.13022064064236],
                  [178.47378054953364, -18.126102515093372],
                  [178.47766356877506, -18.11993924255189],
                  [178.4763336015281, -18.11916523499177],
                  [178.47556154406374, -18.118146504619368],
                  [178.47401506336402, -18.117858671730193],
                  [178.47246990991488, -18.118551971125033],
                  [178.4713968267783, -18.11765454983252],
                  [178.47002343207436, -18.118021665526783],
                  [178.46903627601804, -18.117573022802194],
                  [178.46822044619594, -18.11753243360895],
                  [178.46622525612565, -18.119474754854004],
                  [178.46553892362925, -18.123354444090186],
                  [178.4666979147124, -18.12754000492679],
                  [178.46513179394643, -18.131174963881026],
                  [178.46124275559325, -18.13356338710728],
                  [178.45941355836624, -18.137012147734026],
                  [178.45545467452254, -18.137588152743895],
                  [178.45381322030445, -18.134860648256403],
                  [178.45182846901602, -18.135110316211495],
                  [178.4511417128655, -18.138597205118366],
                  [178.45281535501883, -18.140177523725683],
                  [178.4524290754346, -18.142614268504317],
                  [178.44942467065434, -18.149465695016602],
                  [178.44079724302617, -18.161781914005577],
                  [178.44490690466196, -18.15882561957116]]]),
            {
              "system:index": "0"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51572819827194, -18.191706673967452],
                  [178.51782030782888, -18.19279728765488],
                  [178.520599089843, -18.19319478419986],
                  [178.52066346285935, -18.18557063324616],
                  [178.5165650474846, -18.18589680695692]]]),
            {
              "system:index": "1"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51352878687973, -18.175560880966515],
                  [178.5195583927452, -18.175010429866333],
                  [178.51947256205673, -18.171850399219384],
                  [178.50749918101425, -18.16856798423277],
                  [178.50849690237865, -18.172798402905183]]]),
            {
              "system:index": "2"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-169.9668935964081, -19.153966779328503],
                  [-169.73137418722843, -19.161101588113382],
                  [-169.69292203879093, -18.93588258037662],
                  [-169.9559072682831, -18.93393411148554]]]),
            {
              "system:index": "3"
            })]),
    Water = 
    /* color: #ffc82d */
    /* displayProperties: [
      {
        "type": "rectangle"
      },
      {
        "type": "rectangle"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      }
    ] */
    ee.FeatureCollection(
        [ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45946055771734, -18.14024940244376],
                  [178.45946055771734, -18.1419214701311],
                  [178.46014720322515, -18.1419214701311],
                  [178.46014720322515, -18.14024940244376]]], null, false),
            {
              "landcover": 0,
              "system:index": "0"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51167331687887, -18.171721986737307],
                  [178.51167331687887, -18.174209243289525],
                  [178.51313243858297, -18.174209243289525],
                  [178.51313243858297, -18.171721986737307]]], null, false),
            {
              "landcover": 0,
              "system:index": "1"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[176.98223148357476, -17.66817224634502],
                  [176.99064289104547, -17.67896715613972],
                  [176.98291812908258, -17.657867409893967]]]),
            {
              "landcover": 0,
              "system:index": "2"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[176.99253116619195, -17.64118213816592],
                  [176.99768100750055, -17.65246940292988],
                  [177.0000842667779, -17.649197804788333]]]),
            {
              "landcover": 0,
              "system:index": "3"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.0356181718072, -17.709058632453957],
                  [177.0414546586236, -17.71985108506209],
                  [177.04093967449273, -17.713473805186815]]]),
            {
              "landcover": 0,
              "system:index": "4"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[176.98583637249078, -17.629321569109344],
                  [176.98849712383355, -17.63512925561634],
                  [176.98961292278375, -17.632184536545765]]]),
            {
              "landcover": 0,
              "system:index": "5"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[176.95167685935525, -17.61083697331754],
                  [176.95450927207497, -17.61689067127004],
                  [176.95614005515603, -17.61378204092966]]]),
            {
              "landcover": 0,
              "system:index": "6"
            })]),
    CoralReef = 
    /* color: #d63000 */
    /* displayProperties: [
      {
        "type": "rectangle"
      },
      {
        "type": "rectangle"
      },
      {
        "type": "rectangle"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "marker"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      }
    ] */
    ee.FeatureCollection(
        [ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51937575523817, -18.186962569312442],
                  [178.51937575523817, -18.18805320156292],
                  [178.52004094307387, -18.18805320156292],
                  [178.52004094307387, -18.186962569312442]]], null, false),
            {
              "landcover": 1,
              "system:index": "0"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51681692783797, -18.187431440087725],
                  [178.51681692783797, -18.18846600918673],
                  [178.5174713868376, -18.18846600918673],
                  [178.5174713868376, -18.187431440087725]]], null, false),
            {
              "landcover": 1,
              "system:index": "1"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51973517124617, -18.185826062093724],
                  [178.51973517124617, -18.186478407949856],
                  [178.52055592720473, -18.186478407949856],
                  [178.52055592720473, -18.185826062093724]]], null, false),
            {
              "landcover": 1,
              "system:index": "2"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.5159259727102, -18.191362576810135],
                  [178.51641949916893, -18.191994516830988],
                  [178.51709541584069, -18.1923614486952],
                  [178.51827558780724, -18.192784438636],
                  [178.51888713146263, -18.190552263270703]]]),
            {
              "landcover": 1,
              "system:index": "3"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51609763408715, -18.19174479887066],
                  [178.5158347776037, -18.191433924991703],
                  [178.51577040458733, -18.19174479887066]]]),
            {
              "landcover": 1,
              "system:index": "4"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51494649970618, -18.17230156574618],
                  [178.51544002616492, -18.17255640838955],
                  [178.5171995552787, -18.172138466259174],
                  [178.51747850501624, -18.171588004365358]]]),
            {
              "landcover": 1,
              "system:index": "5"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51966848399815, -18.19291427534797],
                  [178.52036585834202, -18.193056970411913],
                  [178.52053751971897, -18.193077355411507],
                  [178.52054824855503, -18.192109065298183],
                  [178.5204731467026, -18.191670784847336],
                  [178.51986160304722, -18.191691170009015]]]),
            {
              "landcover": 1,
              "system:index": "6"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.76692712786047, -20.683507480991466],
                  [-178.7693303871378, -20.684872549945474],
                  [-178.76804292681066, -20.68599671515808],
                  [-178.76563966753332, -20.688245020598128],
                  [-178.76048982622473, -20.69105535554726]]]),
            {
              "landcover": 1,
              "system:index": "7"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.86008634548926, -21.03602088526002],
                  [-178.86021509152198, -21.03830401693714],
                  [-178.85982885342383, -21.04022662699794],
                  [-178.8593567846372, -21.041428245683868],
                  [-178.8581122396543, -21.040907545443524]]]),
            {
              "landcover": 1,
              "system:index": "8"
            }),
        ee.Feature(
            ee.Geometry.Point([176.9214062576305, -17.52742537081536]),
            {
              "landcover": 1,
              "system:index": "9"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[176.91762970733754, -17.533318174059946],
                  [176.92792938995473, -17.534300289320043],
                  [176.92003296661488, -17.528080136185178]]]),
            {
              "landcover": 1,
              "system:index": "10"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[176.96260498809926, -17.616778956784746],
                  [176.97187470245473, -17.61743339885858],
                  [176.9646649246227, -17.603689616938496]]]),
            {
              "landcover": 1,
              "system:index": "11"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[176.99350403595082, -17.66716404525899],
                  [177.00655030059926, -17.66716404525899],
                  [176.98801087188832, -17.649498151126316]]]),
            {
              "landcover": 1,
              "system:index": "12"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[176.91350983429066, -17.51252882037471],
                  [176.91865967559926, -17.51383845623948],
                  [176.91453980255238, -17.506635342101397]]]),
            {
              "landcover": 1,
              "system:index": "13"
            })]),
    Urban = /* color: #98ff00 */ee.FeatureCollection(
        [ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45997184094549, -18.138171853668013],
                  [178.4615758019364, -18.13710639929882],
                  [178.45985382374883, -18.13663229550734]]]),
            {
              "landcover": 2,
              "system:index": "0"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.4491518097794, -18.150875240061172],
                  [178.44933419999242, -18.150714669308243],
                  [178.44906329688192, -18.150803875300298]]]),
            {
              "landcover": 2,
              "system:index": "1"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.44950049695134, -18.150204919908695],
                  [178.4495836454308, -18.15008257983159],
                  [178.44950317916036, -18.150077482326516]]]),
            {
              "landcover": 2,
              "system:index": "2"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45109104689718, -18.15058723209773],
                  [178.45116078433156, -18.15041646608998],
                  [178.45086842354894, -18.150408819846923]]]),
            {
              "landcover": 2,
              "system:index": "3"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45053046521306, -18.148619589774594],
                  [178.4514531451142, -18.14876232105092],
                  [178.45149606045842, -18.14838000487019],
                  [178.4508416014588, -18.148033370810207]]]),
            {
              "landcover": 2,
              "system:index": "4"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.4527325588143, -18.150077482326516],
                  [178.4526252704537, -18.15057193962624],
                  [178.45354795035482, -18.150750351710524],
                  [178.45356940802694, -18.15058723209773],
                  [178.45281838950277, -18.150439404817142]]]),
            {
              "landcover": 2,
              "system:index": "5"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45149606045842, -18.146468411423726],
                  [178.45198690470815, -18.146534680346296],
                  [178.451780374614, -18.147250892868243],
                  [178.45216661271215, -18.14727892959768],
                  [178.45226317223668, -18.14677936538131],
                  [178.45259844836355, -18.14683289018706],
                  [178.4527137833512, -18.14631293423752],
                  [178.45237046059728, -18.146261958080792],
                  [178.45119297083974, -18.146050406871538],
                  [178.4510749536431, -18.146567814798193],
                  [178.4514370518601, -18.14663663248584]]]),
            {
              "landcover": 2,
              "system:index": "6"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.4515282469666, -18.148795455080542],
                  [178.4526467281258, -18.1490095517359],
                  [178.45268427905202, -18.148889759592002],
                  [178.45152020033956, -18.148698601745533]]]),
            {
              "landcover": 2,
              "system:index": "7"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45301955517888, -18.148354517095072],
                  [178.4532555895722, -18.148430980409266],
                  [178.4533816533959, -18.148395297533465],
                  [178.4535291748917, -18.148362163427997],
                  [178.45304369506002, -18.148250017178245]]]),
            {
              "landcover": 2,
              "system:index": "8"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.72281551840155, -20.66357708016902],
                  [-178.72314274790136, -20.663024954271908],
                  [-178.7233465957865, -20.663306036797536],
                  [-178.7230756930027, -20.66361472494251]]]),
            {
              "landcover": 2,
              "system:index": "9"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.7207019376978, -20.662091355008233],
                  [-178.72051418306677, -20.66200100638794],
                  [-178.72060001375525, -20.661840386485878],
                  [-178.7208628702387, -20.66183034773635]]]),
            {
              "landcover": 2,
              "system:index": "10"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.42325731357283, -18.142875004120153],
                  [178.42406734069533, -18.1429769583938],
                  [178.42422827323622, -18.14128960751552],
                  [178.42325731357283, -18.14130490079902]]]),
            {
              "landcover": 2,
              "system:index": "11"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45481754684528, -18.148265397273175],
                  [178.45506431007465, -18.148461653145745],
                  [178.45515282297214, -18.148349506959836],
                  [178.45513136530002, -18.148308726510752],
                  [178.45493824625095, -18.14814560461928],
                  [178.4548416867264, -18.148216970465523]]]),
            {
              "landcover": 2,
              "system:index": "12"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45495165729602, -18.148117568028848],
                  [178.4551796450623, -18.148278141167687],
                  [178.455305708886, -18.148094628996976],
                  [178.4551796450623, -18.14802581188329],
                  [178.455128683091, -18.148089531433936],
                  [178.45506162786563, -18.14800287283937]]]),
            {
              "landcover": 2,
              "system:index": "13"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45518500948032, -18.14772250650382],
                  [178.45542104387363, -18.1478626897278],
                  [178.45547468805393, -18.147814262808595],
                  [178.4552279247865, -18.147663884443162]]]),
            {
              "landcover": 2,
              "system:index": "14"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45331819200595, -18.147865238512665],
                  [178.45345766687473, -18.14795954352598],
                  [178.45354617977222, -18.147895823927943],
                  [178.45377148532947, -18.148064043616454],
                  [178.4538492693909, -18.147974836226048],
                  [178.45347912454685, -18.147704664994517]]]),
            {
              "landcover": 2,
              "system:index": "15"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.4539753332146, -18.14785249458803],
                  [178.454045070649, -18.14776073830333],
                  [178.45392705345233, -18.147674079545766],
                  [178.45391096019824, -18.147625652574305],
                  [178.45392705345233, -18.1476001646892],
                  [178.45402629518588, -18.147656238031526],
                  [178.45411480808338, -18.147589969534106],
                  [178.4541013970383, -18.14752624980138],
                  [178.45414699459155, -18.147485469160234],
                  [178.45402361297687, -18.14738861509937],
                  [178.45392973566135, -18.147437042136506],
                  [178.4540102019318, -18.147493115531155],
                  [178.45395119333347, -18.14755428648663],
                  [178.45387340927203, -18.147513505852036],
                  [178.45384122276386, -18.147765835875962]]]),
            {
              "landcover": 2,
              "system:index": "16"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.74548775078193, -20.664793965046556],
                  [-178.74527585626976, -20.66472369515217],
                  [-178.7450854194297, -20.664522923846192],
                  [-178.74518466116325, -20.664427557382893],
                  [-178.74563795448677, -20.664686050552504]]]),
            {
              "landcover": 2,
              "system:index": "17"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.74609124781028, -20.66398083999468],
                  [-178.7460241925849, -20.66409126463984],
                  [-178.74571305633918, -20.663915589030264],
                  [-178.74578279377357, -20.66379763529272]]]),
            {
              "landcover": 2,
              "system:index": "18"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.74501836420433, -20.66504743832413],
                  [-178.74454093099968, -20.664766359020504],
                  [-178.7446374905242, -20.664658444506827],
                  [-178.7450693261756, -20.664919446920003]]]),
            {
              "landcover": 2,
              "system:index": "19"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.74547165752784, -20.66400091720886],
                  [-178.74580156923668, -20.664246862867216],
                  [-178.7456969630851, -20.66430709439634],
                  [-178.74539655567543, -20.664113851489184]]]),
            {
              "landcover": 2,
              "system:index": "20"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.72111280078627, -20.663850870441898],
                  [-178.7211181652043, -20.66437789688565],
                  [-178.72069974059798, -20.664362839012625],
                  [-178.7206836473439, -20.663906082821747],
                  [-178.7208123933766, -20.663760522867953]]]),
            {
              "landcover": 2,
              "system:index": "21"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.72215349788405, -20.663183300988475],
                  [-178.72235198135115, -20.663509557102877],
                  [-178.72213204021193, -20.6640616812388],
                  [-178.7216385137532, -20.664036584730695],
                  [-178.7220354806874, -20.663524615060478]]]),
            {
              "landcover": 2,
              "system:index": "22"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45369693874605, -18.14455082656668],
                  [178.45424947380312, -18.14332738237195],
                  [178.45334288715608, -18.14325601452967],
                  [178.45296201347597, -18.14425006399574]]]),
            {
              "landcover": 2,
              "system:index": "23"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.4534257021365, -18.147127135793937],
                  [178.45354371933317, -18.146953817607653],
                  [178.4535329904971, -18.146749913639088],
                  [178.45367246536588, -18.14667344958955],
                  [178.4536778297839, -18.14640837395892],
                  [178.45382266907072, -18.146362495443583],
                  [178.4539460506854, -18.145832342837355],
                  [178.45391922859525, -18.145414336764304],
                  [178.45366710094785, -18.145378653272736],
                  [178.45353567270612, -18.14599801569949],
                  [178.4537851181445, -18.146097419341434],
                  [178.45377170709943, -18.14623760386861],
                  [178.45361345676756, -18.146217213398916],
                  [178.4535517659602, -18.146337007374292],
                  [178.45345520643568, -18.14638798350914],
                  [178.4533130493579, -18.1471143918155]]]),
            {
              "landcover": 2,
              "system:index": "24"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45435805375038, -18.143333184749366],
                  [178.454470706529, -18.14340455256015],
                  [178.45462091023384, -18.14322613297857],
                  [178.45491595322548, -18.142787729804844],
                  [178.45471210534035, -18.142772436651036],
                  [178.45442242676674, -18.14266538453692],
                  [178.4541273837751, -18.142425791472622],
                  [178.45386452729164, -18.142604211870864],
                  [178.45405764634071, -18.142787729804844],
                  [178.45398790890633, -18.142833609258226],
                  [178.45372505242287, -18.14308849489123]]]),
            {
              "landcover": 2,
              "system:index": "25"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.457528424806, -18.14385314956091],
                  [178.45763571316658, -18.14353199500695],
                  [178.45568306500374, -18.142160709401903],
                  [178.45527536923348, -18.142466573293987],
                  [178.45521099621712, -18.142736752620273],
                  [178.45498569065987, -18.14303751779437],
                  [178.45481402928291, -18.14321083986311],
                  [178.45483548695503, -18.14346572494621],
                  [178.45461018139778, -18.143481018039385],
                  [178.45444388443886, -18.143690023511983],
                  [178.45440096909462, -18.143899028734747],
                  [178.45451362187325, -18.143995884728877],
                  [178.45467991883217, -18.14403156850265],
                  [178.454830122537, -18.144138619780307],
                  [178.45456726605354, -18.144439382542956],
                  [178.45487303788124, -18.144729949466345],
                  [178.45518417412697, -18.144893074545013],
                  [178.45559723431526, -18.144592312562835],
                  [178.45603175217568, -18.144663679859665],
                  [178.45606930310188, -18.14484209797428],
                  [178.45678277069985, -18.14490326985738],
                  [178.45700271183907, -18.144653484533325],
                  [178.45678277069985, -18.144235475641718],
                  [178.4570456271833, -18.14373590272863],
                  [178.4567720418638, -18.143603362736588],
                  [178.4564340835279, -18.143695121203315],
                  [178.45644481236397, -18.143929614843945],
                  [178.4562141423887, -18.14424567099242],
                  [178.45629997307717, -18.144368015154566],
                  [178.45617122704445, -18.14457701956686],
                  [178.45584399754463, -18.144541335904407],
                  [178.4559298282331, -18.143791977310403],
                  [178.4561444049543, -18.1435370927029],
                  [178.45629460865914, -18.143307696238445],
                  [178.45633752400337, -18.143215937568417],
                  [178.4566111093229, -18.143287305427048],
                  [178.45734067017494, -18.143766488866365],
                  [178.45743186528145, -18.143827661125805]]]),
            {
              "landcover": 2,
              "system:index": "26"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[179.19917491668403, -8.51468969254331],
                  [179.19895765775382, -8.515477527455577],
                  [179.19941631549537, -8.51561546471445],
                  [179.19953701490104, -8.515374074478789],
                  [179.19951555722892, -8.51466847137846]]]),
            {
              "landcover": 2,
              "system:index": "27"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[179.19997421497047, -8.513140544417025],
                  [179.2000144481057, -8.51327582986285],
                  [179.1998374223107, -8.513328882965805],
                  [179.19988570207298, -8.513424378532603],
                  [179.20007882112205, -8.513387241370578],
                  [179.20011100763023, -8.513286440484034],
                  [179.20016733401954, -8.513323577655843],
                  [179.20029608005225, -8.513294398449714],
                  [179.20026389354408, -8.513175028946968],
                  [179.2001619696015, -8.51316176566658],
                  [179.20010027879417, -8.513137891760746],
                  [179.20019952052772, -8.512994648294406],
                  [179.20007613891303, -8.512811614898292],
                  [179.19988570207298, -8.512782435653133],
                  [179.19984010451972, -8.51291241590987],
                  [179.1997274517411, -8.512944247802757],
                  [179.1998025535935, -8.513188292226886]]]),
            {
              "landcover": 2,
              "system:index": "28"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[179.1992768406266, -8.513320925000835],
                  [179.19937608236015, -8.513278482518174],
                  [179.19930634492576, -8.513071575347812],
                  [179.19912931913078, -8.513124628479073],
                  [179.1990595816964, -8.513055659406984],
                  [179.1989710687989, -8.51309544925778],
                  [179.19905153506934, -8.513164418322706],
                  [179.19912395471275, -8.513243997997543]]]),
            {
              "landcover": 2,
              "system:index": "29"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[179.1966826669576, -8.505690370090754],
                  [179.19688651484273, -8.50565057947053],
                  [179.19685432833455, -8.505528554876065],
                  [179.19691333693288, -8.505475500692448],
                  [179.1968650571706, -8.505385308563422],
                  [179.19676849764608, -8.505377350433369],
                  [179.196755086601, -8.505276547438422],
                  [179.19663170498632, -8.505311032676524],
                  [179.19663170498632, -8.50555242925631]]]),
            {
              "landcover": 2,
              "system:index": "30"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[179.19619986933492, -8.504236683416279],
                  [179.19647613686345, -8.504186281772062],
                  [179.19631788653157, -8.503714097625393],
                  [179.19612744969152, -8.50358411424694],
                  [179.19603893679403, -8.503432909036931],
                  [179.19611672085546, -8.503361285495583],
                  [179.19604966563008, -8.503154372967575],
                  [179.19587532204412, -8.503162331143809],
                  [179.1957653514745, -8.503188858396662],
                  [179.1958377711179, -8.503337410978816],
                  [179.19592091959737, -8.503427603589886],
                  [179.19598797482274, -8.5035071852879],
                  [179.19608185213826, -8.503769804774095],
                  [179.196033572376, -8.503873260885904],
                  [179.19611403864644, -8.503984675128933]]]),
            {
              "landcover": 2,
              "system:index": "31"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[176.31807998181355, -6.287028909580621],
                  [176.31816849471105, -6.287028909580621],
                  [176.3181765413381, -6.286847616272231],
                  [176.31829187632573, -6.286860946664701],
                  [176.31829455853475, -6.286690317615222],
                  [176.31835356713307, -6.2866769872183745],
                  [176.318340156088, -6.286602336989707],
                  [176.3181068039037, -6.286543683231093],
                  [176.3178922271825, -6.2864690329833035],
                  [176.3174952602483, -6.286421043532632],
                  [176.31747916699422, -6.2865196885097525],
                  [176.31756231547368, -6.286527686750327],
                  [176.31756231547368, -6.2864770312246545],
                  [176.31765082837117, -6.286490361626616],
                  [176.31764814616216, -6.286620999547878],
                  [176.3175864553548, -6.286863612743163],
                  [176.31762668849004, -6.286999582726198],
                  [176.31770862388748, -6.287002503309057],
                  [176.31768716621536, -6.28691185665817],
                  [176.31770594167847, -6.286783884888827],
                  [176.3177649502768, -6.286570598536617],
                  [176.31811631965775, -6.286634584451464]]]),
            {
              "landcover": 2,
              "system:index": "32"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.31413460331586, -8.027349665817098],
                  [178.3142848070207, -8.027222181279445],
                  [178.3140702302995, -8.02693534092348],
                  [178.31420434075025, -8.026797232531704],
                  [178.3139951284471, -8.02652632747304],
                  [178.31374836521772, -8.026690995275327],
                  [178.31368399220136, -8.026802544393801],
                  [178.3137108142915, -8.02693534092348],
                  [178.31381273823408, -8.027052201833682],
                  [178.31375372963575, -8.027131879707737],
                  [178.31382883148817, -8.027243428705157],
                  [178.31393075543073, -8.027163750852983]]]),
            {
              "landcover": 2,
              "system:index": "33"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.31461203652051, -8.026568822396152],
                  [178.3148427064958, -8.026802544393801],
                  [178.31495535927442, -8.026717554592027],
                  [178.31478369789747, -8.026366971471793]]]),
            {
              "landcover": 2,
              "system:index": "34"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.314338451201, -8.026382907074733],
                  [178.31447256165174, -8.026207615408069],
                  [178.3144296463075, -8.026149184835708],
                  [178.31427944260267, -8.026266045972033]]]),
            {
              "landcover": 2,
              "system:index": "35"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.42383141803708, -18.13585071615805],
                  [178.42437858867612, -18.136391094834686],
                  [178.42487211513486, -18.13697225494454],
                  [178.4258698968884, -18.137849089350595],
                  [178.42680330562558, -18.138991006638197],
                  [178.4283375291821, -18.138929832686806],
                  [178.42755432414975, -18.137849089350595],
                  [178.42695350933042, -18.136737751974692],
                  [178.4265887289044, -18.135065634728953],
                  [178.426159575462, -18.13404604319823]]]),
            {
              "landcover": 2,
              "system:index": "36"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.42359538364377, -18.134759757893693],
                  [178.42450733470884, -18.134800541502628],
                  [178.4255802183148, -18.132282135817796],
                  [178.4247540979382, -18.131976254115934]]]),
            {
              "landcover": 2,
              "system:index": "37"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.42301602649655, -18.14481262975835],
                  [178.4241961984631, -18.14475145784354],
                  [178.42377777385678, -18.14349742887327],
                  [178.4226297883984, -18.143405670302784]]]),
            {
              "landcover": 2,
              "system:index": "38"
            })]),
    Vegetation = 
    /* color: #0b4a8b */
    /* displayProperties: [
      {
        "type": "rectangle"
      },
      {
        "type": "rectangle"
      },
      {
        "type": "rectangle"
      },
      {
        "type": "rectangle"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      },
      {
        "type": "polygon"
      }
    ] */
    ee.FeatureCollection(
        [ee.Feature(
            ee.Geometry.Polygon(
                [[[178.45473514540254, -18.140313351523663],
                  [178.45473514540254, -18.14156740332945],
                  [178.45522867186128, -18.14156740332945],
                  [178.45522867186128, -18.140313351523663]]], null, false),
            {
              "landcover": 3,
              "system:index": "0"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.4616981600052, -18.13873302953109],
                  [178.4616981600052, -18.13924281237926],
                  [178.46222387297212, -18.13924281237926],
                  [178.46222387297212, -18.13873302953109]]], null, false),
            {
              "landcover": 3,
              "system:index": "1"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.44566262127813, -18.152658287945886],
                  [178.44566262127813, -18.15302020509191],
                  [178.44584501149114, -18.15302020509191],
                  [178.44584501149114, -18.152658287945886]]], null, false),
            {
              "landcover": 3,
              "system:index": "2"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.5179756421324, -18.186488600834497],
                  [178.5179756421324, -18.18756904322816],
                  [178.51859255020582, -18.18756904322816],
                  [178.51859255020582, -18.186488600834497]]], null, false),
            {
              "landcover": 3,
              "system:index": "3"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.7344639564899, -20.63740791490803],
                  [-178.73405894292864, -20.638090658032304],
                  [-178.73326500906023, -20.63802037579335],
                  [-178.73338839067492, -20.63718702677061],
                  [-178.7334581281093, -20.636986219094666],
                  [-178.73420914663348, -20.637202087335634]]]),
            {
              "landcover": 3,
              "system:index": "4"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.73479386819872, -20.638557532081162],
                  [-178.73407235397372, -20.6383215851603],
                  [-178.73484751237902, -20.637897381797032],
                  [-178.73505672468218, -20.638176000707396]]]),
            {
              "landcover": 3,
              "system:index": "5"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.74302761868864, -20.64043253311953],
                  [-178.74099986867338, -20.64115040307717],
                  [-178.73927789048582, -20.639052004447784],
                  [-178.74009328202635, -20.638324084285212],
                  [-178.7434138567868, -20.63992048394668]]]),
            {
              "landcover": 3,
              "system:index": "6"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.7364079268399, -20.64168252827181],
                  [-178.7327011139813, -20.640307031030744],
                  [-178.73516338185698, -20.63825882230738],
                  [-178.73787777738005, -20.638580111773972]]]),
            {
              "landcover": 3,
              "system:index": "7"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.74354656014327, -20.64315982143467],
                  [-178.74451215538863, -20.647276158724186],
                  [-178.74240930352096, -20.650267963214368],
                  [-178.73556430611495, -20.64542884064943],
                  [-178.7352853563774, -20.64332046107668],
                  [-178.73970563683395, -20.642477101063182]]]),
            {
              "landcover": 3,
              "system:index": "8"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.48506573914364, -18.16518913720325],
                  [177.48699692963436, -18.16747260704814],
                  [177.48846678862807, -18.16863473497781],
                  [177.49259738205745, -18.16602505373636],
                  [177.4862459111102, -18.16243670082761]]]),
            {
              "landcover": 3,
              "system:index": "9"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.48013802929762, -18.1673739821305],
                  [177.48028823300245, -18.167944845248158],
                  [177.48174735470656, -18.16842396106695],
                  [177.48153277798536, -18.16903559637184],
                  [177.48200484677199, -18.169524903072947],
                  [177.48538443013075, -18.166752146967767],
                  [177.4841452287637, -18.164315755347644],
                  [177.4814469472969, -18.16306186650306]]]),
            {
              "landcover": 3,
              "system:index": "10"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.50933967101338, -18.17305120584843],
                  [177.5100263165212, -18.17304101217873],
                  [177.50995121466877, -18.172821848136316],
                  [177.50920556056263, -18.172877913382663]]]),
            {
              "landcover": 3,
              "system:index": "11"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.49665282237294, -18.17022243966974],
                  [177.49820313918354, -18.170941103529074],
                  [177.49833188521626, -18.170808584174644],
                  [177.49775789248707, -18.17036005638247],
                  [177.49663136470082, -18.169794298091954]]]),
            {
              "landcover": 3,
              "system:index": "12"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.5003757284856, -18.171384532945005],
                  [177.5009336279607, -18.17151195498262],
                  [177.50122867095234, -18.171348854757813],
                  [177.5003757284856, -18.171206141936093]]]),
            {
              "landcover": 3,
              "system:index": "13"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.50671110617878, -18.172480359428747],
                  [177.50777326094868, -18.17265365228922],
                  [177.50675402152302, -18.17228667896934]]]),
            {
              "landcover": 3,
              "system:index": "14"
            })]),
    Suva = ee.Image("users/kishan2196/Suva/NDSV/2020_S2_NDSV_Suva_clipped"),
    table = ee.FeatureCollection("users/kishan2196/Suvapoints"),
    Sand = /* color: #d63000 */ee.FeatureCollection(
        [ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51897107516305, -18.189650408820107],
                  [178.51924466048257, -18.189647860645188],
                  [178.51920174513833, -18.189459295598997],
                  [178.5188530579664, -18.189500066437077]]]),
            {
              "landcover": 4,
              "system:index": "0"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51880183739956, -18.187190469700646],
                  [178.51887693925198, -18.186948389503723],
                  [178.51875087542828, -18.18694584128933]]]),
            {
              "landcover": 4,
              "system:index": "1"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51447603146653, -18.173094338347838],
                  [178.5147227946959, -18.17288536805175],
                  [178.51433655659775, -18.172864980692385]]]),
            {
              "landcover": 4,
              "system:index": "2"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.8447570243594, -21.04205941478167],
                  [-178.84393895060987, -21.04143107221583],
                  [-178.84419912488437, -21.041205769454937],
                  [-178.8449903765437, -21.04186164987749]]]),
            {
              "landcover": 4,
              "system:index": "3"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.84776646287412, -21.04108560784309],
                  [-178.8463958540675, -21.04185413981267],
                  [-178.8463583031413, -21.04169142164841],
                  [-178.84740704686612, -21.04111815162253]]]),
            {
              "landcover": 4,
              "system:index": "4"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.85175157316473, -21.036074082407534],
                  [-178.8516013694599, -21.036915239953938],
                  [-178.85129023321417, -21.036835129916046],
                  [-178.85152090318945, -21.035943902839076]]]),
            {
              "landcover": 4,
              "system:index": "5"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.844946809328, -21.03604404005809],
                  [-178.84472954996374, -21.036719971661114],
                  [-178.84441841371802, -21.03669493724608],
                  [-178.84458471067694, -21.035588411899578],
                  [-178.84497094877509, -21.035648494934836]]]),
            {
              "landcover": 4,
              "system:index": "6"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.84461335053754, -21.035281506114547],
                  [-178.84525708070112, -21.033363839227878],
                  [-178.84597054829908, -21.03229233836566],
                  [-178.84630314221693, -21.03237245084616],
                  [-178.84599737038923, -21.032707921390465],
                  [-178.84556821694684, -21.03321362929101],
                  [-178.84541264882398, -21.033839503028965],
                  [-178.84491912236524, -21.03552684552803]]]),
            {
              "landcover": 4,
              "system:index": "7"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.85095409135516, -21.033929627256033],
                  [-178.85054103246046, -21.03327872028224],
                  [-178.85079316010786, -21.033198608288984],
                  [-178.85134569516492, -21.034595554998273],
                  [-178.85109893193555, -21.03466064538581]]]),
            {
              "landcover": 4,
              "system:index": "8"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.85061076989484, -21.03302837016031],
                  [-178.8503854643376, -21.03312851025956],
                  [-178.84989730229688, -21.03240750004283],
                  [-178.85021916737867, -21.032327387581176]]]),
            {
              "landcover": 4,
              "system:index": "9"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.74393250451897, -20.637514730754617],
                  [-178.74243046747063, -20.63669141937225],
                  [-178.7415077875695, -20.636570934406],
                  [-178.7409284304223, -20.63604883178292],
                  [-178.74041344629143, -20.635767698859137],
                  [-178.74097134576652, -20.635386160488817],
                  [-178.74245192514275, -20.636309883318415],
                  [-178.74367501245354, -20.63695246980516]]]),
            {
              "landcover": 4,
              "system:index": "10"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.74493331787892, -20.6388992260919],
                  [-178.74490113137074, -20.638788783145987],
                  [-178.74477238533802, -20.638587977584343],
                  [-178.74455780861683, -20.638447413533463],
                  [-178.74442906258412, -20.638256647828115],
                  [-178.74447197792836, -20.638015680279402],
                  [-178.74468655464955, -20.638367091160365],
                  [-178.74495477555104, -20.63837713145932],
                  [-178.7450727927477, -20.63891930661892],
                  [-178.74487967369862, -20.639079950739532]]]),
            {
              "landcover": 4,
              "system:index": "11"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.73950452683277, -20.63554574089891],
                  [-178.73921484825917, -20.635445336010914],
                  [-178.73878569481678, -20.635375052549882],
                  [-178.73843164322682, -20.63525456654115],
                  [-178.7385925757677, -20.6350637968322],
                  [-178.7390753733904, -20.63529472855465],
                  [-178.73974056122609, -20.635475457484276],
                  [-178.73969764588185, -20.63563610524146]]]),
            {
              "landcover": 4,
              "system:index": "12"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.73042793152638, -20.636790756005258],
                  [-178.7304493891985, -20.636369058481222],
                  [-178.73052449105091, -20.636248573259678],
                  [-178.73059959290333, -20.636118047495305],
                  [-178.73064250824757, -20.63604776434504],
                  [-178.73078198311634, -20.636017642984996],
                  [-178.73077125428028, -20.636168249725632],
                  [-178.73069615242787, -20.636218451939364],
                  [-178.73067469475575, -20.636389139342196],
                  [-178.73056740639515, -20.636630109467255]]]),
            {
              "landcover": 4,
              "system:index": "13"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.73215527413197, -20.63503367527735],
                  [-178.73239130852528, -20.63502363475774],
                  [-178.73287410614796, -20.634953351101892],
                  [-178.7328097331316, -20.63503367527735],
                  [-178.73269171593495, -20.635113999410382]]]),
            {
              "landcover": 4,
              "system:index": "14"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.73245568154164, -20.634973432149746],
                  [-178.7330457675249, -20.63488306741362],
                  [-178.7329170214922, -20.635043715796293],
                  [-178.73244495270558, -20.635134080437002]]]),
            {
              "landcover": 4,
              "system:index": "15"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.73694033501454, -20.634963391626158],
                  [-178.73647899506398, -20.634963391626158],
                  [-178.73665065644093, -20.63473245940079],
                  [-178.73747677681752, -20.634923229525175]]]),
            {
              "landcover": 4,
              "system:index": "16"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.5189184843037, -18.186984450245347],
                  [178.5190445481274, -18.186609862408446],
                  [178.5191867052052, -18.186423842299316],
                  [178.51929131135677, -18.186191953940234],
                  [178.5192242561314, -18.186056897940677],
                  [178.51905795917247, -18.185957517044002],
                  [178.51863417014812, -18.1859014560004],
                  [178.5184142290089, -18.18594222767026],
                  [178.5180548130009, -18.186067090849946],
                  [178.5178751049969, -18.186158827006622],
                  [178.51760688409541, -18.1863193651647],
                  [178.51748350248073, -18.186416197633044],
                  [178.5174405871365, -18.186688857189303],
                  [178.51749154910777, -18.186948775247945],
                  [178.51759615525935, -18.186971709175676],
                  [178.51756128654216, -18.186846846643654],
                  [178.51759883746837, -18.18670669471551],
                  [178.51764443502162, -18.186551253354374],
                  [178.5178375540707, -18.18633975095198],
                  [178.51808163509105, -18.186146085876633],
                  [178.5182345210049, -18.18609002489364],
                  [178.51863953456615, -18.186005933385353],
                  [178.51900431499217, -18.18602122275351],
                  [178.5190820990536, -18.18611295893432],
                  [178.51900163278316, -18.186372877851216],
                  [178.51878705606197, -18.186885069877242]]]),
            {
              "landcover": 4,
              "system:index": "17"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.51873126771238, -18.187380523988022],
                  [178.5186239793518, -18.18812459915237],
                  [178.51887074258116, -18.188889056631197],
                  [178.5190424039581, -18.189194838684376],
                  [178.5188600137451, -18.189531198323554],
                  [178.5180875375488, -18.188797321910684],
                  [178.51799097802427, -18.18846096085593],
                  [178.5179587915161, -18.18801247844036],
                  [178.51791587617186, -18.187676115871763],
                  [178.51784077431944, -18.1872480171189],
                  [178.51773348595884, -18.1870747387535],
                  [178.5175618245819, -18.187125702996493],
                  [178.51771202828672, -18.187400909651274],
                  [178.51776567246702, -18.187584380513194],
                  [178.51776567246702, -18.187961514456518],
                  [178.5177978589752, -18.188328454807735],
                  [178.51787296082762, -18.188573081279603],
                  [178.51794806268003, -18.18885847839639],
                  [178.51805535104063, -18.189164260503194],
                  [178.51818409707334, -18.189357922226833],
                  [178.51832357194212, -18.189551583735387],
                  [178.518484504483, -18.189622932658015],
                  [178.51877418305662, -18.189775823108178],
                  [178.51916042115477, -18.18982678656177],
                  [178.51934281136778, -18.189633125358853],
                  [178.51935354020384, -18.189255995030674],
                  [178.51919260766294, -18.18883809290354],
                  [178.51897803094175, -18.18824691257419],
                  [178.51887074258116, -18.187890164854107],
                  [178.5189673021057, -18.18742129531213]]]),
            {
              "landcover": 4,
              "system:index": "18"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.8497657291341, -21.03226592981771],
                  [-178.84893960875752, -21.031785253945014],
                  [-178.84842462462666, -21.031705141148805],
                  [-178.84755558890583, -21.031715155250687],
                  [-178.84736246985676, -21.031454788383083],
                  [-178.84884304923298, -21.03152488719989],
                  [-178.8498730174947, -21.0319955498301]]]),
            {
              "landcover": 4,
              "system:index": "19"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-179.15839595368624, -16.087563466865785],
                  [-179.1582269745183, -16.08762016463949],
                  [-179.15814382603884, -16.08756862120952],
                  [-179.15839595368624, -16.08744749409647]]]),
            {
              "landcover": 4,
              "system:index": "20"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-179.15425902655187, -16.091238528111692],
                  [-179.154280484224, -16.0911148261106],
                  [-179.15434485724035, -16.091132865990556],
                  [-179.15430998852315, -16.0912462594842]]]),
            {
              "landcover": 4,
              "system:index": "21"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-179.16095405142983, -16.085583390034966],
                  [-179.15952175181587, -16.08631531280072],
                  [-179.15926694195946, -16.086116869218678],
                  [-179.16087358515938, -16.085506074092415]]]),
            {
              "landcover": 4,
              "system:index": "22"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-179.16213154118736, -16.085106607909903],
                  [-179.16335462849815, -16.084307673135644],
                  [-179.16344850581368, -16.08441076167412],
                  [-179.16328220885475, -16.08452158179333],
                  [-179.16222273629387, -16.085194232818697]]]),
            {
              "landcover": 4,
              "system:index": "23"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-179.15411958652848, -16.090114923025627],
                  [-179.1551978343566, -16.088877892775148],
                  [-179.15554652148765, -16.088908818622485],
                  [-179.15511066411673, -16.089184573835773],
                  [-179.1547606361394, -16.089563414352874],
                  [-179.1545634932972, -16.089924214419902],
                  [-179.15431270555413, -16.090202545724512]]]),
            {
              "landcover": 4,
              "system:index": "24"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.31518066483167, -8.026239486625817],
                  [178.31587803917554, -8.027514333285696],
                  [178.31595314102796, -8.027163750852983],
                  [178.3151860292497, -8.02568705183101]]]),
            {
              "landcover": 4,
              "system:index": "25"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.31157497925832, -8.02918341540704],
                  [178.31268004937246, -8.029895200127143],
                  [178.31231526894643, -8.0295233724439],
                  [178.31178955597952, -8.028822210743611]]]),
            {
              "landcover": 4,
              "system:index": "26"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[178.31225894255712, -8.026700126859847],
                  [178.31344984335973, -8.025815700642179],
                  [178.3133881525524, -8.025653688423716],
                  [178.31247083831946, -8.026210107738883],
                  [178.31159643693044, -8.02678777259639]]]),
            {
              "landcover": 4,
              "system:index": "27"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.4956281453696, -18.167466271318638],
                  [177.49493077146425, -18.167810317883394],
                  [177.4954672125625, -18.168411760804688],
                  [177.4983210691403, -18.170481128268662],
                  [177.49947978773545, -18.170975527453738],
                  [177.50106766498558, -18.171143715433796],
                  [177.50278427821758, -18.17138327635243],
                  [177.50441506183617, -18.171255846986497],
                  [177.50064383871498, -18.169038719283765],
                  [177.49852222362753, -18.168358268788623],
                  [177.49680831706993, -18.167657407591605]]]),
            {
              "landcover": 4,
              "system:index": "28"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.51083640335872, -18.17350734787359],
                  [177.51590041397884, -18.173405411416265],
                  [177.51120118378475, -18.17210061950442]]]),
            {
              "landcover": 4,
              "system:index": "29"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.5053217794337, -18.171754033602728],
                  [177.50703839539364, -18.17218216928467],
                  [177.50778941568214, -18.172202557366695],
                  [177.50841168640926, -18.17128511960656],
                  [177.50710276841, -18.170755042629928],
                  [177.50540761231258, -18.17124434451165]]]),
            {
              "landcover": 4,
              "system:index": "30"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.49346641777828, -18.169643914513465],
                  [177.4941825668169, -18.170375323704217],
                  [177.49497918366268, -18.17045432516798],
                  [177.4945466771749, -18.16990640557947],
                  [177.4947900878199, -18.169603138332118],
                  [177.49452052587014, -18.16938397003388],
                  [177.49399749516323, -18.168685688256826],
                  [177.4924283972677, -18.168135220567034],
                  [177.49089149712398, -18.167421637064933]]]),
            {
              "landcover": 4,
              "system:index": "31"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.4864958474188, -18.1678262175309],
                  [177.4869625517874, -18.16811674556257],
                  [177.4873219677954, -18.168524503388088],
                  [177.48728978128722, -18.16865702447634],
                  [177.4869625517874, -18.16863153965948],
                  [177.48737024755766, -18.168871096791012],
                  [177.48803543539336, -18.169029102378715],
                  [177.48830365629485, -18.168998520663223],
                  [177.48790668936064, -18.168794642423045],
                  [177.48688744993498, -18.167907769307895],
                  [177.48600232096007, -18.1674694280587]]]),
            {
              "landcover": 4,
              "system:index": "32"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.48497771711638, -18.169533699910904],
                  [177.48533176870635, -18.16966112329949],
                  [177.48581456632903, -18.169523506035805],
                  [177.48469876737883, -18.16905458713752]]]),
            {
              "landcover": 4,
              "system:index": "33"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.48886692018797, -18.169166720032134],
                  [177.4891083189993, -18.169447051953536],
                  [177.48976277799895, -18.169447051953536],
                  [177.4897359559088, -18.169105556644002],
                  [177.4892853447943, -18.168769157626272]]]),
            {
              "landcover": 4,
              "system:index": "34"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.49078738184264, -18.169931260575705],
                  [177.49345349760344, -18.170354304904233],
                  [177.49351250620177, -18.16998222982623],
                  [177.49206411333373, -18.169865000527786],
                  [177.49016510935118, -18.16965092943182]]]),
            {
              "landcover": 4,
              "system:index": "35"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.49173151941588, -18.167795635604715],
                  [177.49383437128355, -18.168437854929845],
                  [177.49371098966887, -18.16827475183296],
                  [177.49343203993132, -18.168050484825972],
                  [177.49308871717741, -18.16797912708152],
                  [177.49184953661253, -18.167596852953878]]]),
            {
              "landcover": 4,
              "system:index": "36"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.48636441911347, -18.169368049079175],
                  [177.4865655848532, -18.16959996008457],
                  [177.4869250008612, -18.169737577287993],
                  [177.48912977667143, -18.169533699910904],
                  [177.48885082693388, -18.169579572341505],
                  [177.48807298631957, -18.16948273052942],
                  [177.48565363378813, -18.16861624876758]]]),
            {
              "landcover": 4,
              "system:index": "37"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.4839745709448, -18.17009945904597],
                  [177.48583065958312, -18.169992423674547],
                  [177.48587893934538, -18.16975796501261],
                  [177.484907979682, -18.16975796501261],
                  [177.48344349355986, -18.16958976621334]]]),
            {
              "landcover": 4,
              "system:index": "38"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.51437210115424, -18.17255053354078],
                  [177.51518749293535, -18.17293789486761],
                  [177.51675390300005, -18.174324229299813],
                  [177.5182773977205, -18.174161131702856],
                  [177.51505874690264, -18.171633099460955]]]),
            {
              "landcover": 4,
              "system:index": "39"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.51000010046374, -18.170618817668018],
                  [177.51001619395464, -18.171210058231456],
                  [177.51193129119127, -18.17151077439152],
                  [177.51102470454424, -18.169482204492763],
                  [177.51005910929888, -18.17014480529325]]]),
            {
              "landcover": 4,
              "system:index": "40"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.4903786130507, -18.16841489590782],
                  [177.49095260577988, -18.16870032599823],
                  [177.49133884387803, -18.168572901908384],
                  [177.49031424003434, -18.168236501864328]]]),
            {
              "landcover": 4,
              "system:index": "41"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.50814727317007, -18.171702083824513],
                  [177.5096117592922, -18.172507388073527],
                  [177.51092067729147, -18.172772424102373],
                  [177.51060954104574, -18.171987508540347],
                  [177.50904849539907, -18.171049682721495]]]),
            {
              "landcover": 4,
              "system:index": "42"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.50416687499197, -18.171003810677195],
                  [177.50523975859792, -18.171279042762198],
                  [177.50522902976186, -18.17081522547953],
                  [177.5041025019756, -18.170331019309184]]]),
            {
              "landcover": 4,
              "system:index": "43"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[177.49879864929989, -18.171403225011716],
                  [177.50088004349544, -18.172076012247278],
                  [177.50545052765682, -18.172667247374193],
                  [177.506201546181, -18.17230027408289],
                  [177.5047424244769, -18.171749812699336],
                  [177.50203875778988, -18.171484775117772],
                  [177.50047234772518, -18.171586712696758],
                  [177.49888447998836, -18.171117799340504]]]),
            {
              "landcover": 4,
              "system:index": "44"
            })]),
    Suva1 = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[178.41037442793356, -18.16178121073107],
          [178.4431617509316, -18.170915089195802],
          [178.48951032270895, -18.1190413657963],
          [178.42324903120505, -18.11643096063503]]]),
    Fijivec = ee.FeatureCollection("users/kishan2196/Fiji"),
    image = ee.Image("users/kishan2196/Fiji_islands/landcover/2019_FijiAll_split_urb_veg_class"),
    geometry2 = /* color: #d63000 */ee.Geometry.MultiPoint(),
    Bareland = /* color: #98ff00 */ee.FeatureCollection(
        [ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.7451428346762, -20.642181708710954],
                  [-178.74494971562711, -20.642914630412296],
                  [-178.74459566403715, -20.643165630183027],
                  [-178.74439181615202, -20.643015030370304],
                  [-178.74439181615202, -20.642442749721877],
                  [-178.74478878308622, -20.642151588565483]]]),
            {
              "landcover": 6,
              "system:index": "0"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.73453201581327, -20.643085310301462],
                  [-178.7342316084036, -20.643506989206454],
                  [-178.73394192983, -20.643386509638614],
                  [-178.7341565065512, -20.642994950383994]]]),
            {
              "landcover": 6,
              "system:index": "1"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.73625003744667, -20.637814526777724],
                  [-178.7356814091355, -20.637764325090654],
                  [-178.73581551958625, -20.637538417293836],
                  [-178.73612665583198, -20.637568538352784]]]),
            {
              "landcover": 6,
              "system:index": "2"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.73884641577308, -20.63785970828192],
                  [-178.73866402556007, -20.637970151902366],
                  [-178.73853527952735, -20.637824567113164],
                  [-178.73870157648628, -20.63779946627335],
                  [-178.738830322519, -20.637638820800504],
                  [-178.73903417040412, -20.637638820800504],
                  [-178.73904489924018, -20.637904889772713]]]),
            {
              "landcover": 6,
              "system:index": "3"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.73325669218605, -20.636805469687676],
                  [-178.7329777424485, -20.63715688336417],
                  [-178.73262369085853, -20.636880772686684],
                  [-178.73306357313697, -20.63668498481173]]]),
            {
              "landcover": 6,
              "system:index": "4"
            }),
        ee.Feature(
            ee.Geometry.Polygon(
                [[[-178.73676502157753, -20.63785970828192],
                  [-178.73659336020057, -20.637759304921044],
                  [-178.73638951231544, -20.637824567113164],
                  [-178.73627149511879, -20.63771914355813],
                  [-178.73635732580726, -20.63766392166682],
                  [-178.73676502157753, -20.637638820800504],
                  [-178.73684548784797, -20.637478175158023],
                  [-178.7370868866593, -20.637488215515663]]]),
            {
              "landcover": 6,
              "system:index": "5"
            })]),
    table2 = ee.FeatureCollection("users/kishan2196/ROI/Fiji"),
    Qelelevu = ee.Image("users/kishan2196/Fiji_islands/NDSV/2020_S2_NDSV_Qelelevu_clipped"),
    FijiLeftOut = ee.FeatureCollection("users/kishan2196/FijiLeftOut"),
    table3 = ee.FeatureCollection("users/kishan2196/_Reprojected"),
    V2020 = ee.Image("users/kishan2196/Fiji_islands/NDSV/2020_S2_NDSV_Vetaua_clipped"),
    image2 = ee.Image("users/kishan2196/Tuvalu/landcover/2020TuvaluLand"),
    Sigatoka = ee.Image("users/kishan2196/Fiji_islands/NDSV/2020_S2_NDSV_Sigatoka_clippedModified"),
    Siga = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[177.4731567311175, -18.15516857119123],
          [177.4760535168536, -18.178207103590474],
          [177.5193336415179, -18.175577194766703],
          [177.50695256470516, -18.156086092977727]]]),
    OnoSar = ee.Image("users/kishan2196/Fiji_islands/NDSV/2020opt_sarono_clipped"),
    QelelevuSar = ee.Image("users/kishan2196/Fiji_islands/NDSV/2020opt_sarqelelvu_clipped"),
    SigatokaSar = ee.Image("users/kishan2196/Fiji_islands/NDSV/2020opt_sarSigatoka_clipped"),
    SuvaSar = ee.Image("users/kishan2196/Fiji_islands/NDSV/2020opt_sarSuvaall_clipped"),
    TuvanaIRaSar = ee.Image("users/kishan2196/Fiji_islands/NDSV/2020opt_sarTuvanaIRa_clipped"),
    FijiSar2019 = ee.Image("users/kishan2196/Fiji_islands/NDSV/2019opt_sarFiji_clipped"),
    image6 = ee.Image("projects/ee-kishan21962/assets/Tonga/NDSV/2021opt_sarTonga_clipped");
    ```

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
Now that a training dataset has been prepared, it is time to choose a classifier to use and train it. GEE has a many different classifiers available, each with a range of user adjustable parameters. In the code snippet below the classifier variable is defined as the GEE Support Vector Machine (SVM) classifier. The full list (as well as the adjustable parameters for each) can be viewed in the ee.Classifier section of the Docs tab to the left of the map view.
```javascript

// Combine required training samples into one combined sample set
//var join = training.merge(training2)

var classifier = ee.Classifier.libsvm();
```
To train the selected classifier, call the train() function on it. The required arguments are the input training features (the training data generated by the code snippet above), the class property (whatever property contains the class code) and the input properties (in this case the bands of the NDSV image being classified).


```javascript

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
```
With the classifier trained, it is simply a case of calling the classify() function on the image to be classified, using the trained classifier object as the argument.
```javascript

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
  
