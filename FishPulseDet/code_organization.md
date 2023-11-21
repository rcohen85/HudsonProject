## Description of code to be translated to python

* [de_detector.m](https://github.com/ScrippsWhaleAcoustics/DE_Detector/blob/master/de_detector.m)

  - Sets input/output locations
  - Loads settings from matlab scripts:
    - general low resolution settings: [dLoad_STsettings_broad.m](https://github.com/ScrippsWhaleAcoustics/DE_Detector/blob/master/dLoad_HRsettings_broad.m)
    
    - general high resolution settings: [dLoad_HRsettings_broad.m](https://github.com/ScrippsWhaleAcoustics/DE_Detector/blob/master/dLoad_STsettings_broad.m)

  Calls: [detection/dtSTbatch.m](https://github.com/ScrippsWhaleAcoustics/DE_Detector/blob/master/detection/dtST_batch.m) and  [detection/dHighres_click_batch.m](https://github.com/ScrippsWhaleAcoustics/DE_Detector/blob/master/detection/dHighres_click_batch.m) 
  

* [detection/dtSTbatch.m](https://github.com/ScrippsWhaleAcoustics/DE_Detector/blob/master/detection/dtST_batch.m) 
  - Loads 75 second chunks of data into memory and uses simple estimate of detection threshold to identify times exceeding that threshold.
  - Saves times to text file to be revisited later in high res step.
  
* [detection/dHighres_click_batch.m](https://github.com/ScrippsWhaleAcoustics/DE_Detector/blob/master/detection/dHighres_click_batch.m) 
  - Loads text file from dtSTbatch step, and revisits times that were flagged as exceeding the amplitude threshold.
  
  
  
**NOTE: This detector is no longer under active developement**  
A newer version without the 2 stage design is being maintained, revised, and further developed at:  
https://github.com/kfrasier/triton1.93.20160524_remoras/tree/master/Remoras/SPICE-Detector  
as part of a NOAA toolbox project.
