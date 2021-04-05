# ETBRASIL - geeSEBAL
<img src="https://github.com/et-brasil/EESEBAL/blob/master/Images/geeSEBAL_logo_update_cut.png?raw=true" width="200">


### geeSEBAL is a open-source implementation of Surface Energy Balance Algorithm for Land (SEBAL) using Google Earth Engine (GEE).

Input Collections
=================

* LANDSAT/LC08/C01/T1_SR
* LANDSAT/LE07/C01/T1_SR
* LANDSAT/LT05/C01/T1_SR
* LANDSAT/LT04/C01/T1_SR
* LANDSAT/LC08/C02/T1_L2
* LANDSAT/LE07/C02/T1_L2
* LANDSAT/LT05/C02/T1_L2
* LANDSAT/LT04/C02/T1_L2

Model Design
============
Working.


| SPACECRAFT_ID   |    Band Names                               |
| --------------- | ------------------------------------------- |
| LANDSAT_4       |    B1, B2, B3, B4, B5, B7, B6, pixel_qa     |           
| LANDSAT_5       |    B1, B2, B3, B4, B5, B7, B6, pixel_qa     | 
| LANDSAT_7       |    B1, B2, B3, B4, B5, B7, B6, pixel_qa     | 
| LANDSAT_8       |    B1, B2, B3, B4, B5, B6, B7, B10, pixel_qa| 


|Property       |    Description|
| --- | ---|
| system:index      |  - Landsat Scene ID |
|                   |  - Must be in the Earth Engine format (e.g. LC08_044033_20170716) |
|                   |  - Used to lookup the scene specific c-factor |
| system:time_start | - Image datetime in milliseconds since 1970 |
| SPACECRAFT_ID     | - Used to determine which Landsat type |
|                  | - Must be: LANDSAT_4, LANDSAT_5, LANDSAT_7, or LANDSAT_8 |
| SUN_ELEVATION     | - Used to correct Correct declivity and aspect effects from LST |
|                  | - Used to estimate instantaneous net radiation |


Model Output
------------
Working.

References
==========

#### [Bastiaanssen et al. 1998] [Bastiaanssen, W.G.M., Menenti, M., Feddes, R.A., Holtslag, A.A.M., (1998). A remote sensing surface energy balance algorithm for land (SEBAL): 1. Formulation. J. Hydrol. 212–213, 198–212.](https://doi.org/10.1016/S0022-1694(98)00253-4)
#### [Bastiaanssen et al. 1998] [Bastiaanssen, W.G.M., Pelgrum, H., Wang, J., Ma, Y., Moreno, J.F., Roerink, G.J., van der Wal, T., (1998). A remote sensing surface energy balance algorithm for land (SEBAL): 2. Validation. J. Hydrol. 212–213, 213–229.](https://doi.org/10.1016/S0022-1694(98)00254-6)
#### [Laipelt et al. 2020] [Laipelt, L.; Ruhoff, A.L.; Fleischmann, A.S.; Kayser, R.H.B.; Kich, E.M.; da Rocha, H.R.; Neale, C.M.U. Assessment of an Automated Calibration of the SEBAL Algorithm to Estimate Dry-Season Surface-Energy Partitioning in a Forest–Savanna Transition in Brazil. Remote Sens. 2020, 12, 1108.](https://doi.org/10.3390/rs12071108)
#### [Laipelt et al. submitted] [Laipelt, L.; Kayser, R.H.B.; Fleischmann, A.S; Ruhoff, A.L; Bastiaanssen, W.; Erickson, T.; Melton, F. Long-term monitoring of evapotranspiration using the SEBAL algorithm and Google Earth Engine cloud computing.]
 
