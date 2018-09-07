geoConfig.pilot, specifies DC telescope structure
Used in root DC telescope shadowing ray tracing 

VERSION3.0 
24Jan2012

E. Pueschel

GEOST <standard number, for later reference>
      <standard identification, see enum in VGA_Definition.h, 0 refers to veritas-like
      structure>
This is the only implemented DC geometry
* GEOST 1 STANDARD 0

* GEOST 1 TOPVOL 0.0

focus box 
dimensions (3 values, 3rd is not important) and rotation1/2/3
* GEOST 1 FOCBOX 0.47 0.36 0.45 30.0 0.0 0.0

edge boxes, one for each corner of focus box, all the same size 
trapesoid dimensions (3 values), offset, rotation1/2/3
* GEOST 1 EDGEBOX 1 0.07 0.2 5. 0.06 -120.0 90.0 0.0 
* GEOST 1 EDGEBOX 2 0.07 0.2 5. 0.06 -60.0 90.0 0.0 
* GEOST 1 EDGEBOX 3 0.07 0.2 5. 0.06 120.0 90.0 0.0 
* GEOST 1 EDGEBOX 4 0.07 0.2 5. 0.06 60.0 90.0 0.0 

shutter 
long edge, thickness, rotation1/2/3
* GEOST 1 SHUTTER 0.52 0.01 90.0 70.0 0.0

quad arms diameter, unused value, unused value
* GEOST 1 QUADARMSIZE  0.06 0.06 0.0 

individual quad arms 
quad arm number, bottom location x and y (meters), top location fixed by
location and size of focus box
* GEOST 1 QUADARM 1 -0.95 1.8
* GEOST 1 QUADARM 2 0.95 -1.8
* GEOST 1 QUADARM 3 0.95 1.8
* GEOST 1 QUADARM 4 -0.95 -1.8

quadrupod supports
diameter, unused value, distance below focus
* GEOST 1 CROSSBAR 1 0.05 0.05 3.2
* GEOST 1 CROSSBAR 2 0.05 0.05 3.2

camera radius in meters, NOT USED since camera radius defined in veritas.cfg.
ray tracing with ROOT Geometry only for shadowing.
Leave in for now.
*  GEOST 1 CAMERARADIUS .3


facet from GrISU configuration file
