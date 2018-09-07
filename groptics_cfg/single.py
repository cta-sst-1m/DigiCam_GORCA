# 'array' is the geometry distrubution of the telescopes: will be converted direcly in 'telescopes.corsikaIO'
# == Usually this is the one you need to edit! ==
array="""Telescope location in the field (TLLOC) the parameters are:
   -the telescope identification number 
   -the X (oriented towrd east) location cordinate in meters. 
   -the Y (oriented toward north)location cordinate in meters. 
   -the Z (oriented toward up) location cordinate in meters. 
   -the telescope rotation offset in meters.
        The telescope rotation offset is the distance below the focal point
        of the mirror along the optical axis about which the telescope
        rotates. (intersection of the elevation and azimuth axis)
	#AP: Offset to locate the position of the rotator, moving along the
	#    position axis with the 0 in the focal point. >0 toward the mirror.
	#    i.e.: rotation point in the center of the mirror => focal length
   -the pointing offset in the x direction
   -the pointing offset in the y direction
#AP: These are coordinates taken from CORSIKA card, in which the
#    coords are X(North), Y(West), Z(down). To be consistent I exchanged
#    X <-> Y and changed sign to this new X. Z is set as 0.

* TLLOC   1       0.0       0.0       0.0  5.6  0.0  0.0
"""

# 'arraystruct' is a list of telescope factories (it index+1 will be the standard number!) used to create
# 'telescopes.arraytel' and 'telescopes.mirrors'. Each factory is a tuple
# (<N. of telescopes>, <factory type>, <printing mode (see 'arraytel_format')>)
# Use 'len(telescopes)' to set 'all telescopes'.
# == Edit in case of a multi-factory array! ==
arraystruct = [('len(telescopes)', 'DC', 0)]

# 'mirrors' is a tuple to build 'telescopes.mirrors':
# (<help mirror design>, <mirror designs>, <help mirror element>, ,mirror elements>)
# <mirror designs/elements> are lists corresponding one to one with arraystruct: each factory will build
# its corresponding design/element for each telescope ID
mirrors=("""Mirror design (MIROR). The parameters are 
   -the telescope identification number
   -the mirror radius in meters
   -focal length in meters
   -the focusing error in meters (>0 when camera is too 
	far away from the dish)
   -the mirror type (1=DavisCotton) 
   -the number of mirror elements. 
#AP: Set according to the TDR (23Mar2015)""",
["""
* MIROR   {telID:>2} 2.0 5.6   0.0 1 18
"""],
"""Mirror element (MIREL) characterisation. For each we have 
   -the telescope identification number
   -the mirror element identification number
   -the element shape (1=circular, 2=hexagonal)
   -the rotation angle of facet dir. in tel.plane 
   -the external radius 
   -the curvature radius 
   -the x position of the element on the dish (in meters)
   -the y position of the element on the dish (in meters)
   -the maximum mis-alignement in degrees
   -the maximum blur radius in degrees
   -the degradation factor (1.0=perfect 0.0=missing)
   -the reflectivity curve identifier

     Note that the mirror misalignment and blur radius will together influence 
     the PSF 

#AP: taken from mirror_DC4m_78cm_s2cm_f5.600_nocent.dat in sim_telarray.
#    TO BE REVISED! It is assumed without mis-align., no blur and perfect (no degrad.)""",
["""
* MIREL   {telID:>2}  1 2   90.0   0.45   11.2    0.0   -1.6    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2}  2 2   90.0   0.45   11.2  0.693   -1.2    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2}  3 2   90.0   0.45   11.2  1.386   -0.8    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2}  4 2   90.0   0.45   11.2 -0.693   -1.2    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2}  5 2   90.0   0.45   11.2    0.0   -0.8    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2}  6 2   90.0   0.45   11.2  0.693   -0.4    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2}  7 2   90.0   0.45   11.2  1.386    0.0    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2}  8 2   90.0   0.45   11.2 -1.386   -0.8    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2}  9 2   90.0   0.45   11.2 -0.693   -0.4    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2} 10 2   90.0   0.45   11.2  0.693    0.4    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2} 11 2   90.0   0.45   11.2  1.386    0.8    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2} 12 2   90.0   0.45   11.2 -1.386    0.0    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2} 13 2   90.0   0.45   11.2 -0.693    0.4    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2} 14 2   90.0   0.45   11.2    0.0    0.8    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2} 15 2   90.0   0.45   11.2  0.693    1.2    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2} 16 2   90.0   0.45   11.2 -1.386    0.8    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2} 17 2   90.0   0.45   11.2 -0.693    1.2    0.0    0.0    1.0 {factoryID}
* MIREL   {telID:>2} 18 2   90.0   0.45   11.2    0.0    1.6    0.0    0.0    1.0 {factoryID}
"""]
)

# 'reflectivity' is a tuple to build 'telescopes.mirrors': (<help reflectivity>, <reflectivities>). 
# <reflectivities> is a list corresponding one to one with arraystruct: each factory will build
# its corresponding reflectivity
reflectivities=("""Mirror reflectivity curve (RFCRV). The flag is followed by 
   -the the reflectivity curve identification index 
   -the number of points given on the quantum efficiency spectral curve. 
This is followed by a serie of 
   -wave length(in nano-meter)
   -reflectivity(fraction of unity) 

You can add additional curves and assign any curve to any individual facet.  
(Hence, each standard telescope might have different mirror reflectivities). 
Mirror reflectivity

#AP: Taken from ref_AR100.dat in sim_telarray""",
["""
* RFCRV {factoryID} 100
300.0    0.810078
305.825  0.814127
310.995  0.818131
316.396  0.822058
322.222  0.825922
328.456  0.82967 
335.086  0.833246
342.363  0.836703
350.633  0.840155
359.959  0.843649
369.965  0.847095
380.005  0.850326
389.478  0.853187
398.036  0.855602
405.617  0.857577
412.357  0.859178
418.468  0.860489
424.155  0.86159 
429.585  0.86254 
434.88   0.863384
440.141  0.864147
445.468  0.864846
450.98   0.865491
456.823  0.866093
463.161  0.866661
470.149  0.867203
477.895  0.867726
486.422  0.868234
495.636  0.868726
505.336  0.869194
515.238  0.869627
525.035  0.87001 
534.453  0.870331
543.291  0.870581
551.451  0.870756
558.925  0.870857
565.778  0.870889
572.117  0.87086 
578.06   0.870775
583.718  0.870643
589.181  0.870467
594.518  0.870252
599.778  0.870001
604.993  0.869716
610.185  0.869398
615.365  0.869048
620.542  0.868668
625.72   0.868258
630.898  0.867819
636.077  0.867351
641.253  0.866855
646.427  0.866331
651.596  0.865781
656.762  0.865203
661.923  0.864599
667.081  0.863968
672.237  0.863307
677.389  0.862617
682.538  0.861895
687.684  0.861141
692.826  0.860354
697.965  0.859532
703.1    0.858677
708.234  0.85779 
713.365  0.856872
718.493  0.855925
723.617  0.85495 
728.736  0.853949
733.849  0.852923
738.956  0.851875
744.055  0.850808
749.149  0.849729
754.237  0.848649
759.322  0.847581
764.405  0.84654 
769.488  0.845543
774.576  0.844608
779.673  0.843753
784.785  0.842994
789.915  0.842351
795.067  0.841845
800.242  0.841501
805.433  0.841345
810.635  0.841405
815.836  0.841705
821.027  0.842263
826.197  0.84309 
831.342  0.844187
836.461  0.845545
841.556  0.847145
846.638  0.84896 
851.726  0.850959
856.861  0.853118
862.122  0.855425
867.673  0.857883
873.803  0.860462
880.895  0.863016
888.965  0.86514 
896.451  0.866203
900.0    0.866279
"""]
)

# 'arraytel_format' is the tuple to build 'telescopes'.arraytel: (<help>,<format>)
# == THIS SHOULD NEVER BE TOUCHED!
arraytel_format=("""Using grisudet coor. system: ground system X(East), Y(North), Z(up)
  - number of telescope in the array (don't use 0)
    Need not be sequential, can be a subset of the 
    array used to create the photon file.
  - standard number (see factory configuration file, e.g. veritas.cfg and
    stdSCTelescopes.cfg
  - x telescope location (meters)
  - y telescope location (meters)
  - z telescope location (meters)
  - pointing offset x (Left on tangent plane) in degrees
  - pointing offset y (Down on tangent plane) in degrees
  - telescope print mode: fully implemented for DC Telescopes
    options 1,2,3 identical for SCTelescopes
     0:  no printing
     1:  print summary information
     2:  add geometry details 
     3:  add facet details
  Tangent plane axes are parallel to the ground system when 
  telescope at zenith, i.e. X(East), Y(North), Z(up). When telescope
  is at stow position, telescope coor. are X(East), Y(down), Z(North) 

    here's a 5 telescope array with a mixture of DC and SC telescopes
 ARRAYTEL 1 1 SC   91.84  -6.4   -0.25 0.0 0.0 0
 ARRAYTEL 2 1 DC    1.87 -49.1   -2.59 0.0 0.0 1
 ARRAYTEL 3 2 DC  -14.36  60.7    2.86 0.0 0.0 0
 ARRAYTEL 4 1 DC  -79.36  11.7   -0.02 0.0 0.0 0
 ARRAYTEL 5 1 SC  0.0  0.0   0.0 0.0 0.0 3

Use the following array for the Config/photon.cph test file from
the svn repository, can replace a DC telescope with an SC telescope;
Be sure to activate the SC factory with a new TELFAC line above.

#AP: These are coordinates taken from CORSIKA card, in which the
#    coords are X(North), Y(West), Z(down). To be consistent I exchanged
#    X <-> Y and changed sign to this new X. Z is set as 0.""",
"""
* ARRAYTEL {telID:>2} {factoryID} {factoryType:>6}   {x:>11.3f}  {y:>11.3f}  {z:>11.3f}  {xoffset:>3.2f}  {yoffset:>3.2f}  {print_mode}
""")

if __name__=='__main__':
    import os, sys
    output = sys.argv[1] if len(sys.argv)>1 else '.'
    telescopes = [tel_line for tel_line in array.splitlines() if tel_line.startswith('*')]
    corsikaIO_help = '\n'.join([tel_line for tel_line in array.splitlines() if not tel_line.startswith('*')])+'\n'
    start, mirror_reflectivities, mirror_desings, mirror_elements, arraytels = 0, '', '\n', '', '\n'
    for ifactory,factory_def in enumerate(arraystruct):
        factory = dict(factoryID=ifactory+1, factoryType=factory_def[1], print_mode=factory_def[2])
        mirror_reflectivities += reflectivities[1][ifactory].format(**factory)
        ntels = eval(factory_def[0]) if isinstance(factory_def[0],str) else factory_def[0]
        for tel in telescopes[start:start+ntels]:
            _,_,telID, x, y, z, _, xoffset, yoffset = tel.split()
            tel = dict(telID=int(telID), x=float(x), y=float(x), z=float(z), xoffset=float(xoffset), yoffset=float(yoffset), **factory)
            mirror_desings += mirrors[1][ifactory].format(**tel).strip()+'\n'
            mirror_elements += '\n'+mirrors[3][ifactory].format(**tel).strip()+'\n'
            arraytels += arraytel_format[1].format(**tel).strip()+'\n'
        start += ntels
    with open(os.path.join(output,'telescopes.corsikaIO'),'w') as f_corsikaIO:
        f_corsikaIO.write(array)
        #f_corsikaIO.write(corsikaIO_help+'\n'.join(sorted(telescopes, key=lambda x: int(x.split()[2])))) # have a sorted telescope. Is it better?
    with open(os.path.join(output,'telescopes.mirrors'),'w') as f_mirrors:
        f_mirrors.write('\n'.join([mirrors[0],mirror_desings,reflectivities[0],mirror_reflectivities,mirrors[2],mirror_elements]))
    with open(os.path.join(output,'telescopes.arraytel'),'w') as f_arraytel:
        f_arraytel.write('\n'.join([arraytel_format[0],arraytels]))
    
