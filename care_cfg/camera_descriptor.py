#!/usr/bin/env python

import sys,os,docopt

if __name__=='__main__':
    version = "Version: 0.1 (A. Porcelli, ale_led@yahoo.it)" \
              "Version: 0.1.1 (I. Al Samarai, imen.alsamarai@gmail.com)"
    doc = """
Tool to build/modify the camera descriptor.

Usage:
  %(prog)s make <output> [-a ROT -c CELLS... -d SPACE... -g GROWTH -l LINE -m MIN -n NPIXELS -r R -s -t TELTYPE -x CROSSTALK...]
  %(prog)s update <camera_descriptor> [<output> -c CELLS... -d SPACE... -p PTYPE -r R -s -t TELTYPE -x CROSSTALK...]

Info options:
  -h --help                  Print this help
  -v --version               Print the version

Commands:
  make                       Create a camera descriptor from the scratch
  update                     Take a camera descriptor and update the values, rebuilding neighbors and groups

Arguments:
  <camera_descriptor>        File containing the pixels description (lines starting with '* PMPIX'). The
                             results will be appended here, unless <output> is defined or the -s option is
                             selected
  <output>                   Output name where the camera descriptor is saved. In 'update' mode, if defined,
                             the resulting lines will be appended to this file, otherwise into <camera_descriptor> 
General options:
  -a --angle-rotation ROT    Rotation of the angle [default: 90.0]
  -c --cells CELLS           Number of cells per pixels [default: 36840]
  -d --dead-space SPACE      Coronal dead space surrounding the pixel. If multiple values are set, the first
                             will be the x dead space, the second the y (the remainings are ignored),
                             relatively to pixel type=0 and rotation=90.0. [default: 0.]
  -g --growth GROWTH         Pixel incrementation line by line (=face-to-face adjacent pixels) [default: 3]
  -l --line-length LINE      Number of maximum face-to-face adjacent pixels (=line) [default: 36]
  -m --minimum MIN           How many pixels are in the shortest line (=face-to-face adjacent pixels) [default: 1]
  -n --n-pixels NPIXELS      Number of pixels [default: 1296]
  -p --pixel-type PTYPE      Pixel type (0: hexagonal; 1: squared) [default: 0]
  -r --radius R              External radius of a pixel [default: 1.0]
  -s --silenced              If this option is set, the results is only printed, not saved
  -t --tel-type TELTYPE      Telescope type [default: 0]
  -x --x-talk CROSSTALK      Optical cross talk [default: 0.18]

%(version)s
Enjoy!

""" % dict(prog=os.path.basename(sys.argv[0]), version=version)
    args = docopt.docopt(doc, help=True, version=version)

import numpy as np

def slice_info(line,splitted_slice=[0,1],columns_skipped=0):
    """
    It will split a first part of the line (use splitted_slice to define it)
    by " " and the last elemnts of the list is the remaining of the line.

    Parameters
    ----------
    line: string
        line to be splitted
    splitted_slice: int or list or tuple (default [0,2])
        indeces of the slice that must be splitted. if int is given or
        len(splitted_slice)<2, splitted_slice=[0,splitted_slice]
    columns_skipped: int (default: 0)
        after the end of the slice, the skipping of <columns_skipped> colums
        before get the remaing part of the string is possible

    Returns
    ----------
    splitted_line: list of string
        It will be built as elements[splitted_slice[0]:splitted_slice[1]]+
        " ".join(elements[splitted_slice[1]+columns_skipped:])
    """
    splitted = line.split()
    if isinstance(splitted_slice,int):
        splitted_slice = [0,splitted_slice]
    elif len(splitted_slice)<2:
        splitted_slice = [0,splitted_slice[0]]
    return splitted[splitted_slice[0]:splitted_slice[1]]+[" ".join(splitted[splitted_slice[1]+columns_skipped:])]

def read_pixels_coordinates(camera_descriptor, tel_type='0'):
    """
    It reads the camera configuration file in CARE format and returns
    the relevant informations.
    
    Parameters
    ----------
    camera_descriptor: string
        file name of the configuration file
    tel_type: string (defualt: '0')
        type of telescioe to be selected
    
    Returns
    ----------
    pixel_type: string
        pixel type
    pixels: ndarray with record acces ['x','y','r','rotation']
        pixels coordinates
    """
    with file(camera_descriptor) as f:
        lines_array = np.array([slice_info(l,[1,3]) for l in f.readlines() if l[0]=='*'])
    pixels = np.array([tuple(l[0].split()) for l in lines_array[(lines_array[:,0]=="PMPIX")&(lines_array[:,1]==tel_type)][:,2:]],
                      dtype=np.dtype({'names':['type','id','x','y','r','rotation'],'formats':[np.int]*2+[np.float_]*4}))
    pixel_type = pixels['type'][0]
    if not all(pixel_type==pixels['type']):
        print(" THE PIXELS ARE NOT ALL OF THE SAME TYPE!!!")
        sys.exit(130)
    pixels.sort(order="id")
    return pixels

def set_pixel_coordinates(pixels, side=1., space=0., centering=False, tel_type='0', output=None, rotation=90., cells=36840, crosstalk=0.18):
    step_x = lambda ptype: 1.5 if int(ptype)==0 else 1.
    step_y = lambda ptype: .5*np.sqrt(3) if int(ptype)==0 else 1.
    space = list(space)[:2] if isinstance(space,(list,tuple,np.ndarray)) else [space]*2
    if len(space)==1:
        space += [space[0]]
    for i,p in enumerate(pixels):
        x = step_x(p['type'])*(side+float(space[1])) # type==0 & rotation==90 => x<->y
        y = step_y(p['type'])*(side+float(space[0])) # type==0 & rotation==90 => x<->y
        rot = round(p['rotation']*np.pi/180., 1)
        xproj, yproj = round(np.cos(rot),1), round(np.sin(rot),1)
        pixels['x'][i] *= xproj*x + yproj*y # In positional coordinates, the y axis is negaitve!!! => x*cos - y*sin -> x*cos + y*sin
        pixels['y'][i] *= yproj*x - xproj*y # In positional coordinates, the y axis is negaitve!!! => x*sin + y*cos -> x*sin - y*cos
        pixels['r'][i] = side
    if centering:
        x_center, y_center = pixels['x'].mean(), pixels['y'].mean()
        pixels['x'] -= x_center
        pixels['y'] -= y_center
    plines = "#\n# Pixel definiton from SST-1M_camera_descriptor\n# telescope type | pixel type (0=hex,1=sq) | pixel ID | x[mm] | y[mm] | r[mm] | rotation angle\n#\n"
    for i,p in enumerate(pixels):
        plines += "* PMPIX  {:>1} {:>1} {:>4} {:> 4.2f} {:> 4.2f} {:> 4.1f} {:> 4.1f}".format(tel_type, p['type'], i, p['x'], p['y'], p['r'], p['rotation'])+'\n'
    plines += "# SiPM Parameters, pixel by pixel\n# the number of cells should be the effective number of cells,\n# i.e. number cells - average number of cells in recovery state\n#\n# telescope type | pixel ID | Number of cells | optical cross-talk\n#\n"
    cells = [int(cells)]*len(pixels) if isinstance(cells, (int,float)) else list(cells)
    if len(cells)<len(pixels):
        cells += [cells[-1]]*(len(pixels)-len(cells))
    crosstalk = [int(crosstalk)]*len(pixels) if isinstance(crosstalk, (int,float)) else list(crosstalk)
    if len(crosstalk)<len(pixels):
        crosstalk += [crosstalk[-1]]*(len(pixels)-len(crosstalk))
    for i in xrange(len(pixels)):
        plines += "* SIPMPARS {:>1} {:>4} {:>5} {:>4}".format(tel_type, i, cells[i], crosstalk[i])+'\n'
    print plines
    if output is not None:
        with file(output, 'a') as f:
            f.write(plines)
    return pixels

def coordinates2positions(pixels):
    xmin, ymin = pixels['x'].min(), pixels['y'].min()
    upper_layer = np.where(pixels['y']>pixels['y'][0])[0][0]
    y_upper = pixels['y'][upper_layer]
    y_step = np.abs(y_upper-pixels['y'][0])
    next_colum = np.where((pixels['y']>pixels['y'][0])&(pixels['x']>pixels['x'][0]))[0][0]
    x_next = pixels['x'][next_colum]
    x_step = np.abs(x_next-pixels['x'][0])
    pixels['x'] -= xmin
    pixels['y'] -= ymin
    for i,p in enumerate(pixels):
        pixels['x'][i] = round(pixels['x'][i]/x_step,0)
        pixels['y'][i] = round(pixels['y'][i]/y_step,0)
        pixels['r'][i] = -1
    return pixels

#TO BE FINISHED
def make_pixels_positions(npixels=1296, line_length=36, minimum=1, growth=3, rotation=90.0, ptype=0):
    dtype = np.dtype({'names':['type','id','x','y','r','rotation'],'formats':[np.int]*2+[np.float_]*4})
    rot = rotation*np.pi/360.
    xproj, yproj = round(np.cos(rot),1), round(np.sin(rot),1)
    growing, rest, decrising_line, remaining_pixels, pixels, elements, line, iP = True, 0, int(npixels/line_length), npixels, [], minimum, 0, 0
    
    if ptype==0:
        while remaining_pixels!=0:
            for e in xrange(elements):
                y = 2*(e - int(elements/2)) + (0 if line%2==0 else 1)
                # In positional coordinates, the y axis is inverted: x*cos - y*sin -> x*cos + y*sin; x*sin + y*cos -> x*sin - y*cos  CHECK!!!
                pixels += [(ptype, iP, xproj*line+yproj*y, yproj*line-xproj*y, -1, rotation)]
                iP += 1
                remaining_pixels -= 1
            line += 1
            if elements<line_length:
                elements += (+1 if growing else -1)*growth
                if elements>line_length:
                    growing = False
                    rest = elements-line_length
                    elements = line_length
                    decrising_line = decrising_line #TBD
            elif line==decrising_line:
                elements -= rest if rest!=0 else growth
    else:
        raise SystemExit("Option not yet implemented") #Not yet implemented
    pixels = np.array(pixels,dtype=dtype)
    return pixels


def make_neighbors(baricenters, label, output=None):
    gnlines = "#\n# %s neighbors\n# telescope type | group ID | Number of Neighbor | list with neighbors\n#\n"%('Pixels' if label=='PIXNGHBR' else 'Groups')
    neigbors = []
    for i, barc in enumerate(baricenters):
        distances = np.round(np.sqrt((baricenters['x']-barc['x'])**2 + (baricenters['y']-barc['y'])**2),3)
        closer_dist = np.ceil(distances[distances>0.].min())
        neighbors = np.where((distances>0.)&(distances<=closer_dist))[0]
        gnlines += "* {} {:>2} {:>4} {:>1}   ".format(label, tel_type,i,len(neighbors))+' '.join([str(n) for n in neighbors])+'\n'
    print gnlines
    if output is not None:
        with file(output, 'a') as f:
            f.write(gnlines)

# ONLY FOR type==0 & rotation==90.0
def make_triplets(pixels):
    used_pixels, iGroup, baricenters = [], 0, []
    glines = "#\n# Group definition\n# telescope type | N. pixels | group ID | ID of the grouped pixel\n#\n"
    for i1, pix in enumerate(pixels):
        if i1 in used_pixels:
            continue
        upper_y = pixels['y'][pixels['y']>pix['y']][0]
        i_upper_line = pixels['y']==upper_y
        i2 = np.where((i_upper_line)&(pixels['x']<pix['x']))[0][-1]
        i3 = np.where((i_upper_line)&(pixels['x']>pix['x']))[0][0]
        glines += "* GROUP {:>2} {:>2} {:>4}   ".format(tel_type,3,iGroup)+' '.join([str(i1),str(i2),str(i3)])+'\n'
        baricenters += [(round(pixels[[i1,i2,i3]]['x'].mean(),3), round(pixels[[i1,i2,i3]]['y'].mean(),3))]
        iGroup+= 1
        used_pixels += [i1,i2,i3]
    baricenters = np.array(baricenters, dtype=np.dtype({'names':['x','y'],'formats':[np.float_]*2}))
    print glines
    if output is not None:
        with file(output, 'a') as f:
            f.write(glines)
    return baricenters

if __name__=='__main__':
    output = None if args['--silenced'] else (args['<camera_descriptor>'] if args['<output>'] is None else args['<output>'])
    tel_type = args['--tel-type']
    side = float(args['--radius'])
    space = [float(ds) for ds in args['--dead-space']]
    if args['make']: #THIS WILL CAUSE AN ERROR FOR NOW: make_pixels_positions must be finished
        Pixels = make_pixels_positions(int(args['--n-pixels']), int(args['--line-length']), int(args['--minimum']),
                                       int(args['--growth']), float(args['--angle-rotation']), int(args['--pixel-type']))
    elif args['update']:
        # Read the existing descriptor of the pixels positions
        Pixels = read_pixels_coordinates(args['<camera_descriptor>'], tel_type)
        # Trasnfrom in positional coordinates: (0,0) is  in the bottom-left corner
        Pixels = coordinates2positions(Pixels)
    else:
        raise SystemExit("Option not yet implemented") #Not yet implemented
    # Creat the coordinates according to the specifications through the options
    Pixels = set_pixel_coordinates(Pixels, side, space, tel_type=tel_type, output=output, cells=args['--cells'], crosstalk=args['--x-talk'], centering=True)
    pixels = Pixels[['x','y', 'r']]
    # pixel neighbors
    make_neighbors(pixels, 'PIXNGHBR', output)
    # pixel groups
    baricenters = make_triplets(pixels)
    # group neighbors
    make_neighbors(baricenters, 'GRPNGHBR', output)
