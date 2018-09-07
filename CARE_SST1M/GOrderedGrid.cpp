/*!  GOrderedGrid.cpp

     Version 1.0
     June 17, 2011
     Charlie Duke
     Grinnell College
 */
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <deque>
#include <list>
#include <iterator>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include <limits>

using namespace std;

#include "TMatrixD.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/RotationX.h"
#include "Math/RotationY.h"
#include "Math/RotationZ.h"
#include "Math/Rotation3D.h"

#include "GOrderedGrid.h"

#define DEBUG(x) cout << #x << " = " << x << endl
#define DEBUGS(x) cout << "       "<< #x << " = " << x << endl;

// GridElem constructor
GridElem::GridElem(const int &elemNum, const double &dist) 
  : elemNum(elemNum),dist(dist) {};
/*************** end of GridElem::GridElem *********/

// used in stl sort algorithm
bool gridBinSort(const GridElem &s1, const GridElem &s2) {
  return (s1.dist < s2.dist);
}
/**************** end of gridBinSort ************/

GOrderedGrid::GOrderedGrid( const vector<double> &xe,const vector<double> &ye,
                            const vector<double> &re,
                            const int &nbinsx, const int &nbinsy,
                            const int &option,
                            const string &filename) :vmx(xe),vmy(ye),vmr(re),
                                                     nbinsx(nbinsx),
                                                     nbinsy(nbinsy),
                                                     iGridOption(option),
                                                     sFileName(filename)
 {
  bool debug = false; 

  bGridOption = true;
  if (iGridOption==0 ) bGridOption = false;

  if (debug) {
    cout << "  -- GOrderedGrid::GOrderedGrid " << endl;
    cout << "       nbinsx   " << nbinsx << endl;
    cout << "       nbinsy   " << nbinsy << endl;
    cout << "       iOption  " << iGridOption << endl;
    cout << "       gridfilename " << sFileName << endl;
  }

  if ( (nbinsx == 0) || (nbinsy == 0) ) {
    cout << "     returning from GOrderedGrid constructor: " << endl; 
    cout << "          nbinsx == 0 or nbinsy == 0. " << endl;
    cout << "          full element loop search forced " << endl;
  }
  else {

    if (iGridOption == 0) {
      return;
    }
    else if (iGridOption==1) {
      initialize();
      makeGridParameters();
      makeGrid();
    }
    else if (iGridOption==2) {
      initialize();
      makeGridParameters();
      makeGrid();
      printGrid();
    }
    else if (iGridOption==3) {
      if (!readGrid()) {
        initialize();
        makeGridParameters();
        makeGrid();
        printGrid();
      }
    }
    else {
      cout << "unknown iGridOption in GOrderGrid constructor: "
            << iGridOption << endl;
      exit(0);
    }
  }
};
/********************* end of GOrderedGrid ***********/

GOrderedGrid::~GOrderedGrid() {

  for (unsigned i = 0;i<vGrid.size();i++) {
    delete vGrid[i]; 
  }

};
/********************* end of ~GOrderedGrid ***********/

GOrderedGrid::GOrderedGrid(const GOrderedGrid &orderedGrid) {

};
/********************* end of GOrderedGrid ***********/

void GOrderedGrid::initialize() {

  fXmin =  99999999.9;
  fXmax = -99999999.9;
  fYmin =  99999999.9;
  fYmax = -99999999.9;

  fDelX = 0.0;
  fDelY = 0.0;

};
/********************* end of initialize ***********************/

bool GOrderedGrid::makeGridParameters() {
  bool debug = false;
  if (debug) {
    cout << "  -- GOrderedGrid::makeGridParameters " << endl;
  }

  // find fXmin,fXmax,fYmin,fYmax
  for (unsigned i = 0;i<vmx.size();i++) {
    double diff;
    diff = vmx[i] - vmr[i];
    if ( diff < fXmin ) fXmin = diff;
    diff = vmx[i] + vmr[i];
    if (diff > fXmax) fXmax = diff;
    diff = vmy[i] - vmr[i];
    if ( diff < fYmin ) fYmin = diff;
    diff = vmy[i] + vmr[i];
    if (diff > fYmax) fYmax = diff;
  }
  fDelX = (fXmax-fXmin)/nbinsx;
  fDelY = (fYmax-fYmin)/nbinsy;

  if (debug) {
    DEBUGS(fXmin);DEBUGS(fXmax);
    DEBUGS(fYmin);DEBUGS(fYmax);
    DEBUGS(fDelX);DEBUGS(fDelY);
  }

  return true;
};
/******************** end of makeGridParameters **********************/

bool GOrderedGrid::makeGrid() {

  bool debug = false;
  if (debug) {
    cout << "  -- GOrderedGrid::makeGrid " << endl;
  }

  int numBins = nbinsx*nbinsy;

  // prepare the grid since we know the number of bins
  for (int i = 0;i<numBins;i++) {
    vGrid.push_back(new list<GridElem> );
  }

  if (debug) {
    cout << "vGrid size " << vGrid.size() << endl;
  }  

  double grid_dia = sqrt(fDelX*fDelX + fDelY*fDelY) / 2.0;

  // loop over grid bins
  for (int gb = 0;gb < numBins; ++gb) {
    // get bin center location
    int ybin = gb/nbinsx;
    int xbin = gb % nbinsx;

    double xc = (fXmin + xbin*fDelX) + (fDelX/2);
    double yc = (fYmin + ybin*fDelY) + (fDelY/2);

    // loop over elements
    for (unsigned elem=0;elem<vmx.size();++elem) {
      
      // get x,y,radius for this element
      double xe = vmx[elem];
      double ye = vmy[elem];
      double ra = vmr[elem];

      // get distance from bin center to element location
      double dst = sqrt( ((xe-xc)*(xe-xc)) + ((ye-yc)*(ye-yc)) );

      // test distance
      if (dst < (ra + grid_dia) ) {
        if (debug) {
          cout << "dst " << dst << endl;
        }
        //cout << " gb elem dst  " << gb << " " << elem 
        //    << " " << dst << endl;
        // insert this element into the grid
        GridElem grdf(elem,dst);
        
        vGrid[gb]->push_back(grdf);
      }
    }
  }

  // sort the lists
  for (int gb = 0;gb < numBins; ++gb) {
    vGrid[gb]->sort(gridBinSort);
  }

  return true;
};
/********************* end of makeGrid ************************/

bool GOrderedGrid::getGridBinList(const double &x, const double &y,
                      list<GridElem> *&gridList) {

  bool debug = false;
  if (debug) {
    cout << "  -- GOrderedGrid::getGridBinList" << endl;
  }
  if ( (nbinsx==0) || (nbinsy==0) ||
       (!bGridOption) ) return false;

  bool gridFound = true;

  int xbin,ybin;

  xbin = (int)( floor( (x - fXmin)/fDelX ));
  ybin = (int)( floor( (y - fYmin)/fDelY ));

  int gridkey = (xbin + (nbinsx*ybin) );

  if (debug) {
    DEBUGS(x);DEBUGS(xbin);
    DEBUGS(y);DEBUGS(ybin);
    DEBUGS(gridkey);
  }
 
  if ( (gridkey < 0) || (gridkey > (nbinsx*nbinsy) - 1) ) {
    gridFound = false;
  }
  else {
    gridList = vGrid[gridkey];
    if (debug) {
      cout << "       gridList->size() " << gridList->size() << endl;
    }
  }

  if (debug) {
    cout << "        gridFound " << gridFound << endl;
  }

  return gridFound;
};
/********************* end of getGridBinList ************************/

bool GOrderedGrid::readGrid() {

  bool debug = false;

  if (debug) {
    cout << "  -- GOrderedGrid::readGrid" << endl;
  }

  ifstream inFile(sFileName.c_str(),ios::in);
  if (! inFile) {
    return false;   // no such file
  }
  string fileRec;
  vector<string> tokens;

  // skip the first three documentation records
  getline(inFile,fileRec,'\n');
  getline(inFile,fileRec,'\n');
  getline(inFile,fileRec,'\n');

  // read grid parameters
  getline(inFile,fileRec,'\n');
  tokenizer(fileRec,tokens);

  nbinsx = atoi(tokens[0].c_str());
  fDelX = atof(tokens[1].c_str());
  fXmin = atof(tokens[2].c_str());
  nbinsy = atoi(tokens[3].c_str());
  fDelY = atof(tokens[4].c_str());
  fYmin = atof(tokens[5].c_str());

  if (debug) {
    DEBUGS(nbinsx);DEBUGS(fDelX);DEBUGS(fXmin);
    DEBUGS(nbinsy);DEBUGS(fDelY);DEBUGS(fYmin);
  }

  int numBins = nbinsx*nbinsy;
    
  // start with a new grid if it's already filled
  if (vGrid.size() > 0) {
    for (unsigned i = 0;i<vGrid.size();i++) {
      delete vGrid[i];
    }
    
    vGrid.clear();
  }

  // prepare the grid since we know the number of bins
  for (int i = 0;i<numBins;i++) {
    vGrid.push_back(new list<GridElem> );
  } 
  
  for (int i = 0;i<nbinsx*nbinsy;i++) {
    //for (int i = 0;i<10;i++) {

    string elements;
    string distances;
    vector<string> tokElem;
    vector<string> tokDist;

    // read element list
    getline(inFile,elements,'\n');
    tokenizer(elements,tokElem);

    int gridkey = atoi(tokElem[0].c_str());
    int numElem = atoi(tokElem[3].c_str());

    // read distance list
    getline(inFile,distances,'\n');
    tokenizer(distances,tokDist);

    // loop over elements
    for (int el=0;el<numElem;++el) {

      int elem = atoi(tokElem[4+el].c_str());
      double dist = atof(tokDist[4+el].c_str());
 
      GridElem grdf(elem-1,dist);
      vGrid[gridkey]->push_back(grdf);
    }
  }

  return true;
};
/********************* end of readGrid ************************/

bool GOrderedGrid::printGrid() {

  bool debug = false;
  if (debug) {
    cout << "  -- GOrderedGrid::printGrid " << endl;
    cout << "      opening file " << sFileName << endl;
  }

  if (sFileName=="") {
    return false;
  }

  FILE *fg;
  if ( (fg= fopen(sFileName.c_str(),"w")) == NULL) {
    cerr << "    failed to open element grid file "
         << sFileName << endl;
    exit(0);
  }
  
  fprintf(fg,"first data record: nbinsx delx xmin nbinsy dely ymin\n");
  fprintf(fg,"listing 1: gridKey,xbin,ybin,numElements,list of elements\n");
  fprintf(fg,"  followed by same format, list of distances\n");

  fprintf(fg,"%d %lf %lf %d %lf %lf\n",nbinsx,fDelX,fXmin,
          nbinsy,fDelY,fYmin);

  int numBins = nbinsx*nbinsy;
  
  for (int i = 0;i<numBins;++i) {
    int xb,yb;
    yb = i / nbinsx;
    xb = i % nbinsx;

    // get list for this element
    list<GridElem> *lst = vGrid[i];
    list<GridElem>::iterator iterlst;
    
    unsigned numElem = vGrid[i]->size();
    fprintf(fg,"%4d %3d %3d %2d ",i, xb,yb,numElem);
 
    if (numElem == 0) {
      fprintf(fg,"\n");
    }
    else {
      
      for (iterlst = lst->begin();
           iterlst!= lst->end(); ++iterlst) {

        int elemN = iterlst->elemNum;
        fprintf(fg,"%3d ",elemN + 1);
        
      }
      fprintf(fg,"\n");
    
    }
    fprintf(fg,"%4d %3d %3d %2d ",i, xb,yb,numElem);
    if (numElem == 0) {
      fprintf(fg,"\n");
    }
    else {
      
      for (iterlst = lst->begin();
           iterlst!= lst->end(); ++iterlst) {

        double distg = iterlst->dist;
        fprintf(fg,"%5.2f ",distg);
        
      }
      fprintf(fg,"\n");
     
    }

   }

  fclose(fg);
  return true;
};

/********************* end of printGrid ************************/

void GOrderedGrid::tokenizer(const string& str, 
                             vector<string>& tokens) {

  tokens.clear();
  // tokenizer with multiple delimiters
  char delim[] = " ,";
  string::size_type lgt = 2;
  
  string::size_type idx = 0;
  string::size_type start;

  // find first non-delimiter character for start
  while ( (idx = str.find_first_not_of(delim,idx,lgt) ) !=string::npos ) {
    start = idx; // start of string

    idx = str.find_first_of(delim,idx,lgt);

    tokens.push_back(str.substr(start,idx - start));
  }

};
/********************  end of tokenizer ************/
