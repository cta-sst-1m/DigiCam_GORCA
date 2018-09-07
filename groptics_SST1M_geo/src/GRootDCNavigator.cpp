/*
VERSION3.0
24Jan2012
*/
// standard includes
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <list>
#include <iterator>
#include <algorithm>
#include <bitset>
#include <iomanip>

using namespace std;

#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "TMatrixD.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoShape.h"
#include "TGeoArb8.h"
#include "TGeoNode.h"
#include "TGeoMatrix.h"
#include "TView.h"
#include "TPad.h"
#include "TVector3.h"
#include "TGeoPgon.h"

//#include ".h"
#include "GDefinition.h"
#include "GUtilityFuncts.h"
#include "GGeometryBase.h"
#include "GDCGeometry.h"

#include "GTelescope.h"
#include "GDCTelescope.h"
#include "GRootDCNavigator.h"

ClassImp(GRootDCNavigator);

// debug macro
#define DEBUG(x) *oLog << #x << " = " << x << endl

//****************************************************

GRootDCNavigator:: GRootDCNavigator(GDCTelescope *dcTel,
                                    int makeFacetsFlag, bool debugT) {

  DCTel = dcTel;
  fDebugT = debugT;  
  fMakeFacetsFlag = makeFacetsFlag;

  if (fDebugT) {
    *oLog << "  -- GRootDCNavigator::GRootDCNavigator" << endl;
    *oLog << "       fMakeFacetsFlag: " << fMakeFacetsFlag << endl;
  }
  // set pointers to zero
  gDC = 0;
  fGeom = 0;
  fMatVacuum = 0;
  fMatAl = 0;
  fVacuum = 0;
  fAl = 0;
  fTopPosV = 0;
  fTopVol = 0;

  // initialize single variables
  fFL = 0.0;
  fnodeNum = 0;
  fepsil = 0.0; 
  fFocBoxZTop = 0.0;
  fMoveToTop = false;
  fNodeC = 0;
  fNodeN = 0;

  // initialize arrays
  for (int i=0;i<3;i++) {
    fTopDim[i] = 0.0;
    fFocusBoxDim[i] = 0.0;
    fFocusBoxDim2[i] = 0.0;
    fFocusBoxRot[i] = 0.0;
    fQuadArmPosV[i] = 0;
    fQuadArmR2R1V[i] = 0;
  }
  fQuadArmPosV[3] = 0;
  fQuadArmR2R1V[3] = 0;

  initialize();
  setupGeometry();

  // turn on via setTrackingDebug
  fDebugTr = true;

  // take this out at some point
  if (0) {
  // ============ set up default veritas telescope
    fFL = 12.0;            // focal length
    
    // default focus box dimensions (edge to edge, meters)
    // have to set focus box dimensions prior to setting
    // top dimensions
    // normally 1.5
    fFocusBoxDim[0] = 1.5;    // set focus box x dimension
    fFocusBoxDim[1] = 1.5;    // set focus box y dimension
    fFocusBoxDim[2] = 1.02;    // set focus box z dimension
    
    // focus box rotation (thru 45 degrees about orig. z axis)
    //fFocBoxRot[0] = 45.0;
    //fFocBoxRot[1] =  0.0;
    //fFocBoxRot[2] =  0.0;
    
    fepsil = 1.00;
    // Top volume dimensions (+-x,+-y,+-z) meters
    fTopDim[0] = 10.0;  // 6.5 covers all facets. set at 10.0 for now
    fTopDim[1] = 10.0;
    fTopDim[2] = (fFL + fFocusBoxDim[2] + fepsil)/2.0; 
    

    // location of Top center, in telescope coordinates
    //Double_t topx = 0.0;
    //Double_t topy = 0.0;
    //Double_t topz = (fFL + fFocusBoxDim[2] + fepsil)/2.0;
    //fTopPosV = new TVector3(topx,topy,topz);
  }

  if (fDebugT) {
    *oLog << "       fepsil: " << fepsil << endl;
    *oLog << "       FocusBoxDim, side length: " << fFocusBoxDim[0] << " " 
          << fFocusBoxDim[1] << " " << fFocusBoxDim[2] << endl;
    *oLog << "       Focal Length; " << fFL << endl;
    
    *oLog << "       TopVolume half-side length : " << fTopDim[0]  << " " 
          << fTopDim[1] << " " << fTopDim[2] << endl;
    *oLog << "       fTopPosV: top volume vector, TopCenter tele.coor." 
          << endl;
    fTopPosV->Print();
    double tmp1 = (*fTopPosV)[2] - fepsil - fFocusBoxDim[2];
    *oLog << "       Focal Point to TopVol center: " << tmp1 
          << endl; 
    *oLog << "       tele.base center to TopVol center: " 
          << (*fTopPosV)[2] << endl;
    *oLog << " *************End of constructor printing****************" 
          << endl;
    
  }
  // this creates a new TGeoManager for each telescope. 
  // If facets are added to the Geometry, then individual
  // geometries are required because of differing blur radii.  
  // However, at this time, root doesn't permit multiple
  // geometries. So, I'll leave the code like this for now.
  // but only the last geometry created will be open and will
  // serve for all telescopes. This will work fine for veritas
  // as long as the telescopes are identical.

  ostringstream os;
  os << "geomTel" << DCTel->iTelID;
  //*oLog << "geomTel " << os.str() << endl;
  string nameGeom = os.str();

  gGeoManager = 0;
  //*oLog << endl << " ********** new GeoManager ********" << endl << endl;
  //*oLog << "numGeoMan before " << gROOT->GetListOfGeometries()->GetSize() << endl;
  fGeom = new TGeoManager(nameGeom.c_str(), nameGeom.c_str());
  //*oLog << "numGeoMan after " << gROOT->GetListOfGeometries()->GetSize() << endl;

  gGeoManager = fGeom;

  makeMaterialMedia();
  makeTop();
  makeFocusBox();
  makeEdgeBoxes();
  makeQuadArms();
  makeCrossArms();
  makeSupportWires();
  makeShutter();
  if (fMakeFacetsFlag) {
    makeFacets();
  }
  //fGeom->SetVisLevel(4);
  
  //--- close the geometry
  fGeom->CloseGeometry();
  //fTopVol->Draw("ogl");
   
   //drawTelescope();
};
//****************************************************

GRootDCNavigator::  ~GRootDCNavigator() {

  bool debug = false;
  if (debug) {
    *oLog << "  -- GRootDCNavigator::~GRootDCNavigator" << endl;
  }
  // fGeom owns all the geometry classes (from root)
  if (fGeom) {
    gGeoManager = fGeom;
    SafeDelete(fGeom);
    // gGeoManager is set to zero ok

  }
  
  SafeDelete(fTopPosV);
  
  for (int i = 0; i<3;i++) {
    SafeDelete(fQuadArmR2R1V[i]);
    SafeDelete(fQuadArmPosV[i]);
  }

};
//****************************************************

void GRootDCNavigator::setupGeometry() {

  bool debug = fDebugT;
  if (debug) {
    *oLog << "  -- GRootDCNavigator::setupGeometry" << endl;
  }

  fFL = DCTel->dFocLgt;

  gDC = dynamic_cast<GDCGeometry*>(DCTel->geoStruct);

  // move structural parameters from gDC to here
  fepsil = gDC->epsil;

  // get focus box dimensions
  for (int i = 0;i<3;i++) {
    fFocusBoxDim[i] = gDC->focBoxDim[i];
    fFocusBoxRot[i] = gDC->focBoxRot[i];
    if (debug) *oLog << "        i fFocusBoxDim  " << i << endl;
  } 
  //elisa addition
  fFocusBoxDim2[0]=0.63;
  fFocusBoxDim2[1]=0.13;
  fFocusBoxDim2[2]=0.45;

  for (int i = 0;i<3;i++) {
    if (debug) *oLog << "   i fFocusBoxRot " << i << " " << fFocusBoxRot[i] << endl;
  }
 
  // calculate top volume
  //  add 5.0 meters to x/y dimensions to ensure that can move
  //  all photons to top of topVol, else get seg. fault in tracker
  fTopDim[0] = ( (fFL + fepsil)/2.0 ) + 5.0;
  fTopDim[1] = ( (fFL + fepsil)/2.0 ) + 5.0;
  fTopDim[2] = (fFL + fFocusBoxDim[2] + fepsil)/2.0; 

  for (int i = 0;i<3;i++) {
    if (debug) *oLog << "        i fTopDim  " << i << " " << fTopDim[i] << endl;
  }   

  Double_t topx = 0.0;
  Double_t topy = 0.0;
  Double_t topz = (fFL + fFocusBoxDim[2] + fepsil)/2.0;
  fTopPosV = new TVector3(topx,topy,topz);

  for (int i = 0;i<3;i++) {
    if (debug) *oLog << "        i fTopPosV  " << i << " " << (*fTopPosV)[i] << endl;
  }   
  if (debug) {
    *oLog << "      fFL " << fFL << endl;
    *oLog << "      fepsil " << fepsil << endl;
    *oLog << "  LEAVING SETUP GEOMETRY" << endl;
   }
  
};
//****************************************************

void GRootDCNavigator::initialize() {

  fnodeNum = 1;
  // initialize photon tracking variables
  for (int i = 0;i<3;i++) {
    fPosC[i] = 0.0;
    fDirC[i] = 0.0;
    fPosN[i] = 0.0;
    fDirN[i] = 0.0;
    fNodeC   = 0;
    fNodeN   = 0;
  }

};
//****************************************************

void GRootDCNavigator::makeMaterialMedia() {

   //--- define some materials
   fMatVacuum = new TGeoMaterial("fMatVacuum", 0,0,0);
   fMatVacuum->SetTransparency(80); 
   fMatAl = new TGeoMaterial("fMatAl", 26.98,13,2.7);

   //   //--- define some media
   fVacuum = new TGeoMedium("fVacuum",0, fMatVacuum);
   fAl = new TGeoMedium("fAl",1, fMatAl);

};
//****************************************************

void GRootDCNavigator::makeTop() {

  //--- make the top container volume
  //     +- x and y = (focal Length + 2 meters);
  //     +-z = (focal length + box z)/2.0
  //--- center of top at x=0,y=0,z = (focal length + box z)/2

  Double_t topx = fTopDim[0];
  Double_t topy = fTopDim[1];
  Double_t topz = fTopDim[2];

  // half side length of topbox is topx,topy,topz
  fTopVol = fGeom->MakeBox("fTopVol",fVacuum,topx,topy,topz);
  fGeom->SetTopVolume(fTopVol);
  //  fGeom->SetTopVisible(1);
  fGeom->SetTopVisible(0);

};
//****************************************************

void GRootDCNavigator::makeFocusBox() {

  TGeoCombiTrans *combiFocBox;
  TGeoVolume *fFocBoxVol;

  // dimensions of box, to be created at center to TOP
  Double_t boxR = fFocusBoxDim[0];
  Double_t boxD = fFocusBoxDim[1];
  Double_t boxZ = fFocusBoxDim[2]/2.0;

  // focus box position relative to TOP center
  Double_t boxPosX = 0.0;
  Double_t boxPosY = 0.0;
  Double_t boxPosZ = (*fTopPosV)[2] - fepsil - boxZ;
  fFocBoxZTop = boxPosZ;

  Double_t boxR1 = fFocusBoxRot[0];
  Double_t boxR2 = fFocusBoxRot[1];
  Double_t boxR3 = fFocusBoxRot[2]; 

  if (fDebugT) {
    *oLog << "********** focus box translation & rotation *************" << endl;
    *oLog << "   size+- " 
         << boxR << " " << boxD << " " << boxZ << endl;
    *oLog << "   position wrt TopVol:  "
         << boxPosZ << endl;
    *oLog << "   rot.A  "
         << boxR1 << " " << boxR2 << " " << boxR3 << endl;
    *oLog << "*****************************************" << endl;
  }

  fFocBoxVol = fGeom->MakePgon("fFocBoxVol", fAl, 0, 360, 6, 2);
  TGeoPgon *pgon = (TGeoPgon*)fFocBoxVol->GetShape();
  pgon->DefineSection(0, -boxD, 0.0, boxR);
  pgon->DefineSection(1, boxD, 0.0, boxR);

  fFocBoxVol->SetLineColor(kRed);
  fFocBoxVol->SetVisibility(kTRUE);
  
  combiFocBox = new TGeoCombiTrans("combiFocBox",boxPosX,boxPosY,boxPosZ,
  				   new TGeoRotation("rFocBox",boxR1,boxR2,boxR3)); 

  fTopVol->AddNode(fFocBoxVol, fnodeNum,combiFocBox);
  fnodeNum++;

  TGeoCombiTrans *combiFocBox2;
  TGeoVolume *fFocBoxVol2;

  fFocBoxVol2 = fGeom->MakePgon("fFocBoxVol2", fAl, 0, 360, 6, 2);
  TGeoPgon *pgon2 = (TGeoPgon*)fFocBoxVol2->GetShape();
  pgon2->DefineSection(0, -0.13, 0.55, 0.63);
  pgon2->DefineSection(1, 0.13, 0.55, 0.63);
  
  fFocBoxVol2->SetLineColor(kRed);
  fFocBoxVol2->SetVisibility(kTRUE);
  
  combiFocBox2 = new TGeoCombiTrans("combiFocBox2",boxPosX,boxPosY,boxPosZ,
  				   new TGeoRotation("rFocBox2",boxR1,boxR2,boxR3)); 
  
  fTopVol->AddNode(fFocBoxVol2, fnodeNum,combiFocBox2);
  fnodeNum++;
};
//****************************************************************

void GRootDCNavigator::makeEdgeBoxes() {

  if (fDebugT) {
    *oLog << "  -- GRootDCNavigator::makeEdgeBoxes" << endl;
    for (int i = 0;i<4;i++) {
      *oLog << " i gDC->edgeX/Y/Z " << gDC->edgeX[i] << "  " 
            << gDC->edgeY[i] << " " << gDC->edgeZ[i] << endl;
    }
    for (int i = 0;i<4;i++) {
      *oLog << " i rot X/Y/Z " << gDC->edgeRot1[i] << "  " 
	    << gDC->edgeRot2[i] << " " << gDC->edgeRot3[i] << endl;
    }
  }
  // load edgebox size, offset, and rotation vectors
  vector<double> vsizex;
  vector<double> vsizez;
  vector<double> vangle;
  vector<double> voffset;
  vector<double>  vdiagh;
  vector<double> vrot1;
  vector<double> vrot2;
  vector<double> vrot3;

  for (int i = 0;i<4;i++) {
    vsizex.push_back(gDC->edgeX[i]);
    vsizez.push_back(gDC->edgeY[i]);
    vangle.push_back(gDC->edgeZ[i]);
    voffset.push_back(gDC->edgeOffset[i]);

    vrot1.push_back(gDC->edgeRot1[i]);
    vrot2.push_back(gDC->edgeRot2[i]);
    vrot3.push_back(gDC->edgeRot3[i]);

    vdiagh.push_back(( (fFocusBoxDim2[0]/2.0)*sqrt(2.0) ) + voffset[i]);
   }

  double epsil = 0.03;

  // create edgebox volume vectors
  vector<TGeoVolume *> vEdgeBoxVol;
  for (int i = 0;i<4;i++) { 
    ostringstream os;
    os << "fEdgeBoxVol" << i+1;
    string edgebStr = os.str();
    vEdgeBoxVol.push_back(fGeom->MakeTrap(edgebStr.c_str(),fAl, vsizex[i], 0, 0, vsizez[i], vsizex[i], epsil, vangle[i], vsizez[i], vsizex[i], epsil, vangle[i]));
    vEdgeBoxVol[i]->SetLineColor(kRed);
    vEdgeBoxVol[i]->SetVisibility(kTRUE);
  }

  //  double vloc = vsizez[0]+epsil;
  double vloc = 0.4;
 
  // edgebox locations
  double aXbox[4] = {vloc+voffset[0],-vloc-voffset[0],vloc+voffset[0],-vloc-voffset[0]};
  double aYbox[4] = {vloc+voffset[0],vloc+voffset[0],-vloc-voffset[0],-vloc-voffset[0]};
  double aZbox[4] = {fFocBoxZTop,fFocBoxZTop,fFocBoxZTop,fFocBoxZTop,};

  for (int i = 0;i<4;i++) {
    ostringstream os;
    ostringstream os1;
    os << "combiEdgeBox" << i+1;
    os1 << "rEdgeBox" << i+1;
    string combi = os.str();
    string redge = os1.str();
  
    fTopVol->AddNode(vEdgeBoxVol[1], fnodeNum,new TGeoCombiTrans(combi.c_str(),
								 aXbox[i],aYbox[i],aZbox[i],
								 new TGeoRotation(redge.c_str(),vrot1[i],vrot2[i],vrot3[i])) );
    fnodeNum++;
  }

};
//****************************************************

void GRootDCNavigator::makeShutter() {

  Double_t fShutterDim[3];
  // set dimensions of shutter
  fShutterDim[0] = gDC->shutterX;   
  fShutterDim[1] = fFocusBoxDim[0]/2.;   // width of focus box
  fShutterDim[2] = gDC->shutterZ;   // make it thin

  Double_t xs = fShutterDim[0];
  Double_t ys = fShutterDim[1];
  Double_t zs = fShutterDim[2];


  // Set up translation of shutter from Top center 
  Double_t xt = xs;
  Double_t yt = 0;
  Double_t epsil = 0.01;
  Double_t zt = fFocBoxZTop - fFocusBoxDim[1] - (fFocusBoxDim[2])/2. + epsil;

  TGeoVolume *fShutterVol = fGeom->MakeTrap("fShutterVol",fAl, zs, 0, 0, ys, xs, xs/2., 0, ys, xs, xs/2., 0);

  fShutterVol->SetLineColor(kGreen);
  Double_t rotate1a = gDC->shutterRot1;
  Double_t rotate1b = gDC->shutterRot1 + 180;
  Double_t rotate2 = gDC->shutterRot2;
  Double_t rotate3 = gDC->shutterRot3;
  

  if (fDebugT) {
    *oLog << "************* shutter ************************" << endl;
    *oLog << "   shutter size: " << 2.*xs << " " << 2.*ys << " " << 2.*zs << endl;
    *oLog << "   translation:  " << xt << " " << yt << " " << zt << endl;
    *oLog << "   rotation:     " << rotate1a << "  " << rotate2 << "  " 
	  << rotate3 << endl;
    *oLog << "******************* end shutter ********************" << endl;
  }

  fTopVol->AddNode(fShutterVol,fnodeNum,
                   new TGeoCombiTrans("combiShutter1",xt,yt,zt,
                                      new TGeoRotation("rShutter",
                                                       rotate1a,rotate2,rotate3)) ); 
  fnodeNum++;

  fTopVol->AddNode(fShutterVol,fnodeNum,
                   new TGeoCombiTrans("combiShutter",-xt,-yt,zt,
                                      new TGeoRotation("rShutter",
                                                       rotate1b,rotate2,rotate3)) ); 
  fnodeNum++;

};
//***************************************************

void GRootDCNavigator::makeQuadArms() {

  TGeoVolume *fQuadArmVol[4];

  Double_t fQuadArmDiam = gDC->quadArmX;  // quad arm dimension
  Double_t offset = gDC->quadArmOffset;

  if (fDebugT) {
    *oLog << "*************** quad arms **********************" << endl;
    *oLog << "    quad arm diameter: "
	  << fQuadArmDiam << "  " << offset << endl;
  }

  for (int i = 0;i<4;i++) {
    fQuadArmR2R1V[i] = new TVector3();
    fQuadArmPosV[i]  = new TVector3();
  }
  // locate top of arm, telescope coor.
  //attach at corner of focus box
  //  Double_t diagh = fFocusBoxDim[0];
  Double_t diagh = 0.47;
  // ============== QuadArm0
  TVector3 r2(-diagh,diagh,fFL+0.4);
 
  // locate bottom of arm, telescope coor
  Double_t xArmM = gDC->quadArmBottomX[0]; 
  Double_t yArmM = gDC->quadArmBottomY[0];
  Double_t zArmM = fFL - sqrt(fFL*fFL - xArmM*xArmM - yArmM*yArmM);

  TVector3 r1(xArmM,yArmM,zArmM);

  // create vector parallel to quad arm, get magnitude and direction
  TVector3 R = r2 - r1;
  Double_t Mag = R.Mag();
  TVector3 Unit = R.Unit();
  *fQuadArmR2R1V[0] = R;

  // find location of quad arm in telescope coordinates
  TVector3 Rloc = r1 + r2; 
  Double_t tmp1 = (Rloc.Mag() )/2.0;
  Rloc.SetMag(tmp1);

  // find location of quad arm in TOP coordinates
  TVector3 RlocTop = Rloc - (*fTopPosV);
  Double_t xarmT = RlocTop(0);
  Double_t yarmT = RlocTop(1);
  Double_t zarmT = RlocTop(2);

  *fQuadArmPosV[0] = RlocTop;

  // set up rotation for quadarm1
  Double_t theta = R.Theta()*(TMath::RadToDeg());
  Double_t phi   = R.Phi()*(TMath::RadToDeg());

  Double_t euler0 = -(90.0 - phi);
  Double_t euler1 = -theta;

  Double_t rsize = fQuadArmDiam/2;
  Double_t zsize = Mag/2.;

  if (fDebugT) {
    *oLog << "      Quad Arm 1 " << endl;
    *oLog << "   top and bottom of arms in telescope coordinates" << endl;
    *oLog << "   Top of arm 1" << endl;
    r2.Print();
    *oLog << "   Bottom of arm 1" << endl;
    r1.Print();
    *oLog << "   length of quad arm: " << 2.*zsize << endl;
    *oLog << "   translate (TopVol coor): " << xarmT << " " << yarmT 
         << " " << zarmT << endl;
    *oLog << "   rotation: " << euler0 << " " << euler1 << "  0.0" 
         << endl<< endl;
  }

  fQuadArmVol[0] = fGeom->MakeTube("fQuadArm0Vol",fAl,0.,rsize,zsize);
  fTopVol->AddNode(fQuadArmVol[0],fnodeNum,
                   new TGeoCombiTrans("combiQArm0",xarmT,yarmT,zarmT,
                                      new TGeoRotation("rArm0",
                                                       euler0,euler1,0.0) ) );
  fnodeNum++;

  
  // ============== QuadArm1
  // locate top of arm, telescope coor.
   TVector3 r21(diagh,-diagh,fFL+0.4);
  // locate bottom of arm, telescope coor
  Double_t xArmM1 = gDC->quadArmBottomX[1]; 
  Double_t yArmM1 = gDC->quadArmBottomY[1];
  Double_t zArmM1 = fFL - sqrt(fFL*fFL - xArmM1*xArmM1 - yArmM1*yArmM1);

  TVector3 r11(xArmM1,yArmM1,zArmM1);

  // create vector parallel to quad arm, get magnitude and direction
  TVector3 R1 = r21 - r11;
  Double_t Mag1 = R1.Mag();
  TVector3 Unit1 = R1.Unit();
  *fQuadArmR2R1V[1] = R1;

  // find location of quad arm in telescope coordinates
  TVector3 Rloc1 = r11 + r21; 
  Double_t tmp11 = (Rloc1.Mag() )/2.0;
  Rloc1.SetMag(tmp11);

  // find location of quad arm in TOP coordinates
  TVector3 RlocTop1 = Rloc1 - (*fTopPosV);
  Double_t xarmT1 = RlocTop1(0);
  Double_t yarmT1 = RlocTop1(1);
  Double_t zarmT1 = RlocTop1(2);
  *fQuadArmPosV[1] = RlocTop1;

  // set up rotation for quadarm1
  Double_t theta1 = R1.Theta()*(TMath::RadToDeg());
  Double_t phi1   = R1.Phi()*(TMath::RadToDeg());

  Double_t euler01 = -(90.0 - phi1);
  Double_t euler11 = -theta1;
  Double_t rsize1 = fQuadArmDiam/2.;
  Double_t zsize1 = Mag1/2.;

  if (fDebugT) {
    *oLog << "      Quad Arm 2 " << endl;
    *oLog << "   top and bottom of arms in telescope coordinates" << endl;
    *oLog << "   Top of arm 2" << endl;
    r21.Print();
    *oLog << "   Bottom of arm 2" << endl;
    r11.Print();
    *oLog << "   length of quad arm: " << 2.*zsize1 << endl;
    *oLog << "   translate (TopVol coor): " << xarmT1 << " " << yarmT1 
         << " " << zarmT1 << endl;
    *oLog << "   rotation: " << euler01 << " " << euler11 << "  0.0" 
         << endl<< endl;
  }

  fQuadArmVol[1] = fGeom->MakeTube("fQuadArm1Vol",fAl,0.,rsize1,zsize1);
  fTopVol->AddNode(fQuadArmVol[1],fnodeNum,
                   new TGeoCombiTrans("combiQArm1",xarmT1,yarmT1,zarmT1,
                                      new TGeoRotation("rArm1",
                                                       euler01,euler11,0.0)) );
  fnodeNum++;
  
    
  // ============== QuadArm2
  // locate top of arm, telescope coor.
   TVector3 r22(diagh,diagh,fFL+0.4);
  // locate bottom of arm, telescope coor
  Double_t xArmM2 = gDC->quadArmBottomX[2]; 
  Double_t yArmM2 = gDC->quadArmBottomY[2];
  Double_t zArmM2 = fFL - sqrt(fFL*fFL - xArmM2*xArmM2 - yArmM2*yArmM2);
  TVector3 r12(xArmM2,yArmM2,zArmM2);

  // create vector parallel to quad arm, get magnitude and direction
  TVector3 R2 = r22 - r12;
  Double_t Mag2 = R2.Mag();
  TVector3 Unit2 = R2.Unit();
  *fQuadArmR2R1V[2] = R2;

  // find location of quad arm in telescope coordinates
  TVector3 Rloc2 = r12 + r22; 
  Double_t tmp12 = (Rloc2.Mag() )/2.0;
  Rloc2.SetMag(tmp12);

  // find location of quad arm in TOP coordinates
  TVector3 RlocTop2 = Rloc2 - (*fTopPosV);
  Double_t xarmT2 = RlocTop2(0);
  Double_t yarmT2 = RlocTop2(1);
  Double_t zarmT2 = RlocTop2(2);
  *fQuadArmPosV[2] = RlocTop2;

  // set up rotation for quadarm2
  Double_t theta2 = R2.Theta()*(TMath::RadToDeg());
  Double_t phi2   = R2.Phi()*(TMath::RadToDeg());
  
  Double_t euler02 = -(90.0 - phi2);
  Double_t euler12 = -theta2;

  Double_t rsize2 = fQuadArmDiam/2.;
  Double_t zsize2 = Mag2/2.;

  if (fDebugT) {
    *oLog << "      Quad Arm 3 " << endl;
    *oLog << "   top and bottom of arms in telescope coordinates" << endl;
    *oLog << "   Top of arm 3" << endl;
    r22.Print();
    *oLog << "   Bottom of arm 3" << endl;
    r12.Print();
    *oLog << "   length of quad arm: " << 2.*zsize2 << endl;
    *oLog << "   translate (TopVol coor): " << xarmT2 << " " << yarmT2 
         << " " << zarmT2 << endl;
    *oLog << "   rotation: " << euler02 << " " << euler12 << "  0.0" 
         << endl<< endl;
  }

  fQuadArmVol[2] = fGeom->MakeTube("fQuadArm2Vol",fAl,0.0,rsize2,zsize2);
  fTopVol->AddNode(fQuadArmVol[2],fnodeNum,
                   new TGeoCombiTrans("combiQArm2",xarmT2,yarmT2,zarmT2,
                                      new TGeoRotation("rArm2",
                                                       euler02,euler12,0.0)) );
  fnodeNum++;
  
  
  // ============== QuadArm3
  // locate top of arm, telescope coor.
 
  TVector3 r23(-diagh,-diagh,fFL+0.4);
  // locate bottom of arm, telescope coor
  Double_t xArmM3 = gDC->quadArmBottomX[3]; 
  Double_t yArmM3 = gDC->quadArmBottomY[3];
  Double_t zArmM3 = fFL - sqrt(fFL*fFL - xArmM3*xArmM3 - yArmM3*yArmM3);

  TVector3 r13(xArmM3,yArmM3,zArmM3);

  // create vector parallel to quad arm, get magnitude and direction
  TVector3 R3 = r23 - r13;
  Double_t Mag3 = R3.Mag();
  TVector3 Unit3 = R3.Unit();
  *fQuadArmR2R1V[3] = R3;

  // find location of quad arm in telescope coordinates
  TVector3 Rloc3 = r13 + r23; 
  Double_t tmp13 = (Rloc3.Mag() )/2.0;
  Rloc3.SetMag(tmp13);

  // find location of quad arm in TOP coordinates
  TVector3 RlocTop3 = Rloc3 - (*fTopPosV);
  Double_t xarmT3 = RlocTop3(0);
  Double_t yarmT3 = RlocTop3(1);
  Double_t zarmT3 = RlocTop3(2);
  *fQuadArmPosV[3] = RlocTop3;

  // set up rotation for quadarm3
  Double_t theta3 = R3.Theta()*(TMath::RadToDeg());
  Double_t phi3   = R3.Phi()*(TMath::RadToDeg());
  
  Double_t euler03 = -(90.0 - phi3);
  Double_t euler13 = -theta3;

  Double_t rsize3 = fQuadArmDiam/2.;
  Double_t zsize3 = Mag3/2.;
 
  if (fDebugT) {
    *oLog << "      Quad Arm 4 " << endl;
    *oLog << "   top and bottom of arms in telescope coordinates" << endl;
    *oLog << "   Top of arm 4" << endl;
    r23.Print();
    *oLog << "   Bottom of arm 4" << endl;
    r13.Print();
    *oLog << "   length of quad arm: " << 2.*zsize3 << endl;
    *oLog << "   translate (TopVol coor): " << xarmT3 << " " << yarmT3 
         << " " << zarmT3 << endl;
    *oLog << "   rotation: " << euler03 << " " << euler13 << "  0.0" 
         << endl<< endl;
  }

  fQuadArmVol[3] = fGeom->MakeTube("fQuadArm3Vol",fAl,0.0,rsize3,zsize3);
  fTopVol->AddNode(fQuadArmVol[3],fnodeNum,
                   new TGeoCombiTrans("combiQArm3",xarmT3,yarmT3,zarmT3,
                                      new TGeoRotation("rArm3",
                                                       euler03,euler13,0.0)) );
  fnodeNum++;
  
};
//********************************************************

void GRootDCNavigator::makeFacets() {
  // all vectors in telescope coordinates (origin at base)
  int ifirst = 0;  // flag for first time through facet loop
  TGeoVolume *facetVol = 0;

  TVector3 fpV(0.0,0.0,2*fFL);  // vector to 2FL point

  //open facet file
  ifstream infile("facet_loc_blrad",ios::in);
  if (! infile) {
    cerr << "    --  GRootDCNavigator::makeFacets " << endl;
    cerr << "      could not open file: " << "facet_loc_blrad" << endl;
    exit(0);
  }

  int facNum = 0;
  Double_t facRad = 0.0;
  Double_t fFX = 0.0, fFY=0.0, fFZ=0.0, fUse=0.0;
  
  while (infile >> facNum >> facRad >> fFX >> fFY >> fUse ) {

    if (fUse > 0.0) {

      // create a facet volume if this is first time thru loop
      // replicate this volume after the first time.
      if (ifirst == 0) {
        ifirst = 1;
        Double_t phi = 90.0;
        Double_t dphi = 360.0;
        Int_t nedges = 6;
        Int_t nz = 2;
        facetVol = fGeom->MakePgon("facet0",fAl,phi,dphi,
                                   nedges,nz);
        facetVol->SetLineColor(kBlue);
        Double_t zSect1 = -.1;
        Double_t zSect2 =  0.0;
        Double_t rmin = facRad - 0.001; 
        Double_t rmax = facRad;
        
        TGeoPgon *pgon = (TGeoPgon *)(facetVol->GetShape());
        pgon->DefineSection(0,zSect1,rmin,rmax);
        pgon->DefineSection(1,zSect2,rmin,rmax);
      }      
      
      fFZ = fFL - sqrt( (fFL*fFL) - (fFX*fFX) - (fFY*fFY) );

      TVector3 fLocV(fFX,fFY,fFZ);  // vector to facet center

     // vector from facet location to 2*fFL point
      TVector3 fcV = fpV - fLocV;
  
      Double_t thetaf = fcV.Theta()* (TMath::RadToDeg());
      Double_t phif   = fcV.Phi()* (TMath::RadToDeg());
      
      // Euler angle rotations 1,2,3
      Double_t r1 = -(90.0 - phif);
      Double_t r2 = - thetaf;
      Double_t r3 = - r1;

      r1 = r1 - 90;
      
      // get translation components, facet loc wrt top center
      TVector3 fLocFacTopV = fLocV - (*fTopPosV);
      Double_t xtop = fLocFacTopV[0];
      Double_t ytop = fLocFacTopV[1];
      Double_t ztop = fLocFacTopV[2];

      fTopVol->AddNode(facetVol,fnodeNum,
                       new TGeoCombiTrans("combiFacet",
                                          xtop,ytop,ztop,
                                          new TGeoRotation("rFacet",
                                                           r1,0.,r3)) );
      fnodeNum++;
    }
  }
};
//************************************************

void GRootDCNavigator::makeCrossArms() {

  //fQuadArmR2R1V[3]
  //fQuadArmPosV[4]; // position of quadarms, TOP coordinates
  //fQuadArmR2R1V[4]; // r2-r1 vector for each quad arm

  /*  cross arms attach to the quad arms at a specificied distance
      below the focus box (the arms attach to the mirror at slightly
      different radial distances, so start at the focus box).  
      Assume the length of the quad arm is 12 meters (close), and that 
      the cross arms attach 0.4 *12 meters along the quad arm below the
      focus box, or 4.8 meters below.  
      Set this up in terms of the quad arm length, the position of
      the quad arm in TOP coordinates, and find the position of the cross
      arms on each quad arm, etc. using these vectors.

      fQuadArmR2R1V[3]; 
      fQuadArmPosV[4];  position of quadarms, TOP coordinates

  */
  // cross bar dimensions, length comes later
  // use two arrays since lengths may be different
  Double_t fCrossBar1Dim[3];  // cross bar dimensions
  Double_t fCrossBar2Dim[3];
  fCrossBar1Dim[0] = gDC->crossBarX[0]; 
  fCrossBar2Dim[0] = gDC->crossBarX[1];   // meters 
  fCrossBar1Dim[1] = gDC->crossBarY[0];
  fCrossBar2Dim[1] = gDC->crossBarY[1];

  fCrossBar1Dim[2] = fCrossBar2Dim[2] = 0.0;  // calc later

  Double_t dcross = gDC->crossBarDistBelowFocBox[0];

  // calculate distance of cross arm attachment from fQuadArmPos
  Double_t dc[4];
  for (int i = 0;i<3;i++) {
    dc[i] = 0.0;
  }

  TVector3 crossEndV[4];  // location vector of end of xarm,1,1,2,2
  
  // cross arm end locations, first two indices, arm1; then arm2
  if (fDebugT) {
    *oLog << "************ crossarm construction ********* " << endl;
    *oLog << "    cross arm distances below focus box " 
	  << dcross << endl;
  }

  for (int i = 0;i<4;i++) {
    Double_t armLength = fQuadArmR2R1V[i]->Mag();
    dc[i] = (armLength/2.) - dcross;
    TVector3 dcV = fQuadArmR2R1V[i]->Unit();
    dcV.SetMag(dc[i]);
    crossEndV[i] = (*fQuadArmPosV[i]) + dcV;

    if (fDebugT) {
        *oLog << "       1/2 length of quad arm " << i+1 << " " 
             << armLength/2 << endl;
        *oLog << "       attach above midpoint of quad: " 
             << i+1 << " " << dc[i] << endl;
      
    }
  }

  // vector parallel to crossarm from end to end
  TVector3 R = crossEndV[0] - crossEndV[2];

  // find cross arm location and set up directions
  // *********ARM 1 NEXT.***************
  TVector3 RlocTop = crossEndV[2] - crossEndV[1];
  Double_t tmp1 = (RlocTop.Mag() )/2.0;
  RlocTop.SetMag(tmp1);

  Double_t xarmT = RlocTop(0);
  Double_t yarmT = RlocTop(1);
  Double_t zarmT = RlocTop(2);
 
  // set up rotation for quadarm1
  Double_t theta = RlocTop.Theta()*(TMath::RadToDeg());
  Double_t phi   = RlocTop.Phi()*(TMath::RadToDeg());
  Double_t psi   = 0.0;

  TGeoVolume *fCrossArmVol[2];

  Double_t rsize = fCrossBar1Dim[0]/2.0;
  Double_t zsize = R.Mag()/2.;

  if (fDebugT) {
    *oLog << "      cross arm 1 dimensions: " << 2.*rsize << " " 
         << " " << 2.*zsize << endl;
    *oLog << "      cross arm 1 translation: " << xarmT << " " 
         << yarmT << " " << zarmT << endl;
    *oLog << "      cross arm 1 rotation: " << theta << "  " << phi << "  " << psi << endl;
  }

  fCrossArmVol[0] = fGeom->MakeTube("fCrossArm0Vol",fAl,
                                   0.,rsize,zsize);
 
  fTopVol->AddNode(fCrossArmVol[0],fnodeNum,
                   new TGeoCombiTrans("combiQArm0",xarmT,yarmT,zarmT,
                                      new TGeoRotation("rArm0",
                                                       theta,phi,psi)) );

  fnodeNum++;
  
  //****************** SECOND ARM ***********
  // vector parallel to crossarm from end to end
  R = crossEndV[1] - crossEndV[3];

  RlocTop = crossEndV[3] - crossEndV[0];
  tmp1 = (RlocTop.Mag() )/2.0;
  RlocTop.SetMag(tmp1);

  xarmT = RlocTop(0);
  yarmT = RlocTop(1);
  zarmT = RlocTop(2);
 
  // set up rotation for quadarm1
  theta = RlocTop.Theta()*(TMath::RadToDeg());
  phi   = RlocTop.Phi()*(TMath::RadToDeg());
  psi   = 0.0;

  rsize = fCrossBar2Dim[0]/2.0;
  zsize = R.Mag()/2.;
  
  if (fDebugT) {
    *oLog << "      cross arm 2 dimensions: " << 2.*rsize << " " 
         << " " << 2.*zsize << endl;
    *oLog << "      cross arm 2 translation: " << xarmT << " " 
         << yarmT << " " << zarmT << endl;
    *oLog << "      cross arm 2 rotation: " << theta << "  " << phi << "  " << psi << endl;
    *oLog << " *********** end of cross arm construction *********" << endl;
  }

  fCrossArmVol[1] = fGeom->MakeTube("fCrossArm0Vol",fAl,0.,rsize,zsize);

  fTopVol->AddNode(fCrossArmVol[1],fnodeNum,
                   new TGeoCombiTrans("combiQArm0",xarmT,yarmT,zarmT,
                                      new TGeoRotation("rArm0",
                                                       theta,phi,psi) ) );

  fnodeNum++;

  //****************** THIRD ARM ***********
  // vector parallel to crossarm from end to end
  R = crossEndV[2] - crossEndV[1];

  RlocTop = crossEndV[0] - crossEndV[2];
  tmp1 = (RlocTop.Mag() )/2.0;
  RlocTop.SetMag(tmp1);

  xarmT = RlocTop(0);
  yarmT = RlocTop(1);
  zarmT = RlocTop(2);
 
  // set up rotation for quadarm1
  theta = RlocTop.Theta()*(TMath::RadToDeg()) + 90;
  phi   = RlocTop.Phi()*(TMath::RadToDeg()) + 90;
  psi   = 0.0;

  rsize = fCrossBar2Dim[0]/2.0;
  zsize = R.Mag()/2.;

  if (fDebugT) {
    *oLog << "      cross arm 3 dimensions: " << 2.*rsize << " " 
         << " " << 2.*zsize << endl;
    *oLog << "      cross arm 3 translation: " << xarmT << " " 
         << yarmT << " " << zarmT << endl;
    *oLog << "      cross arm 3 rotation: " << theta << "  " << phi << "  " << psi << endl;
    *oLog << " *********** end of cross arm construction *********" << endl;
  }

  fCrossArmVol[2] = fGeom->MakeTube("fCrossArm1Vol",fAl,0.,rsize,zsize);
  
  fTopVol->AddNode(fCrossArmVol[2],fnodeNum,
                   new TGeoCombiTrans("combiQArm1",xarmT,yarmT,zarmT,
                                      new TGeoRotation("rArm1",
                                                       theta,phi,psi) ) );
  fnodeNum++;

  //****************** FOURTH ARM ***********
  // vector parallel to crossarm from end to end
  R = crossEndV[3] - crossEndV[0];

  RlocTop = crossEndV[1] - crossEndV[3];
  tmp1 = (RlocTop.Mag() )/2.0;
  RlocTop.SetMag(tmp1);

  xarmT = RlocTop(0);
  yarmT = RlocTop(1);
  zarmT = RlocTop(2);
 
  // set up rotation for quadarm1
  theta = RlocTop.Theta()*(TMath::RadToDeg()) + 90;
  phi   = RlocTop.Phi()*(TMath::RadToDeg()) + 90;
  psi   = 0.0;

  rsize = fCrossBar2Dim[0]/2.0;
  zsize = R.Mag()/2.;

  if (fDebugT) {
    *oLog << "      cross arm 4 dimensions: " << 2.*rsize << " " 
         << " " << 2.*zsize << endl;
    *oLog << "      cross arm 4 translation: " << xarmT << " " 
         << yarmT << " " << zarmT << endl;
    *oLog << "      cross arm 4 rotation: " << theta << "  " << phi << "  " << psi << endl;
    *oLog << " *********** end of cross arm construction *********" << endl;
  }

  fCrossArmVol[3] = fGeom->MakeTube("fCrossArm1Vol",fAl,0.,rsize,zsize);
  
  fTopVol->AddNode(fCrossArmVol[3],fnodeNum,
                   new TGeoCombiTrans("combiQArm1",xarmT,yarmT,zarmT,
                                      new TGeoRotation("rArm1",
                                                       theta,phi,psi) ) );
  fnodeNum++;

}; //end make cross arms
//****************************************************

void GRootDCNavigator::makeSupportWires() {
  //quad arm bottoms {-0.95 1.8 -2.6} {0.95 -1.8 -2.6} {0.95 1.8 -2.6} {-0.95 -1.8 -2.6}
  //{-0.95 1.8 -2.6}->{-0.7 -1.15 0}
 
  TGeoVolume *fSupportWireVol[16];

  Double_t rsize = 0.008/2.;

  /////////////////////////////////////////////////////
  // First wire!
  /////////////////////////////////////////////////////
  Double_t xb = 0.9;
  Double_t yb = 1.8;
  Double_t xm1 = 0.7;
  Double_t ym1 = 1.15;
  TVector3 r1(-xb, yb, -2.6);
  TVector3 r2(-xm1, -ym1, 0);

  TVector3 Rloc = r1 + r2; 
  Double_t tmp1 = (Rloc.Mag() )/2.0;
  Rloc.SetMag(tmp1);
  Double_t xarmT = Rloc(0);
  Double_t yarmT = Rloc(1);
  Double_t zarmT = Rloc(2);

  // // set up rotation for wire
  TVector3 R = r1 - r2;
  Double_t Mag = R.Mag();
  Double_t theta = R.Theta()*(TMath::RadToDeg());
  Double_t phi   = R.Phi()*(TMath::RadToDeg());

  Double_t psi = 0;
  Double_t euler0 = phi - 90;
  Double_t euler1 = 180 -theta;

  Double_t zsize = Mag/2.;

  fSupportWireVol[0] = fGeom->MakeTube("fSupportWire0Vol",fAl,
				       0.,rsize,zsize);
 
  fTopVol->AddNode(fSupportWireVol[0],fnodeNum,
		   new TGeoCombiTrans("combiSWire0",xarmT,yarmT,zarmT,
				      new TGeoRotation("rArm0",
						       euler0,euler1,psi)) );
  fnodeNum++;

  /////////////////////////////////////////////////////
  // Second wire!
  /////////////////////////////////////////////////////
  TVector3 r11(-xb, -yb, -2.6);
  TVector3 r21(-xm1, ym1, 0);

  TVector3 Rloc1 = r11 + r21; 
  Double_t tmp2 = (Rloc1.Mag() )/2.0;
  Rloc1.SetMag(tmp2);
  Double_t xarmT1 = Rloc1(0);
  Double_t yarmT1 = Rloc1(1);
  Double_t zarmT1 = Rloc1(2);

  // // set up rotation for wire
  TVector3 R1 = r11 - r21;
  Double_t theta1 = R1.Theta()*(TMath::RadToDeg());
  Double_t phi1   = R1.Phi()*(TMath::RadToDeg());

  Double_t euler01 = phi1 - 90;
  Double_t euler11 = 180 -theta1;

  fSupportWireVol[1] = fGeom->MakeTube("fSupportWire1Vol",fAl,
				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[1],fnodeNum,
		   new TGeoCombiTrans("combiSWire1",xarmT1,yarmT1,zarmT1,
				      new TGeoRotation("rArm1",
						       euler01,euler11,psi)) );

  fnodeNum++;

  /////////////////////////////////////////////////////
  // Third wire!
  /////////////////////////////////////////////////////
  TVector3 r12(xb, yb, -2.6);
  TVector3 r22(xm1, -ym1, 0);

  TVector3 Rloc2 = r12 + r22; 
  Rloc2.SetMag((Rloc2.Mag() )/2.0);
  Double_t xarmT2 = Rloc2(0);
  Double_t yarmT2 = Rloc2(1);
  Double_t zarmT2 = Rloc2(2);

  // // set up rotation for wire
  TVector3 R2 = r12 - r22;
  Double_t theta2 = R2.Theta()*(TMath::RadToDeg());
  Double_t phi2   = R2.Phi()*(TMath::RadToDeg());

  Double_t euler02 = phi2 - 90;
  Double_t euler12 = 180 -theta2;

  fSupportWireVol[2] = fGeom->MakeTube("fSupportWire2Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[2],fnodeNum,
		   new TGeoCombiTrans("combiSWire2",xarmT2,yarmT2,zarmT2,
				      new TGeoRotation("rArm2",
						       euler02,euler12,psi)) );

  fnodeNum++;

  /////////////////////////////////////////////////////
  // 4th wire!
  /////////////////////////////////////////////////////
  TVector3 r13(xb, -yb, -2.6);
  TVector3 r23(xm1, ym1, 0);

  TVector3 Rloc3 = r13 + r23; 
  Rloc3.SetMag((Rloc3.Mag() )/2.0);
  Double_t xarmT3 = Rloc3(0);
  Double_t yarmT3 = Rloc3(1);
  Double_t zarmT3 = Rloc3(2);

  // // set up rotation for wire
  TVector3 R3 = r13 - r23;
  Double_t theta3 = R3.Theta()*(TMath::RadToDeg());
  Double_t phi3   = R3.Phi()*(TMath::RadToDeg());

  Double_t euler03 = phi3 - 90;
  Double_t euler13 = 180 -theta3;

  fSupportWireVol[3] = fGeom->MakeTube("fSupportWire3Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[3],fnodeNum,
		   new TGeoCombiTrans("combiSWire3",xarmT3,yarmT3,zarmT3,
				      new TGeoRotation("rArm3",
						       euler03,euler13,psi)) );

  fnodeNum++;

  /////////////////////////////////////////////////////
  // 5th wire!
  /////////////////////////////////////////////////////
  TVector3 r14(xb, -yb, -2.6);
  TVector3 r24(-xm1, -ym1, 0);

  TVector3 Rloc4 = r14 + r24; 
  Rloc4.SetMag((Rloc4.Mag() )/2.0);
  Double_t xarmT4 = Rloc4(0);
  Double_t yarmT4 = Rloc4(1);
  Double_t zarmT4 = Rloc4(2);

  // // set up rotation for wire
  TVector3 R4 = r14 - r24;
  //recalculate the wire length
  zsize = R4.Mag()/2.;
  Double_t theta4 = R4.Theta()*(TMath::RadToDeg());
  Double_t phi4   = R4.Phi()*(TMath::RadToDeg());

  Double_t euler04 = phi4 - 90;
  Double_t euler14 = 180 -theta4;

  fSupportWireVol[4] = fGeom->MakeTube("fSupportWire4Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[4],fnodeNum,
		   new TGeoCombiTrans("combiSWire4",xarmT4,yarmT4,zarmT4,
				      new TGeoRotation("rArm4",
						       euler04,euler14,psi)) );

  fnodeNum++;

  /////////////////////////////////////////////////////
  // 6th wire!
  /////////////////////////////////////////////////////
  TVector3 r15(-xb, -yb, -2.6);
  TVector3 r25(xm1, -ym1, 0);

  TVector3 Rloc5 = r15 + r25; 
  Rloc5.SetMag((Rloc5.Mag() )/2.0);
  Double_t xarmT5 = Rloc5(0);
  Double_t yarmT5 = Rloc5(1);
  Double_t zarmT5 = Rloc5(2);

  // // set up rotation for wire
  TVector3 R5 = r15 - r25;
  Double_t theta5 = R5.Theta()*(TMath::RadToDeg());
  Double_t phi5   = R5.Phi()*(TMath::RadToDeg());

  Double_t euler05 = phi5 - 90;
  Double_t euler15 = 180 -theta5;

  fSupportWireVol[5] = fGeom->MakeTube("fSupportWire5Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[5],fnodeNum,
		   new TGeoCombiTrans("combiSWire5",xarmT5,yarmT5,zarmT5,
				      new TGeoRotation("rArm5",
						       euler05,euler15,psi)) );

  fnodeNum++;

  /////////////////////////////////////////////////////
  // 7th wire!
  /////////////////////////////////////////////////////
  TVector3 r16(xb, yb, -2.6);
  TVector3 r26(-xm1, ym1, 0);

  TVector3 Rloc6 = r16 + r26; 
  Rloc6.SetMag((Rloc6.Mag() )/2.0);
  Double_t xarmT6 = Rloc6(0);
  Double_t yarmT6 = Rloc6(1);
  Double_t zarmT6 = Rloc6(2);

  // // set up rotation for wire
  TVector3 R6 = r16 - r26;
  Double_t theta6 = R6.Theta()*(TMath::RadToDeg());
  Double_t phi6   = R6.Phi()*(TMath::RadToDeg());

  Double_t euler06 = phi6 - 90;
  Double_t euler16 = 180 -theta6;

  fSupportWireVol[6] = fGeom->MakeTube("fSupportWire6Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[6],fnodeNum,
		   new TGeoCombiTrans("combiSWire6",xarmT6,yarmT6,zarmT6,
				      new TGeoRotation("rArm6",
						       euler06,euler16,psi)) );

  fnodeNum++;

  /////////////////////////////////////////////////////
  // 8th wire!
  /////////////////////////////////////////////////////
  TVector3 r17(-xb, yb, -2.6);
  TVector3 r27(xm1, ym1, 0);

  TVector3 Rloc7 = r17 + r27; 
  Rloc7.SetMag((Rloc7.Mag() )/2.0);
  Double_t xarmT7 = Rloc7(0);
  Double_t yarmT7 = Rloc7(1);
  Double_t zarmT7 = Rloc7(2);

  // // set up rotation for wire
  TVector3 R7 = r17 - r27;
  Double_t theta7 = R7.Theta()*(TMath::RadToDeg());
  Double_t phi7   = R7.Phi()*(TMath::RadToDeg());

  Double_t euler07 = phi7 - 90;
  Double_t euler17 = 180 -theta7;

  fSupportWireVol[7] = fGeom->MakeTube("fSupportWire7Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[7],fnodeNum,
		   new TGeoCombiTrans("combiSWire7",xarmT7,yarmT7,zarmT7,
				      new TGeoRotation("rArm7",
						       euler07,euler17,psi)) );

  fnodeNum++;

  //now the upper wires
  Double_t xm2 = 0.7;
  Double_t ym2 = 1.15;
  Double_t xt = 0.46;
  Double_t yt = 0.55;
  /////////////////////////////////////////////////////
  // 9th wire!
  /////////////////////////////////////////////////////
  TVector3 r19(-xt, yt, 2.6);
  TVector3 r29(xm2, ym2, 0);

  TVector3 Rloc9 = r19 + r29; 
  Rloc9.SetMag((Rloc9.Mag() )/2.0);
  Double_t xarmT9 = Rloc9(0);
  Double_t yarmT9 = Rloc9(1);
  Double_t zarmT9 = Rloc9(2);

  // // set up rotation for wire
  TVector3 R9 = r19 - r29;
  //recalc z size
  zsize = R9.Mag()/2.;
  Double_t theta9 = R9.Theta()*(TMath::RadToDeg());
  Double_t phi9   = R9.Phi()*(TMath::RadToDeg());

  Double_t euler09 = phi9 - 90;
  Double_t euler19 = 180 -theta9;

  fSupportWireVol[9] = fGeom->MakeTube("fSupportWire9Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[9],fnodeNum,
		   new TGeoCombiTrans("combiSWire9",xarmT9,yarmT9,zarmT9,
				      new TGeoRotation("rArm9",
						       euler09,euler19,psi)) );

  fnodeNum++;

  /////////////////////////////////////////////////////
  // 10th wire!
  /////////////////////////////////////////////////////
  TVector3 r110(xt, yt, 2.6);
  TVector3 r210(-xm2, ym2, 0);

  TVector3 Rloc10 = r110 + r210; 
  Rloc10.SetMag((Rloc10.Mag() )/2.0);
  Double_t xarmT10 = Rloc10(0);
  Double_t yarmT10 = Rloc10(1);
  Double_t zarmT10 = Rloc10(2);

  // // set up rotation for wire
  TVector3 R10 = r110 - r210;
  Double_t theta10 = R10.Theta()*(TMath::RadToDeg());
  Double_t phi10   = R10.Phi()*(TMath::RadToDeg());

  Double_t euler010 = phi10 - 90;
  Double_t euler110 = 180 -theta10;

  fSupportWireVol[10] = fGeom->MakeTube("fSupportWire10Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[10],fnodeNum,
		   new TGeoCombiTrans("combiSWire10",xarmT10,yarmT10,zarmT10,
				      new TGeoRotation("rArm10",
						       euler010,euler110,psi)) );

  fnodeNum++;
  /////////////////////////////////////////////////////
  // 11th wire!
  /////////////////////////////////////////////////////
  TVector3 r111(-xt, -yt, 2.6);
  TVector3 r211(xm2, -ym2, 0);

  TVector3 Rloc11 = r111 + r211; 
  Rloc11.SetMag((Rloc11.Mag() )/2.0);
  Double_t xarmT11 = Rloc11(0);
  Double_t yarmT11 = Rloc11(1);
  Double_t zarmT11 = Rloc11(2);

  // // set up rotation for wire
  TVector3 R11 = r111 - r211;
  //recalc z size
  zsize = R11.Mag()/2.;
  Double_t theta11 = R11.Theta()*(TMath::RadToDeg());
  Double_t phi11   = R11.Phi()*(TMath::RadToDeg());

  Double_t euler011 = phi11 - 90;
  Double_t euler111 = 180 -theta11;

  fSupportWireVol[11] = fGeom->MakeTube("fSupportWire11Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[11],fnodeNum,
		   new TGeoCombiTrans("combiSWire11",xarmT11,yarmT11,zarmT11,
				      new TGeoRotation("rArm11",
						       euler011,euler111,psi)) );

  fnodeNum++;

  /////////////////////////////////////////////////////
  // 12th wire!
  /////////////////////////////////////////////////////
  TVector3 r112(xt, -yt, 2.6);
  TVector3 r212(-xm2, -ym2, 0);

  TVector3 Rloc12 = r112 + r212; 
  Rloc12.SetMag((Rloc12.Mag() )/2.0);
  Double_t xarmT12 = Rloc12(0);
  Double_t yarmT12 = Rloc12(1);
  Double_t zarmT12 = Rloc12(2);

  // // set up rotation for wire
  TVector3 R12 = r112 - r212;
  Double_t theta12 = R12.Theta()*(TMath::RadToDeg());
  Double_t phi12   = R12.Phi()*(TMath::RadToDeg());

  Double_t euler012 = phi12 - 90;
  Double_t euler112 = 180 -theta12;

  fSupportWireVol[12] = fGeom->MakeTube("fSupportWire12Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[12],fnodeNum,
		   new TGeoCombiTrans("combiSWire12",xarmT12,yarmT12,zarmT12,
				      new TGeoRotation("rArm12",
						       euler012,euler112,psi)) );

  fnodeNum++;

  /////////////////////////////////////////////////////
  // 13th wire!
  /////////////////////////////////////////////////////
  TVector3 r113(xt, -yt, 2.6);
  TVector3 r213(xm2, ym2, 0);

  TVector3 Rloc13 = r113 + r213; 
  Rloc13.SetMag((Rloc13.Mag() )/2.0);
  Double_t xarmT13 = Rloc13(0);
  Double_t yarmT13 = Rloc13(1);
  Double_t zarmT13 = Rloc13(2);

  // // set up rotation for wire
  TVector3 R13 = r113 - r213;
  //recalc z size
  zsize = R13.Mag()/2.;
  Double_t theta13 = R13.Theta()*(TMath::RadToDeg());
  Double_t phi13   = R13.Phi()*(TMath::RadToDeg());

  Double_t euler013 = phi13 - 90;
  Double_t euler113 = 180 -theta13;

  fSupportWireVol[13] = fGeom->MakeTube("fSupportWire13Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[13],fnodeNum,
		   new TGeoCombiTrans("combiSWire13",xarmT13,yarmT13,zarmT13,
				      new TGeoRotation("rArm13",
						       euler013,euler113,psi)) );

  fnodeNum++;

  /////////////////////////////////////////////////////
  // 14th wire!
  /////////////////////////////////////////////////////
  TVector3 r114(xt, yt, 2.6);
  TVector3 r214(xm2, -ym2, 0);

  TVector3 Rloc14 = r114 + r214; 
  Rloc14.SetMag((Rloc14.Mag() )/2.0);
  Double_t xarmT14 = Rloc14(0);
  Double_t yarmT14 = Rloc14(1);
  Double_t zarmT14 = Rloc14(2);

  // // set up rotation for wire
  TVector3 R14 = r114 - r214;
  Double_t theta14 = R14.Theta()*(TMath::RadToDeg());
  Double_t phi14   = R14.Phi()*(TMath::RadToDeg());

  Double_t euler014 = phi14 - 90;
  Double_t euler114 = 180 -theta14;

  fSupportWireVol[14] = fGeom->MakeTube("fSupportWire14Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[14],fnodeNum,
		   new TGeoCombiTrans("combiSWire14",xarmT14,yarmT14,zarmT14,
				      new TGeoRotation("rArm14",
						       euler014,euler114,psi)) );

  fnodeNum++;

  /////////////////////////////////////////////////////
  // 15th wire!
  /////////////////////////////////////////////////////
  TVector3 r115(-0.46, -0.55, 2.6);
  TVector3 r215(-0.7, 1.15, 0);

  TVector3 Rloc15 = r115 + r215; 
  Rloc15.SetMag((Rloc15.Mag() )/2.0);
  Double_t xarmT15 = Rloc15(0);
  Double_t yarmT15 = Rloc15(1);
  Double_t zarmT15 = Rloc15(2);

  // // set up rotation for wire
  TVector3 R15 = r115 - r215;
  Double_t theta15 = R15.Theta()*(TMath::RadToDeg());
  Double_t phi15   = R15.Phi()*(TMath::RadToDeg());

  Double_t euler015 = phi15 - 90;
  Double_t euler115 = 180 -theta15;

  fSupportWireVol[15] = fGeom->MakeTube("fSupportWire15Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[15],fnodeNum,
		   new TGeoCombiTrans("combiSWire15",xarmT15,yarmT15,zarmT15,
				      new TGeoRotation("rArm15",
						       euler015,euler115,psi)) );

  fnodeNum++;

  /////////////////////////////////////////////////////
  // 16th wire!
  /////////////////////////////////////////////////////
  TVector3 r116(-0.46, 0.55, 2.6);
  TVector3 r216(-0.7, -1.15, 0);

  TVector3 Rloc16 = r116 + r216; 
  Rloc16.SetMag((Rloc16.Mag() )/2.0);
  Double_t xarmT16 = Rloc16(0);
  Double_t yarmT16 = Rloc16(1);
  Double_t zarmT16 = Rloc16(2);

  // // set up rotation for wire
  TVector3 R16 = r116 - r216;
  Double_t theta16 = R16.Theta()*(TMath::RadToDeg());
  Double_t phi16   = R16.Phi()*(TMath::RadToDeg());

  Double_t euler016 = phi16 - 90;
  Double_t euler116 = 180 -theta16;

  fSupportWireVol[16] = fGeom->MakeTube("fSupportWire16Vol",fAl,
   				       0.,rsize,zsize);

  fTopVol->AddNode(fSupportWireVol[16],fnodeNum,
		   new TGeoCombiTrans("combiSWire16",xarmT16,yarmT16,zarmT16,
				      new TGeoRotation("rArm16",
						       euler016,euler116,psi)) );

  fnodeNum++;
}; //end make support wires
//****************************************************

void GRootDCNavigator::movePositionToTopOfTopVol() {

  if (fDebugTr) {
    *oLog << "  -- GRootDCNavigator::movePositionToTopOfTopVol " << endl;
    *oLog << "        position prior to move to top ";
    *oLog << fPosC[0] << "  " << fPosC[1] << "  " << fPosC[2] << endl;
  }

  double rfx = fPosC[0];
  double rfy = fPosC[1];
  double rfz = fPosC[2];

  double Z = fepsil + fFocusBoxDim[2];

  double dl = fDirC[0];
  double dm = fDirC[1];
  double dn = fDirC[2];

  double Rx = rfx - (rfz - Z)*dl/dn;
  double Ry = rfy - (rfz - Z)*dm/dn;
  double Rz = Z;

  fPosC[0] = Rx;
  fPosC[1] = Ry;
  fPosC[2] = Rz;

  if (fDebugTr) {
    *oLog << "        TopVolPos in focal point coor.  ";
    for (int i = 0;i<3;i++) {
      *oLog << fPosC[i] << " ";
    }
    *oLog << endl;
  }
  
};
//****************************************************

const char * GRootDCNavigator::setPositionDirection( double *pos, double *dir) {

  gGeoManager = fGeom;
  // position and direction: in telescope coordinates
  // translate to TOP coordinates
  for (int i = 0;i< 3;i++) {
    fPosC[i] = pos[i];
    fDirC[i] = dir[i];
  }
  if (fDebugTr) {
    *oLog << "  -- GRootDCNavigator::setPositionDirection " << endl;
    *oLog << "       initial position, focal point coor: " << pos[0] 
          << " " << pos[1] << " " << pos[2] << endl;
    *oLog << "       initial direction: " << dir[0] 
          << " " << dir[1] << " " << dir[2] << endl;
  }
  // move position to intersection with topVol top surface
  if (dir[2] < 0) {
    fMoveToTop = true;
    //*oLog << "moving to top of topVol" << endl;
    movePositionToTopOfTopVol();
  } 

  // this is z component of position relative to the topVol Center
  fPosC[2] = fPosC[2] +(*fTopPosV)[2] - fFocusBoxDim[2] - fepsil; 

  fGeom->InitTrack(fPosC,fDirC);
  TGeoNode *cnode = fGeom->GetCurrentNode();
  fNodeC = cnode->GetName();

  if (fDebugTr) {
    *oLog << "           movePosition to top if dir[2] < 0 " << endl;
    *oLog << "       initial pos. TopVol coor " << fPosC[0] << " " 
          << fPosC[1] << " " << fPosC[2] << endl;

    //Double_t zOffset = (*fTopPosV)[2];
    Double_t zOffset = (*fTopPosV)[2] - fFocusBoxDim[2] - fepsil;

    *oLog << "       ftopPosV[2]: " << (*fTopPosV)[2] << endl;
    
    const Double_t *cpos = fGeom->GetCurrentPoint();
    *oLog << "       current position, tele. coor.: " << cpos[0] << " " 
         << cpos[1] << " " << cpos[2] - zOffset << endl;
    
    const Double_t *cdir = fGeom->GetCurrentDirection();
    *oLog << "       current direction: " << cdir[0] << " " 
         << cdir[1] << " " << cdir[2] << endl;
    *oLog << "       InitialNode: " << fNodeC << endl;
  }
  return fNodeC;
};

//*****************************************************
const char * GRootDCNavigator:: getNextNodeName() {

  gGeoManager = fGeom;

  TGeoNode *nextNode = fGeom->FindNextBoundary();
  //nextNode = fGeom->Step();
  fNodeN = nextNode->GetName();
  nextNode = fGeom->Step();

  if (fDebugTr) {
    //Double_t zOffset = (*fTopPosV)[2];
    Double_t zOffset = (*fTopPosV)[2] - fFocusBoxDim[2] - fepsil;
    //Double_t zOffset = 0.0;
    *oLog << "============ DebugT: First Boundary ===========" << endl;
    *oLog << "ftopPosV[2]: " << (*fTopPosV)[2] << endl;
    
    const Double_t *cpos = fGeom->GetCurrentPoint();
    *oLog << "current position, tele. coor.: " << cpos[0] << " " 
         << cpos[1] << " " << cpos[2] - zOffset << endl;
    
    const Double_t *cdir = fGeom->GetCurrentDirection();
    *oLog << "current direction: " << cdir[0] << " " 
         << cdir[1] << " " << cdir[2] << endl;
    *oLog << "InitialNode: " << fNodeN << endl << endl;
  }
  //exit(0);
  return fNodeN;
};

//************************************************************
const char * GRootDCNavigator::getNextNodeName(const double *pos, 
                                              const double *dir) {
  gGeoManager = fGeom;
  TGeoNode *nextNode = fGeom->FindNextBoundary();
  fNodeN = nextNode->GetName();
  fGeom->Step();

  TGeoNode *cnode = fGeom->GetCurrentNode();
  fNodeC = cnode->GetName();

  pos = fGeom->GetCurrentPoint();
  dir = fGeom->GetCurrentDirection();

  return fNodeC;
};

//****************************************************
const char * GRootDCNavigator::getPositionDirection(double *pos1, 
                                                    double *dir1) {
  gGeoManager = fGeom;

  const double *p,*d;
  
  p = fGeom->GetCurrentPoint();
  d = fGeom->GetCurrentDirection();
  for (int i = 0;i<3;i++) {
    pos1[i] = p[i];
    dir1[i] = d[i];
  }

  TGeoNode *cnode = fGeom->GetCurrentNode();
  const char * currnode = cnode->GetName();


  return currnode;
};

//****************************************************
void GRootDCNavigator::drawTelescope(const int &option) {
  gGeoManager = fGeom;
  gGeoManager->GetMasterVolume()->Draw("ogl");

  //   Draw("pad");        // use defaults in Draw if no openGL graphics
};

//****************************************************
void GRootDCNavigator::setTrackingDebug(bool setDebug) {
  fDebugTr = setDebug;

};
double GRootDCNavigator::getFocalBoxZBottomTopVolCoor() {
  double tmp = (*fTopPosV)[2] - fFocusBoxDim[2] - fepsil;
  return tmp;
};

//****************************************************
void GRootDCNavigator::printVariables() {

  *oLog << "============== printVariables ===========" << endl;
  for (int i = 0;i < 3;i++) {
    DEBUG(fPosC[i]);
  }
  *oLog << endl;
  for (int i = 0;i < 3;i++) {
    DEBUG(fDirC[i]);
  }
  *oLog << endl;

  DEBUG(fNodeC);

  for (int i = 0;i < 3;i++) {
    DEBUG(fPosN[i]);
  }
  *oLog << endl;

  for (int i = 0;i < 3;i++) {
    DEBUG(fDirN[i]);
  }
  *oLog << endl;

  DEBUG(fNodeN);
  *oLog << "============= finished printVariables ===" << endl;

};

