#include "ReadConfig.h"
#include <TMath.h>

ReadConfig::ReadConfig(TRandom3 *random)
{

  rand = random;

  fDebug = kTRUE;
  bUseNSB=0;                      //set to true if we want to use NSB in the simulation
  iTelescopeMultiplicity=-1;      //How many telescopes need to be in a cluster for a trigger
  bArrayTriggerRequiresNextNeighbor=0;
  fArrayCoincidence=-1;

  bMakeBiasCurve=0;              //Set to 1 if we want to make a bias curve
  uBiasCurveTrials=-1;
  fBiasCurveStart=-1;
  fBiasCurveStop=-1;
  fBiasCurveStep=-1;

  bLoopOverEvents=0;
  bWriteVFB=0;

  iNumberOfTelescopes = -1;          //The number of Telescopes in the array
  iNumberOfTelescopeTypes = -1; 


}

/*!
   2d vectors [telescope][pixel]
*/
void ReadConfig::resetCamVectors()
{
    cout<<"Resetting the Camera vectors"<<endl;
    vector< double > d_tel;
    vector< float > i_tel;
    vector< int > i_typ;
    for( int i = 0; i < iNumberOfTelescopeTypes; i++ )
    {
       //Quantum/PDE of the Photon detectors
       wl.push_back( i_tel ); 
       qe.push_back( i_tel );

       i_tel.assign( fCNChannels[i], 0. );
       d_tel.assign( fCNChannels[i], 0. );
       i_typ.assign( fCNChannels[i], 0 );
       fSizeTubeMM.push_back( d_tel );
       fRotAngle.push_back( d_tel );
       fXTubeMM.push_back( d_tel );
       fYTubeMM.push_back( d_tel );
       iTubeSides.push_back( i_typ );
       fXTube.push_back( i_tel );                //!< x-position of tube in [deg] (in camera coordinates)
       fYTube.push_back( i_tel );                //!< y-position of tube in [deg] (in camera coordiantes)
       fSizeTube.push_back( i_tel ); 

       vNumCellsPerSIPM.push_back( i_typ );      //the number of cells in one SiPM
       vSiPMOpticalCrosstalk.push_back( i_tel );

    }
    resetNeighbourGroupLists();
}

void ReadConfig::resetNeighbourGroupLists()
{
   vector< int > i_pix;

   vector< vector< int > > i_N;

   for( int i = 0; i < iNumberOfTelescopeTypes; i++ )
   {
       i_N.assign( fCNChannels[i], i_pix );
       fNeighbour.push_back( i_N );
       i_N.assign( iNumberGroups[i], i_pix );
       fPixelInGroup.push_back( i_N );
       fNeighbourGroups.push_back( i_N );
       i_N.clear();
       vPatch.push_back( i_N );        // pattern trigger patches

   }
}


/*! 
*/
void ReadConfig::resetTelVectors()
{
    vector< int > i_neighb;
	vTelescopeNeighbors.assign( iNumberOfTelescopes, i_neighb );
    iTelIDInSuperArray.assign( iNumberOfTelescopes, 0 );
    iTelType.assign( iNumberOfTelescopes, 0 );
    bBlurPSF.assign( iNumberOfTelescopes, 0 ); 
    fBlurSigma.assign( iNumberOfTelescopes, 0 );                 //Additional Bluring of the optical psf
    fSigmaElectronicNoise.assign( iNumberOfTelescopes, 0 );                 //The relative gain of each telescope
    fRelativeTelescopeGain.assign( iNumberOfTelescopes, 0 );                 //The relative gain of each telescope
    fWinstonConeEfficiency.assign( iNumberOfTelescopes, 0 );                 //The efficiency of the Winstoncone
}

//Function to reset all the vectors that hold the variables defining each telescope
void ReadConfig::resetTelTypeVectors()
{

  fCNChannels.assign( iNumberOfTelescopeTypes, 0 );
  iNumberGroups.assign( iNumberOfTelescopeTypes, 0 );  
   //Telescope trigger configuration
  bCrosstalk.assign( iNumberOfTelescopeTypes, 0 );                  //Use crosstalk between pixel
  fCrosstalkValue.assign( iNumberOfTelescopeTypes, 0 );             
  bSiPM.assign( iNumberOfTelescopeTypes, 0 );                       //we use SiPM as photon detectors
  bUseSumTrigger.assign( iNumberOfTelescopeTypes, 0 );              //Sum pixels before discriminator
  fClippingLevel.assign( iNumberOfTelescopeTypes, 0 );              //The level in mV at which the signals are clipped 
  bDoClipping.assign( iNumberOfTelescopeTypes, 0 );                 //Do we clip the signals before summing
    bCameraSnapshot.assign( iNumberOfTelescopeTypes, 0 );             //Use the camera snapshot logic (SST-1M)
  iSnapshotBits.assign( iNumberOfTelescopeTypes, 0 );               //Number of bit resolution for the camera snapshot (SST-1M)
  iSnapshotScalingDivisor.assign( iNumberOfTelescopeTypes, 0 );     //Value of the divisor to scale the group digitization for the camera snapshot trigger (SST-1M)
  iSnapshotFADCOffset.assign( iNumberOfTelescopeTypes, 0 );         //Offset (pedestal) to be added to the group digitization for the camera snapshot trigger (SST-1M)
  iSnapshotCircle.assign( iNumberOfTelescopeTypes, 0 );             //Number of circles to form the hexagonal pattern for the camera snapshot (SST-1M)
  iSnapshotNeighbors.assign( iNumberOfTelescopeTypes, 0 );          //Number of neighbors to be checked to build the pattern for the camera snapshot (SST-1M)
  iSnapshotComboMode.assign( iNumberOfTelescopeTypes, 0 );          //Number of neighbors to be checked to build the pattern for the camera snapshot (SST-1M)
  iSnapshotSamplingWindow.assign( iNumberOfTelescopeTypes, 0 );     //Number of neighbors to be checked to build the pattern for the camera snapshot (SST-1M)
  fDiscThreshold.assign( iNumberOfTelescopeTypes, 0 );              //Discriminator threshold of pixel
  fDiscWidth.assign( iNumberOfTelescopeTypes, 0 );                  //Width of Discriminator output
  fDiscDelay.assign( iNumberOfTelescopeTypes, 0 );                  //Delay of the inverted signal in the CFD
  fDiscAttenuation.assign( iNumberOfTelescopeTypes, 0 );            //Attenuation of the non-inverted signal in the CFD
  fDiscRFBConstant.assign( iNumberOfTelescopeTypes, 0 );            //The constant in the RFB feedback units pe/MHz
  fDiscPEtomVConversion.assign( iNumberOfTelescopeTypes, 0 );       //The conversion factor at the input of the Discriminator mV per pe.
  fFADCdctomVConversion.assign( iNumberOfTelescopeTypes, 0 );       //The conversion factor at the input of the FADC mV per pe.
  fPileUpWindow.assign( iNumberOfTelescopeTypes, 0 );       //The width of the window used to integrate all photons and get the right pulse shape.
  fDiscRFBDynamic.assign( iNumberOfTelescopeTypes, 0 );             //The value in the RFB feedback in units pe applied as offset.
                                               //If the RFB circuit is used this is just a start value 
  bDiscUseCFD.assign( iNumberOfTelescopeTypes, 0 );                 //Do we use the CFD part of the discriminator
  bDiscUseRFBCircuit.assign( iNumberOfTelescopeTypes, 0 );          //Do we use the RFB circuit
  iGroupMultiplicity.assign( iNumberOfTelescopeTypes, 0 );          //How many groups need to be in a cluster for a trigger  

  bUsePatches.assign( iNumberOfTelescopeTypes, 0 );                 //Is the trigger topology divided into patches


  bUseAfterPulsing.assign( iNumberOfTelescopeTypes, 0 );             //Do we simulate Afterpulsing: true yes false else
  bFlatfieldCamera.assign( iNumberOfTelescopeTypes, 0 );             //Flatfield the camera response
  fGainSigma.assign( iNumberOfTelescopeTypes, -1 );                  //sigma of the gain distribution
  fQESigma.assign( iNumberOfTelescopeTypes, -1 );                    //sigma of the QE distribution

  uMinNumPhotonsRequired.assign( iNumberOfTelescopeTypes, 0 );       //The minimum number of photons required in a camera to simulate the event
  


  //Trace related variables
  fFWHMofSinglePEPulse.assign( iNumberOfTelescopeTypes, 0 );        //The full width at hald maximum of the single pe pulse
  sHighGainPulseShapeFile.assign( iNumberOfTelescopeTypes, "" );         //name of the file that stores the single pe pulse shape
  sLowGainPulseShapeFile.assign( iNumberOfTelescopeTypes, "" );  //name of the file that stores the single pe pulse shape
  fSigmaSinglePEPulseHeightDistribution.assign( iNumberOfTelescopeTypes, 0 ); //The sigma of the single PE pulse height distribution
  fSampleWidthAveragePulse.assign( iNumberOfTelescopeTypes, 0 );    //The sample width used in the average single pe pulse
  fNSBRatePerPixel.assign( iNumberOfTelescopeTypes, 0 );            //the NSB rate per pixel in the focal plane;
  fAfterPulsingConstant.assign( iNumberOfTelescopeTypes, 0 );       //Constant of a fit to the rate vs. threshold curve of a single pe
                                       //Fit function exp(a+b*x), where a=constant b=slope
  fAfterPulsingSlope.assign( iNumberOfTelescopeTypes, 0 );          //Slope of a fit to the rate vs. threshold curve of a single pe
                                       //Fit function exp(a+b*x), where a=constant b=slope
  fTransitTimeSpread.assign( iNumberOfTelescopeTypes, 0 );       //Transit time spread

  fSamplingTime.assign( iNumberOfTelescopeTypes, 0 );       //The sampling rate or resolution of the simulated trace
  fTraceLength.assign( iNumberOfTelescopeTypes, 0 );        //the length of the simulated trace per group
  fStartSamplingBeforeAverageTime.assign( iNumberOfTelescopeTypes, 0 );    //Start sampling before the average photon arrival time
  //Cherenkovphoton throughput



  //FADC parameters
  iFADCSamples.assign( iNumberOfTelescopeTypes, 0 );                //The number of FADC samples
  fFADCSamplingWidth.assign( iNumberOfTelescopeTypes, 0 );          //The sampling width of the FADC
  fFADCDCtoPEconversion.assign( iNumberOfTelescopeTypes, 0 );           //The DC to PE conversion
  iFADCDynamicRange.assign( iNumberOfTelescopeTypes, 0 );           //The dynamic range of the FADC
  fFADCTimeOffsetFromTrigger.assign( iNumberOfTelescopeTypes, 0 );  //The offset between the trigger time and the start of the readout window
  fFADCHiLoGainThreshold.assign( iNumberOfTelescopeTypes, 0 );      //The Threshold in FADC counts at which switching to the low gain is activated
  fFADCLowHiGainRatio.assign( iNumberOfTelescopeTypes, 0 );         //Gain ratio logain/higain
  fFADCHighGainPedestal.assign( iNumberOfTelescopeTypes, 0 );               //The FADC high gain pedestal in dc counts
  fFADCLowGainPedestal.assign( iNumberOfTelescopeTypes, 0 );               //The FADC low gain pedestal in dc counts


  fMirFocalLength.assign( iNumberOfTelescopeTypes, 0 );                    //the focal length of the mirror in m



}

void ReadConfig::convertMMtoDeg()
{
    if( fDebug ) cout << "ReadConfig::convertMMtoDeg" << endl;
    for( int i = 0; i < iNumberOfTelescopeTypes; i++ )
    {
       for( unsigned int j = 0; j < fXTube[i].size(); j++ )
       {

	  fXTube[i][j] = atan2( (double)fXTubeMM[i][j] / 1000., (double)fMirFocalLength[i] ) * 45. / atan( 1. );
	  fYTube[i][j] = atan2( (double)fYTubeMM[i][j] / 1000., (double)fMirFocalLength[i] ) * 45. / atan( 1. );
	  fSizeTube[i][j] = atan2( (double)fSizeTubeMM[i][j] / 1000., (double)fMirFocalLength[i] ) * 45. / atan( 1. );
       }
    }
}


//---------------------------------------------------------------------------------------------------------
// Function to create and return the distribution of the relative QE between pixels in the camera
vector< Float_t > ReadConfig::GetRelQE(UInt_t telType)
{

  vector< Float_t > vRelQE;
  if(fQESigma[telType]<=0) //all photon sensors have the same QE, perfect world
    vRelQE.assign(fCNChannels[telType],1.0);
  else
    {  
	   cout<<"The photon sensors get different rel QEs, distributed following a normal distribution with a sigma of "<<fQESigma[telType]<<endl;
       for(UInt_t i = 0; i<fCNChannels[telType]; i++)
        {
      
           Float_t relQE = -1;
           while(relQE<=0)         
             relQE = rand->Gaus(1.0,fQESigma[telType]);

           vRelQE.push_back(relQE);
    
        }
    }

  return vRelQE;

}

//----------------------------------------------------------------------------------------------------
// Function to create and return the distributions of the relative gains between pixels in a camera
vector< Float_t > ReadConfig::GetRelGain(UInt_t telType, UInt_t telID, const vector< Float_t > vRelQE)
{

   vector< Float_t > vRelGain;

   if(fGainSigma[telType]<=0) //no gain distribution, perfect world for pixel. Only fill with relative telescope gain.
      vRelGain.assign(fCNChannels[telType],GetRelativeTelescopeGain(telID));
   else
     {
		cout<<"The photon sensors get different rel gains, distributed following a normal distribution with a sigma of "<<fGainSigma[telType]<<endl;
        for(UInt_t i = 0; i<fCNChannels[telType]; i++)
           {   
               Float_t relGain = -1;
               while(relGain<=0)
                  relGain = rand->Gaus(1.0,fGainSigma[telType]);

               //if the camera response is flatfielded with respect to an external light source
			   //the gain fluctuations are what you get from the flatfielded pixel responses.
               if(bFlatfieldCamera[telType])
				  { 
                    relGain /= vRelQE[i];
				  }  
               vRelGain.push_back(relGain*GetRelativeTelescopeGain(telID));
           }
     }
      
  return vRelGain;
}

void ReadConfig::ReadCommandLine( int argc, char **argv)
{
   cout<<"Parsing the command line"<<endl;
   //Dummy input file stream needed as argument by ReadLine
   std::ifstream *inFileStream = NULL;
   for (int i = 1; i < argc; ++i) {
      string iline="* ";
      iline+=argv[i];
      if (i + 1 < argc) { // Make sure we aren't at the end of argv!
         iline+=" ";
         iline+=argv[++i];
      }else { // Uh-oh, there was no argument to the option.
         std::cerr << "There is an option missing, so far found: "<<iline<< std::endl;
      }
     //cout<<"Have created "<<iline<<endl;
     ReadLine(iline,inFileStream);
   }
   delete inFileStream;

   cout<<"Finished reading the command line"<<endl<<endl; 
}

Bool_t ReadConfig::ReadConfigFile( string iFile )
{
   if( fDebug ) cout << "ReadConfig::ReadConfigFile " << iFile << endl;

   //iFile.insert( 0, fConfigDir );
   std::ifstream *inFileStream = new ifstream( iFile.c_str() );
   if( !*inFileStream )
   {
      cout << "TriggerRead::ReadConfigFile() error: config file not found: " << iFile << endl;
      exit( -1 );
      return kFALSE;
   }

   string iline;

   while( getline( *inFileStream, iline ) )
   {
     // '*' in line 
     if( iline.substr( 0, 1 ) != "*" ) continue;

     // '#' in line
     if( iline.substr( 0, 1 ) == "#" ) continue;

     ReadLine(iline,inFileStream);     

   }

   convertMMtoDeg();

   delete inFileStream;

   if( fDebug ) cout << "END: ReadConfig::ReadConfigFile " << iFile << endl;

   return true;

}

//Read in one line 
void ReadConfig::ReadLine(string iline, ifstream *inFileStream)
{

   string i_char;
   unsigned int i_telType = 0;
   unsigned int i_chan = 0;
   unsigned int i_NN = 0;


     istringstream i_stream( iline );

     cout<<iline.c_str()<<endl;
     //the number of telescopes
     if( iline.find( "NBRTL " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >> iNumberOfTelescopes;
    cout<<" We have "<<iNumberOfTelescopes<<" telescopes in the array"<<endl;

    resetTelVectors();     
      }

     //the number of telescope types
     if( iline.find( "NBRTELTYPES " ) < iline.size() )
      {
	      i_stream >> i_char; i_stream >> i_char; 
	      i_stream >> iNumberOfTelescopeTypes;
          cout<<" We have "<<(Int_t)iNumberOfTelescopeTypes<<" telescope TYPES in the array"<<endl;
          resetTelTypeVectors();
       
      }

    if( iline.find( "CAMRA " ) < iline.size() )
      {     

        //   -telescope type
        //  -the number of phototubes.
        //  -the number of groups of summed pixel
        //  * CAMRA 0 11328 2832
         i_stream >> i_char; i_stream >> i_char;
	 i_stream >> i_telType;                                   
         i_stream >> fCNChannels[i_telType];
         i_stream >> iNumberGroups[i_telType];
         cout << "FPI configuration of telescope type "<<i_telType<<endl; 
	 cout << "total number of channels: " << fCNChannels[i_telType]  << endl;
	 cout << "the number of sum groups: " << iNumberGroups[i_telType]  << endl;

        if(i_telType == (UInt_t)(iNumberOfTelescopeTypes-1))
             resetCamVectors();
      }


     // Focal length of the mirror. This is used to convert the focal plane from mm to deg (needed to write the vbf file)
     if( iline.find( "TELESCOPEFOCALLENGTH " ) < iline.size() )
      {
        i_stream >> i_char; i_stream >> i_char;
	i_stream >> i_telType;
        i_stream >> fMirFocalLength[i_telType];
        cout<<"Telescope type "<<i_telType<<"  Focal Lengths of telescopes in m for conversion of pixel from mm to deg: "<<fMirFocalLength[i_telType]<<endl;
       }
 
       
     


     //The Start of sampling the Trace before the average photon arrival time
     if( iline.find( "STARTSAMPLINGBEFOREAVERAGEPHOTONARRIVALTIME " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
        i_stream >> i_telType;
	i_stream >> fStartSamplingBeforeAverageTime[i_telType];
	cout<<"Telescope type "<<i_telType<<" The sampling of the trace start ns "<<fStartSamplingBeforeAverageTime[i_telType]<<" before the average photon arrival time"<<endl;
      }

     //Length of the trace
     if( iline.find( "TRACELENGTH " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
        i_stream >> i_telType;
	i_stream >> fTraceLength[i_telType];
	cout<<"Telescope type "<<i_telType<<" Trace length in ns set to "<<fTraceLength[i_telType]<<endl;
      }

      //Width of one sample  of the trace
     if( iline.find( "TRACESAMPLEWIDTH " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >> i_telType;
	i_stream >> fSamplingTime[i_telType];
	cout<<"Telescope type "<<i_telType<<" Width of one sample in ns of the Trace set to "<<fSamplingTime[i_telType]<<endl;
      }


     //The FWHM in ns of the single pe pulse shape used in the simulation
     if( iline.find( "SINGLEPEWIDTH " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
        i_stream >> i_telType;
	i_stream >>  fFWHMofSinglePEPulse[i_telType];
	if(fFWHMofSinglePEPulse[i_telType]>0)
	{
      	  
	  cout<<"Telescope type "<<i_telType<<" The simulated single pe pulse shape is a gauss with FWHM "<<fFWHMofSinglePEPulse[i_telType]<<endl;
	}
	else
	{
	  cout<<"will look for a single pe pulse later...I hope you have set a text file!!!"<<endl;
	}
      }

     //The filename of the highgain pulse shapes that is used if singlepewidth is 0
     if( iline.find( "HIGHGAINPULSESHAPE " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 	
        i_stream >> i_telType;
	i_stream >>  sHighGainPulseShapeFile[i_telType];
      }

   //The sigma of the single PE pulse height distribution in units of PE
     if( iline.find( "SINGLEPEAMPLSIGMA " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
        i_stream >> i_telType;
	i_stream >>   fSigmaSinglePEPulseHeightDistribution[i_telType];
	cout<<"Telescope type "<<i_telType<<" The sigma of the single PE pulse height distribution "<< fSigmaSinglePEPulseHeightDistribution[i_telType]<<endl;
      }

     //The sampling width used in the sample single pe pulse
     if( iline.find( "SINGLEPESAMPLING " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
        i_stream >> i_telType;
	i_stream >>  fSampleWidthAveragePulse[i_telType];
	cout<<"Telescope type "<<i_telType<<" The sampling width used in the  single pe pulse in ns is "<<fSampleWidthAveragePulse[i_telType]<<endl;
      }

     //The filename of the file that contains the low gain pulse shapes
     if( iline.find( "LOWGAINPULSESHAPE " ) < iline.size() )
      {               
	i_stream >> i_char; i_stream >> i_char; 	
        i_stream >> i_telType;
	i_stream >>  sLowGainPulseShapeFile[i_telType];
      }



     //If we want to use NSB in the simulation
     if( iline.find( "USENSB " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >>  bUseNSB;
	cout<<"Do we want to use NSB "<<bUseNSB<<endl;
      }

     //The rate of NSB per pixel in the focal plane
     if( iline.find( "NSBRATEPERPIXEL " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >>  fNSBRatePerPixel[i_telType];
        fNSBRatePerPixel[i_telType]*=1000;
	cout<<"Telescope type "<<i_telType<<" The NSB rate in kHz per pixel in the focal plane is set to "<<fNSBRatePerPixel[i_telType]<<endl;
      }

     //If we want to use Afterpulsing in the simulation
     if( iline.find( "USEAFTERPULSING " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char;
	i_stream >> i_telType;
	int tmp;
	i_stream >> tmp;
	bUseAfterPulsing[i_telType] = (Bool_t)tmp;
	cout<<"Telescope type "<<i_telType<<" Do we want to use Afterpulsing "<<bUseAfterPulsing[i_telType]<<endl;
      }
      
	  //If we want to use Afterpulsing in the simulation
     if( iline.find( "AFTERPULSINGCONSTANT " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >>  fAfterPulsingConstant[i_telType];
	if(fAfterPulsingConstant[i_telType] > 0)
         {
           cout<<"Telescope type "<<i_telType<<" TraceGenerator: fAfterPulsingConstant set to a value >0, bad!"<<endl;
           exit(1);
         }	
          cout<<"Telescope type "<<i_telType<<" Constant of a fit of the rate vs. Threshold measurement with exp(a+b*x), where a=constant= "
             <<fAfterPulsingConstant[i_telType]<<endl;
      }


      //If we want to use Afterpulsing in the simulation
     if( iline.find( "AFTERPULSINGSLOPE " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >>  fAfterPulsingSlope[i_telType];
	if(fAfterPulsingSlope[i_telType] > 0)
         {
           cout<<"Telescope type "<<i_telType<<" TraceGenerator: fAfterPulsingSlope set to a value >0, bad!"<<endl;
           exit(1);
         }	
	cout<<"Telescope type "<<i_telType<<" Slope of a fit of the rate vs. Threshold measurement with exp(a+b*x), where b=slope= "
                <<fAfterPulsingSlope[i_telType]<<endl;
      }

    //If we want to use a transit time spread for the photo sensors
     if( iline.find( "TRANSITTIMESPREAD " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >>  fTransitTimeSpread[i_telType];
	if(fTransitTimeSpread[i_telType] < 0)
         {
           cout<<"Telescope type "<<i_telType<<" TraceGenerator: fTransitTimeSpread set to a value <0, bad!"<<endl;
           exit(1);
         }	
	cout<<"Telescope type "<<i_telType<<" Transit time spread (RMS) [ns] = "
                <<fTransitTimeSpread[i_telType]<<endl;
      }

     if( iline.find( "QESIGMA " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >>  fQESigma[i_telType] ;
	
	cout<<"Telescope type "<<i_telType<<" The sigma of the QE distribution is: "<<fQESigma[i_telType]<<endl;
      }

     if( iline.find( "GAINSIGMA " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >>  fGainSigma[i_telType] ;
	
	cout<<"Telescope type "<<i_telType<<" The sigma of the gain distribution is: "<<fGainSigma[i_telType]<<endl;
      }
    
	
	if( iline.find( "FLATFIELDCAMERA " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
    int tmp;
	i_stream >> tmp;
	bFlatfieldCamera[i_telType] = (Bool_t)tmp;
	
	cout<<"Telescope type "<<i_telType<<" The Cameras are flatfielded: "<<bFlatfieldCamera[i_telType]<<endl;
      }                                                                                  


     //Sets the factor by which the Cherenkov Photons need to be scaled down
     //to reflect the proper efficiency of the PMTs
     if( iline.find( "QUEFF " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	Int_t datapoints;
	i_stream >>  datapoints ;
	wl[i_telType].clear();
	qe[i_telType].clear();
	cout<<endl<<"Telescope type "<<i_telType<<" Reading in QE"<<endl;
	cout<<datapoints<<endl;
	for(Int_t i=0;i<datapoints;i++)
	  {
	    getline( *inFileStream, iline );
	    istringstream d_stream( iline );
	    Float_t w;
	    Float_t q;
	    d_stream >> w;
	    d_stream >> q;
	    wl[i_telType].push_back(w);
	    qe[i_telType].push_back(q);
	    cout<<w<<"  "<<q<<endl;
	      
	  }
	cout<<endl;
      }




     //Clipping level at X mV
     if( iline.find( "DISCCLIPPINGLEVEL " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fClippingLevel[i_telType];
	cout<<"Telescope type "<<i_telType<<" Signal level in mV at which signals will be clipped before summing "<<fClippingLevel[i_telType]<<endl;
      }

     
     /* SST-1M DEV PART */
     //Enable/Disable the camera snapshot logic
     if( iline.find( "USETRIGGERCAMERASNAPSHOT " ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char; 
	 i_stream >> i_telType;
	 int tmp;
	 i_stream >> tmp;
	 bCameraSnapshot[i_telType] = (Bool_t)tmp;
	 cout<<"Telescope type "<<i_telType<<" Do we want to use the 'camera snapshot' logic for the telescope trigger: "<<bCameraSnapshot[i_telType]<<endl;
       }
  
     //Resoution in bits of the dynamic range of the digitization of a group for the camera snapshot logic
     if( iline.find( "SNAPSHOTRESOLUTION " ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char; 
	 i_stream >> i_telType;
	 i_stream >> iSnapshotBits[i_telType];
	 cout<<"Telescope type "<<i_telType<<" The resolution in bits for the 'camera snapshot' logic is "<<iSnapshotBits[i_telType]<<endl;
       }
  
     //Scaling divisor to reduce the digitization of a group for the camera snapshot logic (hint: values are usually powers of 2)
     if( iline.find( "SNAPSHOTSCALINGDIVISOR " ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char; 
	 i_stream >> i_telType;
	 i_stream >> iSnapshotScalingDivisor[i_telType];
	 cout<<"Telescope type "<<i_telType<<" The scaling divisor (values are usually powers of 2) for the 'camera snapshot' logic is "<<iSnapshotScalingDivisor[i_telType]<<endl;
       }

     //Offset (pedestal), in DC, to be summed up to the digitization of a group for the camera snapshot logic
     if( iline.find( "SNAPSHOTFADCOFFSET " ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char; 
	 i_stream >> i_telType;
	 i_stream >> iSnapshotFADCOffset[i_telType];
	 cout<<"Telescope type "<<i_telType<<" The offset (pedestal) in DC of the digitization for the 'camera snapshot' logic is "<<iSnapshotFADCOffset[i_telType]<<endl;
       }


     //Number of circles surrounding the centeral group
     if( iline.find( "SNAPSHOTPATCHCIRCLES " ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char; 
	 i_stream >> i_telType;
	 i_stream >> iSnapshotCircle[i_telType];
	 cout<<"Telescope type "<<i_telType<<" The number of circles to find the patterns for the 'camera snapshot' logic is "<<iSnapshotCircle[i_telType]<<endl;
       }
  
     //How many neighbors are surrounding a group to form a pattern?
     if( iline.find( "SNAPSHOTSURROUNDINGNEIGHBORS " ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char;
	 i_stream >> i_telType;
	 i_stream >> iSnapshotNeighbors[i_telType];
	 cout<<"Telescope type "<<i_telType<<" The number of neighbors to let a group to form a pattern for the 'camera snapshot' logic is "<<iSnapshotNeighbors[i_telType]<<endl;
       }

     //Which mode of the combination of the snapshosts is used and its samplings window lenght 
     if( iline.find( "SNAPSHOTSCOMBINATION" ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char; 
	 i_stream >> i_telType;
	 i_stream >> iSnapshotComboMode[i_telType];
	 i_stream >> iSnapshotSamplingWindow[i_telType];
	 cout<<"Telescope type "<<i_telType<<" The snapshots combination mode is "<<iSnapshotComboMode[i_telType];
	 switch ( iSnapshotComboMode[i_telType] )
	   {
	   case (0):
	     cout << "the Blind mode: sequential samplings (";
	     break;
	   case (1):
	     cout << "the Edge mode: contained edges of samplings windows (";
	     break;
	   case (2):
	     cout << "the Level mode: ... (";
	     break;
	   default:
	     cout << "unknown: forced to Blind mode! (";
	     iSnapshotComboMode[i_telType] = 0;
	   }
	 cout <<iSnapshotComboMode[i_telType]<<") with a sampling window of "<<iSnapshotSamplingWindow[i_telType]<<endl;
       }
     /* END SST-1M DEV PART */
  
     //Clipping the signal before summing
     if( iline.find( "USEDISCCLIPPING " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
    int tmp;
	i_stream >> tmp;
    bDoClipping[i_telType] = (Bool_t)tmp;
	cout<<"Telescope type "<<i_telType<<" Signals will be clipped before summing "<<bDoClipping[i_telType]<<endl;
      }

     //The threshold of the discriminator
     if( iline.find( "DISCTHRESHOLD " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fDiscThreshold[i_telType];
	cout<<"Telescope type "<<i_telType<<" The Threshold of the discriminator in mV "<<fDiscThreshold[i_telType]<<endl;
      }

     //The width of the output signal of the discriminator
     if( iline.find( "DISCWIDTH " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fDiscWidth[i_telType];
	cout<<"Telescope type "<<i_telType<<" The width of the discriminator output signal in ns "<<fDiscWidth[i_telType]<<endl;
      }

     //The delay of the inverted signal in the CFD
     if( iline.find( "DISCDELAY " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fDiscDelay[i_telType];
	cout<<"Telescope type "<<i_telType<<" The delay of the inverted signal in the CFD in ns "<<fDiscDelay[i_telType]<<endl;
      }

     //The attenuation of the non-inverted signal in the CFD
     if( iline.find( "DISCATTENUATION " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fDiscAttenuation[i_telType];
	cout<<"Telescope type "<<i_telType<<" The attenuation of the non-inverted signal in the CFD "<<fDiscAttenuation[i_telType]<<endl;
      }

     //Do we use the CFD part of the discriminator?
     if( iline.find( "DISCUSECFD " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType; 
    int tmp;
    i_stream >> tmp;
    bDiscUseCFD[i_telType] = (Bool_t)tmp;
	cout<<"Telescope type "<<i_telType<<" Do we use the CFD part of the discriminator: "<<bDiscUseCFD[i_telType]<<endl;
      }

     //Do we use the RFB circuit?
     if( iline.find( "DISCUSERFBCIRCUIT " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
    int tmp;
    i_stream >> tmp;
    bDiscUseRFBCircuit[i_telType] = (Bool_t)tmp;
	cout<<"Telescope type "<<i_telType<<" Do we use the RFB circuit: "<<bDiscUseRFBCircuit[i_telType]<<endl;
      }

     //The Constant Value in the RFB in the discriminator
     if( iline.find( "DISCRFBCONSTANT " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fDiscRFBConstant[i_telType];
	cout<<"Telescope type "<<i_telType<<" If the RFB circuit is used this is the magnitude of the rate dependend feedback  in mV/MHz "<<fDiscRFBConstant[i_telType]<<endl;
      }

     //The dyamic value in the RFB in the discriminator
     if( iline.find( "DISCRFBDYNAMIC " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fDiscRFBDynamic[i_telType];
	cout<<"Telescope type "<<i_telType<<" The initial dynamic value in the RFB in mV (will be mulitplied with 0.18 in the sims). If RFB circuit is not used this stays constant "<<fDiscRFBDynamic[i_telType]<<endl;
      }
                             
     //The Conversion factor at the input of the discriminator mV per PE
     if( iline.find( "DISCPETOMVCONVERSION " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fDiscPEtomVConversion[i_telType];
	cout<<"Telescope type "<<i_telType<<" The conversion factor at the input of the discriminator mV per pe (Amplitude)"<<fDiscPEtomVConversion[i_telType]<<endl;
      }

     //The Conversion factor at the input of the FADC mV per dc. This is needed to find the right pulse shape when the trace gets assembled
     if( iline.find( "FADCDCTOMVCONVERSION " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fFADCdctomVConversion[i_telType];
	cout<<"Telescope type "<<i_telType<<" The conversion factor at the input of the FADC mV per dc (Amplitude) where the FADC response is linear. This number is used to pick the right pulse shape when assembling the trace. Note that for the trigger the high gain trace is used"<<fFADCdctomVConversion[i_telType]<<endl;
      }

     //The width of the window that is used to determine the number of photons that pile up to produce one output signal. 
     //That number is used to find the right pulse shape 
     if( iline.find( "PILEUPWINDOW " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fPileUpWindow[i_telType];
        fPileUpWindow[i_telType]=fPileUpWindow[i_telType]/2.0;
	cout<<"Telescope type "<<i_telType<<" The half width of the pileup window used to find the right pulse shape for each photon is"<<fPileUpWindow[i_telType]<<" ns"<<endl;
      }

     //How many goups need to be in a cluster for a telescope trigger
     if( iline.find( "GROUPMULTIPLICITY " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> iGroupMultiplicity[i_telType];
	cout<<"Telescope type "<<i_telType<<" The multiplicity requirement for a telescope trigger is "<<iGroupMultiplicity[i_telType]<<endl;
      }

     //How many Telescopes need to trigger to get an array trigger
     if( iline.find( "TELMULTIPLICITY " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >> iTelescopeMultiplicity;
	cout<<"The multiplicity requirement for an array trigger is "<<iTelescopeMultiplicity<<endl;
      }

     //Defines wether the array trigger requires next neighbor requirement 1 or not 0
     if( iline.find( "TELNEXTNEIGHBOR " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >> bArrayTriggerRequiresNextNeighbor;
	cout<<"The array next neighbor requirement is "<<bArrayTriggerRequiresNextNeighbor<<endl;
      }


     //The coincidence window on array level in ns
     if( iline.find( "ARRAYCOINCIDENCE " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >>  fArrayCoincidence;
	cout<<"The coincidence window in ns on array level is "<<fArrayCoincidence<<endl;
      }

     //Set to true if we want to make a bias curve
     if( iline.find( "MAKEBIASCURVE " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >>  bMakeBiasCurve;
	cout<<"We make a bias curve "<<bMakeBiasCurve<<endl;
      }

     //Number of trials for the bias curve
     if( iline.find( "BIASCURVERTELESCOPEID " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >> iBiasCurveTelescopeID ;
	cout<<"the bias curve will be done for array telescope ID: "<<iBiasCurveTelescopeID<<endl;
      }
     
      //Number of trials for the bias curve
     if( iline.find( "BIASCURVETRIALS " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >> uBiasCurveTrials ;
	cout<<"Trials that are done to make the bias curve "<<uBiasCurveTrials<<endl;
      }

     //Start of the scan range for the bias curve in pe
     if( iline.find( "BIASCURVESTART " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >> fBiasCurveStart ;
	cout<<"Start of the scan range for the bias curve in mV: "<<fBiasCurveStart<<endl;
      }

     //Stop of the scan range for the bias curve in pe
     if( iline.find( "BIASCURVESTOP " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >> fBiasCurveStop ;
	cout<<"Stop of the scan range for the bias curve in mV: "<<fBiasCurveStop<<endl;
      }

     //Step in the scan of the discriminator values for the bias curve in pe
     if( iline.find( "BIASCURVESTEP " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >> fBiasCurveStep ;
	cout<<"Step in the scan of the discriminator values for the bias curve in mV: "<<fBiasCurveStep<<endl;
      }

     //Set to true if we want to loop over events
     if( iline.find( "LOOPOVEREVENTS " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >>  bLoopOverEvents ;
	cout<<"We loop over events "<< bLoopOverEvents<<endl;
      }

     //The name of the person executing the simulation
     if( iline.find( "SIMULATORNAME " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
	i_stream >>  sSimulatorName ;
	cout<<"The Person executing the simulations: "<<sSimulatorName<<endl;
      }

     //The atmospheric model used in the simulation
     if( iline.find( "ATMOSPHERICMODEL " ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char; 
	 i_stream >>  iAtmosphericModel ;
	 cout<<"The atmospheric model used: "<<iAtmosphericModel<<endl;
      }

     //The day to which the datums of the simulated events are set
     if( iline.find( "DAYOFSIMULATEDEVENTS " ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char; 
	 i_stream >>  sDayOfSimulatedEvents ;
	 cout<<"The simulated events will have a datum of: "<<sDayOfSimulatedEvents<<endl;
      }

     //The number of pedestal events to simulate and write to disk
	 if( iline.find( "NUMBERPEDESTALEVENTS " ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char; 
	 i_stream >>  iNumberPedestalEvents ;
	 cout<<"The number of pedestal events that will be simulated: "<<iNumberPedestalEvents<<endl;
	 cout<<"Note that this requires that you also set the write pedestal flag when calling CARE"<<endl;
      }

     //The number of pedestal events to simulate in order to stabilize the discriminator
     if( iline.find( "NUMBERPEDESTALEVENTSTOSTABILZE " ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char; 
	 i_stream >>  iNumberPedestalEventsToStabilize ;
	 cout<<"The number of pedestal events that will be simulated to stabilize the discriminator: "<<iNumberPedestalEventsToStabilize<<endl;
      }

     //The day to which the datums of the simulated events are set
     if( iline.find( "WRITEVBF " ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char; 
	 i_stream >>  bWriteVFB ;
	 cout<<"Is a VBF file written?: "<<bWriteVFB<<endl;
      }

     //the number of sample in a FADC trace
     if( iline.find( "FADCSAMPLES " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> iFADCSamples[i_telType];         
    cout<<"Telescope type "<<i_telType<<" The number of FADC samples per recorded trace "<<iFADCSamples[i_telType]<<endl;
      }


     //the sampling width of the FADC
     if( iline.find( "FADCSAMPLINGWIDTH " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fFADCSamplingWidth[i_telType];
	cout<<"Telescope type "<<i_telType<<" The sampling time of the FADC in ns is: "<<fFADCSamplingWidth[i_telType]<<endl;
      }


     //the dynamic range of the FADC
     if( iline.find( "FADCDYNAMICRANGE " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> iFADCDynamicRange[i_telType];
    cout<<"Telescope type "<<i_telType<<" The dynamic range of the FADC is: "<<iFADCDynamicRange[i_telType]<<endl;
      }

   //the Offset of the beginning of the readout window from the trigger decision
     if( iline.find( "FADCOFFSETFROMTRIGGER " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fFADCTimeOffsetFromTrigger[i_telType];
    cout<<"Telescope type "<<i_telType<<"  The offset between the beginning of the readout window from the trigger time in ns is: "<<fFADCTimeOffsetFromTrigger[i_telType]<<endl;
      }
 
     //the threshold when the low gain Channel is activated
     if( iline.find( "FADCHILOGAINTHRESHOLD " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fFADCHiLoGainThreshold[i_telType];
    cout<<"Telescope type "<<i_telType<<"  The lo gain is activated at dc : "<<fFADCHiLoGainThreshold[i_telType]<<endl;

      }

   //the gain ratio of the low gain channel / high gain channel 
     if( iline.find( "FADCLOHIGHGAINRATIO " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fFADCLowHiGainRatio[i_telType];
    cout<<"Telescope type "<<i_telType<<" The gain ratio between the low gain channel /  the high gain channel (amplitude) "<<fFADCLowHiGainRatio[i_telType]<<endl;
      }

   //the FADC pedestal 
     if( iline.find( "FADCHIGHGAINPEDESTAL " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fFADCHighGainPedestal[i_telType];
    cout<<"Telescope type "<<i_telType<<" The FADC high gain pedestal in dc is: "<<fFADCHighGainPedestal[i_telType]<<endl;
      }
     if( iline.find( "FADCLOWGAINPEDESTAL " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fFADCLowGainPedestal[i_telType];
    cout<<"Telescope type "<<i_telType<<" The FADC low gain pedestal in dc is: "<<fFADCLowGainPedestal[i_telType]<<endl;
      }


    //the conversion from DC to PE
     if( iline.find( "FADCPETODCCONVERSION " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fFADCDCtoPEconversion[i_telType];
	cout<<"Telescope type "<<i_telType<<" The FADC conversion factor from DC counts amplitude  to PE is [dc/pe]: "<<fFADCDCtoPEconversion[i_telType]<<endl;
      }

    //the minimum number of Cherenkov photons in a telescope
     if( iline.find( "MINNUMPHOTONSREQUIRED " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> uMinNumPhotonsRequired[i_telType];
	cout<<"Telescope type "<<i_telType<<" If more than "<<uMinNumPhotonsRequired[i_telType]<<" Cherenkov photons are in the focal plane of the camera the event is simulated"<<endl;
      }


    //Do we use Crosstalk between pixel
    if( iline.find( "USECROSSTALK " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
    int tmp;
    i_stream >> tmp;
    bCrosstalk[i_telType] = (Bool_t)tmp;
	cout<<"Telescope type "<<i_telType<<" Will use Crosstalk between pixel : "<<bCrosstalk[i_telType]<<endl;
      }
    
    //What is the crosstalk value 
    if( iline.find( "CROSSTALKVALUE " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
	i_stream >> fCrosstalkValue[i_telType];
    cout<<"Telescope type "<<i_telType<<" The Crosstalk value is: "<<fCrosstalkValue[i_telType]<<endl;
      }


    //Do we use SiPMs
    if( iline.find( "USESIPM " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
    int tmp;
    i_stream >> tmp;
    bSiPM[i_telType] = (Bool_t)tmp;
	cout<<"Telescope type "<<i_telType<<" Will use SiPM : "<<bSiPM[i_telType]<<endl;
      }

    // SiPM parameters 
    if( iline.find( "SIPMPARS " ) < iline.size() )
      {
        i_stream >> i_char; i_stream >> i_char;
	    i_stream >> i_telType;
        i_stream >> i_chan;
        i_stream >> vNumCellsPerSIPM[i_telType][i_chan]; 
	    i_stream >> vSiPMOpticalCrosstalk[i_telType][i_chan]; 
        cout<<"Chan "<<i_chan<<" teltype  "<<i_telType<<" number of cells "<<vNumCellsPerSIPM[i_telType][i_chan]<<" optical X-talk "<<vSiPMOpticalCrosstalk[i_telType][i_chan]<<endl;
       } 

    //Do we use a sum trigger
     if( iline.find( "USESUMTRIGGER " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
    int tmp;
    i_stream >> tmp;
    bUseSumTrigger[i_telType] = (Bool_t)tmp;
	cout<<"Telescope type "<<i_telType<<" Will use Sum trigger : "<<bUseSumTrigger[i_telType]<<endl;
      }

   //Do we want to use the Patches 
     if( iline.find( "USETRIGGERSUBFIELDS " ) < iline.size() )
      {
	i_stream >> i_char; i_stream >> i_char; 
    i_stream >> i_telType;
    int tmp;
    i_stream >> tmp;
    bUsePatches[i_telType]= (Bool_t)tmp;
	cout<<"Telescope type "<<i_telType<<" Will use patches : "<<bUsePatches[i_telType]<<endl;
      }
   
     if( iline.find( "PATCH " ) < iline.size() )
       {
	 i_stream >> i_char; i_stream >> i_char;
     i_stream >> i_telType;
     i_stream >> i_char;
	 vector< int > i_tPatch;
	 cout<<"Telescope type "<<i_telType<<" Reading PATCH "<<vPatch[i_telType].size()<<" : ";
	 while(i_stream.good())
	   {
	     Int_t pix;
	     i_stream >> pix;
	     i_tPatch.push_back(pix);
	     cout<<pix<<" ";
	   }
	 cout<<endl;
	 vPatch[i_telType].push_back( i_tPatch );
       }
     


// telescope positions

     if( iline.find( "TLCFG " ) < iline.size() )
      {
         int id;
         i_stream >> i_char; i_stream >> i_char;
	     i_stream >> id;
	     i_stream >> iTelIDInSuperArray[id];
	     i_stream >> iTelType[id];
	     i_stream >> fSigmaElectronicNoise[id];
	     i_stream >> fRelativeTelescopeGain[id];
	     i_stream >> fWinstonConeEfficiency[id];
	     i_stream >> fBlurSigma[id];

         bBlurPSF[id] = fBlurSigma[id]<=0 ? kFALSE : kTRUE; 
         
         if( fSigmaElectronicNoise[id]<0)
          {
           cout<<"TelescopeID: "<<id<<" the electronic noise for this telescope is set wrong. It has to be larger than 0 "<<endl;
           cout<<"You tried to set it to "<<fSigmaElectronicNoise[id]<<endl;
           exit(1);
          }

         if( fRelativeTelescopeGain[id]<=0 ||  fRelativeTelescopeGain[id]>1)
          {
           cout<<"TelescopeID: "<<id<<" the relative gain of this telescope is set wrong. It has to be between 0 and 1"<<endl;
           cout<<"You tried to set it to "<<fRelativeTelescopeGain[id]<<endl;
           exit(1);
          }
         
         if( fWinstonConeEfficiency[id]<=0 ||  fWinstonConeEfficiency[id]>1)
          {
           cout<<"TelescopeID: "<<id<<" the efficiency of the Winston cone is set wrong. It has to be between 0 and 1"<<endl;
           cout<<"You tried to set it to "<<fWinstonConeEfficiency[id]<<endl;
           exit(1);
          }
         

         cout <<"Telescope :"<<id<<" SuperArrayID: "<<iTelIDInSuperArray[id]<<" telType: "<<iTelType[id]<<" (id=" << id << ")" 
     <<" sigma electronic noise: "<<fSigmaElectronicNoise[id]<<" relative gain: "<<fRelativeTelescopeGain[id]<<";  Winston Cone efficieny : "<<fWinstonConeEfficiency[id]<<" BlurSigma [mm]: "<<fBlurSigma[id]<<endl;
      }


// tube stuff 
     if( iline.find( "PMPIX " ) < iline.size() )
      {
        i_stream >> i_char; i_stream >> i_char;
	    i_stream >> i_telType;
        int tubetype;
	    i_stream >> tubetype;
	    i_stream >> i_chan;
		int sides = 0;
		if(tubetype == 0)//hexagon
			sides = 6;
		if(tubetype == 1)//square
			sides = 4;
		if(tubetype == 2)//circle
			sides = 1;

        iTubeSides[i_telType][i_chan] = sides; 
	    i_stream >> fXTubeMM[i_telType][i_chan]; 
	    i_stream >> fYTubeMM[i_telType][i_chan];
	    i_stream >> fSizeTubeMM[i_telType][i_chan];
		i_stream >> fRotAngle[i_telType][i_chan];
		fRotAngle[i_telType][i_chan]*=(TMath::DegToRad());
        cout<<"Chan "<<i_chan<<" telType  "<<i_telType<<" Tubesides "<<iTubeSides[i_telType][i_chan]<<" x-Pos [mm] "<<fXTubeMM[i_telType][i_chan]<<" y-Pos [mm]  "<<fYTubeMM[i_telType][i_chan]<<" diameter [mm] "<<fSizeTubeMM[i_telType][i_chan]<<" rotAngle [deg]  "<<fRotAngle[i_telType][i_chan]<<endl;
       } 



// neighbour list
      if( iline.find( "PIXNGHBR " ) < iline.size() )
       {
           // * NGHBR 0  1  5  0  2  8  9  10
           i_stream >> i_char; i_stream >> i_char;
	       i_stream >> i_telType;
	       i_stream >> i_chan;
	       i_stream >> i_NN;
           cout<<" Chan "<<i_chan<<" TelType  "<<i_telType<<" NumNGHBRs "<<i_NN<<" List of neighbrs: ";
	      for( unsigned int j = 0; j < i_NN; j++ )
	      {
	         if( !i_stream.eof() )
		       {
                  Int_t value;
		          i_stream >> value;
                  fNeighbour[i_telType][i_chan].push_back(value);
                  cout<<" "<<fNeighbour[i_telType][i_chan][j];
                }
		           else cout<<"Something is wrong with the neighbors of pixel "<<i_chan<<endl;
           }
              cout<<endl;
	    }

// member list of the groups
      if( iline.find( "GROUP " ) < iline.size() )
       {
           //  * GROUP 0 4 4 16  17  24  25
           i_stream >> i_char; i_stream >> i_char;
	       i_stream >> i_telType;
	       i_stream >> i_NN;
	       i_stream >> i_chan;
           cout<<i_chan<<" "<<i_telType<<" "<<i_NN;
	      for( unsigned int j = 0; j < i_NN; j++ )
	      {
	         if( !i_stream.eof() )
		       {  
                  Int_t value;
		          i_stream >> value;
                  fPixelInGroup[i_telType][i_chan].push_back(value);
                  cout<<" "<<fPixelInGroup[i_telType][i_chan][j];
                }
		           else cout<<"Something is wrong with the members of group "<<i_chan<<endl;;
           }
           cout<<endl;
	    }

// neighbour list of groups
      if( iline.find( "GRPNGHBR " ) < iline.size() )
       {
           //  * GRPNGHBR 0  7  8  2  3  6  10  11  112  116  120
           i_stream >> i_char; i_stream >> i_char;
	       i_stream >> i_telType;
	       i_stream >> i_chan;
	       i_stream >> i_NN;
           cout<<i_chan<<" "<<i_telType<<" "<<i_NN;
	      for( unsigned int j = 0; j < i_NN; j++ )
	      {
	         if( !i_stream.eof() )
		       {                                                   
                  Int_t value;
		          i_stream >> value;
                  fNeighbourGroups[i_telType][i_chan].push_back(value);
                  cout<<" "<<fNeighbourGroups[i_telType][i_chan][j];
                }
		           else cout<<"Something is wrong with the neighbors of group "<<i_chan<<endl;;
           }
           cout<<endl;
	    }


    // The neighbors of each telescope considered in an array trigger
      if( iline.find( "TLNEIGHBRS " ) < iline.size() )
       {
           //  * TLNEIGHBRS 0 3 1 2 3
           i_stream >> i_char; i_stream >> i_char;
	       i_stream >> i_telType;
	       i_stream >> i_NN;
           cout<<i_telType<<" "<<i_NN;
	      for( unsigned int j = 0; j < i_NN; j++ )
	      {
	         if( !i_stream.eof() )
		       {                                                   
                  Int_t value;
		          i_stream >> value;
                  vTelescopeNeighbors[i_telType].push_back(value);
                  cout<<" "<<vTelescopeNeighbors[i_telType][j];
                }
		           else cout<<"Something is wrong with the neighbors of telescope "<<i_telType<<endl;;
           }
           cout<<endl;
	    }



}
                
