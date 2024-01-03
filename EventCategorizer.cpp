/**
 *  @copyright Copyright 2021 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file EventCategorizer.cpp
 */

#include "EventCategorizer.h"
#include "CalibrationTools.h"
#include "EventCategorizerTools.h"
#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>

using namespace jpet_options_tools;
using namespace std;

EventCategorizer::EventCategorizer(const char* name) : JPetUserTask(name) {}

EventCategorizer::~EventCategorizer() {}

bool EventCategorizer::init()
{
  INFO("Event categorization started.");


/*   
 * this is and example of how we can read parameter from config file
  if (isOptionSet(fParams.getOptions(), kLORPosZCutParamKey))
  {
    fLORPosZCut = getOptionAsDouble(fParams.getOptions(), kLORPosZCutParamKey);
  }
  else
  {
    WARNING(Form("No value of the %s parameter provided by the user. Using default value of %lf.", kLORPosZCutParamKey.c_str(), fLORPosZCut));
  }
*/
  // Source position

  // Reading file with constants to property tree
  if (isOptionSet(fParams.getOptions(), kConstantsFileParamKey))
  {
    boost::property_tree::read_json(getOptionAsString(fParams.getOptions(), kConstantsFileParamKey), fConstansTree);
  }

  // Input events type
  fOutputEvents = new JPetTimeWindow("JPetEvent");

  // Initialise hisotgrams
  initialiseHistograms();
  std::cout << "Hello Manish --  version(2.0)\n";

  return true;
}

bool EventCategorizer::exec()


{
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent))
  { 
     vector<JPetEvent> events;
      for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++)
      {
      const auto& event = dynamic_cast<const JPetEvent&>(timeWindow->operator[](i));    
      manishdas(event); 
    }
    saveEvents(events);
  }
  else
  {
    return false;
  }
  return true;
}

bool EventCategorizer::terminate()
{
  INFO("Event categorization completed.");
  return true;
}

void EventCategorizer::saveEvents(const vector<JPetEvent>& events)
{
  for (const auto& event : events)
  {
    fOutputEvents->add<JPetEvent>(event);
  }
}

void EventCategorizer::initialiseHistograms()
{
  auto minScinID = getParamBank().getScins().begin()->first;
  auto maxScinID = getParamBank().getScins().rbegin()->first;

  

  getStatistics().createHistogramWithAxes( new TH1D("multiplicity", "Multiplicity of all events inside the file",20,0,20 ),
		  			"Multiplicity", "Counts"
		  );
		  
  getStatistics().createHistogramWithAxes( new TH1D("raw_TOT", "Time over Threshold", 1000, 0, 9000000),
		 			"TOT [ps]", "Counts" 
		  );
  // Histograms for 2 gamama events

  getStatistics().createHistogramWithAxes( new TH1D("TOT_2multi", "Time over Threshold for multiplicity 2", 1000, 0, 9000000),
		 			"TOT [ps]", "Counts" 
		  );

  getStatistics().createHistogramWithAxes( new TH1D("Angle_annhilation_2multi","Angle between annhilation photons", 180, 0, 190),"Angle[degree]", "Counts" );

  getStatistics().createHistogramWithAxes( new TH2D("raw_XY_2multi","Annhilation point_XY_multiplicity_2", 800,-40,40, 800,-40,40 ),
		  			"X [cm]", "Y [cm]"
		  );
  
  getStatistics().createHistogramWithAxes( new TH2D("raw_XZ_2multi","Annhilation point_XZ_multiplicity_2", 800,-40,40, 500, -25, 25 ),
		  			"X [cm]", "Z [cm]"
		  );

  getStatistics().createHistogramWithAxes( new TH1D("TOT_2multi_TOT_cut", "Time over Threshold (6.25<TOT<7.15) ", 1000, 0, 9000000),
		 			"TOT [ps]", "Counts" 
		  );

  getStatistics().createHistogramWithAxes( new TH2D("XY_2multi_TOT_cut","Annhilation point_XY_multiplicity_2 (6.25<TOT<7.15)", 800,-40,40, 800,-40,40 ),
		  			"X [cm]", "Y [cm]"
		  );
  getStatistics().createHistogramWithAxes( new TH2D("XZ_2multi_TOT_cut","Annhilation point_XZ_multiplicity_2 (6.25<TOT<7.15)", 800,-40,40, 500,-25,25 ),
		  			"X [cm]", "Z [cm]"
		  );

  getStatistics().createHistogramWithAxes( new TH1D("Angle_annhilation_2multi_TOT_cut","Angle between annhilation photons(6.25<TOT<7.15)", 180, 0, 190),"Angle[degree]", "Counts" );

  getStatistics().createHistogramWithAxes( new TH2D("XY_2multi_TOT_ST_angle_cut","Annhilation point_XY_multiplicity_2 (6.25<TOT<7.15)(ST)(120<Angle<180)", 800,-40,40, 800,-40,40 ),
		  			"X [cm]", "Y [cm]"
		  );
  getStatistics().createHistogramWithAxes( new TH2D("XZ_2multi_TOT_ST_angle_cut","Annhilation point_XZ_multiplicity_2 (6.25<TOT<7.15)(ST)(120<Angle<180)", 800,-40,40, 500,-25,25 ),
		  			"X [cm]", "Z [cm]"
		  );

  getStatistics().createHistogramWithAxes( new TH1D("Angle_annhilation_2multi_TOT_ST_angle_cut","Angle between annhilation photons(6.25<TOT<7.15)(ST)(120<Angle<180)", 180, 0, 190),"Angle[degree]", "Counts" );

getStatistics().createHistogramWithAxes( new TH2D("XY_2multi_TOT_ST_angle_cut_zpos<2","Annhilation point_XY_multiplicity_2 (6.25<TOT<7.15)(ST)(120<Angle<180)(0<|z1-z2|<2)", 800,-40,40, 800,-40,40 ),
		  			"X [cm]", "Y [cm]"
		  );
  getStatistics().createHistogramWithAxes( new TH2D("XZ_2multi_TOT_ST_angle_cut_zpos<2","Annhilation point_XZ_multiplicity_2 (6.25<TOT<7.15)(ST)(120<Angle<180)(0<|z1-z2|<2)", 800,-40,40, 500,-25,25 ),
		  			"X [cm]", "Z [cm]"
		  );

getStatistics().createHistogramWithAxes(new TH1D("lor_2multi_TOT_ST_angle_cut_zpos<2", "LOR_2multi_TOT_ST_angle_cut_0<zpos<2", 2000, 0, 100), 
                                          "LOR [cm]", "Counts");
getStatistics().createHistogramWithAxes( new TH2D("XY_2multi_TOT_ST_angle_cut_zpos<8","Annhilation point_XY_multiplicity_2 (6.25<TOT<7.15)(ST)(120<Angle<180)(2<|z1-z2|<8)", 800,-40,40, 800,-40,40 ),
		  			"X [cm]", "Y [cm]"
		  );
  getStatistics().createHistogramWithAxes( new TH2D("XZ_2multi_TOT_ST_angle_cut_zpos<8","Annhilation point_XZ_multiplicity_2 (6.25<TOT<7.15)(ST)(120<Angle<180)(2<|z1-z2|<8)", 800,-40,40, 500,-25,25 ),
		  			"X [cm]", "Z [cm]"
		  );
getStatistics().createHistogramWithAxes(new TH1D("lor_2multi_TOT_ST_angle_cut_zpos<8", "LOR_2multi_TOT_ST_angle_cut_2<zpos<8", 2000, 0, 100), 
                                          "LOR [cm]", "Counts");


getStatistics().createHistogramWithAxes( new TH2D("XY_2multi_TOT_ST","Annhilation point_XY_multiplicity_2 (6.25<TOT<7.15)(ST)", 800,-40,40, 800,-40,40 ),
		  			"X [cm]", "Y [cm]"
		  );
  getStatistics().createHistogramWithAxes( new TH2D("XZ_2multi_TOT_ST","Annhilation point_XZ_multiplicity_2 (6.25<TOT<7.15)(ST)", 800,-40,40, 500,-25,25 ),
		  			"X [cm]", "Z [cm]"
		  );

  getStatistics().createHistogramWithAxes( new TH1D("Angle_annhilation_2multi_TOT_ST","Angle between annhilation photons(6.25<TOT<7.15)(ST)", 180, 0, 190),"Angle[degree]", "Counts" );


  
  getStatistics().createHistogramWithAxes( new TH2D("XY_2multi_TOT_and_angle_cut2","Annhilation point_XY_multiplicity_2 (6.25<TOT<7.15)(150<Angle<180)", 800,-40,40, 800,-40,40 ),
		  			"X [cm]", "Y [cm]"
		  );
  getStatistics().createHistogramWithAxes( new TH2D("XZ_2multi_TOT_and_angle_cut2","Annhilation point_XZ_multiplicity_2 (6.25<TOT<7.15)(150<Angle<180)", 800,-40,40, 500,-25,25 ),
		  			"X [cm]", "Z [cm]"
		  );

  getStatistics().createHistogramWithAxes( new TH1D("Angle_annhilation_2multi_TOT_and_angle_cut2","Angle between annhilation photons(6.25<TOT<7.15)(150<Angle<180)", 180, 0, 190),"Angle[degree]", "Counts" );

getStatistics().createHistogramWithAxes( new TH2D("XY_2multi_TOT_and_angle_cut3","Annhilation point_XY_multiplicity_2 (6.25<TOT<7.15)(130<Angle<150)", 800,-40,40, 800,-40,40 ),
		  			"X [cm]", "Y [cm]"
		  );
  getStatistics().createHistogramWithAxes( new TH2D("XZ_2multi_TOT_and_angle_cut3","Annhilation point_XZ_multiplicity_2 (6.25<TOT<7.15)(130<Angle<150)", 800,-40,40, 500,-25,25 ),
		  			"X [cm]", "Z [cm]"
		  );

  getStatistics().createHistogramWithAxes( new TH1D("Angle_annhilation_2multi_TOT_and_angle_cut3","Angle between annhilation photons(6.25<TOT<7.15)(130<Angle<150)", 180, 0, 190),"Angle[degree]", "Counts" );


  getStatistics().createHistogramWithAxes(new TH1D("ScatterTest", "Test if Annihilation Hits come from scattering", 2200, -19.05, 200.95), 
                                          "Scatter Test [ns]", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D("ScatterTestDistance", "Test if Annihilation Hits come from scattering", 2000, -100.05, 99.95), 
                                          "Scatter Test [cm]", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D("ScatterTest_2multi_TOT", "Test if Annihilation Hits come from scattering (6.25<TOT<7.15)", 2200, -19.05, 200.95), 
                                          "Scatter Test [ns]", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D("ScatterTestDistance_2multi_TOT", "Test if Annihilation Hits come from scattering (6.25<TOT<7.15)", 2000, -100.05, 99.95), 
                                          "Scatter Test [cm]", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D("ScatterTest_2multi_TOT_ST_angle_cut", "Test if Annihilation Hits come from scattering (6.25<TOT<7.15)(ST)(120<Angle<180)", 2200, -19.05, 200.95), 
                                          "Scatter Test [ns]", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D("ScatterTestDistance_2multi_TOT_ST_angle_cut", "Test if Annihilation Hits come from scattering (6.25<TOT<7.15)(ST)(120<Angle<180)", 2000, -100.05, 99.95), 
                                          "Scatter Test [cm]", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D("ScatterTest_2multi_TOT_and_angle_cut2", "Test if Annihilation Hits come from scattering (6.25<TOT<7.15)(150<Angle<180)", 2200, -19.05, 200.95), 
                                          "Scatter Test [ns]", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D("ScatterTestDistance_2multi_TOT_and_angle_cut2", "Test if Annihilation Hits come from scattering (6.25<TOT<7.15)(150<Angle<180)", 2000, -100.05, 99.95), 
                                          "Scatter Test [cm]", "Counts");
getStatistics().createHistogramWithAxes(new TH1D("ScatterTest_2multi_TOT_and_angle_cut3", "Test if Annihilation Hits come from scattering (6.25<TOT<7.15)(130<Angle<150)", 2200, -19.05, 200.95), 
                                          "Scatter Test [ns]", "Counts");
  getStatistics().createHistogramWithAxes(new TH1D("ScatterTestDistance_2multi_TOT_and_angle_cut3", "Test if Annihilation Hits come from scattering (6.25<TOT<7.15)(130<Angle<150)", 2000, -100.05, 99.95), 
                                          "Scatter Test [cm]", "Counts");
getStatistics().createHistogramWithAxes(new TH1D("lor_value_2multi", "LOR for 2 multi", 2000, 0, 100), 
                                          "LOR [cm]", "Counts");
getStatistics().createHistogramWithAxes(new TH1D("z1-z2_2multi", "|Z1-Z2| for 2 multi", 2000, -100, 100), 
                                          "LOR [cm]", "Counts");

getStatistics().createHistogramWithAxes( new TH1D("TOT_3multi_TOT_cut", "Time over Threshold (6.25<TOT<7.15) ", 1000, 0, 9000000),
		 			"TOT [ps]", "Counts" 
		  );
getStatistics().createHistogramWithAxes( new TH1D("ScinID_2multi", "Number of hits per scintillator ", 1000, 200,512 ),
		 			"SciID", "Counts" 
		  );
getStatistics().createHistogramWithAxes( new TH1D("ScinID_2multi_TOT", "Number of hits per scintillator ", 1000, 200,512 ),
		 			"SciID", "Counts" 
		  );
getStatistics().createHistogramWithAxes( new TH1D("ScinID_2multi_TOT_ST", "Number of hits per scintillator ", 1000, 200,512 ),
		 			"SciID", "Counts" 
		  );
getStatistics().createHistogramWithAxes( new TH1D("ScinID_2multi_TOT_ST_angle_cut", "Number of hits per scintillator ", 1000, 200,512 ),
		 			"SciID", "Counts" 
		  );
getStatistics().createHistogramWithAxes( new TH1D("ScinID_2multi_TOT_ST_angle_cut_zpos<2", "Number of hits per scintillator ", 1000, 200,512 ),
		 			"SciID", "Counts" 
		  );
getStatistics().createHistogramWithAxes( new TH1D("ScinID_2multi_TOT_ST_angle_cut_zpos<8", "Number of hits per scintillator ", 1000, 200,512 ),
		 			"SciID", "Counts" 
		  );
getStatistics().createHistogramWithAxes( new TH1D("ScinID_2multi_TOT_ST_angle_cut_zpos<14", "Number of hits per scintillator ", 1000, 200,512 ),
		 			"SciID", "Counts" 
		  );
getStatistics().createHistogramWithAxes( new TH1D("ScinID_2multi_TOT_ST_angle_cut_zpos<20", "Number of hits per scintillator ", 1000, 200,512 ),
		 			"SciID", "Counts" 
		  );

getStatistics().createHistogramWithAxes( new TH2D("XY_2multi_TOT_ST_angle_cut_zpos<14","Annhilation point_XY_multiplicity_2 (6.25<TOT<7.15)(ST)(120<Angle<180)(8<|z1-z2|<14)", 800,-40,40, 800,-40,40 ),
		  			"X [cm]", "Y [cm]"
		  );
  getStatistics().createHistogramWithAxes( new TH2D("XZ_2multi_TOT_ST_angle_cut_zpos<14","Annhilation point_XZ_multiplicity_2 (6.25<TOT<7.15)(ST)(120<Angle<180)(8<|z1-z2|<14)", 800,-40,40, 500,-25,25 ),
		  			"X [cm]", "Z [cm]"
		  );

getStatistics().createHistogramWithAxes(new TH1D("lor_2multi_TOT_ST_angle_cut_zpos<14", "LOR_2multi_TOT_ST_angle_cut_8<zpos<14", 2000, 0, 100), 
                                          "LOR [cm]", "Counts");
getStatistics().createHistogramWithAxes( new TH2D("XY_2multi_TOT_ST_angle_cut_zpos<20","Annhilation point_XY_multiplicity_2 (6.25<TOT<7.15)(ST)(120<Angle<180)(14<|z1-z2|<20)", 800,-40,40, 800,-40,40 ),
		  			"X [cm]", "Y [cm]"
		  );
  getStatistics().createHistogramWithAxes( new TH2D("XZ_2multi_TOT_ST_angle_cut_zpos<20","Annhilation point_XZ_multiplicity_2 (6.25<TOT<7.15)(ST)(120<Angle<180)(14<|z1-z2|<20)", 800,-40,40, 500,-25,25 ),
		  			"X [cm]", "Z [cm]"
		  );

getStatistics().createHistogramWithAxes(new TH1D("lor_2multi_TOT_ST_angle_cut_zpos<20", "LOR_2multi_TOT_ST_angle_cut_14<zpos<20", 2000, 0, 100), 
                                          "LOR [cm]", "Counts");



// histogram for 3 hits

 getStatistics().createHistogramWithAxes( new TH2D("raw_XY_Multiplicity_3","Emission points in the XY plane for multiplicity 3", 800,-40,40, 800,-40,40 ),
		  			"X [cm]", "Y [cm]"
		  );

  getStatistics().createHistogramWithAxes( new TH2D("raw_XZ_Multiplicity_3","Emission points in the XZ plane for multiplicity 3", 800,-40,40, 500,-25,25 ),
		  			"X [cm]", "Z [cm]"
		  );
 
   getStatistics().createHistogramWithAxes( new TH2D("XY_3multi_TOT_and_angle_cut","Annhilation point_XY_multiplicity_3 (6.25<TOT<7.15)(120<Angle<180)", 800,-40,40, 500,-40,40 ),
		  			"X [cm]", "Z [cm]"
		  );
  getStatistics().createHistogramWithAxes( new TH2D("XZ_3multi_TOT_and_angle_cut","Annhilation point_XY_multiplicity_3 (6.25<TOT<7.15)(120<Angle<180)", 800,-40,40, 500,-25,25 ),
		  			"X [cm]", "Z [cm]"
		  );

  getStatistics().createHistogramWithAxes( new TH1D("Multiplicity3_raw_TOT","Time over Threshold for multiplicity 3", 1000, 0, 10000000),
		  			"TOT [ps]", "Counts"
		  );
  
  getStatistics().createHistogramWithAxes( new TH1D("Lifetime_naive","Lifetime based on time diff between prompt and first annhilation, no corrections and cuts", 240, -2, 10),		  			"Lifetime [ns]", "Counts"
		  );
  getStatistics().createHistogramWithAxes( new TH1D("Lifetime_totCuts","Lifetime based on time diff between prompt and first annhilation, no corrections, TOT cut 3-5 on both annihilation gammas and larger than 5 on prompt", 240, -2, 10),		  			"Lifetime [ns]", "Counts"
		  );
  
  getStatistics().createHistogramWithAxes( new TH1D("Lifetime_TOFcorr","Lifetime based on time diff between prompt and mean time of emission, TOF corrected, TOT cut 3-5 on both annihilation gammas and larger than 5 on prompt", 240, -2, 10),		  			"Lifetime [ns]", "Counts"
		  ); 

  getStatistics().createHistogramWithAxes( new TH1D("Lifetime_totCuts_angle","Lifetime final no correction", 240, -2, 10),		  			"Lifetime [ns]", "Counts"
		  );
  getStatistics().createHistogramWithAxes( new TH1D("Lifetime_TOFcorr_anglecut","Lifetime final, TOF corrected", 240, -2, 10),		  			"Lifetime [ns]", "Counts"
		  ); 

getStatistics().createHistogramWithAxes( new TH1D("Lifetime","Lifetime final", 240, -2, 10),		  			"Lifetime [ns]", "Counts"
		  ); 

getStatistics().createHistogramWithAxes( new TH3F("task1","multiplicity as a function of scintillator number and z position ", 800,200,512, 50,-25,25,10,0,10 ),
		  			"scinID", "Z [cm]", "multiplicity"
		  );
 

}




  void EventCategorizer::manishdas(const JPetEvent &event) { 
       vector<JPetEvent> events;
       vector<int> ScinID;
       vector<int> counts;
       vector<int> unique_numbers;

      int multiplicity = event.getHits().size(); 
      getStatistics().fillHistogram("multiplicity",multiplicity);

      for (int k = 0; k < multiplicity; ++k) {
           auto reconstructed_hit = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(k));
              
             int scinIDx = reconstructed_hit->getScin().getID();
             double zposii = reconstructed_hit->getPosZ();
             getStatistics().fillHistogram("task1",scinIDx,zposii,multiplicity);

             double HitEnergy = reconstructed_hit->getToT();
             getStatistics().fillHistogram("raw_TOT", HitEnergy);


//Two_Gamma_Event 
      
     
      if( multiplicity == 2 )
        {


         auto hit1 = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(0));
         auto hit2 = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(1));
         
         //Angle_between two photons  
         TVector3 Hit1Pos = hit1->getPos();
         TVector3 Hit2Pos = hit2->getPos();
         double theta1 = TMath::RadToDeg() * Hit1Pos.Angle(Hit2Pos);

         //second method to draw the angle	
         //double theta =  CalcAngle3D(*hit1,*hit2);
          
         getStatistics().fillHistogram("Angle_annhilation_2multi", theta1);        
  

         TVector3 annPoint = EventCategorizerTools::calculateAnnihilationPoint( hit1, hit2);
      
         getStatistics().fillHistogram("raw_XY_2multi", annPoint.X(), annPoint.Y() );
         getStatistics().fillHistogram("raw_XZ_2multi", annPoint.X(), annPoint.Z() );
         int scinID1 = hit1->getScin().getID();
         int scinID2 = hit2->getScin().getID();
         getStatistics().fillHistogram("ScinID_2multi", scinID1);
         getStatistics().fillHistogram("ScinID_2multi", scinID2);
         
         
         ScinID.push_back(scinID1);
         ScinID.push_back(scinID2);
         


         //TOT_FOR_MULTIPLICITY_2
         double hit1Energy = hit1->getToT();
         double hit2Energy = hit2->getToT();
	 getStatistics().fillHistogram("TOT_2multi", hit1Energy);
      	 getStatistics().fillHistogram("TOT_2multi", hit2Energy);

         double ST = ScatterTest(*hit1, *hit2); 
         getStatistics().fillHistogram("ScatterTest", ST); 
         double STD = ScatterTestD(*hit1, *hit2); 
         getStatistics().fillHistogram("ScatterTestDistance", STD);

        double z1 =  hit1->getPosZ();
        double z2 =  hit2->getPosZ();
         
        double modz1 =  fabs(hit1->getPosZ());
        double modz2 = fabs(hit2->getPosZ());
        
        double zpos = fabs(z1- z2);
        getStatistics().fillHistogram("z1-z2_2multi", zpos);
        //DLOR_CALCULATION
        double distance = fabs(sqrt(pow((hit2->getPosX() - hit1->getPosX()),2)+ pow((hit2->getPosY() - hit1->getPosY()),2)+ pow((hit2->getPosZ() - hit1->getPosZ()),2)));
        getStatistics().fillHistogram("lor_value_2multi", distance);




	if ( ( hit1Energy >6250000 && hit1Energy < 7150000) && (hit2Energy >6250000 && hit2Energy < 7150000) )
	   {
                double hit1Energy = hit1->getToT();
                double hit2Energy = hit2->getToT();
            
		TVector3 annPoint = EventCategorizerTools::calculateAnnihilationPoint( hit1, hit2);
		double theta2 = TMath::RadToDeg() * Hit1Pos.Angle(Hit2Pos);
		getStatistics().fillHistogram("XY_2multi_TOT_cut", annPoint.X(), annPoint.Y() );
		getStatistics().fillHistogram("XZ_2multi_TOT_cut", annPoint.X(), annPoint.Z() );
    		events.push_back(event);
                getStatistics().fillHistogram("Angle_annhilation_2multi_TOT_cut", theta2);
                getStatistics().fillHistogram("TOT_2multi_TOT_cut", hit1Energy);
      	        getStatistics().fillHistogram("TOT_2multi_TOT_cut", hit2Energy);

                double ST = ScatterTest(*hit1, *hit2); 
                getStatistics().fillHistogram("ScatterTest_2multi_TOT", ST); 
                double STD = ScatterTestD(*hit1, *hit2); 
                getStatistics().fillHistogram("ScatterTestDistance_2multi_TOT", STD);
                getStatistics().fillHistogram("ScinID_2multi_TOT", scinID1);
                getStatistics().fillHistogram("ScinID_2multi_TOT", scinID2);
       

       if (STD<-20){ 

                    getStatistics().fillHistogram("XY_2multi_TOT_ST", annPoint.X(), annPoint.Y() );
		    getStatistics().fillHistogram("XZ_2multi_TOT_ST", annPoint.X(), annPoint.Z() );
                    getStatistics().fillHistogram("Angle_annhilation_2multi_TOT_ST", theta2);
                    getStatistics().fillHistogram("ScinID_2multi_TOT_ST", scinID1);
                    getStatistics().fillHistogram("ScinID_2multi_TOT_ST", scinID2);

                    

       if (theta2 > 120 && theta2 <180){

                TVector3 annPoint = EventCategorizerTools::calculateAnnihilationPoint( hit1, hit2);
		double theta2 = TMath::RadToDeg() * Hit1Pos.Angle(Hit2Pos);
		getStatistics().fillHistogram("XY_2multi_TOT_ST_angle_cut", annPoint.X(), annPoint.Y() );
		getStatistics().fillHistogram("XZ_2multi_TOT_ST_angle_cut", annPoint.X(), annPoint.Z() );
    		events.push_back(event);
                getStatistics().fillHistogram("Angle_annhilation_2multi_TOT_ST_angle_cut", theta2);
                 
                double ST = ScatterTest(*hit1, *hit2); 
                getStatistics().fillHistogram("ScatterTest_2multi_TOT_ST_angle_cut", ST); 
                double STD = ScatterTestD(*hit1, *hit2); 
                getStatistics().fillHistogram("ScatterTestDistance_2multi_TOT_ST_angle_cut", STD); 
                getStatistics().fillHistogram("ScinID_2multi_TOT_ST_angle_cut", scinID1);
                getStatistics().fillHistogram("ScinID_2multi_TOT_ST_angle_cut", scinID2);

          if(zpos < 2 && zpos >0 )
                      { 
                getStatistics().fillHistogram("XY_2multi_TOT_ST_angle_cut_zpos<2", annPoint.X(), annPoint.Y() );
		getStatistics().fillHistogram("XZ_2multi_TOT_ST_angle_cut_zpos<2", annPoint.X(), annPoint.Z() );
                getStatistics().fillHistogram("lor_2multi_TOT_ST_angle_cut_zpos<2", distance);
                getStatistics().fillHistogram("ScinID_2multi_TOT_ST_angle_cut_zpos<2", scinID1);
                getStatistics().fillHistogram("ScinID_2multi_TOT_ST_angle_cut_zpos<2", scinID2);
                
                       }

         
          if(zpos > 2 && zpos < 8)
                      { 
                getStatistics().fillHistogram("XY_2multi_TOT_ST_angle_cut_zpos<8", annPoint.X(), annPoint.Y() );
		getStatistics().fillHistogram("XZ_2multi_TOT_ST_angle_cut_zpos<8", annPoint.X(), annPoint.Z() );
                getStatistics().fillHistogram("lor_2multi_TOT_ST_angle_cut_zpos<8", distance);
                getStatistics().fillHistogram("ScinID_2multi_TOT_ST_angle_cut_zpos<8", scinID1);
                getStatistics().fillHistogram("ScinID_2multi_TOT_ST_angle_cut_zpos<8", scinID2);
                
                       }     
         if( zpos > 8 && zpos < 14)
                      { 
                getStatistics().fillHistogram("XY_2multi_TOT_ST_angle_cut_zpos<14", annPoint.X(), annPoint.Y() );
		getStatistics().fillHistogram("XZ_2multi_TOT_ST_angle_cut_zpos<14", annPoint.X(), annPoint.Z() );
                getStatistics().fillHistogram("lor_2multi_TOT_ST_angle_cut_zpos<14", distance);
                getStatistics().fillHistogram("ScinID_2multi_TOT_ST_angle_cut_zpos<14", scinID1);
                getStatistics().fillHistogram("ScinID_2multi_TOT_ST_angle_cut_zpos<14", scinID2);
                
                       }     
         if(zpos > 14 && zpos < 20)
                      { 
                getStatistics().fillHistogram("XY_2multi_TOT_ST_angle_cut_zpos<20", annPoint.X(), annPoint.Y() );
		getStatistics().fillHistogram("XZ_2multi_TOT_ST_angle_cut_zpos<20", annPoint.X(), annPoint.Z() );
                getStatistics().fillHistogram("lor_2multi_TOT_ST_angle_cut_zpos<20", distance);
                getStatistics().fillHistogram("ScinID_2multi_TOT_ST_angle_cut_zpos<20", scinID1);
                getStatistics().fillHistogram("ScinID_2multi_TOT_ST_angle_cut_zpos<20", scinID2);
                
                       }          

                         if (theta2 > 150 && theta2 <180){

                             TVector3 annPoint = EventCategorizerTools::calculateAnnihilationPoint( hit1, hit2);
		              double theta2 = TMath::RadToDeg() * Hit1Pos.Angle(Hit2Pos);
		              getStatistics().fillHistogram("XY_2multi_TOT_and_angle_cut2", annPoint.X(), annPoint.Y() );
		              getStatistics().fillHistogram("XZ_2multi_TOT_and_angle_cut2", annPoint.X(), annPoint.Z() );
    		              events.push_back(event);
                              getStatistics().fillHistogram("Angle_annhilation_2multi_TOT_and_angle_cut2", theta2);
                    
                              double ST = ScatterTest(*hit1, *hit2); 
                              getStatistics().fillHistogram("ScatterTest_2multi_TOT_and_angle_cut2", ST); 
                              double STD = ScatterTestD(*hit1, *hit2); 
                              getStatistics().fillHistogram("ScatterTestDistance_2multi_TOT_and_angle_cut2", STD);                  
              }               
                         if (theta2 > 130 && theta2 <150){

                             TVector3 annPoint = EventCategorizerTools::calculateAnnihilationPoint( hit1, hit2);
		              double theta2 = TMath::RadToDeg() * Hit1Pos.Angle(Hit2Pos);
		              getStatistics().fillHistogram("XY_2multi_TOT_and_angle_cut3", annPoint.X(), annPoint.Y() );
		              getStatistics().fillHistogram("XZ_2multi_TOT_and_angle_cut3", annPoint.X(), annPoint.Z() );
    		              events.push_back(event);
                              getStatistics().fillHistogram("Angle_annhilation_2multi_TOT_and_angle_cut3", theta2);
                    
                              double ST = ScatterTest(*hit1, *hit2); 
                              getStatistics().fillHistogram("ScatterTest_2multi_TOT_and_angle_cut3", ST); 
                              double STD = ScatterTestD(*hit1, *hit2); 
                              getStatistics().fillHistogram("ScatterTestDistance_2multi_TOT_and_angle_cut3", STD);                  
              }               
              }   
            }
          

	} }

    





}









		

      if( multiplicity == 3)
    	{
          auto prompt = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(0));
          auto ann1 = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(1));
          auto ann2 = dynamic_cast<const JPetPhysRecoHit*>(event.getHits().at(2));
	
	  TVector3 annPoint = EventCategorizerTools::calculateAnnihilationPoint( ann1, ann2 );

	  getStatistics().fillHistogram("raw_XY_Multiplicity_3", annPoint.X(), annPoint.Y() );
	  getStatistics().fillHistogram("raw_XZ_Multiplicity_3", annPoint.X(), annPoint.Z() );

	  double ann1Energy = ann1->getToT();
	  double ann2Energy = ann2->getToT();

	  getStatistics().fillHistogram("Multiplicity3_raw_TOT", ann1Energy);
	  getStatistics().fillHistogram("Multiplicity3_raw_TOT", ann2Energy);

	  double timeStart = prompt->getTime();
          double timeStop1 = ann1->getTime();
          double timeStop2 = ann2->getTime();
          double timeStop = (timeStop1+timeStop2)/2.0;

	  getStatistics().fillHistogram("Lifetime_naive", (timeStop-timeStart)/1000.0);

	  if( prompt->getToT() > 5000000 )
  	  {
             if(ann1Energy > 3000000 && ann1Energy < 5000000 && ann2Energy >3000000 && ann2Energy < 5000000)

	     {
	        getStatistics().fillHistogram("Lifetime_totCuts", (timeStop-timeStart)/1000.0);
		TVector3 emissionPoint = EventCategorizerTools::calculateAnnihilationPoint( ann1, ann2);
		/*double distanceAnn1 = sqrt( (ann1->getPos() - emissionPoint)*(ann1->getPos() - emissionPoint) );
		double distanceAnn2 = sqrt( (ann2->getPos() - emissionPoint)*(ann2->getPos() - emissionPoint) );
		double distanceprompt = sqrt( (prompt->getPos() - emissionPoint)*(prompt->getPos() - emissionPoint) );*/
                TVector3 vec1(ann1->getPosX() - emissionPoint.X(), ann1->getPosY() - emissionPoint.Y(), ann1->getPosZ() - emissionPoint.Z());
                double distanceAnn1 = vec1.Mag();
                TVector3 vec2(ann2->getPosX() - emissionPoint.X(), ann2->getPosY() - emissionPoint.Y(), ann2->getPosZ() - emissionPoint.Z());
                double distanceAnn2 = vec2.Mag();

TVector3 vec3(prompt->getPosX() - emissionPoint.X(), prompt->getPosY() - emissionPoint.Y(), prompt->getPosZ() - emissionPoint.Z());
                double distanceprompt = vec3.Mag();

		double c = 30 /1000.0; //speed of light in cm / ps
		double timeStartCorr = timeStart - distanceprompt/c;
		
		double hitAnn1Corr = ann1->getTime() - distanceAnn1/c;
		double hitAnn2Corr = ann2->getTime() - distanceAnn2/c;

		double timeStopCorr = (hitAnn1Corr + hitAnn2Corr)/2.0;
		
	        getStatistics().fillHistogram("Lifetime_TOFcorr", (timeStopCorr-timeStartCorr)/1000.0);
    		events.push_back(event);






                //another method for lifetime
                    

                /*TVector3 AnnihilationPosition = EventCategorizerTools::calculateAnnihilationPoint(&AnniCandidates_Broader.at(0), &AnniCandidates_Broader.at(1));
   
    double Lifetime = 0.5*(NormalizeTimeToPoint(AnniCandidates_Broader.at(0), AnnihilationPosition) + 
                            NormalizeTimeToPoint(AnniCandidates_Broader.at(1), AnnihilationPosition)) - 
                            NormalizeTimeToPoint(DeexCandidates.at(0), AnnihilationPosition);*/ 
        

                //for abs
                TVector3 ann1Pos = ann1->getPos();
                TVector3 ann2Pos = ann2->getPos();
                double theta13 = TMath::RadToDeg() * ann1Pos.Angle(ann2Pos);

                getStatistics().fillHistogram("TOT_3multi_TOT_cut", ann1Energy);
      	        getStatistics().fillHistogram("TOT_3multi_TOT_cut", ann2Energy);
                if (theta13>120 && theta13<180)
                     {
                          TVector3 annPoint3 = EventCategorizerTools::calculateAnnihilationPoint( ann1, ann2);
                          getStatistics().fillHistogram("XY_3multi_TOT_and_angle_cut", annPoint3.X(), annPoint3.Y() );
		          getStatistics().fillHistogram("XZ_3multi_TOT_and_angle_cut", annPoint3.X(), annPoint3.Z() );

                          getStatistics().fillHistogram("Lifetime_totCuts_angle", (timeStop-timeStart)/1000.0);
                          getStatistics().fillHistogram("Lifetime_TOFcorr_anglecut", (timeStopCorr-timeStartCorr)/1000.0);
                          //getStatistics().fillHistogram("Lifetime", Lifetime);

                          
                      }

                



	     }
	  }

	}

}













double NormalizeTimeToPoint(JPetPhysRecoHit Hit1, TVector3 Point)
{
  TVector3 vec1(Hit1.getPosX() - Point(0), Hit1.getPosY() - Point(1), Hit1.getPosZ() - Point(2));
  double Length0 = vec1.Mag();
  return Hit1.getTime()/1000 - (Length0)/29.979246;
}

double ScatterTest(JPetPhysRecoHit Hit1, JPetPhysRecoHit Hit2)
{
  double TDiff = fabs(Hit2.getTime()/1000 - Hit1.getTime()/1000);
  TVector3 vec1(Hit2.getPosX() - Hit1.getPosX(), Hit2.getPosY() - Hit1.getPosY(), Hit2.getPosZ() - Hit1.getPosZ());
  double Distance = vec1.Mag();
  double LightVel = 29.979246; // cm/ns
  double ScattTime = Distance/LightVel;
  return TDiff - ScattTime;
}

double ScatterTestD(JPetPhysRecoHit Hit1, JPetPhysRecoHit Hit2)
{
  double TDiff = fabs(Hit2.getTime()/1000 - Hit1.getTime()/1000);
  TVector3 vec1(Hit2.getPosX() - Hit1.getPosX(), Hit2.getPosY() - Hit1.getPosY(), Hit2.getPosZ() - Hit1.getPosZ());
  double Distance = vec1.Mag();
  double LightVel = 29.979246; // cm/ns
  double ScattDistance = TDiff*LightVel;
  return ScattDistance - Distance;
}

double CalcAngle2D(JPetPhysRecoHit Hit1, JPetPhysRecoHit Hit2)
{
  double scalarProd = Hit1.getPosX()*Hit2.getPosX() + Hit1.getPosY()*Hit2.getPosY();
  double magProd = sqrt( ( pow(Hit1.getPosX(),2) + pow(Hit1.getPosY(),2) )*
                ( pow(Hit2.getPosX(),2) + pow(Hit2.getPosY(),2) ) );
  double Angle = acos( scalarProd/magProd )*180/3.14159265;
  return Angle;
}

double CalcAngle3D(JPetPhysRecoHit Hit1, JPetPhysRecoHit Hit2)
{
  double scalarProd = Hit1.getPosX()*Hit2.getPosX() + Hit1.getPosY()*Hit2.getPosY() + Hit1.getPosZ()*Hit2.getPosZ();
  double magProd = sqrt( ( pow(Hit1.getPosX(),2) + pow(Hit1.getPosY(),2) + pow(Hit1.getPosZ(),2) )*
                ( pow(Hit2.getPosX(),2) + pow(Hit2.getPosY(),2)+pow(Hit2.getPosZ(),2)));
  double Angle = acos( scalarProd/magProd )*180/3.14159265;
  return Angle;
}



double CalcPhiAngleOfAnnihilationPoint(TVector3 AnnihilationPoint)
{
  double X = AnnihilationPoint(0);
  double Radius = sqrt(AnnihilationPoint(0)*AnnihilationPoint(0) + AnnihilationPoint(1)*AnnihilationPoint(1));
  double Y = AnnihilationPoint(1);
  double PhiAngle = 0;
  if( Y >=0 )
    PhiAngle = acos(X/Radius)*180/3.14159265;
  else
    PhiAngle = 360 - acos(X/Radius)*180/3.14159265;
  return PhiAngle;
}

double CalcScattAngle(JPetPhysRecoHit Hit1, JPetPhysRecoHit Hit2, TVector3 sourcePos)
{
    double scalarProd = (Hit1.getPosX() - sourcePos(0))*(Hit2.getPosX()-Hit1.getPosX()) + 
                        (Hit1.getPosY() - sourcePos(1))*(Hit2.getPosY()-Hit1.getPosY()) + 
                        (Hit1.getPosZ() - sourcePos(2))*(Hit2.getPosZ()-Hit1.getPosZ());
    double magProd = sqrt( ( pow(Hit1.getPosX() - sourcePos(0),2)	// Pos in cm
                +pow(Hit1.getPosY() - sourcePos(1),2)
                +pow(Hit1.getPosZ() - sourcePos(2),2) )*( pow(Hit2.getPosX()-Hit1.getPosX(),2)
                +pow(Hit2.getPosY()-Hit1.getPosY(),2)
                +pow(Hit2.getPosZ()-Hit1.getPosZ(),2) ) );
    double ScattAngle = acos(scalarProd/magProd)*180/3.14159265;
    return ScattAngle;
}

double CalcScattDistance(JPetPhysRecoHit Hit1, TVector3 AnniPos)
{
    TVector3 HitPos = Hit1.getPos();
    TVector3 DistanceVec = HitPos - AnniPos;
    return DistanceVec.Mag();
}

double CalcScattDistanceXY(JPetPhysRecoHit Hit1, TVector3 AnniPos)
{
    TVector3 HitPos = Hit1.getPos();
    HitPos(2) = 0;
    TVector3 AnniPosTemp = AnniPos;
    AnniPosTemp(2) = 0;
    TVector3 DistanceVec = HitPos - AnniPosTemp;
    return DistanceVec.Mag();
}

double CorrectionOnWalkPlease(double HitEnergy)
{
  double TimeWalkAMax = 11.59826, TimeWalkY0Max = 2.47393, TimeWalkTauMax = 4919.9666;
  double TimeWalkAMin = 8.68292, TimeWalkY0Min = 1.83881, TimeWalkTauMin = 5042.34889;
  double CorrectionMax = TimeWalkAMax*exp(-HitEnergy/TimeWalkTauMax) + TimeWalkY0Max;
  double CorrectionMin = TimeWalkAMin*exp(-HitEnergy/TimeWalkTauMin) + TimeWalkY0Min;
  return (CorrectionMax + CorrectionMin)/2;
}



