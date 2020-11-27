#define WZAnalysis_cxx
// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.

#include "WZAnalysis.h"
#include "WZhistograms.h"
#include <iostream>
#include <cstring>
#include <string>

#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLorentzVector.h>

string name;

void WZAnalysis::Begin(TTree * )
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
}

void WZAnalysis::SlaveBegin(TTree * )
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  printf("Starting analysis with process option: %s \n", option.Data());

  name=option;

  define_histograms();

  FillOutputList();
}

Bool_t WZAnalysis::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fChain->GetTree()->GetEntry(entry);
  // Cutflow index
    int cut1 = 0;
    int cut2 = 0;
    int cut3 = 0;
    int cut4 = 0;
    int cut5 = 0;
    int cut6 = 0;
    int cut7 = 0;
    int cut8 = 0;

  if(fChain->GetTree()->GetEntries()>0)
  {
    //Do analysis

    //SF
    Float_t scaleFactor = scaleFactor_ELE*scaleFactor_MUON*scaleFactor_TRIGGER;
    //EventW
    Float_t eventWeight = mcWeight*scaleFactor_PILEUP*scaleFactor_ZVERTEX;
    //weight = SF * EventW
    Double_t weight = scaleFactor*eventWeight;

    // Make difference between data and MC
    if (weight == 0. && channelNumber != 110090 && channelNumber != 110091) weight = 1.;

    // Missing Et of the event in GeV
    Float_t missingEt = met_et/1000.;
      
   

    // Preselection cut for electron/muon trigger, Good Run List, and good vertex
    if(trigE || trigM)
    {
        cut1++;
        hist_cutflow->Fill(1);
      if(passGRL)
      {
          cut2++;
          hist_cutflow->Fill(2);
        if(hasGoodVertex)
        {
            cut3++;
            hist_cutflow->Fill(3);
          //Find the good leptons
          /*int goodlep_index = 0;
          int goodlep_n = 0;

          for(int i=0; i<lep_n; i++)
          {
            if(lep_pt[i]>25000. && (lep_ptcone30[i]/lep_pt[i]) < 0.15 && (lep_etcone20[i]/lep_pt[i]) < 0.15)
              {
          // electron selection in fiducial region excluding candidates in the transition region between the barrel and endcap     electromagnetic calorimeters
            if ( lep_type[i]==11 && TMath::Abs(lep_eta[i]) < 2.47 && ( TMath::Abs(lep_eta[i]) < 1.37 || TMath::Abs(lep_eta[i]) > 1.52 ) ) {
              goodlep_n++;
              goodlep_index = i;
              }
           
              }
              
          }*/
          cut4++;
          hist_cutflow->Fill(4);
          //std::cout<<goodlep_n<<std::endl;
          int goodlep1 = 0;
          int goodlep2 = 0;
          int goodlep3 = 0;
          float mass_lep1 = 0;
          float mass_lep2 = 0;
          float mass_lep3 = 0;
          TLorentzVector lep1 = TLorentzVector ();
          TLorentzVector lep2 = TLorentzVector ();
          TLorentzVector lep3 = TLorentzVector ();
          int a = 0;
          int b = 0;
          int c = 0;
          
          cut5++;
          hist_cutflow->Fill(5);
          //Zero cut
         // if(goodlep_n==3)
         // {
              cut6++;
             hist_cutflow->Fill(6);
             for(int i = 0; i < lep_n; i++){
                     for(int j = i+1; j < lep_n; j++){
                             for(int k = j+1; k < lep_n; k++){
                                 if(lep_pt[i]>25000. && (lep_ptcone30[i]/lep_pt[i]) < 0.15 && (lep_etcone20[i]/lep_pt[i]) < 0.15){
          // electron selection in fiducial region excluding candidates in the transition region between the barrel and endcap     electromagnetic calorimeters
            if ( lep_type[i]==11 && TMath::Abs(lep_eta[i]) < 2.47 && ( TMath::Abs(lep_eta[i]) < 1.37 || TMath::Abs(lep_eta[i]) > 1.52 ) ) {
                goodlep1 = i;
                lep1.SetPtEtaPhiE(lep_pt[goodlep1], lep_eta[goodlep1], lep_phi[goodlep1], lep_E[goodlep1]);   }   }                 
                                 
                                 if(lep_pt[j]>25000. && (lep_ptcone30[j]/lep_pt[j]) < 0.15 && (lep_etcone20[j]/lep_pt[j]) < 0.15){
          // electron selection in fiducial region excluding candidates in the transition region between the barrel and endcap     electromagnetic calorimeters
            if ( lep_type[j]==11 && TMath::Abs(lep_eta[j]) < 2.47 && ( TMath::Abs(lep_eta[j]) < 1.37 || TMath::Abs(lep_eta[j]) > 1.52 ) ) {
                goodlep2 = j;
                lep2.SetPtEtaPhiE(lep_pt[goodlep2], lep_eta[goodlep2], lep_phi[goodlep2], lep_E[goodlep2]);   }   } 
                                 if(lep_pt[k]>25000. && (lep_ptcone30[k]/lep_pt[k]) < 0.15 && (lep_etcone20[k]/lep_pt[k]) < 0.15){
          // electron selection in fiducial region excluding candidates in the transition region between the barrel and endcap     electromagnetic calorimeters
            if ( lep_type[k]==11 && TMath::Abs(lep_eta[k]) < 2.47 && ( TMath::Abs(lep_eta[k]) < 1.37 || TMath::Abs(lep_eta[k]) > 1.52 ) ) {
                goodlep3 = k;
                lep3.SetPtEtaPhiE(lep_pt[goodlep3], lep_eta[goodlep3], lep_phi[goodlep3], lep_E[goodlep3]);   }   } 
                  
                                  
                                   //GoodLelectrons mass
                                   mass_lep1 = lep1.M()/1000.;
                                   mass_lep2 = lep2.M()/1000.;
                                   mass_lep3 = lep3.M()/1000.;
                                   
                                   
                           if(mass_lep1 == mass_lep2 && lep_charge[goodlep1] == -lep_charge[goodlep2]){a=goodlep1, b=goodlep2, c=goodlep3;}
                           if(mass_lep1 == mass_lep3 && lep_charge[goodlep1] == -lep_charge[goodlep3]){a=goodlep1, b=goodlep3, c=goodlep2;}
                           if(mass_lep2 == mass_lep3 && lep_charge[goodlep2] == -lep_charge[goodlep3]){a=goodlep3, b=goodlep2, c=goodlep1;}
                               }
                      }
             }
          
            cut7++;
            hist_cutflow->Fill(7);
            // Z Selection
                       //TLorentzVector definitions
                       TLorentzVector l1 = TLorentzVector ();
                       l1.SetPtEtaPhiE(lep_pt[a], lep_eta[a], lep_phi[a], lep_E[a]);
                       TLorentzVector l2 = TLorentzVector ();
                       l2.SetPtEtaPhiE(lep_pt[b], lep_eta[b], lep_phi[b], lep_E[b]);
                       
                       float mass_compare = (l1+l2).M()/1000.;
                       
                                                                                                
                        
            // MTW                                                                                 
              
            // TLorentzVector definitions
            TLorentzVector l3 = TLorentzVector();
            TLorentzVector MeT = TLorentzVector();
            //TLorentzVector lep_w_MeT = TLorentzVector();

            lep3.SetPtEtaPhiE(lep_pt[c], lep_eta[c], lep_phi[c],lep_E[c]);
            MeT.SetPtEtaPhiE(met_et, 0, met_phi , met_et);

           
            float MTW = sqrt(2*l3.Pt()*MeT.Et()*(1-cos(l3.DeltaPhi(MeT))));

                                                                                                
                   if(MTW > 30000.)
                      {
                       cut8++;
                       hist_cutflow->Fill(8);
                       
                      // Fill histograms
                                                                                                
                       if(TMath::Abs(mass_compare-91.18)<101.18 && TMath::Abs(mass_compare-91.18)>81.18 ){
                        hist_Zmass->Fill(mass_compare, weight);}
                       
                       hist_Wmass->Fill(MTW, weight);
                       
                  
                         
                     // Plot all the distributions
                    double names_of_global_variable[]={missingEt, vxp_z, (double)pvxp_n, MTW};
                    double names_of_leadlep_variable[]={l3.Pt()/1000., l3.Eta(), l3.E()/1000.,l3.Phi(),lep_charge[a], (double)lep_type[a], lep_ptcone30[a], lep_etcone20[a], lep_z0[a], lep_trackd0pvunbiased[a], lep_charge[b], (double)lep_type[b], lep_ptcone30[b], lep_etcone20[b], lep_z0[b], lep_trackd0pvunbiased[b]};
                    //double names_of_jet_variable[]={(double)jet_n, jet_pt[0]/1000., jet_eta[0], jet_m[0]/1000., jet_jvf[0], jet_MV1[0]};
                    
                    TString histonames_of_global_variable[]={"hist_etmiss","hist_vxp_z","hist_pvxp_n", "hist_mass"};
                    TString histonames_of_leadlep_variable[]={"hist_leadleptpt", "hist_leadlepteta","hist_leadleptE","hist_leadleptphi","hist_leadleptch","hist_leadleptID","hist_leadlept_ptc","hist_leadleptetc","hist_leadlepz0","hist_leadlepd0","hist_leadleptch2","hist_leadleptID2","hist_leadlept_ptc2","hist_leadleptetc2","hist_leadlepz02","hist_leadlepd02"};
                    //TString histonames_of_jet_variable[]={"hist_n_jets","hist_leadjet_pt","hist_leadjet_eta","hist_leadjet_m", "hist_leadjet_jvf", "hist_leadjet_MV1"};
                    
                    int length_global = sizeof(names_of_global_variable)/sizeof(names_of_global_variable[0]);
                    int length_leadlep = sizeof(names_of_leadlep_variable)/sizeof(names_of_leadlep_variable[0]);
                    //int length_leadjet = sizeof(names_of_jet_variable)/sizeof(names_of_jet_variable[0]);
                    
                    for (int i=0; i<length_global; i++)
                      {
                        FillHistogramsGlobal( names_of_global_variable[i], weight, histonames_of_global_variable[i]);
                      }
                    for (int i=0; i<length_leadlep; i++)
                      {
                        FillHistogramsLeadlept( names_of_leadlep_variable[i], weight, histonames_of_leadlep_variable[i]);
                      }
                   /* for (int i=0; i<length_leadjet; i++)
                      {
                        FillHistogramsLeadJet( names_of_jet_variable[i], weight, histonames_of_jet_variable[i]);
                      }*/
                      }
                 
          //}
        }
          }
        }
      }
  //std::cout << "Done!" << std::endl;
  //std::cout << "Cut1:" << cut1 << std::endl;
  //std::cout << "Cut2:" << cut2 << std::endl;
  //std::cout << "Cut3:" << cut3 << std::endl;
  //std::cout << "Cut4:" << cut4 << std::endl;
  //std::cout << "Cut5:" << cut5 << std::endl;
  //std::cout << "Cut6:" << cut6 << std::endl;
  //std::cout << "Cut7:" << cut7 << std::endl;
  return kTRUE;
}


void WZAnalysis::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

void WZAnalysis::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  name="output_WZ/"+name+".root";

  const char* filename = name.c_str();

  TFile physicsoutput_WZ(filename,"recreate");
  WriteHistograms();
  physicsoutput_WZ.Close();


}
