#include "test.h"
#include "cmath"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"
// AUTHOR: test
//  EMAIL: test
void Test::initialize() {
  setAnalysisName("test");          
  setInformation(""
    "# test \n"
    "# test\n"
  "");
  setLuminosity(12.9*units::INVFB);      
  bookSignalRegions("1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29");
  // You can also book cutflow regions with bookCutflowRegions("CR1;CR2;..."). Note that the regions are
  //  always ordered alphabetically in the cutflow output files.

  // You should initialize any declared variables here
}

void Test::analyze() {
  // Your eventwise analysis code goes here
  // The following objects are always defined unless they are 'ignored' above. They form std::vector objects of the respective Delphes class type (except for Etmiss which is a single object)
  // All std::vector members and etmiss have the common properties PT, Eta, Phi and P4() with the latter giving access to the full ROOT TLorentzVector.
  // Within a std::vector, all members are ordered with highest pt coming first.

  // electronsLoose, electronsMedium, electronsTight   are list of electrons that passed respective efficiency and reconstruction cuts
  // muonsCombinedPlus, muonsCombined                  as above for muons
  // photonsMedium                                     as above for photons
  // jets are all reconstructed jets                   as above for jets. Note that electrons are most certainly also reconstructed as a jet -> overlap removal do avoid double counting necessary!
  // tracks, towers                                    calorimeter and tracker information. Usually not needed.
  // missingET                                         rec missing ET EXCLUDING muons.

  
  // Here is a couple of useful functions and lines:  
  //------------Phase Space Cuts (defined for jets, electronsXYZ, muonsXYZ, photonsXYZ)
  // jets = filterPhaseSpace(jets, 20., -2.8, 2.8)  // The vector 'jets' only contains jets with pt >= 20 GeV and -2.8 < eta < 2.8. This function is applicable to other particles too (electronsMedium, ... ).
  // jets = overlapRemoval(jets, electronsLoose, 0.2) Removes all jets for which there exists any electron in 'electronsLoose' with deltaR < 0.2.
  // jets = overlapRemovel(jets, 0.2) If two jets overlap within deltaR < 0.2, only the harder jet is stored.
  
  //------------Isolation Checks (defined for electronsXYZ, muonsXYZ, photonsXYZ
  //------------        For each object, if the user entered N isolation conditions, they can be
  //------------        checked individually be the second argument (running from 0 to N-1).
  // electronsMedium = filterIsolation(electronsMedium, 0)            Removes electrons that do not pass the first isolation condition entered into the AnalysisManager by the user
  // std::vector<int> flags; flags.push_back(0); flags.push_back(2);
  // electronsMedium = filterIsolation(electronsMedium, flags)        Same as above, but both the first and the third condition have to be fulfilled
  // electronsMedium = filterIsolation(electronsMedium)               Same as above, but all conditions have to be fulfilled.
  
  //-----------Flavour Tag Checks (defined for jets only)
  //----------          Tau tags "loose", "medium" or "tight" can be individually checked if the user enabled tau tagging in the AM.
  //----------          For b-tags, if N working points have been defined, the ith condition can be tested by putting i-1 into the second argument (if there is only one, the argument can be omitted)
  // if checkTauTag(jets[0], "tight") leadingJetIsTagged = True;
  // if checkBTag(jets[0], 0) leadingJetIsBTagged = True;


  //-----------Auxiliary Information
  // - Always ensure that you don't access vectors out of bounds. E.g. 'if(jets[1]->PT > 150)' should rather be if (jets.size() > 1 && jets[1]->PT > 150). 
  // - Use rand()/(RAND_MAX+1.) for random numbers between 0 and 1. The random seed is determined from system time or by the RandomSeed parameter in CheckMATE.
  // - The 'return' statement will end this function for the current event and hence should be called whenever the current event is to be vetoed.
  // - Many advanced kinematical functions like mT2 are implemented. Check the manual for more information.
  // - If you need output to be stored in other files than the cutflow/signal files we provide, check the manual for how to do this conveniently.  

  // missingET->addMuons(muonsCombined);  // Adds muons to missing ET. This should almost always be done which is why this line is not commented out.
  
  electrons=filterIsolation(electrons, 0);
  muons=filterIsolation(muons, 1);
  std::vector<Electron *> electrons_veto;
  std::vector<Muon *> muons_veto;
  std::vector<Photon *> photons_veto;
  // construct tau
  std::vector<Jet *> taus;
  std::vector<Jet *> taus_veto;
  for (int i = 0; i < jets.size(); i++)
  {
    // if (checkTauTag(jets[i], "tight")) {
    // if (checkTauTag(jets[i], "medium")) {
    if (checkTauTag(jets[i], "loose") and fabs(jets[i]->Charge) == 1)
    {
      taus.push_back(jets[i]);
    }
  }
  taus_veto = filterPhaseSpace(taus, 18., -2.3, 2.3);
  // b-jet construction
  std::vector<Jet *> bjets;
  for (int i = 0; i < jets.size(); i++)
  {
    if (checkBTag(jets[i]))
    {
      bjets.push_back(jets[i]);
    }
  }
  bjets = filterPhaseSpace(bjets, 15., -2.4, 2.4);

  // construct Ak4 jet

  std::vector<fastjet::PseudoJet> particles;
  for (int i = 0; i < true_particles.size(); i++)
  {
    if (true_particles[i]->Status == 1)
    {
      int pid = true_particles[i]->PID;

      // particles.push_back( fastjet::PseudoJet( true_particles[i]->Px, true_particles[i]->Py, true_particles[i]->Pz, true_particles[i]->E ) );
      particles.push_back(fastjet::PseudoJet(true_particles[i]->P4().Px(), true_particles[i]->P4().Py(), true_particles[i]->P4().Pz(), true_particles[i]->P4().E()));
      particles.back().set_user_index(pid);
    }
  }
  fastjet::JetDefinition jetAK4(fastjet::antikt_algorithm, 0.4);
  fastjet::ClusterSequence AK4jet(particles, jetAK4);
  std::vector<fastjet::PseudoJet> AK4 = sorted_by_pt(AK4jet.inclusive_jets());
  /////////////////////////////////////////////////////

  // another ak4 jet for charge index

  std::vector<fastjet::PseudoJet> ak4_charge;
  for (int i = 0; i < true_particles.size(); i++)
  {
    if (true_particles[i]->Status == 1)
    {
      int charge = true_particles[i]->Charge;

      // particles.push_back( fastjet::PseudoJet( true_particles[i]->Px, true_particles[i]->Py, true_particles[i]->Pz, true_particles[i]->E ) );
      ak4_charge.push_back(fastjet::PseudoJet(true_particles[i]->P4().Px(), true_particles[i]->P4().Py(), true_particles[i]->P4().Pz(), true_particles[i]->P4().E()));
      ak4_charge.back().set_user_index(charge);
    }
  }
  fastjet::JetDefinition jetAK4_charge(fastjet::antikt_algorithm, 0.4);
  fastjet::ClusterSequence AK4jet_charge(ak4_charge, jetAK4_charge);
  std::vector<fastjet::PseudoJet> AK4_charge = sorted_by_pt(AK4jet_charge.inclusive_jets());

  // constructing AK8 jets

  std::vector<fastjet::PseudoJet> particles_AK8;
  for (int i = 0; i < true_particles.size(); i++)
  {
    if (true_particles[i]->Status == 1)
    {
      int pid = true_particles[i]->PID;
      // particles.push_back( fastjet::PseudoJet( true_particles[i]->Px, true_particles[i]->Py, true_particles[i]->Pz, true_particles[i]->E ) );
      particles_AK8.push_back(fastjet::PseudoJet(true_particles[i]->P4().Px(), true_particles[i]->P4().Py(), true_particles[i]->P4().Pz(), true_particles[i]->P4().E()));
      particles_AK8.back().set_user_index(pid);
    }
  }
  fastjet::JetDefinition jetAK8(fastjet::antikt_algorithm, 0.8);
  fastjet::ClusterSequence AK8jet(particles_AK8, jetAK8);
  std::vector<fastjet::PseudoJet> AK8 = sorted_by_pt(AK8jet.inclusive_jets());
  ///////////////////////////////////////////////////////////

  electrons_veto = filterPhaseSpace(electrons, 10., -2.5, 2.5);
  muons_veto = filterPhaseSpace(muons, 10.);
  photons_veto = filterPhaseSpace(photons, 15., -2.5, 2.5);

  // filtering AK4 with pt > 30
  std::vector<fastjet::PseudoJet> AK4_filtered_pt30;
  if (AK4.size() > 0)
  {
    for (int i = 0; i < AK4.size(); i++)
    {
      if (AK4[i].pt() > 30.)
      {
        AK4_filtered_pt30.push_back(AK4[i]);
      }
    }
  }

  // filtering AK4 jets with pt > 20 and eta<5
  std::vector<fastjet::PseudoJet> AK4_filtered_pt20_eta5;

  for (int i = 0; i < AK4.size(); i++)
  {
    if (AK4[i].pt() > 20. and fabs(AK4[i].eta()) < 5.)
    {
      AK4_filtered_pt20_eta5.push_back(AK4[i]);
    }
  }

  // htmiss trigger
  // float jet_trig_x=0., jet_trig_y=0.;
  // for (int i = 0; i < AK4_filtered_pt20_eta5.size(); i++)
  // {
  //   jet_trig_x += AK4_filtered_pt20_eta5[i].px();
  //   jet_trig_y += AK4_filtered_pt20_eta5[i].py();
  // }

  // float htmiss = sqrtf(jet_trig_x*jet_trig_x + jet_trig_y*jet_trig_y);
  // if(htmiss<100.)
  //   return;

  /////////////////////////////////////////////////////////////////////////////////////////////////

  // defining hadronic recoil

  float sum_leptons_ptx = 0;
  float sum_leptons_pty = 0;
  float sum_leptons_ptz = 0;
  float sum_leptons_E = 0;
  float missET_x = 0;
  float missET_y = 0;
  float missET_z = 0;
  float missET_E = 0;

  if (electrons.size() > 0)
  {
    for (int i = 0; i < electrons.size(); i++)
    {
      sum_leptons_ptx += electrons[i]->P4().Px();
      sum_leptons_pty += electrons[i]->P4().Py();
      sum_leptons_ptz += electrons[i]->P4().Pz();
      sum_leptons_E += electrons[i]->P4().E();
    }
  }
  if (photons.size() > 0)
  {
    for (int i = 0; i < photons.size(); i++)
    {
      sum_leptons_ptx += photons[i]->P4().Px();
      sum_leptons_pty += photons[i]->P4().Py();
      sum_leptons_ptz += photons[i]->P4().Pz();
      sum_leptons_E += photons[i]->P4().E();
    }
  }
  if (taus.size() > 0)
  {
    for (int i = 0; i < taus.size(); i++)
    {
      sum_leptons_ptx += taus[i]->P4().Px();
      sum_leptons_pty += taus[i]->P4().Py();
      sum_leptons_ptz += taus[i]->P4().Pz();
      sum_leptons_E += taus[i]->P4().E();
    }
  }
  if (muons.size() > 0)
  {
    for (int i = 0; i < muons.size(); i++)
    {
      sum_leptons_ptx += muons[i]->P4().Px();
      sum_leptons_pty += muons[i]->P4().Py();
      sum_leptons_ptz += muons[i]->P4().Pz();
      sum_leptons_E += muons[i]->P4().E();
    }
  }

  // sum the missing components of true particles if size is greater than 0 and satus is 1
  if (true_particles.size() > 0)
  {
    for (int i = 0; i < true_particles.size(); i++)
    {
      if (true_particles[i]->Status == 1)
      {
        missET_x += true_particles[i]->P4().Px();
        missET_y += true_particles[i]->P4().Py();
        missET_z += true_particles[i]->P4().Pz();
        missET_E += true_particles[i]->P4().E();
      }
    }
  }
  float missET_pt = sqrtf(missET_x * missET_x + missET_y * missET_y);
  // sum the leptons components of true particles if size is greater than 0 and satus is 1 and PID is between 11 and 18 or 22
  // if (true_particles.size() > 0)
  // {
  //   for (int i = 0; i < true_particles.size(); i++)
  //   {
  //     if (true_particles[i]->Status == 1)
  //     {
  //       if (true_particles[i]->PID > 10 && true_particles[i]->PID < 19 || true_particles[i]->PID == 22)
  //       {
  //         sum_leptons_ptx += true_particles[i]->P4().Px();
  //         sum_leptons_pty += true_particles[i]->P4().Py();
  //         sum_leptons_ptz += true_particles[i]->P4().Pz();
  //         sum_leptons_E += true_particles[i]->P4().E();
  //       }
  //     }
  //   }
  // }

  float had_recoil_x = missET_x - sum_leptons_ptx;
  float had_recoil_y = missET_y - sum_leptons_pty;
  float had_recoil_z = missET_z - sum_leptons_ptz;
  float had_recoil_E = missET_E - sum_leptons_E;
  float had_recoil = sqrtf(had_recoil_x * had_recoil_x + had_recoil_y * had_recoil_y);

  // float hadx = -missingET->P4().Px() - sum_leptons_ptx;
  // float hady = -missingET->P4().Py() - sum_leptons_pty;
  // float hadz = -missingET->P4().Pz() - sum_leptons_ptz;
  // float hadE = missingET->P4().E() - sum_leptons_E;
  // float had_recoil = sqrtf(hadx * hadx + hady * hady);

  // defining sum of muons pt components
  float sum_muons_ptx = 0;
  float sum_muons_pty = 0;

  if (muons.size() > 0)
  {
    for (int i = 0; i < muons.size(); i++)
    {
      sum_muons_ptx += muons[i]->P4().Px();
      sum_muons_pty += muons[i]->P4().Py();
    }
  }

  float muons_pt = sqrtf(sum_muons_ptx * sum_muons_ptx + sum_muons_pty * sum_muons_pty);

  ///////////////////////////////////////////////////////////////////////////////////////////////

  // defining a pseudojet for hadronic recoil

  std::vector<fastjet::PseudoJet> hadronic_recoil_jet;

  hadronic_recoil_jet.push_back(fastjet::PseudoJet(-had_recoil_x, -had_recoil_y, -had_recoil_z, had_recoil_E));

  

  /////veto//////

  // trigger

  if (had_recoil < 200.0)
  {
    return;
  }

  if (photons_veto.size() != 0 || electrons_veto.size() != 0 || taus_veto.size() != 0 || bjets.size() != 0)
    return;

  if (muons_veto.size() != 2)
    return;
  if (muons_veto[0]->Charge * muons_veto[1]->Charge > 0)
    return;

  if (muons_veto[0]->PT < 20.0 && muons_veto[1]->PT < 20.0)
    return;

  if ((muons_veto[0]->P4() + muons_veto[1]->P4()).M() < 60.0 || (muons_veto[0]->P4() + muons_veto[1]->P4()).M() > 120.0)
    return;

  if (AK4.size() > 0)
  {
    if (AK4[0].pt() < 100.0 || fabs(AK4[0].eta()) > 2.5)
      return;
  }

  std::vector<fastjet::PseudoJet> constituents = AK4[0].constituents();
  std::vector<fastjet::PseudoJet> constituents_charge = AK4_charge[0].constituents();
  // print user index of constituents

  float photons_E = 0;
  for (int i = 0; i < constituents.size(); i++)
  {
    if (constituents[i].user_index() == 22)
      std::vector<fastjet::PseudoJet> hadronic_recoil_jet;

    // hadronic_recoil_
    {
      photons_E += constituents[i].E();
    }
  }

  float neutral_hadron_E = 0;
  float charged_hadron_E = 0;

  for (int i = 0; i < constituents.size(); i++)
  {
    if (constituents_charge[i].user_index() == 0)
    {
      neutral_hadron_E += constituents_charge[i].E();
    }
    else
    {
      charged_hadron_E += constituents_charge[i].E();
    }
  }

  neutral_hadron_E -= photons_E;
  float leading_jet_E = AK4[0].E();

  float fraction_neutral = neutral_hadron_E / leading_jet_E;
  float fraction_charged = charged_hadron_E / leading_jet_E;

  if (fraction_neutral > 0.8)
    return;

  if (fraction_charged < 0.1)
    return;

  // int j=AK4_filtered_pt30.size();
  // if (AK4_filtered_pt30.size()>4){
  //   j=4;
  // }
  // for (int k=0; k<j; k++){
  //   if(fabs(AK4_filtered_pt30[k].phi()-hadronic_recoil_jet[0].phi())<0.5){
  //     return;
  //   }
  // }

  cout << muons_pt << " " << had_recoil << endl;

  // mono-V selection criteria

  bool AK8_leading = false;
  bool AK8_subjettiness = false;
  bool AK8_mass = false;
  if (AK8.size() > 0 && AK8[0].pt() > 250.0 && fabs(AK8[0].eta()) < 2.4)
  {
    AK8_leading = true;
  }

  fastjet::contrib::NsubjettinessRatio nSub21_beta1(2, 1, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::NormalizedMeasure(1, 0.8));

  double tau21_beta1 = nSub21_beta1(AK8[0]);

  if (tau21_beta1 < 0.6)
  {
    AK8_subjettiness = true;
  }

  if (AK8[0].m() > 65.0 && AK8[0].m() < 105.0)
  {
    AK8_mass = true;
  }

  if (AK8_leading && AK8_subjettiness && AK8_mass)
  {
    if (had_recoil > 200.0 && had_recoil < 230.0)
      countSignalEvent("23");

    if (had_recoil > 300.0 && had_recoil < 350.0)
      countSignalEvent("24");

    if (had_recoil > 350.0 && had_recoil < 400.0)
      countSignalEvent("25");

    if (had_recoil > 400.0 && had_recoil < 500.0)
      countSignalEvent("26");

    if (had_recoil > 500.0 && had_recoil < 600.0)
      countSignalEvent("27");

    if (had_recoil > 600.0 && had_recoil < 750.0)
      countSignalEvent("28");

    if (had_recoil > 750.0)
      countSignalEvent("29");
  }

  else
  {
    //////signal regions //////////
    ///////////////////////////////
    ////////////////////////////////
    if (had_recoil > 200.0 && had_recoil < 230.0)
      countSignalEvent("1");
    if (had_recoil > 230.0 && had_recoil < 260.0)
      countSignalEvent("2");
    if (had_recoil > 260.0 && had_recoil < 290.0)
      countSignalEvent("3");
    if (had_recoil > 290.0 && had_recoil < 320.0)
      countSignalEvent("4");
    if (had_recoil > 320.0 && had_recoil < 350.0)
      countSignalEvent("5");
    if (had_recoil > 350.0 && had_recoil < 390.0)
      countSignalEvent("6");

    if (had_recoil > 390.0 && had_recoil < 430.0)
      countSignalEvent("7");

    if (had_recoil > 430.0 && had_recoil < 470.0)
      countSignalEvent("8");

    if (had_recoil > 470.0 && had_recoil < 510.0)
      countSignalEvent("9");

    if (had_recoil > 510.0 && had_recoil < 550.0)
      countSignalEvent("10");

    if (had_recoil < 590.0 && had_recoil > 550.0)
      countSignalEvent("11");

    if (had_recoil < 640.0 && had_recoil > 590.0)
      countSignalEvent("12");

    if (had_recoil < 690.0 && had_recoil > 640.0)
      countSignalEvent("13");

    if (had_recoil < 740.0 && had_recoil > 690.0)
      countSignalEvent("14");

    if (had_recoil < 790.0 && had_recoil > 740.0)
      countSignalEvent("15");

    if (had_recoil < 840.0 && had_recoil > 790.0)
      countSignalEvent("16");

    if (had_recoil < 900.0 && had_recoil > 840.0)
      countSignalEvent("17");

    if (had_recoil < 960.0 && had_recoil > 900.0)
      countSignalEvent("18");

    if (had_recoil > 960.0 && had_recoil < 1020.0)
      countSignalEvent("19");

    if (had_recoil < 1090.0 && had_recoil > 1020.0)
      countSignalEvent("20");

    if (had_recoil > 1090.0 && had_recoil < 1160.0)
      countSignalEvent("21");

    if (had_recoil > 1160.0)
      countSignalEvent("22");
  }
}
  


void Test::finalize() {
  // Whatever should be done after the run goes here
}       
