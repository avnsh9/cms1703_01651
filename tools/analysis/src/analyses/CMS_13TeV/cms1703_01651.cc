#include "cms1703_01651.h"

#include "cmath"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contribs/Nsubjettiness/Nsubjettiness.hh"
#include "fastjet/tools/Pruner.hh"

// AUTHOR: Avinash Verma
//  EMAIL: avinash.verma@students.iiserpune.ac.in
void Cms1703_01651::initialize() {
  setAnalysisName("cms1703_01651");          
  setInformation(""
    "# Search for dark matter produced with an energetic\n"
    "# jet or ahadronically decaying W or Z boson\n"
  "");
  setLuminosity(12.9*units::INVFB);      
  bookSignalRegions("MONOJ230_260;MONOJ260_290;MONOJ290_320;MONOJ320_350;MONOJ350_390;MONOJ390_430;MONOJ430_470;MONOJ470_510;MONOJ510_550;MONOJ550_590;MONOJ590_640;MONOJ640_690;MONOJ690_740;MONOJ740_790;MONOJ790_840;MONOJ840_900;MONOJ900_960;MONOJ960_1020;MONOJ1020_1090;MONOJ1090_1160;MONOJ1160;MONOV250_300;MONOV300_350;MONOV350_400;MONOV400_500;MONOV500_600;MONOV600_750;MONOV750");
  // You can also book cutflow regions with bookCutflowRegions("CR1;CR2;..."). Note that the regions are
  //  always ordered alphabetically in the cutflow output files.

  // You should initialize any declared variables here
}

void Cms1703_01651::analyze() {
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

  missingET->addMuons(muonsCombined);  // Adds muons to missing ET. This should almost always be done which is why this line is not commented out.
  std::vector<Electron *> electrons_veto;
  std::vector<Muon *> muons_veto;
  std::vector<Photon *> photons_veto;
  // cout<<"missing phi: "<<missingET->Phi<<endl;
  // construct tau
  std::vector<Jet *> taus;
  std::vector<Jet *> taus_veto;
  std::vector<Jet *> jets20;
  jets20 = filterPhaseSpace(jets, 20., -5.0, 5.0);
  float jets20_ptx = 0;
  float jets20_pty = 0;
  for (int i = 0; i < jets20.size(); i++) {
    jets20_ptx += jets20[i]->P4().Px();
    jets20_pty += jets20[i]->P4().Py();
  }
  float jets20_pt = sqrt(jets20_ptx * jets20_ptx + jets20_pty * jets20_pty);

  if (jets.size() > 0)
  {
    for (int i = 0; i < jets.size(); i++)
    {
      // if (checkTauTag(jets[i], "tight")) {
      // if (checkTauTag(jets[i], "medium")) {
      if (checkTauTag(jets[i], "loose") and fabs(jets[i]->Charge) == 1)
      {
        taus.push_back(jets[i]);
      }
    }
  }
  taus_veto = filterPhaseSpace(taus, 18., -2.3, 2.3);
  // b-jet construction
  std::vector<Jet *> bjets;
  std::vector<Jet *> bjets_veto;
  if (jets.size() > 0)
  {
    for (int i = 0; i < jets.size(); i++)
    {
      if (checkBTag(jets[i]))
      {
        bjets.push_back(jets[i]);
      }
    }
  }
  bjets_veto = filterPhaseSpace(bjets, 15., -2.4, 2.4);

  

  

  std::vector<fastjet::PseudoJet> particles;
  for (int i = 0; i < true_particles.size(); i++)
  {
    if (true_particles[i]->Status == 1)
    {
      int pid = true_particles[i]->PID;

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
      int AK8_pid = true_particles[i]->PID;
      particles_AK8.push_back(fastjet::PseudoJet(true_particles[i]->P4().Px(), true_particles[i]->P4().Py(), true_particles[i]->P4().Pz(), true_particles[i]->P4().E()));
      particles_AK8.back().set_user_index(AK8_pid);
    }
  }
  fastjet::JetDefinition jetAK8(fastjet::antikt_algorithm, 0.8);
  fastjet::ClusterSequence AK8jet(particles_AK8, jetAK8);
  std::vector<fastjet::PseudoJet> AK8 = sorted_by_pt(AK8jet.inclusive_jets());

  // AK8_charge

  std::vector<fastjet::PseudoJet> ak8_charge;
  for (int i = 0; i < true_particles.size(); i++)
  {
    if (true_particles[i]->Status == 1)
    {
      int AK8_charge = true_particles[i]->Charge;
      ak8_charge.push_back(fastjet::PseudoJet(true_particles[i]->P4().Px(), true_particles[i]->P4().Py(), true_particles[i]->P4().Pz(), true_particles[i]->P4().E()));
      ak8_charge.back().set_user_index(AK8_charge);
    }
  }
  fastjet::JetDefinition jetAK8_charge(fastjet::antikt_algorithm, 0.8);
  fastjet::ClusterSequence AK8jet_charge(ak8_charge, jetAK8_charge);
  std::vector<fastjet::PseudoJet> AK8_charge = sorted_by_pt(AK8jet_charge.inclusive_jets());



  ///////////////////////////////////////////////////////////

  electrons_veto = filterPhaseSpace(electrons, 10., -2.5, 2.5);
  // electrons_veto = filterPhaseSpace(electronsLoose, 10., -2.5, 2.5);   // Use this for control region
  // muons_veto = filterPhaseSpace(muons, 10.); /// use this for control region
  muons_veto = filterPhaseSpace(muons, 10., -2.4, 2.4);

  photons_veto = filterPhaseSpace(photons, 15., -2.5, 2.5);
  // photons_veto = filterPhaseSpace(photonsLoose, 15., -2.5, 2.5);  //use this for control region
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
      if (true_particles[i]->Status == 1 && true_particles[i]->PID !=52)
      {
        missET_x += true_particles[i]->P4().Px();
        missET_y += true_particles[i]->P4().Py();
        missET_z += true_particles[i]->P4().Pz();
        missET_E += true_particles[i]->P4().E();
      }
    }
  }
  float missET_pt = sqrtf(missET_x * missET_x + missET_y * missET_y);
  

  float had_recoil_x = missET_x - sum_leptons_ptx;
  float had_recoil_y = missET_y - sum_leptons_pty;
  float had_recoil_z = missET_z - sum_leptons_ptz;
  float had_recoil_E = missET_E - sum_leptons_E;
  float had_recoil = sqrtf(had_recoil_x * had_recoil_x + had_recoil_y * had_recoil_y);

  float hadx = -missingET->P4().Px() - sum_leptons_ptx;
  float hady = -missingET->P4().Py() - sum_leptons_pty;
  float hadz = -missingET->P4().Pz() - sum_leptons_ptz;
  float hadE = missingET->P4().E() - sum_leptons_E;
  float hadronic_recoil = sqrtf(hadx * hadx + hady * hady);

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
  // defining sum of all neutrinos
  float sum_neutrinos_ptx = 0;
  float sum_neutrinos_pty = 0;
  for (int i=0; i<true_particles.size(); i++)
  {
    if (true_particles[i]->PID == 12 || true_particles[i]->PID == 14 || true_particles[i]->PID == 16)
    {
      sum_neutrinos_ptx += true_particles[i]->P4().Px();
      sum_neutrinos_pty += true_particles[i]->P4().Py();
    }
  }
  float neutrino_pt=sqrtf(sum_neutrinos_ptx*sum_neutrinos_ptx+sum_neutrinos_pty*sum_neutrinos_pty);
  // cout<<missingET->PT<<" "<<missET_pt<<endl;

  float muons_pt = sqrtf(sum_muons_ptx * sum_muons_ptx + sum_muons_pty * sum_muons_pty);

  /////veto//////

  // trigger
  if(jets20_pt < 120.)
    return;

  if (missingET->PT < 200.0)
  {
    return;
  }

  if (photons_veto.size() != 0 || electrons_veto.size() != 0 || taus_veto.size() != 0 || bjets_veto.size() != 0)
    return;

  // if (muons_veto.size() != 2)
  //   return;

  if (muons_veto.size() != 0)
    return;
  // if (muons_veto[0]->Charge * muons_veto[1]->Charge > 0)
  //   return;

  // if (muons_veto[0]->PT < 20.0 && muons_veto[1]->PT < 20.0)
  //   return;

  // if ((muons_veto[0]->P4() + muons_veto[1]->P4()).M() < 60.0 || (muons_veto[0]->P4() + muons_veto[1]->P4()).M() > 120.0)
  //   return;

  if (AK4.size() > 0)
  {
    if (AK4[0].pt() < 100.0 || fabs(AK4[0].eta()) > 2.5)
      return;
  }

  //angle cut
  

  int i=AK4_filtered_pt30.size();
  if (i>4){
    i=4;
  }
  for (int j=0;j<i;j++){
    if ((fabs(missingET->Phi - AK4_filtered_pt30[j].phi_std()) < 0.5) || (fabs(missingET->Phi - AK4_filtered_pt30[j].phi_std()) > (2*M_PI - 0.5))){
      return;
    }
  }
  cout<<AK4[0].phi()<<" "<<AK4[0].phi_std()<<endl;
  std::vector<fastjet::PseudoJet> constituents = AK4[0].constituents();
  std::vector<fastjet::PseudoJet> constituents_charge = AK4_charge[0].constituents();
  // print user index of constituents

  float photons_E = 0;
  for (int i = 0; i < constituents.size(); i++)
  {
    if (constituents[i].user_index() == 22)

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

  


  bool AK8_mass = false;
  bool AK8_charged_constituents = false;

  fastjet::contrib::NsubjettinessRatio nSub21_beta1(2, 1, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::NormalizedMeasure(1.0, 0.8));

  // leading AK8 mass
  float rfact = 0.5;
  float Zcut = 0.10;
  fastjet::Pruner pruner(fastjet::cambridge_algorithm, Zcut, rfact);
  fastjet::PseudoJet pruned_jet = pruner(AK8[0]);
  float mass_pruned = pruned_jet.m();
  // cout << mass_pruned << endl;
  bool pruned_leading = false;
  bool pruned_subjettiness = false;
  if (mass_pruned > 65.0 && mass_pruned < 105.0)
  {
    AK8_mass = true;
  }
  if (pruned_jet.pt() > 250.0 && fabs(pruned_jet.eta()) < 2.4)
  {
    pruned_leading = true;
  }
  float tau = nSub21_beta1(pruned_jet);
  if (tau < 0.6)
  {
    pruned_subjettiness = true;
  }

  // cout << "pruned jet energery: " << pruned_jet.E() << endl;
  // cout << "pruned jet transverse energy: " << pruned_jet.Et() << endl;
  // std::vector<fastjet::PseudoJet> constituents_pruned = pruned_jet.constituents();
  // int pruned_size = constituents_pruned.size();
  // float tra_mass=0;
  // for (int i = 0; i < constituents_pruned.size(); i++)
  // {
  //   for (int j = 0; j < pruned_size; j++)
  //   {
  //     if (i != j)
  //     {
  //       float e1=constituents_pruned[i].Et();
  //       float e2=constituents_pruned[j].Et();
  //       float pt1=constituents_pruned[i].pt();
  //       float pt2=constituents_pruned[j].pt();
  //       float b1=pt1/e1;
  //       float b2=pt2/e2;
  //       float m1=constituents_pruned[i].m();
  //       float m2=constituents_pruned[j].m();
  //       float rap1=constituents_pruned[i].rap();
  //       float rap2=constituents_pruned[j].rap();
  //       float phi1=constituents_pruned[i].phi();
  //       float phi2=constituents_pruned[j].phi();
  //       tra_mass += m1*m1+m2*m2 + 2*e1*e2*(cosh(rap1-rap2) - b1*b2*cos(phi1-phi2));
  //     }
  //   }
  // }
  // float mass=sqrtf(tra_mass/2);
  // cout<<"mass: "<<mass<<endl;
  // cout<<"mass_pruned: "<<mass_pruned<<endl;

  // check leading AK8 charged constituents ratio

  std::vector<fastjet::PseudoJet> constituents_AK8 = AK8[0].constituents();
  std::vector<fastjet::PseudoJet> constituents_charge_AK8 = AK8_charge[0].constituents();
  float photons_E_AK8 = 0;
  for (int i = 0; i < constituents_AK8.size(); i++)
  {
    if (constituents_AK8[i].user_index() == 22)

    {
      photons_E_AK8 += constituents_AK8[i].E();
    }
  }

  float neutral_hadron_E_AK8 = 0;
  float charged_hadron_E_AK8 = 0;

  for (int i = 0; i < constituents_AK8.size(); i++)
  {
    if (constituents_charge_AK8[i].user_index() == 0)
    {
      neutral_hadron_E_AK8 += constituents_charge_AK8[i].E();
    }
    else
    {
      charged_hadron_E_AK8 += constituents_charge_AK8[i].E();
    }
  }

  neutral_hadron_E_AK8 -= photons_E_AK8;
  float leading_jet_E_AK8 = AK8[0].E();

  float fraction_neutral_AK8 = neutral_hadron_E_AK8 / leading_jet_E_AK8;
  float fraction_charged_AK8 = charged_hadron_E_AK8 / leading_jet_E_AK8;

  if (fraction_neutral_AK8 < 0.8 && fraction_charged_AK8 > 0.1)
  {
    AK8_charged_constituents = true;
  }

  if (pruned_leading && pruned_subjettiness && AK8_mass && AK8_charged_constituents && missingET->PT > 250.0)
  {
    if (missingET->PT > 250.0 && missingET->PT < 300.0)
      countSignalEvent("MONOV250_300");

    if (missingET->PT > 300.0 && missingET->PT < 350.0)
      countSignalEvent("MONOV300_350");

    if (missingET->PT > 350.0 && missingET->PT < 400.0)
      countSignalEvent("MONOV350_400");

    if (missingET->PT > 400.0 && missingET->PT < 500.0)
      countSignalEvent("MONOV400_500");

    if (missingET->PT > 500.0 && missingET->PT < 600.0)
      countSignalEvent("MONOV500_600");

    if (missingET->PT > 600.0 && missingET->PT < 750.0)
      countSignalEvent("MONOV600_750");

    if (missingET->PT > 750.0)
      countSignalEvent("MONOV750");
  }

  else
  {
    //////signal regions //////////
    ///////////////////////////////
    ////////////////////////////////
    
    if (missingET->PT > 230.0 && missingET->PT < 260.0)
      countSignalEvent("MONOJ230_260");
    if (missingET->PT > 260.0 && missingET->PT < 290.0)
      countSignalEvent("MONOJ260_290");
    if (missingET->PT > 290.0 && missingET->PT < 320.0)
      countSignalEvent("MONOJ290_320");
    if (missingET->PT > 320.0 && missingET->PT < 350.0)
      countSignalEvent("MONOJ320_350");
    if (missingET->PT > 350.0 && missingET->PT < 390.0)
      countSignalEvent("MONOJ350_390");

    if (missingET->PT > 390.0 && missingET->PT < 430.0)
      countSignalEvent("MONOJ390_430");

    if (missingET->PT > 430.0 && missingET->PT < 470.0)
      countSignalEvent("MONOJ430_470");

    if (missingET->PT > 470.0 && missingET->PT < 510.0)
      countSignalEvent("MONOJ470_510");

    if (missingET->PT > 510.0 && missingET->PT < 550.0)
      countSignalEvent("MONOJ510_550");

    if (missingET->PT < 590.0 && missingET->PT > 550.0)
      countSignalEvent("MONOJ550_590");

    if (missingET->PT < 640.0 && missingET->PT > 590.0)
      countSignalEvent("MONOJ590_640");

    if (missingET->PT < 690.0 && missingET->PT > 640.0)
      countSignalEvent("MONOJ640_690");

    if (missingET->PT < 740.0 && missingET->PT > 690.0)
      countSignalEvent("MONOJ690_740");

    if (missingET->PT < 790.0 && missingET->PT > 740.0)
      countSignalEvent("MONOJ740_790");

    if (missingET->PT < 840.0 && missingET->PT > 790.0)
      countSignalEvent("MONOJ790_840");

    if (missingET->PT < 900.0 && missingET->PT > 840.0)
      countSignalEvent("MONOJ840_900");

    if (missingET->PT < 960.0 && missingET->PT > 900.0)
      countSignalEvent("MONOJ900_960");

    if (missingET->PT > 960.0 && missingET->PT < 1020.0)
      countSignalEvent("MONOJ960_1020");

    if (missingET->PT < 1090.0 && missingET->PT > 1020.0)
      countSignalEvent("MONOJ1020_1090");

    if (missingET->PT > 1090.0 && missingET->PT < 1160.0)
      countSignalEvent("MONOJ1090_1160");

    if (missingET->PT > 1160.0)
      countSignalEvent("MONOJ1160");
  }
  
}

void Cms1703_01651::finalize() {
  // Whatever should be done after the run goes here
}       
