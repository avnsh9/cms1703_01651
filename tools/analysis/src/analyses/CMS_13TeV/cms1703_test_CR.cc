#include "cms1703_test.h"
#include "cmath"
// AUTHOR: Avinash Verma
//  EMAIL: avinash.verma@students.iiserpune.ac.in
void Cms1703_test_CR::initialize() {
  setInformation(""
    "@#Search for dark matter produced with an energetic\n"
    "@#jet or a\n"
    "@#√\n"
    "@#hadronically decaying W or Z boson at\n"
  "");
  setLuminosity(12.9*units::INVFB);  
  bookSignalRegions("MONOJ200_230;MONOJ230_260;MONOJ260_290;MONOJ290_320;MONOJ320_350;MONOJ350_390;MONOJ390_430;MONOJ430_470;MONOJ470_510;MONOJ510_550;MONOJ550_590;MONOJ590_640;MONOJ640_690;MONOJ690_740;MONOJ740_790;MONOJ790_840;MONOJ840_900;MONOJ900_960;MONOJ960_1020;MONOJ1020_1090;MONOJ1090_1160;MONOJ1160;MONOV250_300;MONOV300_350;MONOV350_400;MONOV400_500;MONOV500_600;MONOV600_750;MONOV750");

  setAnalysisName("cms1703_test_CR");
  // You should initialize any declared variables here
}

void Cms1703_test_CR::analyze() {
  // Your eventwise analysis code goes here
  ++n;
  std::vector<Electron*> electrons_veto;
  std::vector<Muon*> muons_veto;
  std::vector<Photon*> photons_veto;
  std::vector<Jet*> jets_trigger;
  std::vector<Jet*> jets_signal;
  std::vector<GenParticle*> true_particles;
  //construct tau
  std::vector<Jet *> taus;
  for (int i = 0; i < jets.size(); i++)
  {
    //if (checkTauTag(jets[i], "tight")) {
    //if (checkTauTag(jets[i], "medium")) {
    if (checkTauTag(jets[i], "loose") and fabs(jets[i]->Charge) == 1)
    {
      taus.push_back(jets[i]);
    }
  }
  taus = filterPhaseSpace(taus, 18., -2.3, 2.3);
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
  electrons_veto = filterPhaseSpace(electrons, 10., -2.5, 2.5);
  muons_veto = filterPhaseSpace(muons, 10., -2.4, 2.4);
  photons_veto = filterPhaseSpace(photons, 15., -2.5, 2.5);
  jets_trigger= filterPhaseSpace(jets, 20., -5, 5);
  jets_signal = filterPhaseSpace(jets, 20., -2.5, 2.5);

  // hadronic recoil
  float totmomentumx = 0;
  float totmomentumy = 0;
  
  float hadronic_recoil = 0;

  if (true_particles.size()>0){
    for (int i = 0; i < true_particles.size();i++){
      totmomentumx += true_particles[i]->P4().Px();
      totmomentumy += true_particles[i]->P4().Py();

    }
    for (int i = 0; i < true_particles.size();i++){
      if ((true_particles[i]->PID >= 11 && true_particles[i]->PID <= 18) || true_particles[i]->PID ==22){
        totmomentumx -= true_particles[i]->P4().Px();
        totmomentumy -= true_particles[i]->P4().Py();
      }
    }
  }
  float hadronic_recoil_temp=totmomentumx*totmomentumx+totmomentumy*totmomentumy;
  hadronic_recoil=sqrtf(hadronic_recoil_temp);

  cout<<"event no is :"<<n<<endl;
  cout<<"hadronic recoil is :"<<hadronic_recoil<<endl;
  cout<<"missingET is :"<<missingET->PT<<endl;
  
  





}

void Cms1703_test_CR::finalize() {
  // Whatever should be done after the run goes here
}       
