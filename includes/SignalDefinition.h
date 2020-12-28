#ifndef SignalDefinition_H
#define SignalDefinition_H

#include "includes/CVUniverse.h"

enum SignalDefinition{ 
  kOnePi, kOnePiNoW, kNPi, kNPiNoW, kNSignalDefTypes 
};

double GetWCutValue(SignalDefinition signal_definition) {
  switch (signal_definition) {
    case kOnePi:
      return 1400.;
    case kNPi:
      return 2000.;
    case kOnePiNoW:
    case kNPiNoW:
      return 120000.;
    default:
      std::cout << "ERROR GetWCutValue" << std::endl;
      return -1.;
  }
}

// Number of abs(pdg) == 211 true TG4Trajectories which also:
// (1) are pip, (2) satisfy a KE restriction
int NSignalPions(const CVUniverse& univ) {
  int n_signal_pions = 0;
  int n_true_pions = univ.GetNChargedPionsTrue();
  for (TruePionIdx idx = 0; idx<n_true_pions; ++idx) {
    double t_pi = univ.GetTpiTrue(idx);
    double theta_pi = univ.GetThetapiTrue(idx);
    if(univ.GetPiChargeTrue(idx) > 0 
       && t_pi > 35
       && t_pi < 350
       //&& (theta_pi < 1.39626 || 1.74533 < theta_pi))
    )
      ++n_signal_pions;
  }
  return n_signal_pions;
}

int NOtherParticles(const CVUniverse& univ){
  int n_other_particles = 0;
  n_other_particles = univ.GetInt("truth_N_chargedK") +
                      univ.GetInt("truth_N_K0") +
                      univ.GetInt("truth_N_sigma") +
                      univ.GetInt("truth_N_lambda");
  return n_other_particles;
}

bool IsSignal(const CVUniverse& universe, SignalDefinition signal_definition = kOnePi) {
  int n_signal_pions = NSignalPions(universe);
  if( universe.GetInt("mc_current")  == 1 
       && universe.GetBool("truth_isFidVol") 
       && universe.GetInt("mc_incoming") == 14 
       && universe.GetDouble("truth_mu_theta_wrtbeam") < 0.3491 // 20 deg
       && universe.GetWexpTrue() > 0
       && universe.GetWexpTrue() < GetWCutValue(signal_definition)
       && n_signal_pions > 0
       && NOtherParticles(universe) == 0
       // && TODO Muon or neutrino energy cut
       //&& 1500. < universe.GetDouble("mc_incomingE") && universe.GetDouble("mc_incomingE") < 10000.
  ) {}
  
  else {
    return false;
  } 


  switch (signal_definition) {
    case kOnePi:
    case kOnePiNoW:
      if (n_signal_pions == 1
          && universe.GetInt("truth_N_pi0") == 0                        
          && universe.GetInt("truth_N_pim") == 0
         )
          return true;
      else return false;
    case kNPi:
    case kNPiNoW:
      return true;

    default:
      std::cout << "IsSignal Error Unknown Signal Definition!" << std::endl;
      return false;
  }
}

std::string GetSignalName(SignalDefinition signal_definition){
  switch (signal_definition){
    case kOnePi:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X  (W < 1.4 GeV)";
    case kOnePiNoW:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X";
    case kNPi:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X  (W < 1.8 GeV)";
    case kNPiNoW:
      return "#nu_{#mu} Tracker #rightarrow #mu^{-} 1#pi^{+} X";
  }
}

std::string GetSignalFileTag(SignalDefinition signal_definition){
  switch (signal_definition){
    case kOnePi:
      return "1Pi";
    case kOnePiNoW:
      return "1PiNoW";
    case kNPi:
      return "NPi";
    case kNPiNoW:
      return "NPiNoW";
  }
}

#endif // Signal_H