#ifndef Cuts_cxx
#define Cuts_cxx

#include "Cuts.h"

#include "TruthCategories/Sidebands.h"  // sidebands::kSidebandCutVal
#include "utilities.h"                  // ContainerEraser

//==============================================================================
// Helpers
//==============================================================================
// Get the pion candidate indexes from a map of michels
std::vector<int> GetHadIdxsFromMichels(endpoint::MichelMap michels) {
  std::vector<int> ret;
  for (auto m : michels) ret.push_back(m.second.had_idx);
  return ret;
}

//==============================================================================
// Passes ALL Cuts
//==============================================================================
// cut-by-cut, we fill endpoint_michels and vertex_michels.
// if a track fails a cut, we remove the track's michel from the lists.
// then at the end, return the track indices.
bool PassesCuts(CVUniverse& univ, std::vector<int>& pion_candidate_idxs,
                bool is_mc, SignalDefinition signal_definition,
                std::vector<ECuts> cuts) {
  pion_candidate_idxs.clear();
  static endpoint::MichelMap endpoint_michels;
  static endpoint::MichelMap
      vertex_michels;  // Keep track of these, but not used currently
  endpoint_michels.clear();
  vertex_michels.clear();
  bool pass = true;
  for (auto c : cuts) {
    univ.SetPionCandidates(GetHadIdxsFromMichels(
        endpoint_michels));  // Set the pion candidates to the universe
    pass = pass && PassesCut(univ, c, is_mc, signal_definition,
                             endpoint_michels, vertex_michels);
  }

  // Each endpoint michel has an associated hadron track.
  // Our official pion candidates are those tracks.
  pion_candidate_idxs = GetHadIdxsFromMichels(endpoint_michels);

  return pass;
}

// Check all cuts AND whether we are W sideband. Save a lot of time.
// Strategy is:
// 1. check all cuts except for W
// 2. check if W > 1.5 (sideband)
// 3. check if W < 1.4 (signal)
bool PassesCuts(CVUniverse& universe, std::vector<int>& pion_candidate_idxs,
                const bool is_mc, SignalDefinition signal_definition,
                bool& is_w_sideband, std::vector<ECuts> cuts) {
  // is the W cut even in the cuts vector provided?
  bool do_w_cut = std::find(cuts.begin(), cuts.end(), kWexp) != cuts.end();

  // either way, attempt to remove it
  std::vector<ECuts> w_sideband_cuts = kCutsVector;
  w_sideband_cuts.erase(
      std::remove(w_sideband_cuts.begin(), w_sideband_cuts.end(), kWexp),
      w_sideband_cuts.end());

  // check passes all but w cut
  bool passes_all_but_w_cut = PassesCuts(universe, pion_candidate_idxs, is_mc,
                                         signal_definition, w_sideband_cuts);

  // is w sideband = all cuts but W && W > 1.5
  is_w_sideband = passes_all_but_w_cut &&
                  (universe.GetWexp() >= sidebands::kSidebandCutVal);

  // Finally check all cuts == all cuts but W && W
  bool passes_all_cuts = passes_all_but_w_cut;
  if (do_w_cut)
    passes_all_cuts =
        passes_all_but_w_cut && WexpCut(universe, signal_definition);

  return passes_all_cuts;
}

// NEW! Return passes_all_cuts, is_w_sideband, and pion_candidate_indices
std::tuple<bool, bool, std::vector<int>> PassesCuts(
    CVUniverse& universe, const bool is_mc,
    const SignalDefinition signal_definition, std::vector<ECuts> cuts) {
  //============================================================================
  // passes all cuts but w cut
  //============================================================================
  endpoint::MichelMap endpoint_michels;
  endpoint::MichelMap vertex_michels;
  bool passes_all_but_w_cut = true;
  for (auto c : GetWSidebandCuts()) {
    // Set the pion candidates to the universe. The values set in early cuts
    // are used for later cuts, which is why we assign them to the CVU.
    universe.SetPionCandidates(GetHadIdxsFromMichels(endpoint_michels));

    passes_all_but_w_cut =
        passes_all_but_w_cut && PassesCut(universe, c, is_mc, signal_definition,
                                          endpoint_michels, vertex_michels);
  }

  //============================================================================
  // The cuts function returns a container of endpoint michels which are
  // matched to hadron tracks that have passed the pion candidate cuts. From
  // here on out, use the pion candidates to calculate pion quantities for this
  // event-universe.
  //============================================================================
  std::vector<int> pion_candidate_idxs =
      GetHadIdxsFromMichels(endpoint_michels);

  //============================================================================
  // is in the w sideband
  //============================================================================
  bool is_w_sideband = passes_all_but_w_cut &&
                       (universe.GetWexp() >= sidebands::kSidebandCutVal);

  //============================================================================
  // finally: check the w cut
  //============================================================================
  // is the W cut in the cuts vector provided?
  bool do_w_cut = std::find(cuts.begin(), cuts.end(), kWexp) != cuts.end();

  bool passes_all_cuts = passes_all_but_w_cut;
  if (do_w_cut)
    passes_all_cuts =
        passes_all_but_w_cut && WexpCut(universe, signal_definition);

  return {passes_all_cuts, is_w_sideband, pion_candidate_idxs};
}

//==============================================================================
// Fuction to count the number of events that pass the cuts
//==============================================================================
EventCount PassedCuts(const CVUniverse& univ,
                      std::vector<int>& pion_candidate_idxs, bool is_mc,
                      SignalDefinition signal_definition,
                      std::vector<ECuts> cuts) {
  pion_candidate_idxs.clear();
  static endpoint::MichelMap endpoint_michels;
  static endpoint::MichelMap vertex_michels;
  endpoint_michels.clear();
  vertex_michels.clear();
  EventCount Pass;
  bool pass = true;
  for (auto cu : cuts) Pass[cu] = 0;

  for (auto c : cuts) {
    pass = pass && PassesCut(univ, c, is_mc, signal_definition,
                             endpoint_michels, vertex_michels);
    if (pass) {
      Pass[c] = 1.;
    }
  }

  return Pass;
}

//==============================================================================
// Passes INDIVIDUAL Cut
//==============================================================================
// Updates the michel containers
bool PassesCut(const CVUniverse& univ, const ECuts cut, const bool is_mc,
               const SignalDefinition signal_definition,
               endpoint::MichelMap& endpoint_michels, endpoint::MichelMap& vertex_michels) {
  const bool useOVMichels = false;
  if (IsPrecut(cut) && !is_mc) return true;

  switch (cut) {
    case kNoCuts:
      return true;

    case kGoodObjects:
      return univ.IsTruth() ? GoodObjectsCut(univ) : true;

    case kGoodVertex:
      return univ.IsTruth() ? GoodVertexCut(univ) : true;

    case kFiducialVolume:
      return univ.IsTruth() ? FiducialVolumeCut(univ) : true;

    case kMinosActivity:
      return univ.IsTruth() ? MinosActivityCut(univ) : true;

    case kPrecuts:
      return univ.IsTruth() ? GoodObjectsCut(univ) && GoodVertexCut(univ) &&
                                  FiducialVolumeCut(univ)
                            : true;
      // MinosActivityCut(univ) : true;

    case kVtx:
      return vtxCut(univ);

    case kMinosMatch:
      return MinosMatchCut(univ);

    case kMinosCharge:
      return MinosChargeCut(univ);

    case kMinosCoil:
      return MinosCoilCut(univ);

    case kMinosMuon:
      return MinosMatchCut(univ) && MinosChargeCut(univ);
      //&& MinosCoilCut(univ);

    case kThetaMu:
      return ThetamuCut(univ);

    case kDeadTime:
      return DeadTimeCut(univ);

    case kWexp:
      return WexpCut(univ, signal_definition);

    case kIsoBlobs:
      return IsoBlobCut(univ);

    case kIsoProngs:
      return IsoProngCut(univ);

    case kIsoProngSep:
      return IsoProngSepCut(univ);

    case kPmu:
      return PmuCut(univ);

    // ==== At Least One Michel ====
    // For now, we need at least one ENDPOINT michel (any # of vtx michels).
    // This cut fills our michel containers, which we use to ID pion tracks
    // and subsequently make track cuts (LLR, node).
    case kAtLeastOneMichel: {
      endpoint::MichelMap all_michels = endpoint::GetQualityMichels(univ);
      for (auto m : all_michels) {
        if (m.second.had_idx == -1)
          vertex_michels.insert(m);
        else
          endpoint_michels.insert(m);
      }
      vertex::MichelEvent mehreen_michels = vertex::GetQualityMichels(univ);
      return endpoint_michels.size() > 0 /*|| mehreen_michels.size() = 0*/;
    }

    case kAtLeastOnePionCandidateTrack:
      return GetQualityPionCandidateIndices(univ).size() > 0;

    case kAtLeastOneLLRCandidate:
      return AtLeastOneLLRCandidateCut(univ);

    case kAtLeastOneNodeCandidate:
      return AtLeastOneNodeCandidateCut(univ);

    case kAtLeastOneAnchoredProng:
      return AtLeastOneAnchoredProngCut(univ);

    case kAtLeastOneBrandonMichel:
      return AtLeastOneBrandonMichelCut(univ);

    // If a michel's pion fails the LLR cut, remove it from the michels
    case kLLR: {
      ContainerEraser::erase_if(endpoint_michels,
                                [&univ](std::pair<int, endpoint::Michel> mm) {
                                  return !LLRCut(univ, mm.second.had_idx);
                                });
      return endpoint_michels.size() > 0;
    }

    // If a michel's pion fails the node cut, remove it from the michels
    case kNode: {
      ContainerEraser::erase_if(endpoint_michels,
                                [&univ](std::pair<int, endpoint::Michel> mm) {
                                  return !NodeCut(univ, mm.second.had_idx);
                                });
      return endpoint_michels.size() > 0;
    }

    case kPionMult: {
      if (signal_definition == kOnePi || signal_definition == kOnePiNoW)
        return endpoint_michels.size() == 1 && vertex_michels.size() == 0;
      else
        return endpoint_michels.size() >= 1;
    }

    case kAllCuts:
      return true;

    default:
      std::cout << "PassesCut Error Unknown Cut!" << cut << std::endl;
      return false;
  };
}

//==============================================================================
// Cut functions
//==============================================================================
// Truth precuts
bool GoodObjectsCut(const CVUniverse& univ) {
  return univ.GetBool("truth_reco_hasGoodObjects");
}
bool GoodVertexCut(const CVUniverse& univ) {
  return univ.GetBool("truth_reco_isGoodVertex");
}
bool FiducialVolumeCut(const CVUniverse& univ) {
  return univ.GetBool("truth_reco_isFidVol_smeared");
}
bool MinosActivityCut(const CVUniverse& univ) {
  return univ.GetInt("truth_reco_muon_is_minos_match");
}

// Eventwide reco cuts
bool MinosMatchCut(const CVUniverse& univ) {
  return univ.GetBool("isMinosMatchTrack");
}
// Equivalent to Brandon's, but using standard minos branches
bool MinosChargeCut(const CVUniverse& univ) {
  return univ.GetDouble("MasterAnaDev_minos_trk_qp") < 0.0;
}

bool ThetamuCut(const CVUniverse& univ) { return univ.GetThetamu() < 0.3491; }
bool DeadTimeCut(const CVUniverse& univ) { return univ.GetInt("tdead") <= 1; }
bool WexpCut(const CVUniverse& univ, SignalDefinition signal_definition) {
  switch (signal_definition) {
    case kOnePi:
    case kNPi:
      return univ.GetWexp() < GetWCutValue(signal_definition);
    case kOnePiNoW:
    case kNPiNoW:
      return true;
    default:
      std::cout << "WexpCut SIGNAL DEF ERROR";
      return false;
  }
}

// cut on max number of iso prongs
// PrimaryBlobProngTool::makeShowerBlobProngs
bool IsoProngCut(const CVUniverse& univ) {
  return univ.GetNIsoProngs() < CCNuPionIncConsts::kIsoProngCutVal;
}

bool NodeCut(const CVUniverse& univ, const RecoPionIdx pion_candidate_idx) {
  return 6. < univ.GetEnode01(pion_candidate_idx) &&
         univ.GetEnode01(pion_candidate_idx) < 32. &&
         2. < univ.GetEnode2(pion_candidate_idx) &&
         univ.GetEnode2(pion_candidate_idx) < 22. &&
         0. < univ.GetEnode3(pion_candidate_idx) &&
         univ.GetEnode3(pion_candidate_idx) < 19. &&
         0. < univ.GetEnode4(pion_candidate_idx) &&
         univ.GetEnode4(pion_candidate_idx) < 31. &&
         0. < univ.GetEnode5(pion_candidate_idx) &&
         univ.GetEnode5(pion_candidate_idx) < 60.;
}

bool LLRCut(const CVUniverse& univ, const RecoPionIdx pion_candidate_idx) {
  // if (pion_candidate_idx < 0) return false;
  // else return univ.GetLLRScore(pion_candidate_idx) > 0.;
  return univ.GetLLRScore(pion_candidate_idx) > 0.;
}

// Helper -- get vector of good-enough had tracks, id'ed by unique indices.
std::vector<int> GetQualityPionCandidateIndices(const CVUniverse& univ) {
  std::vector<int> pion_candidate_indices;
  int n_hadrons = univ.GetInt("MasterAnaDev_hadron_number");
  for (int i_hadron = 0; i_hadron != n_hadrons; ++i_hadron)
    if (HadronQualityCuts(univ, i_hadron))
      pion_candidate_indices.push_back(i_hadron);
  return pion_candidate_indices;
}

bool HadronQualityCuts(const CVUniverse& univ,
                       const RecoPionIdx pion_candidate_idx) {
  return univ.GetVecElem("MasterAnaDev_hadron_isForked", pion_candidate_idx) ==
             0 &&
         univ.GetVecElem("MasterAnaDev_hadron_isExiting", pion_candidate_idx) ==
             0 &&
         univ.GetVecElem("MasterAnaDev_hadron_isSideECAL",
                         pion_candidate_idx) == 0 &&
         univ.GetVecElem("MasterAnaDev_hadron_isODMatch", pion_candidate_idx) ==
             0 &&
         univ.GetVecElem("MasterAnaDev_hadron_isTracker", pion_candidate_idx) ==
             1;
};
// Vtx cut for detection volume
bool vtxCut(const CVUniverse& univ) {
  bool pass = true;
  pass = pass && zVertexCut(univ, 8340.0, 5990.0);
  pass = pass && XYVertexCut(univ, 850.0);
  return pass;
}

bool zVertexCut(const CVUniverse& univ, const double upZ, const double downZ) {
  double vtxZ = univ.GetVecElem("vtx", 2);
  if (vtxZ > downZ && vtxZ < upZ)
    return true;
  else
    return false;
}

bool XYVertexCut(const CVUniverse& univ, const double a) {
  const double x = univ.GetVecElem("vtx", 0), y = univ.GetVecElem("vtx", 1);
  if (x < 0) {
    if (x > -a && univ.leftlinesCut(a, x, y))
      return true;
    else
      return false;
  } else {
    if (x < a && univ.rightlinesCut(a, x, y))
      return true;
    else
      return false;
  }
}

bool PmuCut(const CVUniverse& univ) {
  if (univ.GetPmu() / 1000 < 1.5 || 20 < univ.GetPmu() / 1000)
    return false;
  else
    return true;
}

//==============================================================================
// Retired
//==============================================================================
// Brandon's Minos Cut
bool BrandonMinosChargeCut(CVUniverse& univ) {
  return univ.GetDouble("MasterAnaDev_muon_qpqpe") < 0.0;
}

// CCInclusive-style minos charge cut
bool CCIncMinosChargeCut(CVUniverse& univ) {
  if (univ.GetBool("MasterAnaDev_minos_used_curvature"))
    return 1. / univ.GetDouble("MasterAnaDev_minos_trk_eqp_qp") < -5.0;
  else if (univ.GetBool("MasterAnaDev_minos_used_range"))
    return univ.GetBool("MasterAnaDev_minos_trk_qp") < 0.0;
  else
    return false;
}

bool MinosCoilCut(const CVUniverse& univ) {
  const double MINOS_COIL_RADIUS = 210;  // mm
  const double MAX_MINOS_RADIUS = 2500;  // mm
  const double coilXPos = 1219.0;
  const double coilYPos = 393.0;
  const double minos_x =
      univ.GetDouble("MasterAnaDev_minos_trk_end_x") + coilXPos;
  const double minos_y =
      univ.GetDouble("MasterAnaDev_minos_trk_end_y") + coilYPos;
  double minosR = sqrt(pow(minos_x, 2) + pow(minos_y, 2));
  // if (!((pow(minos_x,2) + pow(minos_y,2) )>= pow(MINOS_COIL_RADIUS, 2)) )
  //  cout << minos_x << " " << minos_y << " " << MINOS_COIL_RADIUS << endl;
  return (minosR > MINOS_COIL_RADIUS && minosR < MAX_MINOS_RADIUS);
}

// no isolated blobs allowed
bool IsoBlobCut(const CVUniverse& univ) {
  return univ.GetInt("n_iso_blob_prongs") <
         1;  // RecoilUtils::createIsoBlobProngs
}

// IsoProngSeparation
bool IsoProngSepCut(const CVUniverse& univ) {
  return univ.GetLargestIsoProngSep() < 300;
}

bool AtLeastOneBrandonMichelCut(const CVUniverse& univ) {
  // Get quality hadron candidates
  std::vector<int> pion_candidate_indices =
      GetQualityPionCandidateIndices(univ);

  // Loop quality hadron candidates to see if they have a good brandon michel
  for (auto pion_candidate_idx : pion_candidate_indices) {
    int michel_views = univ.GetVecElem("MasterAnaDev_hadron_endMichel_category",
                                       pion_candidate_idx);
    int michel_ndigits = univ.GetVecElem(
        "MasterAnaDev_hadron_endMichel_ndigits", pion_candidate_idx);
    double michel_energy =
        univ.GetVecElem("MasterAnaDev_hadron_endMichel_energy",
                        pion_candidate_idx);  /// TODO sys universe function
    double michel_slice_energy = univ.GetVecElem(
        "MasterAnaDev_hadron_endMichel_slice_energy", pion_candidate_idx);

    if (michel_views < 1)
      continue;                    // no michel
    else if (michel_views == 1) {  // 1 view
      if (michel_energy < 55.0 && michel_ndigits < 35 &&
          michel_slice_energy < 100.0)
        return true;
      else
        continue;
    } else if (michel_views > 1) {  // 2+3 views
      if (michel_energy < 55.0 && michel_ndigits < 35 &&
          michel_ndigits >= michel_views)
        return true;
      else
        continue;
    }
  }
  return false;
}

// ntrk are tracks not including the muon prong
bool AtLeastOneAnchoredProngCut(const CVUniverse& univ) {
  int ntrk = univ.GetInt("n_anchored_long_trk_prongs") +
             univ.GetInt("n_anchored_short_trk_prongs");
  return ntrk > 0;
}

bool ExactlyOneEndpointMichelCut(const CVUniverse& univ,
                                 SignalDefinition signal_definition) {
  if (signal_definition == kNPi || signal_definition == kNPiNoW) {
    return true;
  } else if (signal_definition == kOnePi || signal_definition == kOnePiNoW) {
    endpoint::MichelMap mm = endpoint::GetQualityMichels(univ);
    if (mm.size() == 1) {  // require only one michel
      endpoint::Michel m = (mm.begin())->second;
      if (m.vtx == 0)
        return false;  // so this has a no-vertex michel cut baked in
      // pion_candidate_idx = m.vtx - 1;   // SELECT OUR PION
      return true;
    } else
      return false;
  } else {
    std::cout << "ExactlyOneEndpointMichelcut SIGNAL DEFINITION ERROR"
              << std::endl;
    return false;
  }
}

bool AtLeastOneLLRCandidateCut(const CVUniverse& univ) {
  // Get quality hadron candidates
  std::vector<int> pion_candidate_indices =
      GetQualityPionCandidateIndices(univ);
  for (auto pion_candidate_idx : pion_candidate_indices) {
    if (LLRCut(univ, pion_candidate_idx)) return true;
  }
  return false;
}

bool AtLeastOneNodeCandidateCut(const CVUniverse& univ) {
  // Get quality hadron candidates
  std::vector<int> pion_candidate_indices =
      GetQualityPionCandidateIndices(univ);
  for (auto pion_candidate_idx : pion_candidate_indices) {
    if (NodeCut(univ, pion_candidate_idx)) return true;
  }
  return false;
}

//==============================================================================
// Cut Names
//==============================================================================
std::string GetCutName(ECuts cut) {
  switch (cut) {
    case kNoCuts:
      return "No Cuts";

    case kGoodObjects:
      return "Good Objects";

    case kGoodVertex:
      return "Good Vertex";

    case kFiducialVolume:
      return "Fiducial Volume";

    case kMinosActivity:
      return "MINOS Activity";

    case kPrecuts:
      return "Anatool Precuts";

    case kVtx:
      return "vertex position Cut";

    case kMinosMatch:
      return "MINOS Muon";

    case kMinosCharge:
      return "MINOS Charge";

    case kMinosCoil:
      return "MINOS Coil";

    case kMinosMuon:
      return "MINOS Muon";

    case kThetaMu:
      return "Muon Angle";

    case kDeadTime:
      return "Dead Time";

    case kWexp:
      return "$W_{experimental}$";

    case kIsoBlobs:
      return "$<$1 Isolated Blobs";

    case kIsoProngs:
      return "$<$2 Isolated Prongs";

    case kIsoProngSep:
      return "Iso Prong Sep $<$ 300";

    case kNProngs:
      return "Max 2 Hadr Prongs";

    case kNPionCandidates:
      return "$\\pi$ candidate";

    case kPionCandidateQuality:
      return "Quality $\\pi$ candidate";

    case kAtLeastOneMichel:
      return "$>$= 1 Michel";

    case kAtLeastOneBrandonMichel:
      return "$>$= 1 Brandon Michel";

    case kAtLeastOneAnchoredProng:
      return "$>$= 1 Anchored Prong";

    case kAtLeastOnePionCandidateTrack:
      return "$>$= 1 Hadron Track";

    case kAtLeastOneLLRCandidate:
      return "$>$= 1 Track Pass LLR Cut";

    case kAtLeastOneNodeCandidate:
      return "$>$= 1 Track Pass Node Cut";

    case kExactlyOneEndpointMichel:
      return "== 1 Michel";

    case kNode:
      return "Node";

    case kPionMult:
      return "Pion Multiplicity";

    case kOldMichel:
      return "Old Michel Cut";

    case kdEdx:
      return "dEdx PID";

    case kLLR:
      return "LLR PID";

    case kAllCuts:
      return "Total";

    case kPmu:
      return "1.5 GeV $<$ Pmu $<$ 20 GeV";

    default:
      std::cout << "ERROR: GetCutName unknown cut!" << std::endl;
      return "";
  };
}



#endif  // Cuts_cxx
