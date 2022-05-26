#ifndef Variable2D_cxx
#define Variable2D_cxx

#include "Variable.h"
#include "Variable2D.h"
#include "TruthMatching.h"
#include <stdlib.h> // exit


// CTOR -- default
Variable2D::Variable2D()
  : m_labelX(),
    m_labelY(),
    m_unitsX(),
    m_unitsY(),
    m_pointer_to_GetValueX(&CVUniverse::GetDummyVar),
    m_pointer_to_GetValueY(&CVUniverse::GetDummyVar),
    m_hists2D(),
    m_is_true(false)
{}


// CTOR -- uniform binning
Variable2D::Variable2D(const std::string labelX, const std::string labelY, 
                       const std::string xaxisX, const std::string xaxisY,
                       const std::string unitsX, const std::string unitsY,
                       const int nbinsX, const double xminX, const double xmaxX,
                       const int nbinsY, const double xminY, const double xmaxY,
                       PointerToCVUniverseFunction px,
                       PointerToCVUniverseFunction py,
                       const bool is_true)
  : m_labelX(labelX),
    m_labelY(labelY),
    m_unitsX(unitsX),
    m_unitsY(unitsY),
    m_pointer_to_GetValueX(px),
    m_pointer_to_GetValueY(py),
    m_hists2D(m_labelX, labelY, xaxisX, xaxisY, nbinsX, nbinsY, xminX, xminY, xmaxX, xmaxY),
    m_is_true(is_true)
{}


// CTOR -- variable binning
Variable2D::Variable2D(const std::string labelX, const std::string labelY, 
		       const std::string xaxisX, const std::string xaxisY,
		       const std::string unitsX, const std::string unitsY,
		       const TArrayD& bins_arrayX, const TArrayD& bins_arrayY, 
		       PointerToCVUniverseFunction px,
		       PointerToCVUniverseFunction py,
		       const bool is_true)
  : m_labelX(labelX),
    m_labelY(labelY),
    m_unitsX(unitsX),
    m_unitsY(unitsY),
    m_pointer_to_GetValueX(px),
    m_pointer_to_GetValueY(py),
    m_hists2D(m_labelX, labelY, xaxisX, xaxisY, bins_arrayX, bins_arrayY),
    m_is_true(is_true)
{}

// CTOR -- using other declared histograms

Variable2D::Variable2D(const Variable* x,
               const Variable* y)
  : m_labelX(x->m_label),
    m_labelY(y->m_label),
    m_unitsX(x->m_units),
    m_unitsY(y->m_units),
    m_pointer_to_GetValueX(x->m_aux_pointer_to_GetValue),
    m_pointer_to_GetValueY(y->m_aux_pointer_to_GetValue),
    m_hists2D(x->m_label, y->m_label, x->m_xlabel, y->m_xlabel, x->m_hists.m_bins_array, y->m_hists.m_bins_array),
    m_is_true(x->m_is_true)
{}

Variable2D::Variable2D(const std::string name,
               const Variable* x,
               const Variable* y)

  : m_labelX(name),
    m_labelY(y->m_label),
    m_unitsX(x->m_units),
    m_unitsY(y->m_units),
    m_pointer_to_GetValueX(x->m_aux_pointer_to_GetValue),
    m_pointer_to_GetValueY(y->m_aux_pointer_to_GetValue),
    m_hists2D(name, y->m_label, x->m_xlabel, y->m_xlabel, x->m_hists.m_bins_array, y->m_hists.m_bins_array),
    m_is_true(x->m_is_true)
{}

// GetValue defines this variable
double Variable2D::GetValueX (const CVUniverse& universe, 
                          const int hadron_ID) const { 
  return m_pointer_to_GetValueX(universe); 
}

double Variable2D::GetValueY (const CVUniverse& universe,
                          const int hadron_ID) const {
  return m_pointer_to_GetValueY(universe);
}


// Histogram Initialization
template<typename T>
void Variable2D::InitializeAllHists(T systematic_univs, T systematic_univs_truth) {
  m_hists2D.InitializeAllHists(systematic_univs, systematic_univs_truth);
}


template<typename T>
void Variable2D::InitializeSidebandHists(T systematic_univs) {
  m_hists2D.InitializeSidebandHists(systematic_univs);
}


void Variable2D::InitializeStackedHists(){
  m_hists2D.InitializeStackedHists();
}

void Variable2D::InitializeDataHists(){
  m_hists2D.InitializeDataHists();
}


// Histogram Access
template <typename T>
PlotUtils::MnvH2D* Variable2D::GetStackComponentHist(T type) const {
  std::map<T, PlotUtils::MnvH2D*> stack_map = m_hists2D.GetStackMap(type);
  return stack_map[type];
}


// Save with the object names that hists were initialized with
void Variable2D::WriteMCHists(TFile& fout) const {
  fout.cd();
  m_hists2D.m_selection_mc.hist -> Write();
  m_hists2D.m_bg.hist           -> Write();
  m_hists2D.m_bg_loW.hist       -> Write();
  m_hists2D.m_bg_midW.hist      -> Write();
  m_hists2D.m_bg_hiW.hist       -> Write();
  m_hists2D.m_effnum.hist       -> Write();
  m_hists2D.m_effden.hist       -> Write();
/*if (NameX() == sidebands::kFitVarString) {
    m_hists2D.m_wsidebandfit_sig.hist ->Write();
    m_hists2D.m_wsidebandfit_loW.hist ->Write();
    m_hists2D.m_wsidebandfit_midW.hist->Write();
    m_hists2D.m_wsidebandfit_hiW.hist ->Write();
  }*/
//if (!m_is_true && NameX() != sidebands::kFitVarString)
//  m_hists2D.m_migration.hist->Write();
}


void Variable2D::LoadDataHistsFromFile(TFile& fin) {
  m_hists2D.LoadDataHistsFromFile(fin);
}


void Variable2D::LoadMCHistsFromFile(TFile& fin, UniverseMap& error_bands) {
  m_hists2D.LoadMCHistsFromFile(fin, error_bands, m_is_true);
}


template <typename T>
TObjArray Variable2D::GetStackArray(T type) const {
  if(std::is_same<T, WType>::value)                     return m_hists2D.m_stacked_w.m_hist_array;
  else if(std::is_same<T, SignalBackgroundType>::value) return m_hists2D.m_stacked_sigbg.m_hist_array;
  else if(std::is_same<T, WBackgroundType>::value)      return m_hists2D.m_stacked_wbg.m_hist_array;
  else if(std::is_same<T, MesonBackgroundType>::value)  return m_hists2D.m_stacked_mesonbg.m_hist_array;
  else if(std::is_same<T, HadronType>::value)           return m_hists2D.m_stacked_hadron.m_hist_array;
  else if(std::is_same<T, FSParticleType>::value)       return m_hists2D.m_stacked_fspart.m_hist_array;
  else if(std::is_same<T, ChannelType>::value)          return m_hists2D.m_stacked_channel.m_hist_array;
  else if(std::is_same<T, NPionsType>::value)           return m_hists2D.m_stacked_npi.m_hist_array;
  else if(std::is_same<T, NPi0Type>::value)             return m_hists2D.m_stacked_npi0.m_hist_array;
  else if(std::is_same<T, NPipType>::value)             return m_hists2D.m_stacked_npip.m_hist_array;
  else if(std::is_same<T, WSidebandType>::value)        return m_hists2D.m_stacked_wsideband.m_hist_array;
  else if(std::is_same<T, CoherentType>::value)         return m_hists2D.m_stacked_coherent.m_hist_array;
  else{
    std::cerr << "GetStackArray: Unknown truth type.\n";
    std::exit(1);
  }
}


#endif // Variable_cxx
