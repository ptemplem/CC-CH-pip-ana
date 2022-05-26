//#ifndef plotting_functions2D_h
//#define plotting_functions2D_h

#include <iostream>
#include <sstream>

#include "../includes/myPlotStyle.h"
#include "PlotUtils/MnvColors.h"
#include "PlotUtils/MnvH1D.h"
#ifndef MNVROOT6
#define MNVROOT6
#include "PlotUtils/MnvPlotter.h"
#endif
#include "PlotUtils/MnvVertErrorBand.h"
#include "SignalDefinition.h"
#include "Systematics.h"  // namespace systematics
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TList.h"
#include "TPad.h"
#include "TPaveStats.h"
//#include "TStyle.h"
#include "TText.h"
#include "Variable.h"
#include "Variable2D.h"
#include "PlotUtils/GridCanvas.h"

class Variable;
class Variable2D;

//==============================================================================
// Container class
//==============================================================================
class EventSelectionPlotInfo2D {
 public:

  // Constructor with Var2D
  EventSelectionPlotInfo2D(Variable2D* variable2D, float mc_pot, float data_pot,
                         bool do_frac_unc, bool do_cov_area_norm,
                         bool include_stat, SignalDefinition signal_definition)
      : m_mnv_plotter(kCCNuPionIncStyle),
        m_variable2D(variable2D),
        m_mc_pot(mc_pot),
        m_data_pot(data_pot),
        m_do_frac_unc(do_frac_unc),
        m_do_cov_area_norm(do_cov_area_norm),
        m_include_stat(include_stat),
        m_signal_definition(signal_definition) {
    m_do_frac_unc_str = m_do_frac_unc ? "Frac" : "Abs";
    m_do_cov_area_norm_str = m_do_cov_area_norm ? "CovAreaNorm" : "";
  }

  // Constructor without var
  EventSelectionPlotInfo2D(float mc_pot, float data_pot, bool do_frac_unc,
                         bool do_cov_area_norm, bool include_stat,
                         SignalDefinition signal_definition)
      : m_mnv_plotter(kCCNuPionIncStyle),
        m_variable2D(nullptr),
        m_mc_pot(mc_pot),
        m_data_pot(data_pot),
        m_do_frac_unc(do_frac_unc),
        m_do_cov_area_norm(do_cov_area_norm),
        m_include_stat(include_stat),
        m_signal_definition(signal_definition) {
    m_do_frac_unc_str = m_do_frac_unc ? "Frac" : "Abs";
    m_do_cov_area_norm_str = m_do_cov_area_norm ? "CovAreaNorm" : "";
  }

  // Members
  MnvPlotter m_mnv_plotter;
  Variable2D* m_variable2D;
  float m_mc_pot;
  float m_data_pot;
  bool m_do_frac_unc;
  bool m_do_cov_area_norm;
  bool m_include_stat;
  SignalDefinition m_signal_definition;
  std::string m_do_frac_unc_str;
  std::string m_do_cov_area_norm_str;

  // Add X label
  void SetXLabel(PlotUtils::MnvH2D* hist) {
    std::string label =
        m_variable2D->m_hists2D.m_xlabelX + " (" + m_variable2D->m_unitsX + ")";
    if (hist) hist->GetXaxis()->SetTitle(label.c_str());
  }

  // Add X label
  void SetXLabel(TH2* hist) {
    std::string label =
        m_variable2D->m_hists2D.m_xlabelX + " (" + m_variable2D->m_unitsX + ")";
    if (hist) hist->GetXaxis()->SetTitle(label.c_str());
  }

  // Add title to the MnvPlotter
  void SetTitle() {
    if (!gPad)
      throw std::runtime_error("Need a TCanvas. Please make one first.");
    std::string title = GetSignalName(m_signal_definition);
    m_mnv_plotter.AddHistoTitle(title.c_str());
  }

  void SetTitle(std::string title) {
    if (!gPad)
      throw std::runtime_error("Need a TCanvas. Please make one first.");
    m_mnv_plotter.AddHistoTitle(title.c_str());
  }
};

//-------------------------------------
//------Event Selection Plots----------
//-------------------------------------

void PlotVar_Selection2D(std::vector<TH2D*> hists, TH2D* data, std::string var, std::string xlabelX, std::string xlabelY, std::string unitsX, std::string unitsY){
    

    std::vector<int> bins = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    auto legend = new TLegend(0.69,0.15,0.98,0.33);
    legend->AddEntry(hists[0], "MC", "l");
    legend->AddEntry(hists[1], "BG", "f");
    legend->AddEntry(data, "Data", "p");

    if (false) {
      hists[0]->Scale(1.,"width");
      hists[1]->Scale(1.,"width");
      data->Scale(1.,"width");
    }
    PlotUtils::GridCanvas* StackedPlot =new  PlotUtils::GridCanvas(Form("2D_Sel_%s",var.c_str()), 4, 3, 1400, 600);
    StackedPlot->SetRightMargin(0.01);
    StackedPlot->SetLeftMargin(0.1);
    StackedPlot->ResetPads();
//    data->SetMarkerStyle(8);
//    StackedPlot->DrawStack( hists, "HIST", false, bins, false, NULL);
    StackedPlot->DrawOneHist( hists[0], "HIST E2", false, bins, false, NULL);
    StackedPlot->DrawOneHist( hists[1], "SAME HIST", false, bins, false, NULL);
    StackedPlot->DrawOneHist(data, "SAME E1", false, bins, false, NULL);
    StackedPlot->DrawBinRanges(data, 2, bins, Form("%s (%s)",xlabelY.c_str(), unitsY.c_str()), 0.03, ".2f", 2);
    StackedPlot->SetXTitle(Form("%s (%s)",xlabelX.c_str(), unitsX.c_str()));
    StackedPlot->SetYTitle("Events");
    //StackedPlot->BuildLegend();
    StackedPlot->ResetPads();
    legend->Draw();
    StackedPlot->Draw();
    StackedPlot->Print(Form("2D_Sel_%s.png",var.c_str()));
 
}




