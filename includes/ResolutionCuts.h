#ifndef ResolutionCuts_h
#define ResolutionCuts_h 

#include "CVUniverse.h"
#include "PlotUtils/ChainWrapper.h"

#include <iostream>
#include <utility>

class CVUniverse;
namespace PlotUtils
{
  class ChainWrapper;
}

int GetNResPions(CVUniverse *cv);
int GetNResPions(PlotUtils::ChainWrapper* chw);
std::vector<int> getResolutionPionIndices( CVUniverse *cv );
std::vector<int> getResolutionPionIndices( PlotUtils::ChainWrapper* chw );
bool hasGoodMomentum( CVUniverse* cv, int iProng );
bool hasGoodMomentum( PlotUtils::ChainWrapper* chw, int iProng );
bool passResolutionNumberNodeCuts( CVUniverse* cv, int iMM );
bool passResolutionLLRCuts( CVUniverse* cv, int iMM );
bool passResolutionMichelCuts( CVUniverse* cv, int iMM );
bool passResolutionMichelCuts( PlotUtils::ChainWrapper* chw, int iMM );
bool passResolutionMichelCuts( double avg_dist, double end_dist );
bool passResolutionNodeCuts2( CVUniverse* cv, int iMM );
bool passResolutionNodeCuts3( CVUniverse* cv, int iMM );
bool passResolutionNodeCuts4( CVUniverse* cv, int iMM );
bool passResolutionNodeCuts5( CVUniverse* cv, int iMM );
bool passResolutionNode01Cuts( CVUniverse* cv, int iMM );
bool passResolutionNumberNodeCuts( PlotUtils::ChainWrapper* chw, int iMM );
bool passResolutionLLRCuts( PlotUtils::ChainWrapper* chw, int iMM );
bool passResolutionNodeCuts2( PlotUtils::ChainWrapper* chw, int iMM );
bool passResolutionNodeCuts3( PlotUtils::ChainWrapper* chw, int iMM );
bool passResolutionNodeCuts4( PlotUtils::ChainWrapper* chw, int iMM );
bool passResolutionNodeCuts5( PlotUtils::ChainWrapper* chw, int iMM );
bool passResolutionNode01Cuts( PlotUtils::ChainWrapper* chw, int iMM );
bool passResolutionNumberNodeCuts( int pion_nNodes);
bool passResolutionLLRCuts( double LLR_Score);
bool passResolutionNodeCuts2( double pion_lastnode_Q );
bool passResolutionNodeCuts3( double pion_lastnode_Q );
bool passResolutionNodeCuts4( double pion_lastnode_Q );
bool passResolutionNodeCuts5( double pion_lastnode_Q );
bool passResolutionNode01Cuts( double pion_lastnode_Q);
