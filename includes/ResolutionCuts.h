#ifndef ResolutionCuts_H
#define ResolutionCuts_H

#include "includes/CVUniverse.h"

bool ResolutionCuts::passResolutionCuts( CVUniverse* cv, int iMM )
{
  if( cv->GetVecElemInt("matched_michel_idx", iMM) < 0) return false;

  if( !hasGoodMomentum( cv, iMM ) ) return false;
  if( !passResolutionMichelCuts(     cv, iMM ) ) return false;
  if( !passResolutionLLRCuts(        cv, iMM ) ) return false;
  if( !passResolutionNumberNodeCuts( cv, iMM ) ) return false;
  if( !passResolutionNode01Cuts(     cv, iMM ) ) return false;
  if( !passResolutionNodeCuts2(       cv, iMM, 2 ) ) return false;
  if( !passResolutionNodeCuts3(       cv, iMM, 3 ) ) return false;
  if( !passResolutionNodeCuts4(       cv, iMM, 4 ) ) return false;
  if( !passResolutionNodeCuts5(       cv, iMM, 5 ) ) return false;
  return true;
}

bool ResolutionCuts::passResolutionMichelCuts( CVUniverse *cv, int iMM )
{
  double matched_michel_avg_dist = cv->GetMichelProngAvg(iMM);
  double matched_michel_end_dist = cv->GetMichelProngEnd(iMM);

  if( !passMichelECut( cv, iMM ) ) return false;
  return passNewResolutionMichelCuts( matched_michel_onv_dist, matched_michel_avg_dist, matched_michel_end_dist );
}

bool ResolutionCuts::passResolutionMichelCuts( double avg_dist, double end_dist )
{
  if(end_dist>0 && 0.0 <= end_dist && end_dist <= 9.25 ) return true;
  if(avg_dist>0 && 0.0 <= avg_dist && avg_dist <= 16.25) return true;
  return false;

}
bool ResolutionCuts::passResolutionLLRCuts( CVUniverse* cv, RecoPionIdx hadron )
{
  double LLR_Score = cv->GetLLRScore( hadron);

  return passResolutionLLRCuts( LLR_Score );
}

bool ResolutionCuts::passResolutionLLRCuts( PlotUtils::ChainWrapper* chw, int ievent, RecoPionIdx iMM )
{
  double LLR_Score = chw->GetValue("MasterAnaDev_hadron_piFit_scoreLLR",ievent,iMM);

  return passResolutionLLRCuts( LLR_Score );
}

bool ResolutionCuts::passResolutionLLRCuts( double LLR_Score )
{
  if( 0.0 <= LLR_Score && LLR_Score <= 50.0 ) return true;
  return false;
}

bool ResolutionCuts::passResolutionNumberNodeCuts( CVUniverse* cv, RecoPionIdx hadron )
{
  int pion_nNodes = cv->GetNNodes(hadron);
  return passResolutionNumberNodeCuts( pion_nNodes );
}

bool ResolutionCuts::passResolutionNumberNodeCuts( PlotUtils::ChainWrapper* chw, int ievent, RecoPionIdx iMM )
{
  int pion_nNodes = chw->GetValue("MasterAnaDev_pion_nNodes",ievent,iMM);
  return passResolutionNumberNodeCuts(  pion_nNodes );
}

bool ResolutionCuts::passResolutionNumberNodeCuts( int pion_nNodes )
{
  if( 0.0 <= pion_nNodes && pion_nNodes <= 44 ) return true;
  return false;
}

bool ResolutionCuts::passResolutionNode01Cuts( CVUniverse* cv, RecoPionIdx iMM )
{
  double pion_lastnode_Q = cv->GetEnode01(iMM);

  return passResolutionNode01Cuts( pion_lastnode_Q );
}

bool ResolutionCuts::passResolutionNode01Cuts( PlotUtils::ChainWrapper* chw, int ievent, RecoPionIdx iMM )
{
  double pion_lastnode_Q = chw->GetValue("MasterAnaDev_pion_lastnode_Q0",ievent,iMM)+chw->GetValue("MasterAnaDev_pion_lastnode_Q1",ievent,iMM);

  return passResolutionNode01Cuts( pion_lastnode_Q );
}

bool ResolutionCuts::passResolutionNode01Cuts( double pion_lastnode_Q )
{
  if( pion_lastnode_Q<0 ) return true;//Doesn't exist.  Don't want to cut on nodes that don't exist

  if( 6.0 <= pion_lastnode_Q && pion_lastnode_Q<= 32.0 ) return true;
  return false;
}

bool NukeCCPionCuts::passResolutionNodeCuts2( CVUniverse* cv, int iMM, int node)
{
  double pion_lastnode_Q = cv->GetLastNodeQ(iMM, node);

  return passResolutionNodeCuts2( pion_lastnode_Q);
}

bool NukeCCPionCuts::passResolutionNodeCuts2( PlotUtils::ChainWrapper* chw, int ievent, int iMM, int node)
{
  double pion_lastnode_Q = chw->GetValue(Form("NukeCCPion_pion_lastnode_Q%d",node),ievent,iMM);

  return passResolutionNodeCuts2( pion_lastnode_Q );
}

bool NukeCCPionCuts::passResolutionNodeCuts2( double pion_lastnode_Q )
{
  if( pion_lastnode_Q<0 ) return true;//Doesn't exist.  Don't want to cut on nodes that don't exist

  if( m_cut_list[Form("res_node_%d",node)][tarNum].lower_cut <= pion_lastnode_Q && pion_lastnode_Q<= m_cut_list[Form("res_node_%d",node)][tarNum].upper_cut ) return true;
  return false;
}

bool NukeCCPionCuts::passResolutionNodeCuts3( CVUniverse* cv, int iMM)
{
  double pion_lastnode_Q = cv->GetLastNodeQ(iMM, node);

  return passResolutionNodeCuts3( tarNum, node, pion_lastnode_Q);
}

bool NukeCCPionCuts::passResolutionNodeCuts3( PlotUtils::ChainWrapper* chw, int ievent, int iMM, int node)
{
  double pion_lastnode_Q = chw->GetValue(Form("NukeCCPion_pion_lastnode_Q%d",node),ievent,iMM);

  return passResolutionNodeCuts3( tarNum, node, pion_lastnode_Q );
}

bool NukeCCPionCuts::passResolutionNodeCuts3( int tarNum, int node, double pion_lastnode_Q )
{
  if( pion_lastnode_Q<0 ) return true;//Doesn't exist.  Don't want to cut on nodes that don't exist

  if( m_cut_list[Form("res_node_%d",node)][tarNum].lower_cut <= pion_lastnode_Q && pion_lastnode_Q<= m_cut_list[Form("res_node_%d",node)][tarNum].upper_cut ) return true;
  return false;
}

bool NukeCCPionCuts::passResolutionNodeCuts4( CVUniverse* cv, int iMM, int node)
{
  double pion_lastnode_Q = cv->GetLastNodeQ(iMM, node);

  return passResolutionNodeCuts4( tarNum, node, pion_lastnode_Q);
}

bool NukeCCPionCuts::passResolutionNodeCuts4( PlotUtils::ChainWrapper* chw, int ievent, int iMM, int node)
{
  double pion_lastnode_Q = chw->GetValue(Form("NukeCCPion_pion_lastnode_Q%d",node),ievent,iMM);

  return passResolutionNodeCuts4( tarNum, node, pion_lastnode_Q );
}

bool NukeCCPionCuts::passResolutionNodeCuts4( int tarNum, int node, double pion_lastnode_Q )
{
  if( pion_lastnode_Q<0 ) return true;//Doesn't exist.  Don't want to cut on nodes that don't exist

  if( m_cut_list[Form("res_node_%d",node)][tarNum].lower_cut <= pion_lastnode_Q && pion_lastnode_Q<= m_cut_list[Form("res_node_%d",node)][tarNum].upper_cut ) return true;
  return false;
}

bool NukeCCPionCuts::passResolutionNodeCuts5( CVUniverse* cv, int iMM, int node)
{
  double pion_lastnode_Q = cv->GetLastNodeQ(iMM, node);

  return passResolutionNodeCuts5( tarNum, node, pion_lastnode_Q);
}

bool NukeCCPionCuts::passResolutionNodeCuts5( PlotUtils::ChainWrapper* chw, int ievent, int iMM, int node)
{
  double pion_lastnode_Q = chw->GetValue(Form("NukeCCPion_pion_lastnode_Q%d",node),ievent,iMM);

  return passResolutionNodeCuts5( tarNum, node, pion_lastnode_Q );
}

bool NukeCCPionCuts::passResolutionNodeCuts5( int tarNum, int node, double pion_lastnode_Q )
{
  if( pion_lastnode_Q<0 ) return true;//Doesn't exist.  Don't want to cut on nodes that don't exist

  if( m_cut_list[Form("res_node_%d",node)][tarNum].lower_cut <= pion_lastnode_Q && pion_lastnode_Q<= m_cut_list[Form("res_node_%d",node)][tarNum].upper_cut ) return true;
  return false;
}


