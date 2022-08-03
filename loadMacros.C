void loadMacros() {
    TString path(
      TString::Format("%s/CC-CH-pip-ana/xsec/", gSystem->Getenv("TOPDIR")));
    TString oldpath = gSystem->GetIncludePath();
    oldpath += " -I";
    oldpath += path;
    gSystem->SetIncludePath(oldpath);
    gSystem->CompileMacro("makeCrossSectionMCInputs.C", "k");
    gSystem->CompileMacro("crossSectionDataFromFile.C", "k");
    gSystem->CompileMacro("plotCrossSectionFromFile.C", "k");
}