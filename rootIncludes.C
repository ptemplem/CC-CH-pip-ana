{
    gInterpreter->AddIncludePath( "${TOPDIR}/CC-CH-pip-ana/includes" );
    gInterpreter->AddIncludePath( "${TOPDIR}/CC-CH-pip-ana" );
    gInterpreter->AddIncludePath( "${TOPDIR}/MAT" );
    gInterpreter->AddIncludePath( "${TOPDIR}/MAT/PlotUtils" );
    gInterpreter->AddIncludePath( "${TOPDIR}/MAT-MINERvA" );
    gInterpreter->AddIncludePath( "${TOPDIR}/MAT-MINERvA/calculators" );
    gInterpreter->AddIncludePath( "${TOPDIR}/MAT-MINERvA/universes" );
    gInterpreter->AddIncludePath( "${TOPDIR}/MAT-MINERvA/utilities" );
    gInterpreter->AddIncludePath( "${TOPDIR}/MAT-MINERvA/weighters" );
    gSystem->Load( gSystem->ExpandPathName("$TOPDIR/opt/lib/libMAT.so") );
    gSystem->Load( gSystem->ExpandPathName("$TOPDIR/opt/lib/libMAT-MINERvA.so") );
}