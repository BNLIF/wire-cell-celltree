{
    TString include = ".include ";
    TString load = ".L ";
    
    TString prefix = "../event";
    gROOT->ProcessLine( include + prefix );
    gROOT->ProcessLine( load + prefix + "/CellEvent.cc+" );
}