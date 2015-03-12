{
    gROOT->Reset();
    gROOT->ProcessLine(".x loadClasses.C" );

    CellEvent *ev = new CellEvent("../data/celltree.root");
    // get total number of events stored
    int nEvent = ev->nEvent;

    // for (int i=0; i<nEvent; i++){
    //     ev->GetEntry(i);
    //     ev->PrintInfo(1);
    // }

    // example: get info from the first(0) event:
    ev->GetEntry(0);
    int nChannel = ev->calib_nChannel;  // number of hit channels
    // loop through hit channel ID
    for (int i=0; i<nChannel; i++) {
        cout << ev->calib_channelId->at(i) << " ";
    }
    cout << endl;

    // get the signal from the first(0) hit channel
    TH1F *h = (TH1F*)ev->calib_wf->At(0);
    h->Draw();
    for (int i=0; i<9600; i++) {
        if (fabs(h->GetBinContent(i))>0.1) {  // most 
            cout << h->GetBinContent(i) << " ";
        }
    }
    cout << endl;


}