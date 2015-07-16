{
    gROOT->Reset();
    gROOT->ProcessLine(".x loadClasses.C" );

    CellEvent *ev = new CellEvent("../data/3.0/celltree-e.root");
    // get total number of events stored
    int nEvent = ev->nEvent;

    // example: get info from the first(0) event:
    ev->GetEntry(0);

    ev->PrintInfo(1);

    cout << "raw channels hit: " << ev->raw_nChannel << endl;

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

    // get the signal from the first(0) raw channel
    TH1F *h2 = (TH1F*)ev->raw_wf->At(144);
    TCanvas c2;
    h2->Draw();
    for (int i=0; i<9600; i++) {
        if (fabs(h2->GetBinContent(i)-2048)>0.1) {  // most
            cout << h2->GetBinContent(i) << " ";
        }
    }
    cout << endl;

    // get the simide energy deposit
    // cout << ev->simide_numElectrons->at(1) << endl;
    // cout << ev->simide_size << endl;
    float totalE = 0;
    for (int i=0; i<ev->simide_size; i++) {
        totalE += ev->simide_numElectrons->at(i);
    }
    cout << "total number of electrons: " << totalE << endl;

}
