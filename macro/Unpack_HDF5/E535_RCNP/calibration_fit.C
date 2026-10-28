void calibration_fit(){
        int runNum;
        std::cout << "Calibration File: ";
        std::cin >> runNum;
        
        TString CalFile = TString::Format("kineC_hists_%dto%d_stripsTest.root", runNum, runNum);

        TFile *file = new TFile(CalFile);
        TH2F *histSiStripvE = (TH2F*)file->Get<TH2F>("histSiStripvE");

        std::cout << "Begin fitting for gain strip-by-strip calibrations" << std::endl;

        std::ofstream out_file;
        std::string filename = (Form("parfile_%d.csv", runNum));
        out_file.open(filename);

        //out_file << "Strip,Gain\n";

        TF1 *func = new TF1("func", "[0]*exp(-pow((x-[1]*5.5)/[2],2))+[3]*exp(-pow((x-[1]*5.8)/[2],2))+[4]", 0, 1);

        int binMax, binLow, binHi, stripMax, stripMin = 0;
        float amp1, amp2, g, o, c = 0.0;
        double xmin, xmax = 0.0;

        TCanvas *c1 = new TCanvas("c1", "fits", 800,600);
        gErrorIgnoreLevel = kWarning;
        c1->Print((Form("fits_%d.pdf[", runNum)));
        if (runNum == 5160){
          stripMin = 1;
          stripMax = 256;
        }
        if (runNum == 5152){
          stripMin = 257;
          stripMax = 515;
        }
        for (int i = stripMin; i <= stripMax ; i++){
                TH1D *proj = histSiStripvE->ProjectionY(Form("proj_%d",i), i, i);
                if (proj->GetEntries() < 20) continue;
                //Find maximum value for range
                binMax = proj->GetMaximumBin();
                binLow = binMax - 15;
                binHi = binMax + 15;
                xmin = proj->GetXaxis()->GetBinCenter(binLow);
                xmax = proj->GetXaxis()->GetBinCenter(binHi);
                double xguess = proj->GetXaxis()->GetBinCenter(binMax);
                func->SetParameters(600,xguess/5.5,0.5,400,0);
                proj->Fit(func,"Q", "", xmin, xmax);
                amp1 = func->GetParameter(0);
                amp2 = func->GetParameter(3);
                g = func->GetParameter(1);
                g = 1/g;
                o = func->GetParameter(2);
                c = func->GetParameter(4);
                out_file << i << "," << g << "\n";
                proj->GetXaxis()->SetRangeUser(0,100);
                proj->Draw();
                func->Draw("same");
                c1->Print((Form("fits_%d.pdf", runNum)));
        }
        c1->Print((Form("fits_%d.pdf]", runNum)));
        gErrorIgnoreLevel = kInfo;
        out_file.close();

}
