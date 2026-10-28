#include <fstream>
#include <iomanip>
//#include <stdio.h>
//#include <windows.h>

TGraph* ReadKinematics(TString kineFile);
Double_t omega(Double_t x, Double_t y, Double_t z);
std::tuple<double, double> kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject);

void kine()
{
   FairRunAna *run = new FairRunAna(); // Forcing a dummy run
   //Bypassing cutfiles edit It just creates a giant cut that allows everything through
   auto makeOpenCut = []() {
      auto *c = new TCutG();
      c->SetVarX("x"); 
      c->SetVarY("y");
      c->SetPoint(0, -1e9, -1e9);
      c->SetPoint(1,  1e9, -1e9);
      c->SetPoint(2,  1e9,  1e9);
      c->SetPoint(3, -1e9,  1e9);
      c->SetPoint(4, -1e9, -1e9);
      return c;
   };
   // File Name Input
   int runStart, runEnd;
   std::cout << "Run Start: ";
   std::cin >> runStart;
   std::cout << "Run End: ";
   std::cin >> runEnd;


   // File iteration
   std::vector<int> runNums;
   for (int i=runStart; i<= runEnd; ++i) {
	   runNums.push_back(i); 
   }
   
   double col_boundaries[6] = {0,28,52,76,99,128};
   double row_boundaries[6] = {0,30,54,78,101,128};
	   
   //TString outfname= Form("kine_%dto%d_delEvEattempt_digisplit_countfixmaybe_rejectplt_multi2.root", runStart, runEnd); 
   TString outfname= Form("kine_%dto%d_multifix_veto_10232026_CUTPID.root", runStart, runEnd); 
   TFile *outfile=new TFile(outfname,"recreate");

   // AtMap to check if a hit belong to a big pad or small pad.
   TString scriptfile = "rcnp_map_size.xml";
   TString dir = getenv("VMCWORKDIR");
   TString mapDir = dir + "/scripts/" + scriptfile;
   AtTpcMap *map = new AtTpcMap();
   map->ParseXMLMap(mapDir.Data());
   map->GeneratePadPlane();

   // Punch through filter.
   double punchThroughThreshold = 20;
   AtTools::AtPunchThroughChecker punchThroughChecker = AtTools::AtPunchThroughChecker();
   punchThroughChecker.SetDistanceThreshold(punchThroughThreshold);

   // ELoss model for kinetic energy estimations.
   double density = 1.3378e-3; // 470Torr
   std::vector<std::tuple<int, int, int>> materialComponents;
   materialComponents.push_back(std::make_tuple(12, 6, 3));
   materialComponents.push_back(std::make_tuple(2, 1, 8));

   std::unique_ptr<AtTools::AtELossCATIMA> eLossModelC3D8_p = std::make_unique<AtTools::AtELossCATIMA>(density, "CATima_C3D8_470Torr_p");
   eLossModelC3D8_p->SetMaterial(materialComponents);
   eLossModelC3D8_p->SetProjectile(1, 1, 1.007825031898);
   eLossModelC3D8_p->SetPDGCode("1000010010");

   std::unique_ptr<AtTools::AtELossCATIMA> eLossModelC3D8_d = std::make_unique<AtTools::AtELossCATIMA>(density, "CATima_C3D8_470Torr_d");
   eLossModelC3D8_d->SetMaterial(materialComponents);
   eLossModelC3D8_d->SetProjectile(2, 1, 2.0135532);
   eLossModelC3D8_d->SetPDGCode("1000020010");

   std::unique_ptr<AtTools::AtELossCATIMA> eLossModelC3D8_3He = std::make_unique<AtTools::AtELossCATIMA>(density, "CATima_C3D8_470Torr_3He");
   eLossModelC3D8_3He->SetMaterial(materialComponents);
   eLossModelC3D8_3He->SetProjectile(3, 2, 3.01602932197);
   eLossModelC3D8_3He->SetPDGCode("1000030020");

   // Cut files.
   /*TFile *cutKineFile = new TFile("./cutFiles/cutKine.root", "READ");
   TCutG *cutKineBProton = (TCutG *)cutKineFile->Get("cutKineBProton");
   cutKineFile->Close();*/
   //Here Amber	
   TFile *cutPIDFilep = new TFile("./cutFiles/PID_proton.root", "READ"); 
   TCutG *cutPIDproton = (TCutG *)cutPIDFilep->Get("cutPIDproton");
   cutPIDFilep->Close();

   TFile *cutPIDFile = new TFile("./cutFiles/PID_deuteron.root", "READ"); 
   TCutG *cutPIDdeuteron = (TCutG *)cutPIDFile->Get("cutPIDdeuteron");
   cutPIDFile->Close();

   TFile *cutPIDFile2 = new TFile("./cutFiles/PID_3He.root", "READ"); 
   TCutG *cutPID3He = (TCutG *)cutPIDFile2->Get("cutPID3He");
   cutPIDFile2->Close();
   /*// Original cut files 
   TFile *cutSiPIDFile = new TFile("./cutFiles/SiPID.root", "READ");
   TCutG *cutSiB = (TCutG *)cutSiPIDFile->Get("cutSiPIDB");
   TCutG *cutSiC = (TCutG *)cutSiPIDFile->Get("cutSiPIDC");
   TCutG *cutSi13B = (TCutG *)cutSiPIDFile->Get("cutSiPID13B");
   TCutG *cutSiBe = (TCutG *)cutSiPIDFile->Get("cutSiPIDBe");
   cutSiPIDFile->Close();
   */

   /*TFile *cutGaggPIDFile = new TFile("./cutFiles/GAGGPID.root", "READ");
   TCutG *cutGaggB = (TCutG *)cutGaggPIDFile->Get("cutGAGGPIDB");
   TCutG *cutGaggBe = (TCutG *)cutGaggPIDFile->Get("cutGAGGPIDBe");
   TCutG *cutGaggLi = (TCutG *)cutGaggPIDFile->Get("cutGAGGPIDLi");
   TCutG *cutGaggHe = (TCutG *)cutGaggPIDFile->Get("cutGAGGPIDHe");
   cutGaggPIDFile->Close();

   TFile *cutSi13BPIDFile = new TFile("./cutFiles/cutFile_B13Testv2.root", "READ");
   TCutG *cutSi13B = (TCutG *)cutSi13BPIDFile->Get("B13_Test");
   cutSi13BPIDFile->Close();

   TFile *cutSiBPIDFile = new TFile("./cutFiles/cutFile_BTestv2.root", "READ");
   TCutG *cutSiB = (TCutG *)cutSiBPIDFile->Get("BTest");
   cutSiBPIDFile->Close();

   TFile *cutSiBePIDFile = new TFile("./cutFiles/cutFile_BeTestv2.root", "READ");
   TCutG *cutSiBe = (TCutG *)cutSiBePIDFile->Get("BeTest");
   cutSiBePIDFile->Close();

   TFile *cutSiCPIDFile = new TFile("./cutFiles/cutFile_CTestv2.root", "READ");
   TCutG *cutSiC = (TCutG *)cutSiCPIDFile->Get("CTest");
   cutSiCPIDFile->Close();*/


   //TFile *cutSi13BPIDFile = new TFile("./cutFiles/cut20to20_13B.root", "READ");
   //TCutG *cutSi13B = (TCutG *)cutSi13BPIDFile->Get("cut20to20_13B");
   //cutSi13BPIDFile->Close();

   //Bypass cut files edit
   /*TCutG *cutPIDproton  = makeOpenCut();
   TCutG *cutPIDdeuteron = makeOpenCut();
   TCutG *cutPID3He     = makeOpenCut();*/

   TCutG *cutSiB   = makeOpenCut();
   TCutG *cutSiC   = makeOpenCut();
   TCutG *cutSi13B = makeOpenCut();
   TCutG *cutSiBe  = makeOpenCut();

   TCutG *cutGaggB  = makeOpenCut();
   TCutG *cutGaggBe = makeOpenCut();
   TCutG *cutGaggLi = makeOpenCut();
   TCutG *cutGaggHe = makeOpenCut();
   
   // Kinematic curve
   TGraph* kinecurve_3HeGS = new TGraph("ang_lab_cm_3HeGS.txt","%lg %*s %lg");
   kinecurve_3HeGS->SetLineWidth(2);
   kinecurve_3HeGS->SetLineColor(kRed);
   TGraph* ang_lab_cm_3HeGS = new TGraph("ang_lab_cm_3HeGS.txt","%lg %*s %lg");// 5 deg pitch in theta_cm
   ang_lab_cm_3HeGS->SetMarkerStyle(8);
   ang_lab_cm_3HeGS->SetMarkerSize(1);
   TLegend *legend = new TLegend(0.6,0.2,0.85,0.5);
   legend->AddEntry(kinecurve_3HeGS,"3He G.S.","l");
   legend->AddEntry(ang_lab_cm_3HeGS,"5 deg pitch in #theta_{cm}","p");
   legend->SetFillColor(0);

   outfile->cd();
   
   // Histogram definitions.
   TH2F *histSiPIDTraceIntegralFronts = new TH2F("histSiPIDTraceIntegralFronts", "histSiPIDTraceIntegralFronts;E_{2}^{Front}[MeV];E_{1}^{Front}[MeV]", 5000, 0, 300, 5000, 0, 300);
   TH2F *histSiPIDTraceIntegralBacks = new TH2F("histSiPIDTraceIntegralBacks", "histSiPIDTraceIntegralBacks;E_{2}^{Back}[MeV];E_{1}^{Back}[MeV]", 5000, 0,3000, 5000, 0, 300);
   
   TH2F *histSiPID_GAGGveto = new TH2F("histSiPID_GAGGveto", "histSiPID_GAGGveto;E_{2}^{Back}[MeV];E_{1}^{Back}[MeV]", 5000, 0, 300, 5000, 0, 250);

   TH2F *histSi1PIDFrontvBack = new TH2F("histSi1PIDFrontvBack", "histSi1PIDFrontvBack;E_{1}^{Front}[MeV];E_{1}^{Back}[MeV]", 5000, 0, 400, 5000, 0, 400);
   TH2F *histSi2PIDFrontvBack = new TH2F("histSi2PIDFrontvBack", "histSi2PIDFrontvBack;E_{2}^{Front}[MeV];E_{2}^{Back}[MeV]", 5000, 0, 400, 5000, 0, 400);

   TH1F *histSi1FrontE = new TH1F("histSi1FrontE", "histSi1Front;E1_Front;Counts", 5000, 0, 300);
   TH1F *histSi1BackE = new TH1F("histSi1BackE", "histSi1Back;E1_Back;Counts", 5000, 0, 300);
   TH1F *histSi2FrontE = new TH1F("histSi2FrontE", "histSi2FrontE;E2_Front;Counts", 5000, 0, 300);
   TH1F *histSi2BackE = new TH1F("histSi2BackE", "histSi2BackE;E2_Back;Counts", 5000, 0, 300);
   TH2F *histSi1V2FrontE = new TH2F("histSi1V2FrontE", "histSi1V2FrontE;E1_Front;E2_Front", 5000, 0, 300, 5000, 0, 300);

   TH1F *histSiMultiplicityFront1 = new TH1F("histSiMultiplicityFront1", "histSiMultiplicityFront1", 10, 0, 10);
   TH1F *histSiMultiplicityFront2 = new TH1F("histSiMultiplicityFront2", "histSiMultiplicityFront2", 10, 0, 10);

   TH1F *histSiMultiplicityBack1 = new TH1F("histSiMultiplicityBack1", "histSiMultiplicityBack1", 10, 0, 10);
   TH1F *histSiMultiplicityBack2 = new TH1F("histSiMultiplicityBack2", "histSiMultiplicityBack2", 10, 0, 10);

   TH1F *histGaggMultiplicity1 = new TH1F("histGaggMultiplicity1", "histGaggMultiplicity1", 25, 0, 25);
   TH1F *histGaggMultiplicity2 = new TH1F("histGaggMultiplicity2", "histGaggMultiplicity2", 16, 0, 16);
   TH2F *histGaggPIDADCMax = new TH2F("histGaggPIDADCMax", "histGaggPIDADCMax", 5000, 0, 10000, 4000, 0, 4000);
   TH2F *histGaggPIDADCMax3HeSiB = new TH2F("histGaggPIDADCMax3HeSiB", "histGaggPIDADCMax3HeSiB", 5000, 0, 10000, 4000, 0, 4000);
   TH2F *histGaggPIDADCMaxSiB = new TH2F("histGaggPIDADCMaxSiB", "histGaggPIDADCMaxSiB", 5000, 0, 10000, 4000, 0, 4000);

   TH1F *histSiStrip = new TH1F("histSiStrip", "histSiStrip;Strip;Counts", 550, 0, 550);
   TH2F *histSiStripvE = new TH2F("histSiStripvE", "histSiStripvE;Strip;E", 550, 0, 550, 5000, 0, 3500);
   TH2F *histSiStripvMaxADC = new TH2F("histSiStripvMaxADC", "histSiStripvMaxADC;Strip;ADC", 550, 0, 550, 5000, 0, 5000);
   TH2F *Si1Position = new TH2F("Si1Position", "Si1Position; X (strip); Y (strip)", 128, 0, 128, 128, 0, 128);
   TH2F *Si2Position = new TH2F("Si2Position", "Si2Position; X (strip); Y (strip)", 128, 0, 128, 128, 0, 128);
   
   TH2F *HitPattern1 = new TH2F("HitPattern1", "HitPattern1; X (strip); Y (strip)", 128, 0, 128, 128, 0, 128);
   TH2F *HitPattern2 = new TH2F("HitPattern2", "HitPattern2; X (strip); Y (strip)", 128, 0, 128, 128, 0, 128);

   TH1F *histSmallPad = new TH1F("histSmallPad", "histSmallPad;Hit Charge;Counts", 4100, 0, 4100);
   TH1F *histBigPad = new TH1F("histBigPad", "histBigPad;Hit Charge;Counts", 600, 0, 600);

   TH1F *histTrackMultiplicity = new TH1F("histTrackMultiplicity", "Track Multiplicty;Tracks Detected;Counts", 11, -0.5, 10.5);

   TH1F *histVertexCal = new TH1F("histVertexCal", "Reconstructed Vertex;Vertex Z;Counts", 1000, 0, 1000);

   //TH2F *GAGG_IDvE = new TH2F("GAGG_IDvE", "GAGG Energy vs ID;Energy;ID", 8000, 0, 500, 26, 0, 26);
   TH2F *GAGG_IDvE = new TH2F("GAGG_IDvE", "GAGG ID vs Energy;ID;Energy", 26, 0, 26, 300, 0, 500);
   //TH2F *GAGG_IDvE = new TH2F("GAGG_IDvE", "GAGG ID vs Energy;ID;Energy", 26, 0, 26, 8000, 0, 40000);

   //TH3F *oddpos = new TH3F("oddpos", "3D Plot", 100, -150, 150, 100, -150, 150, 100, 0, 1000); 

   TH1F *histSiTDiff = new TH1F("histSiTDiff", "Si1v2 TDiff;Time Difference;Counts", 1024, -512, 512);
   TH1F *histSi_delX = new TH1F("histSi_delX", "Si_{1}-S_{2} x;#Delta x;Counts", 5, 0, 5);
   TH1F *histSi_delY = new TH1F("histSi_delY", "Si_{1}-S_{2} y;#Delta y;Counts", 5, 0, 5);
   TH2F *histSiIDTDiff_1 = new TH2F("histSiIDTDiff_1", "Si1 ID TDiff;Time Difference;Strip ID", 1024, -512, 512, 128, 0, 128);
   TH2F *histSiIDTDiff_2 = new TH2F("histSiIDTDiff_2", "Si2 ID TDiff;Time Difference;Strip ID", 1024, -512, 512, 128, 0, 128);
   TH2F *histSiGAGGIDTDiff = new TH2F("histSiGAGGIDTDiff", "Si1 vs GAGG ID TDiff;Time Difference;GAGG ID", 1024, -512, 512, 24, 0, 24);
   TH2F *histSiGAGGIDTDiff_GATE = new TH2F("histSiGAGGIDTDiff_GATE", "Si1 vs GAGG ID TDiff GATE;Time Difference;GAGG ID", 1024, -512, 512, 24, 0, 24);
   
   TH2F *histPIDSivGAGG_FULL = new TH2F("histPIDSivGAGG_FULL","histPIDSivGAGG_FULL;E_{GAGG}[MeV];E_{1}^{Back}+E_{2}^{Back}[MeV]", 200, 0, 300, 4000, 0, 400);
   TH2F *histPIDSivGAGG = new TH2F("histPIDSivGAGG","histPIDSivGAGG;E_{GAGG}[MeV];E_{1}^{Back}+E_{2}^{Back}[MeV]", 200, 0, 300, 4000, 0, 400);
   TH2F *histPIDSivGAGG_reject = new TH2F("histPIDSivGAGG_reject","histPIDSivGAGG_reject;E_{GAGG}[MeV];E_{1}^{Back}+E_{2}^{Back}[MeV]", 200, 0, 400, 2000, 0, 350);
   TH2F *histPIDSivGAGG_timegate = new TH2F("histPIDSivGAGG_timegate","histPIDSivGAGG_timegate;E_{GAGG}[MeV];E_{1}^{Back}+E_{2}^{Back}[MeV]", 400, 0, 70, 2000, 0, 80);
   //TH2F *histPIDSivGAGG_timegate = new TH2F("histPIDSivGAGG_timegate","histPIDSivGAGG_timegate;E_{GAGG}[MeV];E_{1}^{Back}+E_{2}^{Back}[MeV]", 400, 0, 400, 2000, 0, 350);
   TH2F *histPIDSivGAGG_1906 = new TH2F("histPIDSivGAGG_1906","histPIDSivGAGG_1906;E_{GAGG}[MeV];E_{1}^{Back}+E_{2}^{Back}[MeV]", 400, 0, 400, 2000, 0, 350);
   TH2F *histPIDSivGAGG_1905 = new TH2F("histPIDSivGAGG_1905","histPIDSivGAGG_1905;E_{GAGG}[MeV];E_{1}^{Back}+E_{2}^{Back}[MeV]", 400, 0, 400, 2000, 0, 350);
   TH1F *histGMULT = new TH1F("histGMULT", "histGMULT", 11,-0.5, 10.5);

   //Kinematic plots
   TH1F *histTrackTheta = new TH1F("histTrackTheta", "histTrackTheta", 180, 0, 180);
   TH2F *histEstimatedKinEVThetaLABTotal = new TH2F("histEstimatedKinEVThetaLABTotal", "histEstimatedKinEVThetaLABTotal", 360, 0, 180, 500, 0, 20);
   TH2F *histdEdxVTotalRange = new TH2F("histdEdxVTotalRange", "histdEdxVTotalRange", 500, 0, 400, 2000, 0, 3000);

   /*TH2F *Si1PosGAGG[40];
   TH2F *Si2PosGAGG[40];
   for (int i=0; i<40; ++i) {
	   Si1PosGAGG[i] = new TH2F(Form("Si1PosGAGG%i", i), "Si1 Position;X(strip);Y(strip)", 128, 0, 128, 128, 0, 128);
	   Si2PosGAGG[i] = new TH2F(Form("Si2PosGAGG%i", i), "Si2 Position;X(strip);Y(strip)", 128, 0, 128, 128, 0, 128);
   }*/

    //Turning on Calibration Mode
    bool calibrationMode = false;
    if (runNums[0] > 5000) {
	    calibrationMode = true;
	    std::cout << "CALIBRATION MODE ON" << std::endl;
    }

   for (int runNum: runNums) {
      if (runNum == 3012) { //junk run 
	      continue;
      }
      if (runNum == 3095) { //junk run 
	      continue;
      }
      if (runNum >= 3096 && runNum <= 3108) {
	      std::cout << "Beam tuning run! Run Number: " << runNum << std::endl;
	      continue;
      }
      // Open the digitalization file and get the TTree.
      TString unpackFileName = TString::Format("/home/astinson/e535rawdata/UnpackerTestOutput/run_%04d_FULL_CAL.root", runNum);
      TFile *unpackFile = new TFile(unpackFileName, "READ");
      TTree *unpackTree = (TTree *)unpackFile->Get("cbmsim");
      int nUnpackEvents = unpackTree->GetEntries();
      std::cout << " Number of unpacked events in run " << runNum << ": " << nUnpackEvents << std::endl;
      int nEventsWith2Tracks = 0;
      int nEventsWith3He = 0;
      // Creare the TTreeReader to read the AtTrackingEvents and simulation.
      TTreeReader unpackReader("cbmsim", unpackFile);
      TTreeReaderValue<TClonesArray> siArray(unpackReader, "AtSiEvent");
      TTreeReaderValue<TClonesArray> gaggArray(unpackReader, "AtGaggEvent");
      TTreeReaderValue<TClonesArray> patternArray(unpackReader, "AtPatternEvent");

      // Loop over events.
      for (int i = 0; i < nUnpackEvents; i++) {
      //for (int i = 0; i < nGaggEvents; i++) {
         unpackReader.Next();
         // Check the Si data first.
         AtSiEvent *siEvent = (AtSiEvent *)siArray->At(0);

         Int_t multiplicityFront1 = siEvent->GetMultiplicityFront1();
         Int_t multiplicityFront2 = siEvent->GetMultiplicityFront2();

         Int_t multiplicityBack1 = siEvent->GetMultiplicityBack1();
         Int_t multiplicityBack2 = siEvent->GetMultiplicityBack2();

         histSiMultiplicityFront1->Fill(multiplicityFront1);
         histSiMultiplicityFront2->Fill(multiplicityFront2);

         histSiMultiplicityBack1->Fill(multiplicityBack1);
         histSiMultiplicityBack2->Fill(multiplicityBack2);
	 
	 Int_t strip1X = siEvent->GetXStrip1();
	 Int_t strip1Y = siEvent->GetYStrip1();
	 Float_t front1E = siEvent->GetEnergyFront1();
	 Float_t back1E = siEvent->GetEnergyBack1();

         //std::cout << " front1E: " << front1E <<std::endl; //calibration debug
         //std::cout << " back1E: " << back1E <<std::endl;


	 Int_t strip2X = siEvent->GetXStrip2();
	 Int_t strip2Y = siEvent->GetYStrip2();
	 Float_t front2E = siEvent->GetEnergyFront2();
	 Float_t back2E = siEvent->GetEnergyBack2();

	 Int_t time1 = siEvent->GetTimeFront1();
	 Int_t time2 = siEvent->GetTimeFront2();

	 double SiSum = 0;

	 bool beamSiGate = false;
	 bool goodSiPID = false;
	 
	 if (front1E > 20 && back1E > 20) {
		 Si1Position->Fill(strip1X, strip1Y);
	 }

	 if (front2E > 20 && back2E > 20) {
		 Si2Position->Fill(strip2X, strip2Y);
	 }

	 if (front1E > 1 && back1E > 1 && front2E > 1 && back2E > 1) {
		 int tdiff = time1 - time2;
		 //HitPattern1->Fill(strip1X, strip1Y);
		 //HitPattern2->Fill(strip2X, strip2Y);
		 //histSiPIDTraceIntegralFronts->Fill(front2E, front1E);
		 histSi1PIDFrontvBack->Fill(front1E, back1E);
		 histSi2PIDFrontvBack->Fill(front2E, back2E);
		 histSiIDTDiff_1->Fill(time1, strip1Y);
		 histSiIDTDiff_2->Fill(time2, strip2Y);
		 histSi1BackE->Fill(back1E);
		 histSi2BackE->Fill(back2E);
		 histSi1FrontE->Fill(front1E);
		 histSi2FrontE->Fill(front2E);
		 if (abs(tdiff)<10){
			 SiSum = back1E + back2E;
		 	 histSiPIDTraceIntegralFronts->Fill(front2E, front1E);
		 	 histSiPIDTraceIntegralBacks->Fill(back2E, back1E);
			 goodSiPID = true;
		 }
		 histSi_delX->Fill(abs(strip1X-strip2X));
		 histSi_delY->Fill(abs(strip1Y-strip2Y));
		 /*if (cutSi13B->IsInside(back2E, back1E)) {
			beamSiGate = true;
		 	histSiPIDTraceIntegralBacks->Fill(back2E, back1E);
		 	histSiPIDTraceIntegralBacksCAL->Fill(back2E, back1E);
		 }*/
		 //SiSum = back1E + back2E; //AMBER TIME
	 }

	 histSiTDiff->Fill(time1-time2);

	 if (calibrationMode) { //these are specific to the calibration files
		 if (runNum > 5155) { //DSSD1
			 if (multiplicityFront1 != 1 || multiplicityBack1 != 1) continue;
		 } 
		 if (runNum >= 5150 && runNum <= 5154) { //DSSD2
			 if (multiplicityFront2 != 1 || multiplicityBack2 != 1) continue;
		 }
	 }
	 else {
		 if ((multiplicityFront1 != 1 || multiplicityBack1 != 1) || (multiplicityFront2 != 1 || multiplicityBack2 != 1)) continue;
	 }

	 /*for (int j=0; j<4; ++j) {
		 Int_t stripSi1 = siEvent->GetStripFront1(j);
		 Int_t stripSi2 = siEvent->GetStripFront2(j);
		 if (stripSi1 >= 0) {
			 histSiIDTDiff_1->Fill(time1, stripSi1);
		 }
		 if (stripSi2 >= 0) {
			 histSiIDTDiff_2->Fill(time2, stripSi2);
		 }

	 }*/

	 if (calibrationMode) { //these are specific to the calibration files
		 if (runNum > 5155) { //DSSD1
			 if (multiplicityFront1 != 1 || multiplicityBack1 != 1) continue;
		 } 
		 if (runNum >= 5150 && runNum <= 5154) { //DSSD2
			 if (multiplicityFront2 != 1 || multiplicityBack2 != 1) continue;
		 }
	 }

         // Check the GAGG data next. For calibration need to add if statement to allow energies past gate
         AtGaggEvent *gaggEvent = (AtGaggEvent *)gaggArray->At(0);

         Int_t multiGagg1 = gaggEvent->GetMultiplicity1();
         Int_t multiGagg2 = gaggEvent->GetMultiplicity2();

         histGaggMultiplicity1->Fill(multiGagg1);
         histGaggMultiplicity2->Fill(multiGagg2);

         //if (multiGagg1 != 1) continue;
	 int g_mult = 0;
	 double GAGGSum = 0;
	 int GAGGTime = -1;
	 int GAGGID = -1;
	 int GAGGX, GAGGY;
	 bool CleanGAGG = false;

         Double_t maxADCGagg1, maxADCGagg2;
         Double_t gaggADC1{-1}, gaggADC2{-1};
         Double_t totalMaxADCGagg{};

	 //std::cout << "Using Gagg Multiplicity: " << multiGagg1 << std::endl;//testing multiplicity

	 //std::cout << "GAGG event" << std::endl;
         for (int imul = 0; imul < multiGagg1; imul++) {
	    int ID = (int)gaggEvent->GetID1(imul);
            maxADCGagg1 = gaggEvent->GetADCMax1(imul);
	    if (ID == 36) { ID = 24; } //crystal missing. Use crystal behind instead
	    if (ID < 0 || ID > 24) { 
		    continue; 
	    }
	    if (maxADCGagg1 < 40) {
		    CleanGAGG = true;
		    continue; 
	    }
	    GAGGID = ID;
	    GAGGX = 4-GAGGID/5;
	    GAGGY = 4-GAGGID%5;
	    if (beamSiGate) {
		    //GAGG_IDvE->Fill(gaggEvent->GetE1(imul), ID);
		    GAGG_IDvE->Fill(ID, gaggEvent->GetE1(imul));
	    }
	    g_mult++;
	    //std::cout << "g_mult = " << g_mult << std::endl; //testing multiplicity
	    GAGGSum += gaggEvent->GetE1(imul);
	    GAGGTime = gaggEvent->GetTS1(imul);
            gaggADC1 += maxADCGagg1;
            totalMaxADCGagg += maxADCGagg1;
            //std::cout << imul << "  " << ID << " ADC entry = " << maxADCGagg1 << std::endl;
         }
	 //std::cout << "g_mult = " << g_mult << std::endl;
	 histGMULT->Fill(g_mult);

	 if ((g_mult==0) && goodSiPID && CleanGAGG) {
		 //if ((back2E > 700) || (back1E > 700)) {continue; }
		 histSiPID_GAGGveto->Fill(back2E, back1E);
		 HitPattern1->Fill(strip1X, strip1Y);
		 HitPattern2->Fill(strip2X, strip2Y);
	    	 //std::cout << "g_mult:" << g_mult << "; X strip: " << strip1X << "; Y strip: " << strip1Y << "; Energy (back 2): " << back2E << "; Energy (front 2): " << front2E << std::endl; //testing multiplicity
	 }

	 if ((g_mult==1) && goodSiPID) {
		 int tdiff = time1 - GAGGTime;
		 histSiGAGGIDTDiff->Fill(tdiff, GAGGID);
		 histPIDSivGAGG->Fill(GAGGSum, SiSum);
		 //bool pos_gate = true; //test
		 bool pos_gate = false;
		 if ((col_boundaries[GAGGX]-2) < strip1X && strip1X < (col_boundaries[GAGGX+1]+2)) {
			if ((row_boundaries[GAGGY]-2) < strip1Y && strip1Y < (row_boundaries[GAGGY+1]+2)) {
				pos_gate = true;
			}
		 }
		 if (tdiff >= 25 && tdiff <= 50) {
			 //Si1PosGAGG[GAGGID]->Fill(strip1X, strip1Y);
			 //Si2PosGAGG[GAGGID]->Fill(strip2X, strip2Y);
		 	 if (pos_gate) {
		 		histSiGAGGIDTDiff_GATE->Fill(tdiff, GAGGID);
			 	histPIDSivGAGG_timegate->Fill(GAGGSum, SiSum);
				if (GAGGID < 16 && GAGGID !=6) {
					histPIDSivGAGG_1906->Fill(GAGGSum, SiSum);
				}
				else if (GAGGID > 15) {
					histPIDSivGAGG_1905->Fill(GAGGSum, SiSum);
				}
			 }
			 //else {
			 //	 histPIDSivGAGG_reject->Fill(GAGGSum, SiSum);
			 //}
		 }
	 }
	 else {
	 	histPIDSivGAGG_reject->Fill(GAGGSum, SiSum);
	 }
	 histPIDSivGAGG_FULL->Fill(GAGGSum, SiSum); 

         // First, we obtain some rough kinematics just by using the AtPatternEvent.
         AtPatternEvent *patternEvent = (AtPatternEvent *)patternArray->At(0);
         if (!patternEvent) continue;

         // We want to focus on events with 2 or less tracks for now.
         auto &tracks = patternEvent->GetTrackCand();
	 histTrackMultiplicity->Fill(tracks.size());
	 //if (tracks.size() > 4) continue;
         if (tracks.size() > 2) continue;
         //if (tracks.size() > 1) continue;
	 if (tracks.size()==2) nEventsWith2Tracks++;
	 int trackIndex = 0;
	 double thetaLABArray[2] = {0.0,0.0};
	 double estimatedKineEArray[2] = {0.0,0.0};
	 bool isHe3Array[2] = {false,false};
	 // Iterate over AtTracks and extract their kinematics.                                      
	 for (auto &track: tracks) {                                                                    
		 bool isPunchThrough = punchThroughChecker.IsPunchThrough(&track);
		 if (isPunchThrough) continue;
		 
		 auto *pattern = track.GetPattern();

		 auto firstPoint = track.GetFirstPoint();
		 auto lastPoint = track.GetLastPoint();
		 //cut penetrate particles
		 //if (lastPoint.Z()>900) continue;
		 //std::cout<<lastPoint.Z()<<std::endl;
		 double roughRangeEstimation = pattern->DistanceAlongPattern(lastPoint, firstPoint);
		 double trackThetaLAB = track.GetGeoTheta() * 180 / TMath::Pi();
            	 double trackPhi = track.GetGeoPhi() * 180 / TMath::Pi();

		 histTrackTheta->Fill(trackThetaLAB);

            	 double smallPadCharge{};
            	 double bigPadCharge{};
            	 auto braggCurvePairs = track.GetBraggCurveValues();
            	 auto &hits = track.GetHitArray();
            	 double rangeInSmallPads{};
		 for (auto &hit:hits) {
			 int padNum = hit->GetPadNum();
			 int sizeID = map->GetPadSize(padNum);
			 if (sizeID == 1) {
				 bigPadCharge += hit->GetCharge();
				 continue;
			 }
			 smallPadCharge += hit->GetCharge();
			 double currentRangeInSmallPads = pattern->DistanceAlongPattern(hit->GetPosition(), firstPoint);
			 auto pos = hit->GetPosition();
			 if (currentRangeInSmallPads > rangeInSmallPads)
				 rangeInSmallPads = currentRangeInSmallPads;
		 }
		 double dEdx = smallPadCharge / rangeInSmallPads;

		 double rangeInBigPads = roughRangeEstimation - rangeInSmallPads;
		 bool reachedBigPads = true;
		 if (rangeInBigPads / roughRangeEstimation < 0.05)
			 reachedBigPads = false;

		 double estimatedKinE{0.1};
	 
		 if (cutPID3He->IsInside(roughRangeEstimation, dEdx)) {
			 while (eLossModelC3D8_3He->GetRange(estimatedKinE) < roughRangeEstimation)
				 estimatedKinE += 0.01;
		 } else if (cutPIDdeuteron->IsInside(roughRangeEstimation, dEdx)) {
			 while (eLossModelC3D8_d->GetRange(estimatedKinE) < roughRangeEstimation)
				 estimatedKinE += 0.01;
		 } else if (cutPIDproton->IsInside(roughRangeEstimation, dEdx)) {
			 while (eLossModelC3D8_p->GetRange(estimatedKinE) < roughRangeEstimation)
				 estimatedKinE += 0.01;
		 } else {
			 estimatedKinE = 0;
		 }
		 histEstimatedKinEVThetaLABTotal->Fill(trackThetaLAB, estimatedKinE);
		 histdEdxVTotalRange->Fill(roughRangeEstimation, dEdx);
      	}
      }
      // Close files.
      unpackFile->Close();
      //gaggFile->Close();
   }

   outfile->Write();

   // Kinematic lines.
   TGraph *kine_dd = ReadKinematics("./kineFiles/17N_dd_gs.txt");
   TGraph *kine_d3He = ReadKinematics("./kineFiles/13B_d3He_gs.txt");
   TGraph *kine_d3HeEx2_2 = ReadKinematics("./kineFiles/13B_d3He_Ex2_2.txt");
   TGraph *kine_d3HeEx2_7 = ReadKinematics("./kineFiles/13B_d3He_Ex2_7.txt");
   TGraph *kine_d3He_tt = ReadKinematics("./kineFiles/13B_d3He_tt.txt");
   TGraph *kine_12C12C = ReadKinematics("./kineFiles/17N_12C12C_gs.txt");

   
   outfile->Close();
}

TGraph* ReadKinematics(TString kineFile)
{
   Double_t *ThetaCMS = new Double_t[20000];
   Double_t *ThetaLabRec = new Double_t[20000];
   Double_t *EnerLabRec = new Double_t[20000];
   Double_t *ThetaLabSca = new Double_t[20000];
   Double_t *EnerLabSca = new Double_t[20000];
   Double_t *MomLabRec = new Double_t[20000];

   std::ifstream *kineStr = new std::ifstream(kineFile.Data());
   Int_t numKin = 0;

   if (!kineStr->fail()){
      while (!kineStr->eof()){
	//         *kineStr >> ThetaCMS[numKin] >> ThetaLabRec[numKin] >> EnerLabRec[numKin] >>
	//                     ThetaLabSca[numKin] >> EnerLabSca[numKin];
         *kineStr >> ThetaLabRec[numKin] >> EnerLabRec[numKin];
         numKin++;
      }
   } else if (kineStr->fail())
      std::cout << " Warning : No Kinematics file found for this reaction!" << std::endl;

   TGraph *kine = new TGraph(numKin, ThetaLabRec, EnerLabRec);
   return kine;
}

Double_t omega(Double_t x, Double_t y, Double_t z)
{
   return sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

std::tuple<double, double>
kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject)
{

   // in this definition: m1(projectile); m2(target); m3(ejectile); and m4(recoil);
   double Et1 = K_proj + m1;
   double Et2 = m2;
   double Et3 = K_eject + m3;
   double Et4 = Et1 + Et2 - Et3;
   double m4_ex, Ex, theta_cm;
   double s, t, u; //---Mandelstam variables

   s = pow(m1, 2) + pow(m2, 2) + 2 * m2 * Et1;
   u = pow(m2, 2) + pow(m3, 2) - 2 * m2 * Et3;

   m4_ex = sqrt((cos(thetalab) * omega(s, pow(m1, 2), pow(m2, 2)) * omega(u, pow(m2, 2), pow(m3, 2)) -
                 (s - pow(m1, 2) - pow(m2, 2)) * (pow(m2, 2) + pow(m3, 2) - u)) /
                   (2 * pow(m2, 2)) +
                s + u - pow(m2, 2));
   Ex = m4_ex - m4;

   t = pow(m2, 2) + pow(m4_ex, 2) - 2 * m2 * Et4;

   // for inverse kinematics Note: this angle corresponds to the recoil
    theta_cm = TMath::Pi() - acos((pow(s, 2) + s * (2 * t - pow(m1, 2) - pow(m2, 2) - pow(m3, 2) - pow(m4_ex, 2)) +
                                  (pow(m1, 2) - pow(m2, 2)) * (pow(m3, 2) - pow(m4_ex, 2))) /
                                 (omega(s, pow(m1, 2), pow(m2, 2)) * omega(s, pow(m3, 2), pow(m4_ex, 2))));

   /*theta_cm = acos((pow(s, 2) + s * (2 * u - pow(m1, 2) - pow(m2, 2) - pow(m3, 2) - pow(m4_ex, 2)) +
                                  (pow(m1, 2) - pow(m2, 2)) * (pow(m4_ex, 2) - pow(m3, 2))) /
                                 (omega(s, pow(m1, 2), pow(m2, 2)) * omega(s, pow(m4_ex, 2), pow(m3, 2))));*/

   theta_cm = theta_cm * TMath::RadToDeg();
   return std::make_tuple(Ex, theta_cm);
}
