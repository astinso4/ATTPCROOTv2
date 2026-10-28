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

   TString outfname= Form("kine_%dto%d_delEvEattempt.root", runStart, runEnd); 
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
   /*TFile *cutPIDFilep = new TFile("./cutFiles/PID_proton.root", "READ"); 
   TCutG *cutPIDproton = (TCutG *)cutPIDFilep->Get("cutPIDproton");
   cutPIDFilep->Close();

   TFile *cutPIDFile = new TFile("./cutFiles/PID_deuteron.root", "READ"); 
   TCutG *cutPIDdeuteron = (TCutG *)cutPIDFile->Get("cutPIDdeuteron");
   cutPIDFile->Close();

   TFile *cutPIDFile2 = new TFile("./cutFiles/PID_3He.root", "READ"); 
   TCutG *cutPID3He = (TCutG *)cutPIDFile2->Get("cutPID3He");
   cutPIDFile2->Close();
   // Original cut files 
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
   TCutG *cutPIDproton  = makeOpenCut();
   TCutG *cutPIDdeuteron = makeOpenCut();
   TCutG *cutPID3He     = makeOpenCut();

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
   TH2F *histRangeVThetaLAB = new TH2F("histRangeVThetaLAB", "histRangeVThetaLAB", 180, 0, 180, 1030, 0, 1030);
   TH2F *histTrackThetaLABRange3He = new TH2F("histTrackThetaLABRange3He", "histTrackThetaLABRange3He", 1030, 0, 1030, 180, 0, 180);
   TH2F *histTrackThetaLABRange3HeSiGaggBe = new TH2F("histTrackThetaLABRange3HeSiGaggBe", "histTrackThetaLABRange3HeSiGaggBe", 506, 0, 1030, 90, 0, 180);
   TH2F *histEstimatedKinEVThetaLABTotal = new TH2F("histEstimatedKinEVThetaLABTotal", "histEstimatedKinEVThetaLABTotal", 360, 0, 180, 500, 0, 20);
   TH2F *histEstimatedKinEVThetaLAB3He = new TH2F("histEstimatedKinEVThetaLAB3He", "histEstimatedKinEVThetaLAB3He", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLAB3HeSiBe = new TH2F("histEstimatedKinEVThetaLAB3HeSiBe", "histEstimatedKinEVThetaLAB3HeSiBe", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLAB3HeSiGaggBe = new TH2F("histEstimatedKinEVThetaLAB3HeSiGaggBe", "histEstimatedKinEVThetaLAB3HeSiGaggBe", 360, 0, 180, 300, 0, 30);
   TH2F *histEstimatedKinEVThetaLAB3HeSiGaggBetrack1 = new TH2F("histEstimatedKinEVThetaLAB3HeSiGaggBetrack1", "histEstimatedKinEVThetaLAB3HeSiGaggBetrack1", 180, 0, 180, 300, 0, 30);
   TH2F *histEstimatedKinEVThetaLAB3HeSiGaggBetrack2 = new TH2F("histEstimatedKinEVThetaLAB3HeSiGaggBetrack2", "histEstimatedKinEVThetaLAB3HeSiGaggBetrack2", 180, 0, 180, 300, 0, 30);
   TH2F *histEstimatedKinEVThetaLABkinett = new TH2F("histEstimatedKinEVThetaLABkinett", "histEstimatedKinEVThetaLABkinett", 180, 0, 180, 300, 0, 30);
   TH2F *histThetaLABThetaLAB = new TH2F("histThetaLABThetaLAB", "histThetaLABThetaLAB", 360, 0, 180, 360, 0, 180);
   TH2F *histThetaLABThetaLABKine = new TH2F("histThetaLABThetaLABKine", "histThetaLABThetaLABKine", 500, 0, 50, 150, 0, 15);
   TH2F *histThetaLABThetaLAB3HeSiGaggBe = new TH2F("histThetaLABThetaLAB3HeSiGaggBe", "histThetaLABThetaLAB3HeSiGaggBe", 360, 0, 180, 360, 0, 180);
   TH2F *histEstimatedKinEVThetaLAB2H  = new TH2F("histEstimatedKinEVThetaLAB2H", "histEstimatedKinEVThetaLAB2H", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLAB1H  = new TH2F("histEstimatedKinEVThetaLAB1H", "histEstimatedKinEVThetaLAB1H", 180, 0, 180, 250, 0, 20);      
   TH2F *histEstimatedKinEVThetaLABSiB  = new TH2F("histEstimatedKinEVThetaLABSiB", "histEstimatedKinEVThetaLABSiB", 180, 0, 180, 250, 0, 20);      
   TH2F *histEstimatedKinEVThetaLABSiC  = new TH2F("histEstimatedKinEVThetaLABSiC", "histEstimatedKinEVThetaLABSiC", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLABSi13B  = new TH2F("histEstimatedKinEVThetaLABSi13B", "histEstimatedKinEVThetaLABSi13B", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLABSiBe  = new TH2F("histEstimatedKinEVThetaLABSiBe", "histEstimatedKinEVThetaLABSiBe", 180, 0, 180, 250, 0, 20);
   TH2F *histEstimatedKinEVThetaLABSiGaggBe  = new TH2F("histEstimatedKinEVThetaLABSiGaggBe", "histEstimatedKinEVThetaLABSiGaggBe", 180, 0, 180, 250, 0, 20);
   TH2F *histdEdxVTotalRange = new TH2F("histdEdxVTotalRange", "histdEdxVTotalRange", 500, 0, 1030, 1600, 0, 4000);
   TH2F *histESmallVTotalRange = new TH2F("histESmallVTotalRange", "histESmallVTotalRange", 500, 0, 1030, 1600, 0, 160000);
   TH2F *histEBigVBigRange = new TH2F("histEBigVBigRange", "histEBigVBigRange", 500, 0, 1030, 1600, 0, 160000);

   TH2F *histSiPIDTraceIntegralFronts = new TH2F("histSiPIDTraceIntegralFronts", "histSiPIDTraceIntegralFronts;E_{2}^{Front}[MeV];E_{1}^{Front}[MeV]", 5000, 0, 4000, 5000, 0, 4000);
   TH2F *histSiPIDTraceIntegralBacks = new TH2F("histSiPIDTraceIntegralBacks", "histSiPIDTraceIntegralBacks;E_{2}^{Back}[MeV];E_{1}^{Back}[MeV]", 5000, 0, 4000, 5000, 0, 4000);
   TH2F *histSiPIDTraceIntegralFrontsCAL = new TH2F("histSiPIDTraceIntegralFrontsCAL", "histSiPIDTraceIntegralFrontsCAL;E_{2}^{Front}[MeV];E_{1}^{Front}[MeV]", 5000, 0, 300, 5000, 0, 300);
   TH2F *histSiPIDTraceIntegralBacksCAL = new TH2F("histSiPIDTraceIntegralBacksCAL", "histSiPIDTraceIntegralBacksCAL;E_{2}^{Back}[MeV];E_{1}^{Back}[MeV]", 5000, 0, 300, 5000, 0, 300);
   TH2F *histSiPIDADCMaxFront = new TH2F("histSiPIDADCMaxFront", "histSiPIDADCMaxFront;ADC_{2}^{max}[ADC];ADC_{1}^{max}[ADC]", 4000, 0, 4000, 4000, 0, 4000);
   TH2F *histSiPIDADCMaxBack = new TH2F("histSiPIDADCMaxBack", "histSiPIDADCMaxBack;ADC_{2}^{max}[ADC];ADC_{1}^{max}[ADC]", 4000, 0, 4000, 4000, 0, 4000);
   TH2F *histSiPIDADCMax3He = new TH2F("histSiPIDADCMax3He", "histSiPIDADCMax3He", 4000, 0, 4000, 4000, 0, 4000);
   TH2F *histSiPIDADCMax3HeGaggBe = new TH2F("histSiPIDADCMax3HeGaggBe", "histSiPIDADCMax3HeGaggBe", 4000, 0, 4000, 4000, 0, 4000);
   TH2F *histSiPIDADCMax2H = new TH2F("histSiPIDADCMax2H", "histSiPIDADCMax2H", 1000, 0, 4000, 1000, 0, 4000);
   TH2F *histSiPIDADCMax1H = new TH2F("histSiPIDADCMax1H", "histSiPIDADCMax1H", 1000, 0, 4000, 1000, 0, 4000);
   TH2F *histSiPIDADCMaxGaggBe = new TH2F("histSiPIDADCMaxGaggBe", "histSiPIDADCMaxGaggBe", 1000, 0, 4000, 1000, 0, 4000);
   TH2F *histSiPIDADCMaxGaggLi = new TH2F("histSiPIDADCMaxGaggLi", "histSiPIDADCMaxGaggLi", 1000, 0, 4000, 1000, 0, 4000);

   TH2F *histSi1PIDFrontvBack = new TH2F("histSi1PIDFrontvBack", "histSi1PIDFrontvBack;E_{1}^{Front}[MeV];E_{1}^{Back}[MeV]", 5000, 0, 400, 5000, 0, 400);
   TH2F *histSi2PIDFrontvBack = new TH2F("histSi2PIDFrontvBack", "histSi2PIDFrontvBack;E_{2}^{Front}[MeV];E_{2}^{Back}[MeV]", 5000, 0, 400, 5000, 0, 400);

   TH1F *histSi1FrontE = new TH1F("histSi1FrontE", "histSi1Front;E1_Front;Counts", 5000, 0, 50000);
   TH1F *histSi1BackE = new TH1F("histSi1BackE", "histSi1Back;E1_Back;Counts", 5000, 0, 25000);
   TH1F *histSi2FrontE = new TH1F("histSi2FrontE", "histSi2FrontE;E2_Front;Counts", 5000, 0, 50000);
   TH1F *histSi2BackE = new TH1F("histSi2BackE", "histSi2BackE;E2_Back;Counts", 5000, 0, 25000);
   TH2F *histSi1V2FrontE = new TH2F("histSi1V2FrontE", "histSi1V2FrontE;E1_Front;E2_Front", 5000, 0, 25000, 5000, 0, 25000);

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

   TH1F *histSmallPad = new TH1F("histSmallPad", "histSmallPad;Hit Charge;Counts", 4100, 0, 4100);
   TH1F *histBigPad = new TH1F("histBigPad", "histBigPad;Hit Charge;Counts", 600, 0, 600);

   TH1F *histTrackMultiplicity = new TH1F("histTrackMultiplicity", "Track Multiplicty;Tracks Detected;Counts", 11, -0.5, 10.5);

   TH1F *histVertexCal = new TH1F("histVertexCal", "Reconstructed Vertex;Vertex Z;Counts", 1000, 0, 1000);

   //TH2F *GAGG_IDvE = new TH2F("GAGG_IDvE", "GAGG Energy vs ID;Energy;ID", 8000, 0, 500, 26, 0, 26);
   TH2F *GAGG_IDvE = new TH2F("GAGG_IDvE", "GAGG ID vs Energy;ID;Energy", 26, 0, 26, 300, 0, 500);
   //TH2F *GAGG_IDvE = new TH2F("GAGG_IDvE", "GAGG ID vs Energy;ID;Energy", 26, 0, 26, 8000, 0, 40000);

   //TH3F *oddpos = new TH3F("oddpos", "3D Plot", 100, -150, 150, 100, -150, 150, 100, 0, 1000); 

   TH2F *histPIDSivGAGG = new TH2F("histPIDSivGAGG","histPIDSivGAGG;E_{GAGG}[MeV];E_{1}^{Back}+E_{2}^{Back}[MeV]", 5000, 0, 600, 5000, 0, 2000);
    //Turning on Calibration Mode
    bool calibrationMode = false;
    if (runNums[0] > 5000) {
	    calibrationMode = true;
	    std::cout << "CALIBRATION MODE ON" << std::endl;
    }

   for (int runNum: runNums) {
      if (runNum == 3012) { 
	      continue;
      }
      // Open the digitalization file and get the TTree.
      TString unpackFileName = TString::Format("/home/astinson/e535rawdata/UnpackerTestOutput/run_%04d_GAGG_CAL.root", runNum);
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

      // Open the GAGG file of this run.
      /*TString gaggFileName = TString::Format("/data/sustech/user/public/frib/frib-decode/data/hit%04d.root", runNum);
      TFile *gaggFile = new TFile(gaggFileName, "READ");
      TTree *gaggTree = (TTree *)gaggFile->Get("tree");
      int nGaggEvents = gaggTree->GetEntries();
      std::cout << " Number of unpacked GAGG events in run " << runNum << ": " << nGaggEvents << std::endl;

      Int_t idGagg1[25];
      Int_t idGagg2[16];
      Double_t gaggADCMax1[25];
      Double_t gaggADCMax2[16];

      gaggTree->SetBranchAddress("id_g1", idGagg1);
      gaggTree->SetBranchAddress("id_g2", idGagg2);
      gaggTree->SetBranchAddress("ADC_max_g1", gaggADCMax1);
      gaggTree->SetBranchAddress("ADC_max_g2", gaggADCMax2);*/

      // Loop over events.
      for (int i = 0; i < nUnpackEvents; i++) {
      //for (int i = 0; i < nGaggEvents; i++) {
         unpackReader.Next();
         //gaggTree->GetEntry(i);

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
		 histSiPIDTraceIntegralFronts->Fill(front2E, front1E);
		 histSiPIDTraceIntegralFrontsCAL->Fill(front2E, front1E);
		 if (abs(tdiff)<10){
			 SiSum = back1E + back2E;
			 goodSiPID = true;
		 }
		 if (cutSi13B->IsInside(back2E, back1E)) {
			beamSiGate = true;
		 	histSiPIDTraceIntegralBacks->Fill(back2E, back1E);
		 	histSiPIDTraceIntegralBacksCAL->Fill(back2E, back1E);
		 }
		 histSi1PIDFrontvBack->Fill(front1E, back1E);
		 histSi2PIDFrontvBack->Fill(front2E, back2E);
	 }

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

	 for (int j=0; j<4; ++j) { //Testing for looping over Fronts.
	 	Int_t stripSi1 = siEvent->GetStripFront1(j);
		Int_t stripSi2 = siEvent->GetStripFront2(j);
	 	if (stripSi1 >= 0 ) {
			Double_t maxADCFront1 = siEvent->GetADCMaxFront1(j);
			Double_t EFront1 = siEvent->GetEFront1(j);
			histSi1FrontE->Fill(EFront1);
			histSiStrip->Fill(stripSi1);
			if (cutSi13B->IsInside(back2E, back1E)) {
				histSiStripvE->Fill(stripSi1, EFront1);
			}
			histSiStripvMaxADC->Fill(stripSi1, maxADCFront1);
		}
	 	if (stripSi2 >= 0 ) {
			Double_t maxADCFront2 = siEvent->GetADCMaxFront2(j);
			Double_t EFront2 = siEvent->GetEFront2(j);
			histSi2FrontE->Fill(EFront2);
			histSiStrip->Fill(stripSi2);
			if (cutSi13B->IsInside(back2E, back1E)) {
				histSiStripvE->Fill(stripSi2, EFront2);
			}
			histSiStripvMaxADC->Fill(stripSi2, maxADCFront2);
		}
		if (stripSi1 >= 0 && stripSi2 >= 0) {
			Double_t maxADCFront1 = siEvent->GetADCMaxFront1(j);
			Double_t EFront1 = siEvent->GetEFront1(j);
			Double_t maxADCFront2 = siEvent->GetADCMaxFront2(j);
			Double_t EFront2 = siEvent->GetEFront2(j);
			histSiPIDADCMaxFront->Fill(maxADCFront2, maxADCFront1);
		}
	 }

	 for (int j=0; j<4; ++j) { //Testing for looping over Backs.
	 	Int_t stripSi1 = siEvent->GetStripBack1(j);
		Int_t stripSi2 = siEvent->GetStripBack2(j);
	 	if (stripSi1 >= 0 ) {
			Double_t maxADCBack1 = siEvent->GetADCMaxBack1(j);
			Double_t EBack1 = siEvent->GetEBack1(j);
			histSi1BackE->Fill(EBack1);
			histSiStrip->Fill(stripSi1);
			if (cutSi13B->IsInside(back2E, back1E)) {
				histSiStripvE->Fill(stripSi1, EBack1);
			}
			histSiStripvMaxADC->Fill(stripSi1, maxADCBack1);
		}
	 	if (stripSi2 >= 0 ) {
			Double_t maxADCBack2 = siEvent->GetADCMaxBack2(j);
			Double_t EBack2 = siEvent->GetEBack2(j);
			histSi2BackE->Fill(EBack2);
			histSiStrip->Fill(stripSi2);
			if (cutSi13B->IsInside(back2E, back1E)) {
				histSiStripvE->Fill(stripSi2, EBack2);
			}
			histSiStripvMaxADC->Fill(stripSi2, maxADCBack2);
		}
		if (stripSi1 >= 0 && stripSi2 >= 0) {
			Double_t maxADCBack1 = siEvent->GetADCMaxBack1(j);
			Double_t EBack1 = siEvent->GetEBack1(j);
			Double_t maxADCBack2 = siEvent->GetADCMaxBack2(j);
			Double_t EBack2 = siEvent->GetEBack2(j);
			histSiPIDADCMaxBack->Fill(maxADCBack2, maxADCBack1);
		}
	 }

         //if (!cutSiN->IsInside(maxADCFront2, maxADCFront1)) continue;
         //if (!cutSiC->IsInside(maxADCFront2, maxADCFront1)) continue;

         // Check the GAGG data directly from the GAGG decoder.
         /*Double_t totalMaxADCGagg1{};
         for (int j = 0; j < 25; j++) {
            if (idGagg1[j])
               totalMaxADCGagg1 += gaggADCMax1[j];
         }

         Double_t totalMaxADCGagg2{};
         for (int j = 0; j < 16; j++) {
            if (idGagg1[j])
               totalMaxADCGagg2 += gaggADCMax2[j];
         }

         Double_t totalMaxADCGagg = totalMaxADCGagg1 + totalMaxADCGagg2;
         if (totalMaxADCGagg)
            histGaggPIDADCMax->Fill(totalMaxADCGagg, maxADCFront2);

         //if (!cutGaggC->IsInside(totalMaxADCGagg, maxADCFront2)) continue;*/

         // Check the GAGG data next. For calibration need to add if statement to allow energies past gate
         AtGaggEvent *gaggEvent = (AtGaggEvent *)gaggArray->At(0);

         Int_t multiGagg1 = gaggEvent->GetMultiplicity1();
         Int_t multiGagg2 = gaggEvent->GetMultiplicity2();

         histGaggMultiplicity1->Fill(multiGagg1);
         histGaggMultiplicity2->Fill(multiGagg2);

	 //Create ADC plot for calibration 


         //std::cout << " Multiplicity GAGG1: " << multiGagg1 <<std::endl;

         //if (multiGagg1 != 1) continue;
	 int g_mult = 0;
	 double GAGGSum = 0;
	 int GAGGTime = -1;

         Double_t maxADCGagg1, maxADCGagg2;
         Double_t gaggADC1{-1}, gaggADC2{-1};
         Double_t totalMaxADCGagg{};
         for(int imul = 0; imul < multiGagg1; imul++) {
	    int ID = (int)gaggEvent->GetID1(imul);
	    if (ID == 36) { ID = 24; } //crystal missing. Use crystal behind instead
	    //std::cout << "ID = " << ID << std::endl;
	    if (ID < 0 || ID > 24) { continue; }
	    if (beamSiGate) {
		    //GAGG_IDvE->Fill(gaggEvent->GetE1(imul), ID);
		    GAGG_IDvE->Fill(ID, gaggEvent->GetE1(imul));
	    }
	    GAGGSum += gaggEvent->GetE1(imul);
	    GAGGTime = gaggEvent->GetTS1(imul);
	    g_mult++;
            maxADCGagg1 = gaggEvent->GetADCMax1(imul);
            gaggADC1 += maxADCGagg1;
            totalMaxADCGagg += maxADCGagg1;
            //std::cout << " ADC entry = " << maxADCGagg1 << std::endl;
         }
	 if ((g_mult==1) && goodSiPID) {
		 int tdiff = time1 - GAGGTime;
		 if (tdiff >= 25 && tdiff <= 50){
			 histPIDSivGAGG->Fill(GAGGSum, SiSum);
		 }
	 }

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

            double smallPadCharge{};
            double bigPadCharge{};
            auto braggCurvePairs = track.GetBraggCurveValues();
            auto &hits = track.GetHitArray();
            double rangeInSmallPads{};
	         for (auto &hit: hits) {
               int padNum = hit->GetPadNum();
               int sizeID = map->GetPadSize(padNum);
               if (sizeID == 1) {
		  histSmallPad->Fill(hit->GetCharge());
                  bigPadCharge += hit->GetCharge();
                  continue;
               }
	       histSmallPad->Fill(hit->GetCharge());
               smallPadCharge += hit->GetCharge();
               double currentRangeInSmallPads = pattern->DistanceAlongPattern(hit->GetPosition(), firstPoint);
	       auto pos = hit->GetPosition(); //for 3D plot
	       //oddpos->Fill(pos.X(), pos.Y(), pos.Z());
	       histVertexCal->Fill(pos.Z());
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

            histRangeVThetaLAB->Fill(trackThetaLAB, roughRangeEstimation);
            histdEdxVTotalRange->Fill(roughRangeEstimation, dEdx);
            histEstimatedKinEVThetaLABTotal->Fill(trackThetaLAB, estimatedKinE);

            if (cutPID3He->IsInside(roughRangeEstimation, dEdx))
               histTrackThetaLABRange3He->Fill(roughRangeEstimation, trackThetaLAB);

            if (cutPID3He->IsInside(roughRangeEstimation, dEdx))
               histEstimatedKinEVThetaLAB3He->Fill(trackThetaLAB, estimatedKinE);

            if (cutPIDdeuteron->IsInside(roughRangeEstimation, dEdx))
               histEstimatedKinEVThetaLAB2H->Fill(trackThetaLAB, estimatedKinE);

            if (cutPIDproton->IsInside(roughRangeEstimation, dEdx))
               histEstimatedKinEVThetaLAB1H->Fill(trackThetaLAB, estimatedKinE);

            if (!reachedBigPads)
               histESmallVTotalRange->Fill(roughRangeEstimation, smallPadCharge);
            else
               histEBigVBigRange->Fill(rangeInBigPads, bigPadCharge);

	    if(tracks.size()==2 && trackIndex < 2){
	      thetaLABArray[trackIndex] = trackThetaLAB;
	      estimatedKineEArray[trackIndex] = estimatedKinE;
	      isHe3Array[trackIndex] = cutPID3He->IsInside(roughRangeEstimation, dEdx);
	      trackIndex++;
	    }
         }
         if(tracks.size()==2 && trackIndex==2){
	   histThetaLABThetaLAB->Fill(thetaLABArray[0], thetaLABArray[1]);
	   if(thetaLABArray[0]<=thetaLABArray[1]){
	     if(thetaLABArray[0]<15 && thetaLABArray[1]<50){
	       histThetaLABThetaLABKine->Fill(thetaLABArray[1], thetaLABArray[0]);
	       histEstimatedKinEVThetaLABkinett->Fill(thetaLABArray[0],estimatedKineEArray[0]);
	     }
	   }
	   else if(thetaLABArray[0]>thetaLABArray[1]){
	     if(thetaLABArray[0]<50 && thetaLABArray[1]<15){
	       histThetaLABThetaLABKine->Fill(thetaLABArray[0], thetaLABArray[1]);
	       histEstimatedKinEVThetaLABkinett->Fill(thetaLABArray[1],estimatedKineEArray[1]);
	     }
	   }
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
