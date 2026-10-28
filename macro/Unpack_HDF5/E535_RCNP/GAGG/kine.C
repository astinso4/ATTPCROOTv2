#include <fstream>
#include <iomanip>

TGraph* ReadKinematics(TString kineFile);
Double_t omega(Double_t x, Double_t y, Double_t z);
std::tuple<double, double> kine_2b(Double_t m1, Double_t m2, Double_t m3, Double_t m4, Double_t K_proj, Double_t thetalab, Double_t K_eject);

void kine(int r=-1)
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
  //   FILE NAME
  int runNo;
  if (r==-1) {
    std::cout << "Run: ";
    std::cin >> runNo;
  }
  else {
    runNo = r;
  }

  TString outfname= Form("run%d.root", runNo); 
  //read cutfiles
  
  TFile *beamcutfile = new TFile("BeamCut.root");
  TCutG *beamcut = beamcutfile->Get<TCutG>("BeamCut");

  int ind;
  double offset, gain;
  double gagg_cal[40]={1.0};
  for (int i=0; i<40; ++i) { 
    gagg_cal[i] = 1.0;
  }
  std::ifstream gaggcalfile("gagg_cal.dat");
  while (gaggcalfile >> ind >> offset >> gain) {
    gagg_cal[ind] = gain;
  }

  double si_cal[512]={1.0};
  for (int i=0; i<512; ++i) { 
    si_cal[i] = 1.0;
  }
  std::ifstream sicalfile("si_cal.dat");
  while (sicalfile >> ind >> offset >> gain) {
    si_cal[ind] = gain;
  }

  double col_boundaries[6]={0,28,52,76,99,128};
  double row_boundaries[6]={0,30,54,78,101,128};

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

  // Histogram definitions.

  double SiMax = 500;
  double GAGGMax = 8000;
  TH2F *SiPIDTraceIntegralFronts = new TH2F("SiPIDTraceIntegralFronts", "SiPIDTraceIntegralFronts;E_{2}^{Front}[MeV];E_{1}^{Front}[MeV]", 5000, 0, SiMax, 5000, 0, SiMax);
  TH2F *SiPIDTraceIntegralBacks = new TH2F("SiPIDTraceIntegralBacks", "SiPIDTraceIntegralBacks;E_{2}^{Back}[MeV];E_{1}^{Back}[MeV]", 5000, 0, SiMax, 5000, 0, SiMax);
  TH2F *SiPIDADCMaxFront = new TH2F("SiPIDADCMaxFront", "SiPIDADCMaxFront;ADC_{2}^{max}[ADC];ADC_{1}^{max}[ADC]", 4000, 0, 4000, 4000, 0, 4000);
  TH2F *SiPIDADCMaxBack = new TH2F("SiPIDADCMaxBack", "SiPIDADCMaxBack;ADC_{2}^{max}[ADC];ADC_{1}^{max}[ADC]", 4000, 0, 4000, 4000, 0, 4000);
  TH2F *SiPIDNoGAGG = new TH2F("SiPIDNoGAGG", "SiPIDNoGAGG;E_{2}^{Back}[MeV];E_{1}^{Back}[MeV]", 5000, 0, SiMax, 5000, 0, SiMax);
  TH2F *SiGAGGPID = new TH2F("SiGAGGPID", "SiGAGGPID;GAGG sum;E_{1}^{Back} + E_{2}^{Back}", 5000, 0, GAGGMax, 5000, 0, 2*SiMax);
  TH2F *SiGAGGPID_prompt= new TH2F("SiGAGGPID_prompt", "SiGAGGPID_prompt;GAGG sum;E_{1}^{Back} + E_{2}^{Back}", 5000, 0, GAGGMax, 5000, 0, 2*SiMax);
  TH2F *SiGAGGPID_prompt_1906 = new TH2F("SiGAGGPID_prompt_1906", "SiGAGGPID_prompt;GAGG sum;E_{1}^{Back} + E_{2}^{Back}", 5000, 0, GAGGMax, 5000, 0, 2*SiMax);
  TH2F *SiGAGGPID_prompt_1905 = new TH2F("SiGAGGPID_prompt_1905", "SiGAGGPID_prompt;GAGG sum;E_{1}^{Back} + E_{2}^{Back}", 5000, 0, GAGGMax, 5000, 0, 2*SiMax);
  TH2F *SiPID_gated = new TH2F("SiPID_gated", "SiPID_gated;E_{2}^{Back}[MeV];E_{1}^{Back}[MeV]", 5000, 0, SiMax, 5000, 0, SiMax);
  TH2F *SiE_gated = new TH2F("SiE_gated", "Si Energy (gated); ID; E",512,0,512,5000,0,SiMax);

  TH2F *Si1PIDFrontvBack = new TH2F("Si1PIDFrontvBack", "Si1PIDFrontvBack;E_{1}^{Front}[MeV];E_{1}^{Back}[MeV]", 5000, 0, SiMax, 5000, 0, SiMax);
  TH2F *Si2PIDFrontvBack = new TH2F("Si2PIDFrontvBack", "Si2PIDFrontvBack;E_{2}^{Front}[MeV];E_{2}^{Back}[MeV]", 5000, 0, SiMax, 5000, 0, SiMax);

  TH1F *SiMultiplicityFront1 = new TH1F("SiMultiplicityFront1", "SiMultiplicityFront1", 10, 0, 10);
  TH1F *SiMultiplicityFront2 = new TH1F("SiMultiplicityFront2", "SiMultiplicityFront2", 10, 0, 10);

  TH1F *SiMultiplicityBack1 = new TH1F("SiMultiplicityBack1", "SiMultiplicityBack1", 10, 0, 10);
  TH1F *SiMultiplicityBack2 = new TH1F("SiMultiplicityBack2", "SiMultiplicityBack2", 10, 0, 10);

  TH1F *GAGGMultiplicity1 = new TH1F("GAGGMultiplicity1", "GAGGMultiplicity1", 25, 0, 25);
  TH1F *GAGGMultiplicity2 = new TH1F("GAGGMultiplicity2", "GAGGMultiplicity2", 16, 0, 16);
  TH1F *GAGGADCMax = new TH1F("GAGGADCMax", "GAGGADCMax", 5000, 0, 10000);
  TH2F *GAGGEvsID = new TH2F("GAGGEvsID", "GAGGEvsID; E; ID", 8192, 0, 8*8192, 64, 0, 64);
  TH2F *GAGGEvsID_gated = new TH2F("GAGGEvsID_gated", "GAGGEvsID (gated on beam cut); E; ID", 8192, 0, 8*8192, 64, 0, 64);
  TH2F *GAGGMaxvsID = new TH2F("GAGGMaxvsID", "GAGGMaxvsID; E; ID", 8192, 0, 8192, 64, 0, 64);
  TH2F *GAGGHitPattern = new TH2F("GAGGHitPattern", "GAGGHitPattern; ID; ID", 64, 0, 64, 64, 0, 64);
  TH2F *GAGGIDvsMult = new TH2F("GAGGIDvsMult", "GAGGIDvsMult; ID, Multiplicity", 64, 0, 64, 10, 0, 10);

  TH2F *GAGGPosition = new TH2F("GAGGPosition", "GAGGPosition; X; Y", 5,0,5,5,0,5);

  TH2F *Si1Position = new TH2F("Si1Position", "Si1Position; X (strip); Y (strip)", 128, 0, 128, 128, 0, 128);
  TH2F *Si2Position = new TH2F("Si2Position", "Si2Position; X (strip); Y (strip)", 128, 0, 128, 128, 0, 128);
  TH2F *Si1PosGated = new TH2F("Si1PosGated", "Si1Position; X (strip); Y (strip)", 128, 0, 128, 128, 0, 128);

  TH1F *Si1FBTDiff = new TH1F("Si1FBTDiff", "Si1 Front vs Back time difference; Time difference; Counts", 1024, -512, 512);
  TH1F *Si2FBTDiff = new TH1F("Si2FBTDiff", "Si2 Front vs Back time difference; Time difference; Counts", 1024, -512, 512);

  TH1F *SiSiTDiff = new TH1F("SiSiTDiff", "Si1 vs Si2 time difference; Time difference; Counts", 1024, -512, 512);

  TH1F *SiGAGGTDiff = new TH1F("SiGAGGTDiff", "Si1 vs GAGG time difference; Time difference; Counts", 1024, -512, 512);
  TH2F *SiGAGGIDTDiff = new TH2F("SiGAGGIDTDiff", "Si1 vs GAGG time difference; Time difference; GAGG ID", 1024, -512, 512, 24, 0, 24);
  TH2F *SiGAGGIDTDiff_gated = new TH2F("SiGAGGIDTDiff_gated", "Si1 vs GAGG time difference (position gated); Time difference; GAGG ID", 1024, -512, 512, 24, 0, 24);

  TH2F *Si1PosGAGG[40];
  TH2F *Si2PosGAGG[40];
  for (int i=0; i<40; ++i) {
  Si1PosGAGG[i] = new TH2F(Form("Si1PosGAGG%i", i), "Si1Position; X (strip); Y (strip)", 128, 0, 128, 128, 0, 128);
  Si2PosGAGG[i] = new TH2F(Form("Si2PosGAGG%i", i), "Si2Position; X (strip); Y (strip)", 128, 0, 128, 128, 0, 128);
  }
  // Open the digitalization file and get the TTree.
  TString unpackFileName = TString::Format("/home/tgray/e535/rawdata/UnpackerTestOutput/run_%04d.root", runNo);
  TFile *unpackFile = new TFile(unpackFileName, "READ");
  TTree *unpackTree = (TTree *)unpackFile->Get("cbmsim");
  int nUnpackEvents = unpackTree->GetEntries();
  std::cout << " Number of unpacked events in run " << runNo << ": " << nUnpackEvents << std::endl;
  // Creare the TTreeReader to read the AtTrackingEvents and simulation.
  TTreeReader unpackReader("cbmsim", unpackFile);
  TTreeReaderValue<TClonesArray> siArray(unpackReader, "AtSiEvent");
  TTreeReaderValue<TClonesArray> gaggArray(unpackReader, "AtGaggEvent");
  TTreeReaderValue<TClonesArray> patternArray(unpackReader, "AtPatternEvent");

  // Loop over events.
  for (int i = 0; i < nUnpackEvents; i++) {
    unpackReader.Next();

    // Check the Si data first.
    AtSiEvent *siEvent = (AtSiEvent *)siArray->At(0);

    Int_t multiplicityFront1 = siEvent->GetMultiplicityFront1();
    Int_t multiplicityFront2 = siEvent->GetMultiplicityFront2();

    Int_t multiplicityBack1 = siEvent->GetMultiplicityBack1();
    Int_t multiplicityBack2 = siEvent->GetMultiplicityBack2();

    SiMultiplicityFront1->Fill(multiplicityFront1);
    SiMultiplicityFront2->Fill(multiplicityFront2);

    SiMultiplicityBack1->Fill(multiplicityBack1);
    SiMultiplicityBack2->Fill(multiplicityBack2);

    Int_t strip1X = siEvent->GetXStrip1();
    Int_t strip1Y = siEvent->GetYStrip1();

    ind = strip1X; 
    Float_t front1E = si_cal[ind]*siEvent->GetEnergyFront1();
    Float_t front1raw = siEvent->GetEnergyFront1();

    ind = strip1Y+128; 
    Float_t back1E = si_cal[ind]*siEvent->GetEnergyBack1();
    Float_t back1raw = siEvent->GetEnergyBack1();
    Int_t time1 = siEvent->GetTimeFront1();

    Int_t strip2X = siEvent->GetXStrip2();
    Int_t strip2Y = siEvent->GetYStrip2();

    ind = strip2X+256; 
    Float_t front2E = si_cal[ind]*siEvent->GetEnergyFront2();
    Float_t front2raw = siEvent->GetEnergyFront2();
    ind = strip2Y+256+128; 
    Float_t back2E = si_cal[ind]*siEvent->GetEnergyBack2();
    Float_t back2raw = siEvent->GetEnergyBack2();

    Int_t time2 = siEvent->GetTimeFront2();

    bool goodSiPID = false;
    bool beamSiPID = false;
    double SiSum = 0;

    if (front1E > 20 && back1E > 20) {
      Si1Position->Fill(strip1X, strip1Y);
    }

    if (front2E > 20 && back2E > 20) {
      Si2Position->Fill(strip2X, strip2Y);
    }

    if ((front1E > 1) && (back1E > 1) && (front2E > 1) && (back2E > 1)) {
      int tdiff = time1-time2;
      Si1PIDFrontvBack->Fill(front1E, back1E);
      Si2PIDFrontvBack->Fill(front2E, back2E);
      if (abs(tdiff) < 10) {
        SiSum = back1E + back2E;
        SiPIDTraceIntegralFronts->Fill(front2E, front1E);
        SiPIDTraceIntegralBacks->Fill(back2E, back1E);
        goodSiPID = true;
          beamSiPID = true;
          SiE_gated->Fill(strip1X,front1E);
          SiE_gated->Fill(strip1Y+128,back1E);
          if (beamcut->IsInside(back2raw, back1raw)) {
          SiE_gated->Fill(strip2X+2*128,front2E);
          SiE_gated->Fill(strip2Y+3*128,back2E);
        }
      }
    }

    for (int i = 0; i<multiplicityFront1; ++i) {
      for (int j = 0; j<multiplicityBack1; ++j) {
        Si1FBTDiff->Fill(siEvent->GetTSFront1(i) - siEvent->GetTSBack1(j));
      }
    }
    for (int i = 0; i<multiplicityFront2; ++i) {
      for (int j = 0; j<multiplicityBack2; ++j) {
        Si2FBTDiff->Fill(siEvent->GetTSFront2(i) - siEvent->GetTSBack2(j));
      }
    }

    SiSiTDiff->Fill(time1-time2);

    if ((multiplicityFront1 != 1 && multiplicityBack1 != 1) || (multiplicityFront2 != 1 && multiplicityBack2 != 1)) continue;

    // Check the GAGG data next.
    AtGaggEvent *gaggEvent = (AtGaggEvent *)gaggArray->At(0);

    Int_t multiGagg1 = gaggEvent->GetMultiplicity1();
    Int_t multiGagg2 = gaggEvent->GetMultiplicity2();

    GAGGMultiplicity1->Fill(multiGagg1);
    GAGGMultiplicity2->Fill(multiGagg2);

    int g_mult = 0;
    double GAGGSum = 0;
    int GAGGID = -1;
    int GAGGTime = -1;
    int GAGGX, GAGGY;
    for (int i=0; i<multiGagg1; ++i) {
      int ID = (int)gaggEvent->GetID1(i);
      if (ID == 36) { ID = 24; } //crystal 24 is missing but 36 (second layer) is behind
      if (ID < 0 || ID > 24) { continue; }
      GAGGID = ID;
      GAGGX = 4-GAGGID/5;
      GAGGY = 4-GAGGID%5;
      GAGGPosition->Fill(GAGGX, GAGGY);
      GAGGADCMax->Fill(gaggEvent->GetADCMax1(i));
      GAGGMaxvsID->Fill(gaggEvent->GetADCMax1(i), GAGGID);
      if (beamSiPID) {
        GAGGEvsID_gated->Fill(gaggEvent->GetE1(i), GAGGID);
      }
      GAGGIDvsMult->Fill(gaggEvent->GetID1(i), multiGagg1);
      g_mult++; 
      double energy = gagg_cal[GAGGID]*gaggEvent->GetE1(i);
      //double energy = gaggEvent->GetE1(i);
      GAGGTime = gaggEvent->GetTS1(i);
      GAGGEvsID->Fill(energy, GAGGID);
      GAGGSum += energy;
      for (int j=i+1; j<multiGagg1; ++j) {
        GAGGHitPattern->Fill(gaggEvent->GetID1(i), gaggEvent->GetID1(j));
        GAGGHitPattern->Fill(gaggEvent->GetID1(j), gaggEvent->GetID1(i));
      }
    }

    if ((g_mult == 0) && goodSiPID) {
      SiPIDNoGAGG->Fill(back2E, back1E);
      if (back2E > 1200 && back2E < 1800) {
        if (back1E > 550 && back1E < 700) {
          Si1PosGated->Fill(strip1X, strip1Y);
        }
      }
    }
    else if ((g_mult==1) && goodSiPID) {
      int tdiff = time1-GAGGTime;
      SiGAGGTDiff->Fill(tdiff);
      SiGAGGIDTDiff->Fill(tdiff, GAGGID);
      SiGAGGPID->Fill(GAGGSum, SiSum);
      //SiGAGGPID->Fill(GAGGSum, SiSum);
      bool pos_gate = false;
      if ((col_boundaries[GAGGX]-2) < strip1X && strip1X < (col_boundaries[GAGGX+1]+2)) {
        if ((row_boundaries[GAGGY]-2) < strip1Y && strip1Y < (row_boundaries[GAGGY+1]+2)) {
          pos_gate = true;
          SiGAGGIDTDiff_gated->Fill(tdiff, GAGGID);
        }
      }
      if (tdiff >= 25 && tdiff <=50) {
        Si1PosGAGG[GAGGID]->Fill(strip1X, strip1Y);
        Si2PosGAGG[GAGGID]->Fill(strip2X, strip2Y);
        if (pos_gate) {
          SiGAGGPID_prompt->Fill(GAGGSum, SiSum);
          if (GAGGID < 16 && GAGGID != 6) {
            SiGAGGPID_prompt_1906->Fill(GAGGSum, SiSum);
          }
          else if (GAGGID >15) {
            SiGAGGPID_prompt_1905->Fill(GAGGSum, SiSum);
          }

          if (GAGGSum>2450 && GAGGSum < 2550) {
            SiPID_gated->Fill(back2E, back1E);
          }
        }
      }
    }
  }
  outfile->Write();
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
