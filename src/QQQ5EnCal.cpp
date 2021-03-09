/*************************************************************************
 * QQQ5EnCal.cpp script can be used to alpha calibrate energies for      *
 * different strips of 8 QQQ5 detectors.		                         *
 *                                                                       *
 * Energy calibration (for 4 dets*32 rings ) is done by gaussian	     *
 * fit of each histograms. The channel numbers are then saved in         *
 * QQQ5EnCalChannels.dat file.                				             *
 * 									                                     *
 * To see if the gaussian fit is correct, EnCal.root file is created in  *
 * the Output directory.		   				                         *
 *************************************************************************/

#define Analysis_cxx
#include "EnergyLoss.h"
#include "EffectiveThickness.h"
#include "Analysis.h"


#include <string>
#include <fstream>

TChain* MakeChain();

int main() {
    TChain *chain = MakeChain();

    Analysis t(chain);
    t.Loop();
}

TChain* MakeChain() {
    auto *chain = new TChain("data");
    TString PathToFiles = "/mnt/e/goddessSort-master/Output/Run";

    //Choose your detector (Find QQQ5 file) (All Upstreamstream)
    chain->Add(PathToFiles + "0455.root"); // All upstream QQQ5
    



    return chain;
}

void Analysis::Loop() {
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    int prevRunNumber = -1;


	//Define Histograms here
	//Histograms without the gains applied
	TH1F *QQQ5_EnCal[4][32]; // Energy Histogram for the Calibration of the Energy
	for (Int_t i=0; i<4; i++){ // Loop over detectors
		for (Int_t j=0; j<32; j++){ //Loop over the strips
			std::string nameQQQ5_EnCal = Form("QQQ5_%i_Strip_%i_EnCal",i,j);
			QQQ5_EnCal[i][j] = new TH1F(nameQQQ5_EnCal.c_str(), "QQQ5 detector Energy Spectrum", 1000,1000,4000);
		}//End of Loop over Strips
	}//End of Loop over Detectors

	EnergyLoss* Source = new EnergyLoss("AlphaInAu.dat");// SRIM File for Source casing energy loss
	EnergyLoss* DeadLayer = new EnergyLoss("AlphaInSi.dat");// SRIM File for DeadLayer energy loss


	//Create Output File   
	TFile* outputFile = new TFile("/mnt/e/Analysis/QQQ5Calibration/Output/EnCal.root", "recreate");

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

		
	   	
		//Loop over the Multiplicity
		for(Int_t j=0; j<QQQ5Mul; j++){
						
			//Pedestals Substraction
			Float_t RawEnergy = QQQ5RingADC[j]; //- Pedestals[(BB10Det[j]*8)+(BB10Strip[j])];  
		

			//Filling the histograms here
			QQQ5_EnCal[QQQ5Det[j]][QQQ5Ring[j]]->Fill(RawEnergy);
			
	
		}// End of the multiplicity Loop	
    }// End of event by event analysis

    //Cleaning the dat file to save the gains
    std::ofstream pfile;
    pfile.open("/mnt/e/Analysis/QQQ5Calibration/DatFile/QQQ5EnCalChannels.dat", std::ofstream::out | std::ofstream::trunc);
    pfile.close();

	
    
    //Define your fit function here
    TF1* fit_gaus = new TF1 ("fit_gaus","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp(-pow((x-[8]),2)/(2*[2]*[2]))+[9]*TMath::Exp(-pow((x-[10]),2)/(2*[2]*[2]))",1000,4000);
    
    Int_t npeaks = 5;
    TSpectrum *s = new TSpectrum(npeaks);
    Float_t xp[5] = {0};
	
    outputFile->cd();

    for (Int_t i=0; i<4; i++){
   		for (Int_t j=0; j<32; j++){

			Int_t nfound = s->Search(QQQ5_EnCal[i][j],2.5,"",0.3);
			Double_t *xpeaks = s->GetPositionX();
			for (Int_t p=0; p<nfound; p++){
				xp[p] = xpeaks[p];
			}
			std::sort(xp, xp+5);
			Float_t xmin = std::min(std::min(std::min(std::min(xp[0],xp[1]),xp[2]),xp[3]),xp[4]);
			Float_t xmax = std::max(std::max(std::max(std::max(xp[0],xp[1]),xp[2]),xp[3]),xp[4]);
			fit_gaus->SetParLimits (1,xp[0]-5,xp[0]+5);
			fit_gaus->SetParLimits(4,xp[1]-5,xp[1]+5);
			fit_gaus->SetParLimits(6,xp[2]-5,xp[2]+5);
			fit_gaus->SetParLimits(8,xp[3]-5,xp[3]+5);
			fit_gaus->SetParLimits(10,xp[4]-5,xp[4]+5);
			fit_gaus->SetParLimits(2,1,50);	
			
			QQQ5_EnCal[i][j]->Fit("fit_gaus","Q","",xmin-50,xmax+50);
			//std::cout << fit_gaus->GetParameter (1) << std::endl;
			std::ofstream pfile;
			pfile.open("/mnt/e/Analysis/QQQ5Calibration/DatFile/QQQ5EnCalChannels.dat", std::ofstream::app);
			pfile << i << '\t' << j << '\t' << fit_gaus->GetParameter (1) << '\t' << fit_gaus->GetParameter (4) << '\t' << fit_gaus->GetParameter (6) << '\t' << fit_gaus->GetParameter (8) << '\t' << fit_gaus->GetParameter (10) << std::endl;
			pfile.close();
			QQQ5_EnCal[i][j]->Write();   
		}
    }
    outputFile->Close();


	//Now check the "EnCal.root" file. Is the histogram fitted correctly?
	//If yes, the energy loss for 5 alpha peaks are calculated and "QQQ5EnCalChannels.dat" file is opened. Then Slopes and Intercepts are calculated


	//Open QQQ5 Geometry File for Ring Angle for Effective thickness
	std::ifstream geofile;
	geofile.open("DatFile/QQQ5Geometry.dat");
	Double_t QQQ5RingDistance[32] = {0};
	Double_t QQQ5RingAngle[32]={0};
	for (Int_t i=0; i<32; i++){
		geofile >> QQQ5RingDistance[i] >> QQQ5RingAngle[i];
		//std::cout << QQQ5RingDistance[i] << std::endl;
	}
	geofile.close();
	

	//Open the QQQ5EnCalChannels.dat file
	std::ifstream file;
	file.open("DatFile/QQQ5EnCalChannels.dat");
	Double_t Det[128] = {0};
	Double_t Ring[128] = {0};
	Double_t FirstChannel[128] = {0};
	Double_t SecondChannel[128] = {0};
	Double_t ThirdChannel[128] = {0};
	Double_t FourthChannel[128] = {0};
	Double_t FifthChannel[128] = {0};
	for (Int_t i=0; i<128; i++){
		file >> Det[i] >> Ring[i] >> FirstChannel[i] >> SecondChannel[i] >> ThirdChannel[i] >> FourthChannel[i] >> FifthChannel[i];
		std::cout << Det[i] << '\t' << Ring[i] << '\t' << FirstChannel[i] << '\t' << SecondChannel[i] << '\t' << ThirdChannel[i] << '\t' << FourthChannel[i] << '\t' << FifthChannel[i] << std::endl;
	}
	file.close();

	//Cleaning the dat file to save the gains
    std::ofstream cfile;
    cfile.open("/mnt/e/Analysis/QQQ5Calibration/Output/QQQ5EnCal.dat", std::ofstream::out | std::ofstream::trunc);
    cfile.close();

	//Finding the source casing thickness and deadlayer thickness
	for (Int_t i=0; i<128; i++){
		Int_t RingNum = Ring[i];
		Float_t DeadLayerThickness;
    	Float_t Q5_Phi = 45 + (Det[i]*90.);
    	targ_thick(QQQ5RingAngle[RingNum]*M_PI/180., Q5_Phi*M_PI/180., 0.0001, 0.0001, 0*M_PI/180., &DeadLayerThickness);//Dead layer thickness = 0.0001mm and Depth = thickness
    	//std::cout << DeadLayerThickness << std::endl;

    	Float_t SourceCasingThickness;
    	targ_thick(QQQ5RingAngle[RingNum]*M_PI/180., Q5_Phi*M_PI/180., 0.022940293, 0.022940293, 27*M_PI/180., &SourceCasingThickness); //Effective Source casing thickness = 0.022940293 um (See Analysis Meeting Dead Layer Geometry.ppt)
    	//std::cout << TargetThickness << std::endl;

		//First, Calculate the remainder energy after the alpha passes through source casing 
		double Alpha1 = Source->CalcRemainder(5.42308,(SourceCasingThickness/1000.)); 
		double Alpha2 = Source->CalcRemainder(5.68537,(SourceCasingThickness/1000.));
		double Alpha3 = Source->CalcRemainder(6.288,(SourceCasingThickness/1000.));
		double Alpha4 = Source->CalcRemainder(6.7783,(SourceCasingThickness/1000.));
		double Alpha5 = Source->CalcRemainder(8.78412,(SourceCasingThickness/1000.));

		//Then, Calculate the remainder energy after the alpha passes through the Deadlayer
		double EffAlpha1 = DeadLayer->CalcRemainder(Alpha1,(DeadLayerThickness/1000.)); 
		double EffAlpha2 = DeadLayer->CalcRemainder(Alpha2,(DeadLayerThickness/1000.));
		double EffAlpha3 = DeadLayer->CalcRemainder(Alpha3,(DeadLayerThickness/1000.));
		double EffAlpha4 = DeadLayer->CalcRemainder(Alpha4,(DeadLayerThickness/1000.));
		double EffAlpha5 = DeadLayer->CalcRemainder(Alpha5,(DeadLayerThickness/1000.));

		std::ofstream cfile;
		cfile.open("/mnt/e/Analysis/QQQ5Calibration/Output/QQQ5EnCal.dat", std::ofstream::app);
		cfile << EffAlpha1 << '\t' << EffAlpha2 << '\t' << EffAlpha3 << '\t' << EffAlpha4 << '\t' << EffAlpha5 << '\t' << FirstChannel[i] << '\t' << SecondChannel[i] << '\t' << ThirdChannel[i] << '\t' << FourthChannel[i] << '\t' << FifthChannel[i] << std::endl;
		cfile.close();
	}
}














