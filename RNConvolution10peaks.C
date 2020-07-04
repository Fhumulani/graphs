#include <TROOT.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TCanvas.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

TF1 *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8, *f9, *f10;
 vector<TH1F*> h;

double FitFunction(double *x, double *pars);

ofstream output;


void RNConvolution10peaks()
{



char* Name1 = NULL;
ofstream results_text;
Char_t infoFileName[200];



TFile* cut_run65 = new TFile("CutRun65.root");//cut
  TCutG* incut65;
  incut65 = (TCutG*)gROOT->FindObjectAny("cutRun65");
   TIter next(cut_run65->GetListOfKeys());
   TKey *key;
   while((key=(TKey*)next())){
        if(strcmp(key->GetClassName(),"TCutG") == 0){
            Name1 = (char*)key->GetName();
            cout << "cut found "<<Name1 << endl;
            incut65 = (TCutG*)cut_run65->Get(Name1); 
               
        }
    }
    cut_run65->Close();



 vector<TChain*> Chains;
 vector<double> Angles;
 vector<double> Offset;
 vector<int> K600_ang;


 
 Angles.push_back(-1.8);
 Angles.push_back(-0.9);
 Angles.push_back(0.0);
 Angles.push_back(0.9); 
 Angles.push_back(1.8);
 



 K600_ang.push_back(4);
 Offset.push_back(-0.15);                       // e f g are EX energy calibration parameters 
 Chains.push_back(new TChain("DATA"));

 K600_ang.push_back(7);
 Offset.push_back(-0.25 );
 Chains.push_back(new TChain("DATA"));

 K600_ang.push_back(10); 
 Offset.push_back(-0.2);
 Chains.push_back(new TChain("DATA"));

 K600_ang.push_back(15);
 Offset.push_back(-0.1);  
 Chains.push_back(new TChain("DATA"));

 K600_ang.push_back(21);
 Offset.push_back(-0.1);
 Chains.push_back(new TChain("DATA"));

 K600_ang.push_back(28);
 Offset.push_back(-0.15);
 Chains.push_back(new TChain("DATA"));



   // Set the Chains;  
   Chains[0]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00250.root");  //4 deg
   Chains[0]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00252.root");
   Chains[0]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00251.root"); 
 
  
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00065.root"); //7 deg
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00067.root");
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00068.root");
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00069.root");
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00070.root");
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00071.root");
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00072.root");
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00073.root");
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00074.root");
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00075.root");
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00076.root");
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00077.root");
   Chains[1]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00078.root");
  
   Chains[2]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00079.root");  //10deg
   Chains[2]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00080.root");

   
   Chains[3]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00105.root");//15deg
   Chains[3]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00106.root");
   Chains[3]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00107.root");
   Chains[3]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00108.root");
   Chains[3]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00110.root");
    
   Chains[4]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00120.root");  //21 deg
   Chains[4]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00121.root");

   Chains[5]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00140.root");  //28 deg
   Chains[5]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00141.root");
   Chains[5]->Add("/media/fhumulani/Transcend/tiny_phd/sorted/sorted00142.root");
    

//------------------------------------------------histos from chains-----------------------------------------------//

 
  for(int s=0; s< K600_ang.size();s++){//defining histos 
    for(int c=0; c< Angles.size()-1; c++){
    h.push_back(new TH1F(Form("h%d%d",K600_ang[s],c), Form("h%d%d",K600_ang[s],c), 2400, -1., 11.));
   }
 } 
 
cout <<"h size: " << h.size() <<endl;

 
 
 for(UShort_t t=0; t< K600_ang.size(); t++){//for each k600 angle...
    cout << "K600 angle" << K600_ang[t]<<endl;
    int count = 0;
   
    for(UShort_t i=0; i<Angles.size()-1; i++){//there is angle bite
    
    Chains[t]->Draw(Form("Ex >>h%d%d",K600_ang[t],count),Form("cutRun65 &&!X1flag && !U1flag && !X2flag &&!pulser==1 &&X1pos>0 &&(thetaSCAT+%f)>%f && (thetaSCAT+%f)<%f",Offset[t],Angles[i],Offset[t],Angles[i+1]),"");
    
  count = count + 1;
   }
   
   
 }


//---------------------- fitting doublets wth funtion to get area--------------------------------------------------------------//


 float cxn[100];
 float error[100];
 float CII;
 float tAD = 1.6e19; 
 float eff;
 float qe = 1.602e-19;
 float Q_value = -7.6645 ; //in MeV
 float If[100];
   TGraphErrors *gr[100];
 double Angles_plot[] = {2.5,3.5,4.5,5.5,5.5,6.5,7.5,8.5,8.5,9.5,10.5,11.5,13.5,14.5,15.5,16.5,19.5,20.5,21.5,22.5,26.5,27.5,28.5,29.5};
 double theta_cm[100];
 double gamma = 0.04;
   float SA[] = {
0.7827,1.05678,1.05678,0.7827,
0.7827,1.05678,1.05678,0.7827,
0.7827,1.05678,1.05678,0.7827,
0.7827,1.05678,1.05678,0.7827,
0.7827,1.05678,1.05678,0.7827,
0.7827,1.05678,1.05678,0.7827
};  		// msr   +- 2 degrees
float Erro_x[100] ={0.00000};


   

   /* 
    double peak1=1.57; 
    double peak2=1.95;
    double peak3=2.35; 
    double peak4=2.71; 
    double peak5=2.92; 
    double peak6=3.4; 
    double peak7=3.6; 
    double peak8=3.8;
    double peak9=3.9;
    double peak10=4.0;
    double fitmin=1.;
    double fitmax=4.5;
    double displaymin=1.;
    double displaymax=4.5;
*/
    
    double peak1=4.48; 
    double peak2=4.58;
    double peak3=4.62; 
    double peak4=4.95; 
    double peak5=5.64; 
    double peak6=5.77; 
    double peak7=5.85; 
    double peak8=6.;
    double peak9=6.04;
    double peak10=6.13;
    double fitmin=4.;
    double fitmax=6.5;
    double displaymin=4.;
    double displaymax=6.5;

    double gauswidth=0.015;
    double landwidth=0.005;

    double peakarea[100];
    double peak[100];

    TCanvas *c5 ;//= new TCanvas("c5","c5",800,900);
    //c5->Divide(2,5);
  


    gPad->SetLogy(1);

    int nb_c = (int)h.size()/4;
    if (h.size()%4!=0) {
        nb_c++;
    }
    TCanvas *c2[nb_c];
    int can_id = -1;

   


   
    TF1 *fGaus = new TF1("fGaus","[0] * TMath::Gaus(x,[1],[2],0)",displaymin,displaymax);
    fGaus->SetParameter(0,200);
    fGaus->SetParameter(1,peak1);
    fGaus->SetParameter(2,gauswidth);
    fGaus->SetNpx(2e4);
    fGaus->SetLineColor(6);
    fGaus->Draw();
    
    TF1 *fLandau = new TF1("fLandau","[0] * TMath::Landau(x,[1],[2])",displaymin,displaymax);
    fLandau->SetParameter(0,500);
    fLandau->SetParameter(1,peak1);
    fLandau->SetParameter(2,landwidth);
    fLandau->SetLineColor(6);
    fLandau->Draw("same");
    
    TF1Convolution *fConv = new TF1Convolution(fGaus,fLandau,true);
    fConv->SetRange(displaymin,displaymax);
    fConv->SetNofPointsFFT(1000);
    fConv->SetParameters(20,peak1,gauswidth,50,peak1,landwidth);
    
    f1 = new TF1("fConvFn1",*fConv,displaymin,displaymax,fConv->GetNpar());
    f2 = new TF1("fConvFn2",*fConv,displaymin,displaymax,fConv->GetNpar());
    f3 = new TF1("fConvFn3",*fConv,displaymin,displaymax,fConv->GetNpar());
    f4 = new TF1("fConvFn4",*fConv,displaymin,displaymax,fConv->GetNpar());
    f5 = new TF1("fConvFn5",*fConv,displaymin,displaymax,fConv->GetNpar());
    f6 = new TF1("fConvFn6",*fConv,displaymin,displaymax,fConv->GetNpar());
    f7 = new TF1("fConvFn7",*fConv,displaymin,displaymax,fConv->GetNpar());
    f8 = new TF1("fConvFn8",*fConv,displaymin,displaymax,fConv->GetNpar());
    f9 = new TF1("fConvFn9",*fConv,displaymin,displaymax,fConv->GetNpar());
    f10 = new TF1("fConvFn10",*fConv,displaymin,displaymax,fConv->GetNpar());

    //Parameters for convolved functions
    //0: Gaussian size
    //1: Gaussian position
    //2: Gaussian width
    //3: Landau size
    //4: Landau position
    //5: Landau width
    
   
    for(int i=0;i<fGaus->GetNpar();i++)f1->SetParameter(i,fGaus->GetParameter(i));
    for(int i=0;i<fLandau->GetNpar();i++)f1->SetParameter(i+fGaus->GetNpar(),fLandau->GetParameter(i));
    f1->SetParameter(0,5);//to make the convolution plot nicely in scale with the Gaussian

    for(int i=0;i<fGaus->GetNpar();i++)f2->SetParameter(i,fGaus->GetParameter(i));
    for(int i=0;i<fLandau->GetNpar();i++)f2->SetParameter(i+fGaus->GetNpar(),fLandau->GetParameter(i));
    f2->SetParameter(1,peak2);
    f2->SetParameter(4,peak2);

    for(int i=0;i<fGaus->GetNpar();i++)f3->SetParameter(i,fGaus->GetParameter(i));
    for(int i=0;i<fLandau->GetNpar();i++)f3->SetParameter(i+fGaus->GetNpar(),fLandau->GetParameter(i));
    f3->SetParameter(1,peak3);
    f3->SetParameter(4,peak3);

    for(int i=0;i<fGaus->GetNpar();i++)f4->SetParameter(i,fGaus->GetParameter(i));
    for(int i=0;i<fLandau->GetNpar();i++)f4->SetParameter(i+fGaus->GetNpar(),fLandau->GetParameter(i));
    f4->SetParameter(1,peak4);
    f4->SetParameter(4,peak4);

    for(int i=0;i<fGaus->GetNpar();i++)f5->SetParameter(i,fGaus->GetParameter(i));
    for(int i=0;i<fLandau->GetNpar();i++)f5->SetParameter(i+fGaus->GetNpar(),fLandau->GetParameter(i));
    f5->SetParameter(1,peak5);
    f5->SetParameter(4,peak5);

    for(int i=0;i<fGaus->GetNpar();i++)f6->SetParameter(i,fGaus->GetParameter(i));
    for(int i=0;i<fLandau->GetNpar();i++)f6->SetParameter(i+fGaus->GetNpar(),fLandau->GetParameter(i));
    f6->SetParameter(1,peak6);
    f6->SetParameter(4,peak6);

    for(int i=0;i<fGaus->GetNpar();i++)f7->SetParameter(i,fGaus->GetParameter(i));
    for(int i=0;i<fLandau->GetNpar();i++)f7->SetParameter(i+fGaus->GetNpar(),fLandau->GetParameter(i));
    f7->SetParameter(1,peak7);
    f7->SetParameter(4,peak7);

    for(int i=0;i<fGaus->GetNpar();i++)f8->SetParameter(i,fGaus->GetParameter(i));
    for(int i=0;i<fLandau->GetNpar();i++)f8->SetParameter(i+fGaus->GetNpar(),fLandau->GetParameter(i));
    f8->SetParameter(1,peak8);
    f8->SetParameter(4,peak8);

    for(int i=0;i<fGaus->GetNpar();i++)f9->SetParameter(i,fGaus->GetParameter(i));
    for(int i=0;i<fLandau->GetNpar();i++)f9->SetParameter(i+fGaus->GetNpar(),fLandau->GetParameter(i));
    f9->SetParameter(1,peak9);
    f9->SetParameter(4,peak9);

    for(int i=0;i<fGaus->GetNpar();i++)f10->SetParameter(i,fGaus->GetParameter(i));
    for(int i=0;i<fLandau->GetNpar();i++)f10->SetParameter(i+fGaus->GetNpar(),fLandau->GetParameter(i));
    f10->SetParameter(1,peak10);
    f10->SetParameter(4,peak10);

    f1->SetNpx(1000);
    f2->SetNpx(1000);
    f3->SetNpx(1000);
    f4->SetNpx(1000);
    f5->SetNpx(1000);
    f6->SetNpx(1000);
    f7->SetNpx(1000);
    f8->SetNpx(1000);
    f9->SetNpx(1000);
    f10->SetNpx(1000);

    f1->Draw("same");
    f2->Draw("same");
    f3->Draw("same");
    f4->Draw("same");
    f5->Draw("same");
    f6->Draw("same");
    f7->Draw("same");
    f8->Draw("same");
    f9->Draw("same");
    f10->Draw("same");
   

    for(UShort_t d=0; d<24 ;d++){ 

 
    c5 = new TCanvas();
    gPad->SetLogy(1);
    TF1 *fitF = new TF1("fitF",FitFunction,fitmin,fitmax,22);
    fitF->SetParameter(0,70000);//amplitude first Gaussian
    fitF->SetParameter(1,20000);//amplitude 
    fitF->SetParameter(2,20000);//amplitude
    fitF->SetParameter(3,20000);//amplitude
    fitF->SetParameter(4,20000);//amplitude
    fitF->SetParameter(5,20000);//amplitude
    fitF->SetParameter(6,20000);//amplitude
    fitF->SetParameter(7,20000);//amplitude
    fitF->SetParameter(8,20000);//amplitude
    fitF->SetParameter(9,20000);//amplitude
    fitF->SetParameter(10,peak1);//position first Gaussian
    fitF->SetParameter(11,peak2);//position 
    fitF->SetParameter(12,peak3);//position 
    fitF->SetParameter(13,peak4);//position 
    fitF->SetParameter(14,peak5);//position 
    fitF->SetParameter(15,peak6);//position 
    fitF->SetParameter(16,peak7);//position 
    fitF->SetParameter(17,peak8);//position 
    fitF->SetParameter(18,peak9);//position 
    fitF->SetParameter(19,peak10);//position 
    fitF->SetParameter(20,gauswidth);//Gaussian width
    fitF->SetParameter(21,landwidth);//Landau spreading
    //fitF->FixParameter(9,landwidth);//Landau spreading
  
    fitF->SetNpx(10000);
    fitF->SetLineColor(1);
    fitF->Draw("same");


     h[d]->Fit("fitF","RQ","",4.3,6.5);
   
     f1->Draw("same");
     f2->Draw("same");
     f3->Draw("same");
     f4->Draw("same");
     f5->Draw("same");
     f6->Draw("same");
     f7->Draw("same");
     f8->Draw("same");
     f9->Draw("same");
     f10->Draw("same");


  
     peakarea[0] = 200*f1->Integral(f1->GetParameter(1)-0.05,f1->GetParameter(1)+0.15);
     peakarea[1] = 200*f2->Integral(f2->GetParameter(1)-0.05,f2->GetParameter(1)+0.15);
     peakarea[2] = 200*f3->Integral(f3->GetParameter(1)-0.05,f3->GetParameter(1)+0.15);
     peakarea[3] = 200*f4->Integral(f4->GetParameter(1)-0.05,f4->GetParameter(1)+0.15);
     peakarea[4] = 200*f5->Integral(f5->GetParameter(1)-0.05,f5->GetParameter(1)+0.15);
     peakarea[5] = 200*f6->Integral(f6->GetParameter(1)-0.05,f6->GetParameter(1)+0.15);
     peakarea[6] = 200*f7->Integral(f7->GetParameter(1)-0.05,f7->GetParameter(1)+0.15);
     peakarea[7] = 200*f8->Integral(f8->GetParameter(1)-0.05,f8->GetParameter(1)+0.15);
     peakarea[8] = 200*f9->Integral(f9->GetParameter(1)-0.05,f9->GetParameter(1)+0.15);
     peakarea[9] = 200*f10->Integral(f10->GetParameter(1)-0.05,f10->GetParameter(1)+0.15);

     peak[0] = f1->GetParameter(1);
     peak[1] = f2->GetParameter(1);
     peak[2] = f3->GetParameter(1);
     peak[3] = f4->GetParameter(1);
     peak[4] = f5->GetParameter(1);
     peak[5] = f6->GetParameter(1);
     peak[6] = f7->GetParameter(1);
     peak[7] = f8->GetParameter(1);
     peak[8] = f9->GetParameter(1);
     peak[9] = f10->GetParameter(1);
     
     h[d]->GetXaxis()->SetRangeUser(4., 7.);
     
     //if you want save fits (only for i == 0 as the reapeat for other peaks)
     /*
       if (i == 0){
     c5->Update();
     c5->Print(Form("c5_peak%d_%d.pdf",i,d));
     }
       
     */
    
     
     delete c5;
    
     //gr[i] = new TGraphErrors();
    

     //cout<<"Peak centered at "<< peak[i] << ", area = " << peakarea[i] <<endl;
     //int number;
     //cout << "Enter : ";
     //cin >> number;
     //cout << h[d] << endl;
   

        //gPad->SetLogy(1);
         //if (d%5==0){
          //can_id++;
          //c2[can_id] = new TCanvas();
          //c2[can_id]->Divide(1,5);
        //}
        //c2[can_id]->cd(d%5+1);
	//h[d]->SaveAs(Form("c2%d.pdf",d));

     for(int i = 0; i < 10; i++ ){
         output.open("out1.txt");
	
 
     theta_cm[d] = asin(gamma*sin(Angles_plot[d])) + Angles_plot[d] ; 
 
         if(d<4) {CII = 9585027*6E-9; eff = 0.9;}//9.59e+06*6E-9;//4 degrees runs 250-251 DIFFERENT CIIs
         if(d>3&& d<8){ CII =8304488*20E-9; eff =0.9;}// 10.0e+06*20E-9;//7 degrees runs 65-78
         if(d>7 && d<12) {CII =1498575*20E-9; eff = 0.9;}//1.49e+06*20E-9; //10 degrees runs 79-80
         if(d>11 && d<16){ CII =4588264 *20E-9; eff = 0.9;}//; 4.66e+06*20E-9; //15degress runs 105-110
         if(d>15&& d<20){ CII =2960057*20E-9; eff = 0.9;}// 2.96e+06*20E-9; //21 degrees runs 120-121
         if(d>19 ) {CII = 4463435*20E-9; eff = 0.9;}//4.46e+06*20E-9; //28 degrees runs 140-142 

      If[d] =  CII/(1000*qe);//incident flux
     if(d>3&& d<8) cxn[d]= 2.25*10e27*peakarea[i]/(tAD*If[d]*eff*SA[d]*10E-03);//scale 7 deg down
     else if(d>7 && d<12) cxn[d]= 2.26*10e27*peakarea[i]/(tAD*If[d]*eff*SA[d]*10E-03);
     else if(d>11 && d<16) cxn[d]=  2.68*10e27*peakarea[i]/(tAD*If[d]*eff*SA[d]*10E-03); // 15 deg
     else if(d>15&& d<20) cxn[d]= 2.82*10e27*peakarea[i]/(tAD*If[d]*eff*SA[d]*10E-03);
     else if(d>19 ) cxn[d]= 2.97*10e27*peakarea[i]/(tAD*If[d]*eff*SA[d]*10E-03);
     else cxn[d]= 10e27*peakarea[i]/(tAD*If[d]*eff*SA[d]*10E-03);
    
     error[d]=cxn[d]*sqrt(1/peakarea[i] );
  
    output<< d << " " << theta_cm[d] <<" " << peak[i] << " " << cxn[d] << " " << error[d] << " " << peakarea[i] <<endl;
    if(i == 0){gROOT->ProcessLine(Form(".! head out1.txt >> peak_%d.txt",i));}
    else if(i == 1){gROOT->ProcessLine(Form(".! head out1.txt >> peak_%d.txt",i));}
    else if(i == 2){gROOT->ProcessLine(Form(".! head out1.txt >> peak_%d.txt",i));}
    else if(i == 3){gROOT->ProcessLine(Form(".! head out1.txt >> peak_%d.txt",i));}
    else if(i == 4){gROOT->ProcessLine(Form(".! head out1.txt >> peak_%d.txt",i));}
    else if(i == 5){gROOT->ProcessLine(Form(".! head out1.txt >> peak_%d.txt",i));}
    else if(i == 6){gROOT->ProcessLine(Form(".! head out1.txt >> peak_%d.txt",i));}
    else if(i == 7){gROOT->ProcessLine(Form(".! head out1.txt >> peak_%d.txt",i));}
    else if(i == 8){gROOT->ProcessLine(Form(".! head out1.txt >> peak_%d.txt",i));}
    else if(i == 9){gROOT->ProcessLine(Form(".! head out1.txt >> peak_%d.txt",i));}
    gROOT->ProcessLine(".!rm out1.txt");
    output.close();
      
   
 }
// gROOT->ProcessLine(Form(".!mv out1.txt peak_%d.txt",i));
 
  
     
}
output.close();
//you can read those files with the for loop and draw here with the (5,2) canvas
/*
for (int peak = 0;peak < 10;peak++){
   
}
*/
}

double FitFunction(double *x, double *pars)
{
    double result = 0;
    
    int size = sizeof(pars)/sizeof(pars[0]);
    
    f1->SetParameter(0,pars[0]); //size of first Gaussian
    f2->SetParameter(0,pars[1]); //size of second Gaussian
    f3->SetParameter(0,pars[2]); //size of third Gaussian
    f4->SetParameter(0,pars[3]); //size of third Gaussian
    f5->SetParameter(0,pars[4]); //size of third Gaussian
    f6->SetParameter(0,pars[5]); //size of third Gaussian
    f7->SetParameter(0,pars[6]); //size of third Gaussian
    f8->SetParameter(0,pars[7]); //size of third Gaussian
    f9->SetParameter(0,pars[8]); //size of third Gaussian
    f10->SetParameter(0,pars[9]); //size of third Gaussian

    f1->SetParameter(1,pars[10]); //position of one Gaussian
    f2->SetParameter(1,pars[11]); //position of second Gaussian
    f3->SetParameter(1,pars[12]); //position of third Gaussian
    f4->SetParameter(1,pars[13]); //position of third Gaussian
    f5->SetParameter(1,pars[14]); //position of third Gaussian
    f6->SetParameter(1,pars[15]); //position of third Gaussian
    f7->SetParameter(1,pars[16]); //position of third Gaussian
    f8->SetParameter(1,pars[17]); //position of third Gaussian
    f9->SetParameter(1,pars[18]); //position of third Gaussian
    f10->SetParameter(1,pars[19]); //position of third Gaussian

    f1->SetParameter(2,pars[20]); //Gaussian width
    f2->SetParameter(2,pars[20]); //also Gaussian width, and the same one
    f3->SetParameter(2,pars[20]); //also Gaussian width, and the same one
    f4->SetParameter(2,pars[20]); //also Gaussian width, and the same one
    f5->SetParameter(2,pars[20]); //also Gaussian width, and the same one
    f6->SetParameter(2,pars[20]); //also Gaussian width, and the same one
    f7->SetParameter(2,pars[20]); //also Gaussian width, and the same one
    f8->SetParameter(2,pars[20]); //also Gaussian width, and the same one
    f9->SetParameter(2,pars[20]); //also Gaussian width, and the same one
    f10->SetParameter(2,pars[20]); //also Gaussian width, and the same one

    f1->FixParameter(3,1); //size of the Landau - set to 1 because otherwise there's two amplitudes
    f2->FixParameter(3,1); //same as above
    f3->FixParameter(3,1); //same as above
    f4->FixParameter(3,1); //same as above
    f5->FixParameter(3,1); //same as above
    f6->FixParameter(3,1); //same as above
    f7->FixParameter(3,1); //same as above
    f8->FixParameter(3,1); //same as above
    f9->FixParameter(3,1); //same as above
    f10->FixParameter(3,1); //same as above

    f1->SetParameter(4,0); //to avoid shifting?
    f2->SetParameter(4,0);
    f3->SetParameter(4,0);
    f4->SetParameter(4,0);
    f5->SetParameter(4,0);
    f6->SetParameter(4,0);
    f7->SetParameter(4,0);
    f8->SetParameter(4,0);
    f9->SetParameter(4,0);
    f10->SetParameter(4,0);

    f1->SetParameter(5,pars[21]); //Landau spreading
    f2->SetParameter(5,pars[21]); //Landau spreading
    f3->SetParameter(5,pars[21]); //Landau spreading
    f4->SetParameter(5,pars[21]); //Landau spreading
    f5->SetParameter(5,pars[21]); //Landau spreading
    f6->SetParameter(5,pars[21]); //Landau spreading
    f7->SetParameter(5,pars[21]); //Landau spreading
    f8->SetParameter(5,pars[21]); //Landau spreading
    f9->SetParameter(5,pars[21]); //Landau spreading
    f10->SetParameter(5,pars[21]); //Landau spreading

    result = f1->Eval(x[0]) + f2->Eval(x[0]) + f3->Eval(x[0]) + f4->Eval(x[0]) + f5->Eval(x[0])   + f6->Eval(x[0])   + f7->Eval(x[0])   + f8->Eval(x[0])   + f9->Eval(x[0])   + f10->Eval(x[0])     ;
    


    return result;
}
