/*
Authors: S.Yahia-Cherif, F.Dournac.
Last Update 03/03/2021.
This is the main script of XSAF. XSAF computes the photometric Cl and the photometric Fisher matrix
*/

#include "XSAF_C_gnu.h"

using namespace std;

int main(){

    // Standard or Common bias
    string strModel;
    // Subcase : opt or pess
    string strSubCase;
    // Curvature
    string strCurvature;
    // Gamma/No-GAMMA
    string strGamma;
    // Zcut/No-Zcut
    string strZcut;
    // Compute Pz
    string strPz;
    // Compute Cl's
    string compute_Cl;

    // Number of bins
    int numBins;

    // Fix last issue - GitHub 03/03/21
    int use_GCb = 1;

    //Suppress the sub Fishers before the total Fisher matrix is built.
    system("rm -f output/Fisher_*");

    // Read parameters run of EUCLID survey
    ifstream modelFile("./XCPP_params.txt");
    // Read multiple lines
      modelFile >> strModel;
      modelFile >> strSubCase;
      modelFile >> strCurvature;
      modelFile >> strGamma;
      modelFile >> strZcut;
      modelFile >> strPz;
      modelFile >> compute_Cl;
    modelFile.close();

    // Check different parameters
    cout << "model= " << strModel << endl;
    cout << "subCase = " << strSubCase << endl;
    cout << "curvature = " << strCurvature << endl;
    cout << "gamma = " << strGamma << endl;
    cout << "zcut = " << strZcut << endl;
    cout << "Photoz = " << strPz << endl;

/*
    //Load the parameters files and stock them in dictionaries.
    int count=0;
    ifstream ifile("../QTLauncher/Parameters_W.txt");
    while(!ifile.fail()){
        elements_0.push_back("");
        elements_1.push_back(0.0);
        ifile>>elements_0[count];
        ifile>>elements_1[count];
        count++;
    }
    ifile.close();

    vector<int> PAR_X;
    int PAR_IND = 0;
    for(int i=0; i<elements_0.size(); i++){
        if(elements_0[i] != "Usesp_ch" && elements_0[i] != "Usesv_ch" && elements_0[i] != "UseGCspecbias_ch" && elements_0[i] != "UseGCspecPS_ch"){
            if(elements_0[i].find("Use") != string::npos && elements_1[i] == 0){
                PAR_X.push_back(PAR_IND);
            }
            PAR_IND++;
        }
    }

    ifstream ifile2("../QTLauncher/Codes_W.txt");
    while(!ifile2.fail()){
        elements_0.push_back("");
        elements_1.push_back(0.0);
        ifile2>>elements_0[count];
        ifile2>>elements_1[count];
        count++;
    }
    ifile2.close();
    ifstream ifile3("../QTLauncher/XSAF_W.txt");
    while(!ifile3.fail()){
        elements_0.push_back("");
        elements_1.push_back(0.0);
        ifile3>>elements_0[count];
        ifile3>>elements_1[count];
        count++;
    }
    ifile3.close();

    ifstream ifile4("../QTLauncher/Extra_W.txt");
    while(!ifile4.fail()){
        elements_0.push_back("");
        elements_1.push_back(0.0);
        ifile4>>elements_0[count];
        ifile4>>elements_1[count];
        count++;
    }
    ifile4.close();
*/
    // Reading main parameters file
    vector<string> elements_0; vector<double> elements_1;
    //Load the parameters files and stock them in dictionaries.
    int count=0;
    ifstream ifile("./params.txt");
    while(!ifile.fail()){
        elements_0.push_back("");
        elements_1.push_back(0.0);
        ifile>>elements_0[count];
        ifile>>elements_1[count];
        count++;
     }
     ifile.close();

    map<string, double> XSAF_elts;
    for(int i=0; i<elements_0.size(); i++){
        XSAF_elts[elements_0[i]] = double(elements_1[i]);
    }

    // Include 11th bin for common bias but needs fiducial values for zrange
    if (strModel == "S") {
      XSAF_elts["VRS_bins"] = 10;
      }
    else if (strModel == "E") {
      // Change number of bins : 11 with common bias
      XSAF_elts["VRS_bins"] = 11;
      // Test for rescaling : only 5 bias considered
      //XSAF_elts["VRS_bins"] = 5;
      }
    // Not Recommended since we only process 5 bias
    else if (strModel == "C") {
      XSAF_elts["VRS_bins"] = 5;
      }
    else if (strModel == "O") {
      cout << "Version Equidistant redshift" << endl;
      }
    else { cout << "Error on model for main_intel.cpp" << endl; }

    cout << "1) Number of photo bias : " << XSAF_elts["VRS_bins"] << endl;

    vector<int> PAR_X;
    int PAR_IND = 0;
    // Upperbound : test on 17 first parameters into params.txt
    int limit_params = 17;
    //for(int i=0; i<elements_0.size(); i++){
    for(int i=0; i<limit_params; i++){
      if(elements_0[i] != "Usesp_ch" && elements_0[i] != "Usesv_ch" && elements_0[i] != "UseGCspecbias_ch" && elements_0[i] != "UseGCspecPS_ch"){
        if(elements_0[i].find("Use") != string::npos && elements_1[i] == 0){
          PAR_X.push_back(PAR_IND);
        }
        PAR_IND++;
      }
      // Fix last issue - GitHub 03/03/21
      if(elements_0[i] == "UseGCphotbias_ch" && elements_1[i] == 0){
        use_GCb=0;
      }
    }

    // Fix last issue - GitHub 03/03/21
    if(use_GCb == 0){
      PAR_IND = PAR_X[PAR_X.size()-1]+1;
      for(int i=PAR_IND; i<PAR_IND + XSAF_elts["VRS_bins"]-1; i++){
        PAR_X.push_back(i);
      }
    }

// Parameter subCase : PESSIMISTIC or OPTIMISTIC cases = OPT/PESS
  if (strSubCase == "PESS") {
    // Original version
    XSAF_elts["Vlmax_GC"] = 750;
    XSAF_elts["Vlmax_WL"] = 1500;
  }
  else if (strSubCase == "OPT") {
    // Original version
    XSAF_elts["Vlmax_GC"] = 3000;
    XSAF_elts["Vlmax_WL"] = 5000;
  }
  else {
    cout << "Error on PESS/OPT cases" << endl;
    exit(0);
  }

  if (strCurvature == "F") {
    XSAF_elts["FNF_ch"] == 0;
    XSAF_elts["UseOmegaDE_ch"] == 1;
    }
  else if (strCurvature == "NF") {
    XSAF_elts["FNF_ch"] == 1;
    XSAF_elts["UseOmegaDE_ch"] == 0;
    }
  else { 
    cout << "Error on Curvature Y/N" << endl;
    exit(0);
  }

  if (strGamma == "N") {
    XSAF_elts["Usegamma_ch"] = 1;
  }  
  else if (strGamma == "Y") 
    XSAF_elts["Usegamma_ch"] = 0;
  else { 
    cout << "Error on Gamma Y/N" << endl;
    exit(0);
  }

  if (strZcut== "N") {
    XSAF_elts["zcut_ch"] = 1;
  }  
  else if (strZcut == "Y") 
    XSAF_elts["zcut_ch"] = 0;
  else { 
    cout << "Error on Zcut Y/N" << endl;
    exit(0);
  }

/*
    if(XSAF_elts["zcut_ch"] == 1){
        zcut = "N";
    }
    else{
        zcut = "Y";
    }
    if(XSAF_elts["Usegamma_ch"] == 1){
        gma = "N";
    }
    else{
        gma = "Y";
    }
*/

  //Call default constructor.
  XC XC(XSAF_elts["VRS_bins"], 5, XSAF_elts["Vlnum"], XSAF_elts["Vlmin_GCWL"], XSAF_elts["Vlmax_WL"], XSAF_elts["Vlmax_GC"], 60, XSAF_elts["Vprec_Int_z"], 3, PAR_X.size(), strCurvature, strGamma, strZcut, strModel);

    //Call the methods to compute the Cl.
	for(int X_ind=0; X_ind<PAR_X.size(); X_ind++){
        XC.Initializers_G(PAR_X[X_ind]);
        XC.Initializers_Pk();
        XC.background();
        // DEV - remark : UsePZ_ch not used into main.cpp
        //if(X_ind==0 && XSAF_elts["UsePZ_ch"] == 0){
        if(X_ind==0 && strPz == "Y") {
            //XSAF_elts["UsePZ_ch"] == 0
            cout << "2) Number of photo bias : " << XSAF_elts["VRS_bins"] << endl;
            XC.photoz();
        }
    	XC.photoz_load();
    	XC.windows();
        if (compute_Cl == "Y")
    	  XC.C_l_computing();
    }

    //Call the method to compute the Fisher matrix.
    int run_index = 0;
    // Original : loop on 60 multipoles
    for(int i=0; i<XSAF_elts["Cutting_l_V"]; i++){
    // Only once
        if(i == XSAF_elts["Cutting_l_V"]-1){
            run_index = 1;
        }
        // Comment if we want the FoM of all combinations
        if(XSAF_elts["UseXC_ch"] == 0){
            XC.Fisher("XC", "Fisher_"+strSubCase+"_GCph_WL_XC_"+strCurvature+"_"+strGamma+"_XSAF", XSAF_elts["Vlnum"]/XSAF_elts["Cutting_l_V"]*i, XSAF_elts["Vlnum"]/XSAF_elts["Cutting_l_V"]*(i+1), run_index, i, XSAF_elts["Cutting_l_V"]);
        }
        if(XSAF_elts["UseWL_ch"] == 0){
            XC.Fisher("WL", "Fisher_"+strSubCase+"_WL_"+strCurvature+"_"+strGamma+"_XSAF", XSAF_elts["Vlnum"]/XSAF_elts["Cutting_l_V"]*i, XSAF_elts["Vlnum"]/XSAF_elts["Cutting_l_V"]*(i+1), run_index, i, XSAF_elts["Cutting_l_V"]);
        }
        if(XSAF_elts["UseGC_ch"] == 0){
            XC.Fisher("GCp", "Fisher_"+strSubCase+"_GCph_"+strCurvature+"_"+strGamma+"_XSAF", XSAF_elts["Vlnum"]/XSAF_elts["Cutting_l_V"]*i, XSAF_elts["Vlnum"]/XSAF_elts["Cutting_l_V"]*(i+1), run_index, i, XSAF_elts["Cutting_l_V"]);
        }
    }

    cout<<endl<<endl;
    return 0;
}
