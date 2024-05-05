// josie 23 dec 2023
// run python script to calculate baselines, OPDs, uv coords, WITHOUT a root file to read in - just read a text file for CTA simulations

void rootRunPySim(TString filename){

  // =================== get necessary info from headers and run python script to calculate delays ========================================
  TString runpy;
  TString py = "python3 delays/CalcMvtUVonly.py";
  TString length = "0";
  TString loctime = "0";
  TString nframes = "0";
  TString frameSize = "0";

  cout << "run taken on " << *(source) << "    version " << version << endl;

  if(version != -4){ cout << "wrong version, bye!" << endl; return;}

  double longestRun(0.0);
  int mostFrames(0);
  double frameWidth(0.0);
  TDatime* ltime = new TDatime;  ltime->Set(2030,1,1,23,59,59); // init to some time far in the future, so that all times in runs would be less
  TDatime* timetemp = new TDatime;


  cout << "checking python parameters " << length << "  " << nframes << "  " << frameSize << "  " << *(source) << "  " << loctime << endl << endl;

  // run python script to calculate delays/baselines via bash from within this macro
  cout << "================> Running python script to calculate delays <=====================" << endl;
  runpy = py + " " + length + " " + nframes + " " + frameSize + " " + *(source) + " " + loctime; // note: source is a pointer so must be dereferenced to add to other strings
  gSystem->Exec(runpy);
  cout << "===============================> Python is done! <================================" << endl;

}
