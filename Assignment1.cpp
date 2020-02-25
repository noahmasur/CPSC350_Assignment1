#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cctype>

using namespace std;

//Reads file content in
string readFile(string fName){
  ifstream inFS;
  string fCont;
  string out;

  inFS.open(fName);

//Testing if file was able to be opened
  if(!inFS.is_open()){
    return "N/A";
  }

  //Reading file line by line and adding each line to an output string
  while(!inFS.eof()){
    inFS >> fCont;

    if(!inFS.fail()){
      out += fCont + "\n";
    }
  }

  inFS.close();
  return out;
}

//Calculates summary statistics
string calcStats(string fCont, double &mean, double &stDev){
  string out;

  int lCount = 0; //line count
  int sum = 0; //total char count

  int cCount = 0; //character count per line
  double vSum = 0.0;
  double variance;

  //Finds count of lines and total count of characters for mean
  for(char c : fCont){
    if(c == '\n'){
      lCount += 1;
    } else{
      sum += 1;
    }
  }
  mean = sum/(double)lCount;

  //Calculates variance
  for(char c : fCont){
    if(c == '\n'){
      vSum += pow((cCount - mean), 2);
      cCount = 0;
    } else{
      cCount += 1;
    }
  }
  variance = vSum / lCount;
  stDev = sqrt(variance);

  out = "Sum: " + to_string(sum) + "\nMean: " + to_string(mean) + "\nVariance: " + to_string(variance) + "\nStandard Deviation: " + to_string(stDev);

  return out;
}

//Calculates probabilities of each nucleotide and nucleotide bigram
string calcProb(string fCont, double &a_prob, double &c_prob, double &t_prob, double &g_prob){
  string out;

  int char_count = 0;

  int a_count = 0;
  int c_count = 0;
  int t_count = 0;
  int g_count = 0;

  int aa_count = 0;
  int ac_count = 0;
  int at_count = 0;
  int ag_count = 0;

  int ca_count = 0;
  int cc_count = 0;
  int ct_count = 0;
  int cg_count = 0;

  int ta_count = 0;
  int tc_count = 0;
  int tt_count = 0;
  int tg_count = 0;

  int ga_count = 0;
  int gc_count = 0;
  int gt_count = 0;
  int gg_count = 0;

  char c_before;

  //Iterating through string to find counts
  for(char c : fCont){
    //Resets charcter line count if it reaches delimiter (\n)
    if(c == '\n'){
      c_before = '\0';
    }else{
      c = toupper(c);
      if(c == 'A'){
        a_count += 1;

        if(c_before == 'A'){
          aa_count += 1;
        }else if(c_before == 'C'){
          ca_count += 1;
        }else if(c_before == 'T'){
          ta_count += 1;
        }else if(c_before == 'G'){
          ga_count += 1;
        }

      }else if(c == 'C'){
        c_count += 1;

        if(c_before == 'A'){
          ac_count += 1;
        }else if(c_before == 'C'){
          cc_count += 1;
        }else if(c_before == 'T'){
          tc_count += 1;
        }else if(c_before == 'G'){
          gc_count += 1;
        }

      }else if(c == 'T'){
        t_count += 1;

        if(c_before == 'A'){
          at_count += 1;
        }else if(c_before == 'C'){
          ct_count += 1;
        }else if(c_before == 'T'){
          tt_count += 1;
        }else if(c_before == 'G'){
          gt_count += 1;
        }
      }else if(c == 'G'){
        g_count += 1;

        if(c_before == 'A'){
          ag_count += 1;
        }else if(c_before == 'C'){
          cg_count += 1;
        }else if(c_before == 'T'){
          tg_count += 1;
        }else if(c_before == 'G'){
          gg_count += 1;
        }
      }
      c_before = c;
      char_count += 1;
    }
  }

  //Calculating probablities from charater counts
  a_prob = a_count/(double)char_count;
  c_prob = c_count/(double)char_count;
  t_prob = t_count/(double)char_count;
  g_prob = g_count/(double)char_count;

  double aa_prob = aa_count/(double)char_count;
  double ac_prob = ac_count/(double)char_count;
  double at_prob = at_count/(double)char_count;
  double ag_prob = ag_count/(double)char_count;

  double ca_prob = ca_count/(double)char_count;
  double cc_prob = cc_count/(double)char_count;
  double ct_prob = ct_count/(double)char_count;
  double cg_prob = cg_count/(double)char_count;

  double ta_prob = ta_count/(double)char_count;
  double tc_prob = tc_count/(double)char_count;
  double tt_prob = tt_count/(double)char_count;
  double tg_prob = tg_count/(double)char_count;

  double ga_prob = ga_count/(double)char_count;
  double gc_prob = gc_count/(double)char_count;
  double gt_prob = gt_count/(double)char_count;
  double gg_prob = gg_count/(double)char_count;

  //Adding each variable to output string seperated by a new line
  out = "\nA count: " + to_string(a_count) + "\nC count: " + to_string(c_count) + "\nT count: " + to_string(t_count) + "\nG count: " + to_string(g_count) +
  "\nAA count: " + to_string(aa_count) + "\nAC count: " + to_string(ac_count) + "\nAT count: " + to_string(at_count) + "\nAG count: " + to_string(ag_count) +
  "\nCA count: " + to_string(ca_count) + "\nCC count: " + to_string(cc_count) + "\nCT count: " + to_string(ct_count) + "\nCG count: " + to_string(cg_count) +
  "\nTA count: " + to_string(ta_count) + "\nTC count: " + to_string(tc_count) + "\nTT count: " + to_string(tt_count) + "\nTG count: " + to_string(tg_count) +
  "\nGA count: " + to_string(ga_count) + "\nGC count: " + to_string(gc_count) + "\nGT count: " + to_string(gt_count) + "\nGG count: " + to_string(gg_count) +
  "\n \nA Probability: " + to_string(a_prob) + "\nC probability: " + to_string(c_prob) + "\nT probability: " + to_string(t_prob) + "\nG probability: " + to_string(g_prob) +
  "\nAA Probability: " + to_string(aa_prob) + "\nAC probability: " + to_string(ac_prob) + "\nAT probability: " + to_string(at_prob) + "\nAG probability: " + to_string(ag_prob) +
  "\nCA Probability: " + to_string(ca_prob) + "\nCC probability: " + to_string(cc_prob) + "\nCT probability: " + to_string(ct_prob) + "\nCG probability: " + to_string(cg_prob) +
  "\nTA Probability: " + to_string(ta_prob) + "\nTC probability: " + to_string(tc_prob) + "\nTT probability: " + to_string(tt_prob) + "\nTG probability: " + to_string(tg_prob) +
  "\nGA Probability: " + to_string(ga_prob) + "\nGC probability: " + to_string(gc_prob) + "\nGT probability: " + to_string(gt_prob) + "\nGG probability: " + to_string(gg_prob) + "\n";

  return out;
}

//Outputs 1000 DNA strings that follow a Gaussian distribution
string gaussian(double mean, double stDev, double a_prob, double c_prob, double t_prob, double g_prob){
  string out = "";
  double a;
  double b;
  double c;
  double d;

  //Converting probablities into percentiles
  c_prob += a_prob;
  t_prob += c_prob;
  g_prob += t_prob;

  for(int i = 0; i < 1000; i++){
    //Finding random variables and converting them using gaussian formulas
    a = ((double)rand() / (RAND_MAX));
    b = ((double)rand() / (RAND_MAX));
    c = sqrt(-2 * log(a)) * cos(2 * M_PI * b);
    d = (int)(stDev * c + mean); //length of DNA string rounded

    //Adding letters in DNA string of length d based on probability of letters
    for(int x = 0; x < d; x++){
      double letter_prob = ((double)rand() / (RAND_MAX));

      if(letter_prob < a_prob){
        out += "A";
      } else if(letter_prob < c_prob){
        out += "C";
      } else if(letter_prob < t_prob){
        out += "T";
      } else if(letter_prob < g_prob){
        out += "G";
      }
    }

    out += "\n";
  }

  return out;
}

//Writes statitstics to output file
string writeFile(string stats, string prob, string gaus){
  ofstream outFS;

  outFS.open("NoahMasur.txt", ios::app);
  outFS << stats << endl << prob << endl << gaus << endl;
  outFS.close();

  return "Done";
}

//Writes the header at the top of the file
void writeHeader(){
  ofstream outFS;

  outFS.open("NoahMasur.txt");
  outFS << "Noah Masur \n2327080 \nmasur@chapman.edu \n" << endl;
  outFS.close();
}

int main(int argc, char const *argv[]) {
  string fName; //file name
  string fCont; //file content

  double mean;
  double stDev;

  double a_prob;
  double c_prob;
  double t_prob;
  double g_prob;

  string stats;
  string prob;
  string gaus;

  char choice = 'y';

  //Checks for command line parameter
  if( argc < 2){
    cout << "Invalid command line parameters" << endl;
    return 1;
  }

  writeHeader();

  fName = argv[1];

  while(choice != 'n'){
    fCont = readFile(fName);

    if(choice != 'y'){
      //If choice input was not 'y' or 'n', inform user and skip to question again
      cout << "Please enter a valid input. \n" << endl;
    }else if(fCont == "N/A"){
      //If read file was unable to open
      cout << "Could not open " << fName<< endl;
    } else{
      stats = calcStats(fCont, mean, stDev);
      prob = calcProb(fCont, a_prob, c_prob, t_prob, g_prob);
      gaus = gaussian(mean, stDev, a_prob, c_prob, t_prob, g_prob);
      cout << writeFile(stats, prob, gaus) << endl;
    }

    cout << "Do you want to input from another file? (y/n)" << endl;
    cin >> choice;

    if(choice == 'y'){
      cout << "Enter the name of the file you want to input from: " << endl;
      cin >> fName;
    }
  }

  return 0;
}
