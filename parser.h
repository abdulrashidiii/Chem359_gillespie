#ifndef PARSER_H_
#define PARSER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

// custom data structure for reaction; represented as transition vectors
struct reaction{
  vector<int> rxn_vec;
};

class Parser{
protected:
  ifstream infile; // file stream
  string input_file; // input file name
  int n_species, n_rxn; // number of species and number of reactions
  double volume; // system volume
  int n_iter; // number of iterations
  vector<string> species_vec; // vector for species names
  vector<double> amounts_vec; // vector for species amounts
  vector<double> k_vec; // vector for reaction rate constants
  vector<reaction> rxns_vec; // vector for reaction transition vectors
public:
  Parser(string); // constructor that checks if filename is valid, and then performs the parsing
  void ExtractContents(void); // reads all contents of the input file, and then saves them to appropriate data structures
  void ReadSpecies(void); // read and save species
  void ReadReactions(void); // read and save reactions
  void ReadRateConstant(string*); // extracts rate constant
  void ReactionToTransitionVector(string*, vector<int>*); // transforms a reaction to a transition vector
  void ScalarMultiplication(vector<int>*, int); // scalar multiplication of a vector
  void UpdateReactionVectors(string, vector<int>*); // creates the transition vector of a reaction
  void UpdateSpeciesinReactionVector(int, string, vector<int>*); // updates one of the species in the reaction vector
  bool ContainsDigits(string); // test function checking for numerical values in a string
  int GetIndex(string); // get index of species in species vector

  // returning data
  vector<string> GetNamesVec(){ return species_vec; }; // returns name vector
  vector<double> GetStateVec(){ return amounts_vec; }; // returns state vector
  vector<reaction> GetReactionsVector(){ return rxns_vec; }; // returns reactions vector
  vector<double> GetRateConstantVec(){ return k_vec; }; // returns rate constant vector
  double GetVolume(){ return volume; }; // returns system volume
  int GetNIterations(){ return n_iter; }; // returns number of iterations
};

// constructor
Parser::Parser(string fname) : input_file(fname){
  infile.open(input_file);
  // throws an error if filename is not found
  if(!infile){
    cout << "Unable to open file" << endl;
    exit(1);
  }

  // main parser function
  ExtractContents();
}

void Parser::ExtractContents(void){
  infile >> volume >> n_iter; // first line of input file should be "volume - number of iterations"
  ReadSpecies();
  ReadReactions();
}

// names and amounts of the species are saved into separate vectors
void Parser::ReadSpecies(void){
  string buff; // buffer string variable for unnecessary text

  infile >> buff >> n_species; // number of species are specified beforehand
  for(unsigned i=0; i<n_species; i++){
    string species;
    double amount;
    infile >> species >> amount;
    species_vec.push_back(species);
    amounts_vec.push_back(amount);
  }
}

// reactions are represented as transition vectors
void Parser::ReadReactions(void){
  string buff; // buffer string variable

  infile >> buff >> n_rxn; // number of reactions 
  getline(infile, buff);
  for(unsigned i=0; i<n_rxn; i++){
    vector<int> rxn_i_vec(n_species, 0); // a blank vector with length n_species is initialized

    getline(infile, buff); // each reaction is temporarily saved to buff before being parsed
    ReadRateConstant(&buff); // rate constant is extracted and added to vector of k's
    ReactionToTransitionVector(&buff, &rxn_i_vec); // transforms each reaction to a transition vector

    reaction rxn_i = {rxn_i_vec};
    rxns_vec.push_back(rxn_i);
  }
}

void Parser::ReadRateConstant(string *str_ptr){
  double k; // rate constant saved as double precision float
  auto k_pos = (*str_ptr).find('|');
  k = stof((*str_ptr).substr(k_pos+1));
  k_vec.push_back(k);
  (*str_ptr) = (*str_ptr).substr(0,k_pos-1);
}

void Parser::ReactionToTransitionVector(string *str_ptr, vector<int> *rxn_i_vec_ptr){
  auto delim_pos = (*str_ptr).find('>');
  auto reactants = (*str_ptr).substr(0,delim_pos-1);
  UpdateReactionVectors(reactants, rxn_i_vec_ptr);
  ScalarMultiplication(rxn_i_vec_ptr, -1);
  auto products = (*str_ptr).substr(delim_pos+2);
  UpdateReactionVectors(products, rxn_i_vec_ptr);
}

void Parser::UpdateReactionVectors(string rxnstring, vector<int> *rxn_i_vec_ptr){
  if(rxnstring.compare("0")==0){}
  else{
    size_t delim_ind;
    string token, delimiter=" + ";
    do{
      int coeff;
      delim_ind = rxnstring.find_last_of(delimiter);
      token = rxnstring.substr(delim_ind+1);
      rxnstring = rxnstring.substr(0,delim_ind-2);

      string species;
      if(ContainsDigits(token)){
        istringstream ss(token);
        ss >> coeff >> species;
        UpdateSpeciesinReactionVector(coeff, species, rxn_i_vec_ptr);
      }
      else{
        coeff = 1;
        UpdateSpeciesinReactionVector(coeff, token, rxn_i_vec_ptr);
      }

    } while(delim_ind != string::npos);
  }
}

// checks if a string contains numerical values
bool Parser::ContainsDigits(string s){
  return string::npos != s.find_first_of("0123456789");
}

// updates the value of the species in the reaction vector
void Parser::UpdateSpeciesinReactionVector(int coeff, string species, vector<int> *rxn_i_vec_ptr){
  int ind;
  ind = GetIndex(species);
  rxn_i_vec_ptr->at(ind) += coeff;
}

// locates a species from the species vector
int Parser::GetIndex(string s){
  for(unsigned i=0; i<n_species; i++){
    if(species_vec[i].compare(s)==0){
      return i;
      break;
    }
  }
  // if species is not found in the species vector, the parser ends and throws an error message
  cout<< "Species " << s << " was not initialized!" << endl;
  exit(1);
}

// performs scalar multiplication on a vector
void Parser::ScalarMultiplication(vector<int> *vec_ptr, int multiplier){
  for(size_t i=0; i<(*vec_ptr).size(); i++){
    (*vec_ptr)[i] *= multiplier;
  }
}

#endif /* __PARSER_H__ */
