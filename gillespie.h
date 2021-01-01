#ifndef GILLESPIE_H_
#define GILLESPIE_H_

#include "parser.h"
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>
#include <memory>
#include <cmath>
#include <iomanip>

using namespace std;

// Base Class for Propensity Functions
class Propensity{
protected:
  double *c, *x1, *x2, a1;
  double volume;
public:
  virtual double Get_Propensity(void){ return 0; };
};

// Derived Class for Propensity Function of Spontaenous Formation Reactions
class Form: public Propensity{
public:
  Form(double *k){ c=k; };
  double Get_Propensity(void){ return (double) 1 * (*c); };
};

// Derived Class for Propensity Function of Unimolecular Reactions
class Unimol: public Propensity{
public:
  Unimol(double *k, double *x){ c=k; x1=x; };
  double Get_Propensity(void){ return (*c) * (*x1); };
};

// Derived Class for Propensity Function of Bimolecular Reactions (Same Species)
class Bimol1: public Propensity{
public:
  Bimol1(double *k, double *x, double z){ c=k; x1=x; volume=z;};
  double Get_Propensity(void){ return (*c * 2 / volume) * (*x1) * ((*x1)-1) / 2; };
};

// Derived Class for Propensity Function of Bimolecular Reactions (Different Species)
class Bimol2: public Propensity{
public:
  Bimol2(double *k, double *x, double *y, double z){ c=k; x1=x; x2=y; volume=z;};
  double Get_Propensity(void){ return (*c/volume) * (*x1) * (*x2); };
};

class Gillespie{
protected:
  double t=0; // time
  double a0; // total a_0
  double volume; // system volume
  int n_iter; // number of iterations
  vector<string> names_vec; // vector for names of species
  vector<double> state_vec; // vector for amounts of species
  vector<reaction> rxn_vecs; // reactions are represented as transition vectors
  vector<double> propensity_vec; // vector for propensities
  vector<double> prop_cumul_sum; // cumulative sum of propensities, to get the range in (0,1)
  vector<unique_ptr<Propensity>> propensities; // vector for classes of propensity functions
  vector<double> k_vec; // vector for rate constants of reactions
public:
  Gillespie(vector<string>, vector<double>, vector<reaction>, vector<double>, double, int); // constructor that takes data from the Parser class
  void Start(void); // main function of the class; does the simulation
  void SetPropensityFunctions(void); // creates propensity function class for each reaction
  vector<int> CountReactants(reaction); // counts the number of reactants; determines if reaction is uni or bimolecular
  vector<double> SetPropensities(void); // populates the propensity vector, ai/a0
  double Random(void); // Mersenne Twister random number generator
  vector<double> CumulSum(vector<double>); // takees the cumulative sum of a vector
  int WillReact(void); // determines which reaction will occur
  void React(int); // process of reaction; updates the state vector
  bool IsValid(vector<double>);
};

// takes data from parser and saves them to this class's own data
Gillespie::Gillespie(vector<string> a, vector<double> b, vector<reaction> c, vector<double> d, double e, int f){
  names_vec = a;
  state_vec = b;
  rxn_vecs = c;
  k_vec = d;
  volume = e;
  n_iter = f;
}

void Gillespie::Start(void){
  SetPropensityFunctions();

  // outputs data as Time-Species1-Species2-...Speciesn
  cout << "Time\t";
  for(auto i : names_vec){
    cout << i << '\t';
  }
  cout << endl;

  unsigned int counter=0;
  do{
    // sends data to stdout (doesn't save to file YET)
    printf("%.10f\t", t);
    for(auto i : state_vec){
      printf("%.10f\t", i);
    }
    cout << endl;

    // Step 1. Calculate/Update propensities.
    propensity_vec = SetPropensities();

    // Step 2. Figure out which reaction will react next.
    int will_react = WillReact();

    // Step 3. Perform the chosen reaction.
    React(will_react);

    // Step 4. Compute the time step.
    double rand = Random();
    t += (1 / a0) * log(1 / rand);

    counter++;
  }while(counter < n_iter && IsValid(propensity_vec)); // number of iterations is user-specified
  // prematurely ends simulation if nan or inf values are observed
}

// assigns propensity functions to each reaction
void Gillespie::SetPropensityFunctions(){
  for(size_t i=0; i<k_vec.size(); i++){
    vector<int> reactants = CountReactants(rxn_vecs[i]);
    if(reactants[reactants.size()-1] < 1){
      propensities.emplace_back(new Form(&k_vec[i]));
    }
    // Unimolecular
    else if(reactants[reactants.size()-1] < 2){
      propensities.emplace_back(new Unimol(&k_vec[i], &state_vec[reactants[0]]));
    }
    // Bimolecular, Same Species
    else if(reactants.size() < 3){
      propensities.emplace_back(new Bimol1(&k_vec[i], &state_vec[reactants[0]], volume));
    }
    // Bimolecular, Different Species
    else{
      propensities.emplace_back(new Bimol2(&k_vec[i], &state_vec[reactants[0]], &state_vec[reactants[1]], volume));
    }
  }
}  

// counts the reactants in each reaction
vector<int> Gillespie::CountReactants(reaction x){
  vector<int> n;
  int number=0;
  for(size_t i=0; i<x.rxn_vec.size(); i++){
    // lists indices of reactants
    if(x.rxn_vec[i] < 0){
      n.push_back(i);
      number = number - x.rxn_vec[i];
    }
  }
  // also takes note of number of each species involved
  n.push_back(number);
  return n;
}

// populates the propensity vector, a1/a0
vector<double> Gillespie::SetPropensities(void){
  vector<double> output;
  // calculates propensity of each reaction
  for(size_t i=0; i<propensities.size(); i++){
    output.push_back(propensities[i]->Get_Propensity());
  }
  // takes the total sum
  a0 = accumulate(output.begin(), output.end(), 0);
  // scalar division
  transform(output.begin(), output.end(), output.begin(), bind1st(multiplies<double>(),(1 / a0)));
  return output;
}
  
// draws a random number from the interval (0,1)
double Gillespie::Random(void){
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  mt19937 generator(seed);
  uniform_real_distribution<double> dis(0.0, 1.0);
  return dis(generator);
}

// returns the cumulative sum of a vector
vector<double> Gillespie::CumulSum(vector<double> input){
  vector<double> output;
  double sum=0;
  for(auto i : input){
    sum = sum + i;
    output.push_back(sum);
  }
  return output;
}

// determines which reaction will occur next 
int Gillespie::WillReact(void){
  // takes the cumulative sum of the propensity vector; analogous to laying the fractions on the number line
  prop_cumul_sum = CumulSum(propensity_vec);

  double rand = Random();
  // subtracting the random number from each element
  transform(prop_cumul_sum.begin(), prop_cumul_sum.end(), prop_cumul_sum.begin(), bind1st(plus<double>(),-rand));

  // finding the smallest positive integer, which corresponds to the reaction chosen
  for(unsigned int i=0; i<prop_cumul_sum.size(); i++){
    if(prop_cumul_sum[i] > 0){
      return i;
      break;
    }
  }
  exit(1);
}

// performs the reaction
void Gillespie::React(int ind){
  // pointers are necessary for inplace update of the amounts
  double *state_ptr = &state_vec[0];
  for(size_t i=0;i<state_vec.size(); i++){
    *(state_ptr+i) += rxn_vecs[ind].rxn_vec[i];
  }
}

// ensures amounts don't diverge to inf or nan
// prevents negative values on amounts
// ends the simulation if such divergence is reached
bool Gillespie::IsValid(vector<double> input){
  for(auto i : input){
    if(!isfinite(i)){
      return false;
      break;
    }
  }
  return true;
}

#endif
