#include "parser.h"
#include "gillespie.h"
#include <vector>
#include <string>


using namespace std;

int main(){
  Parser test("input");
  vector<string> names_vec = test.GetNamesVec();
  vector<double> state_vec = test.GetStateVec();
  vector<reaction> rxn_vecs = test.GetReactionsVector();
  vector<double> k_vec = test.GetRateConstantVec();
  double volume =  test.GetVolume();
  int n_iter = test.GetNIterations();

  Gillespie simulator(names_vec, state_vec, rxn_vecs, k_vec, volume, n_iter);
  simulator.Start();
}
