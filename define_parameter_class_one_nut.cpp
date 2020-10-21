#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
class parameters{
public:
  // number of species
  int nb_s;
  // number of basal species
  int nb_b;
  // number of nutrients
  int nb_n = 1;
  
  int n_tot = nb_s + nb_n;

  //used in calculating change in nutrient concentration
  // global nutrient turn over rate (rate of replenishment)
  double D = 0.25;
  double v1;
  // half saturation density of nutrient, or nutrient uptake efficiency
  double K;
  double q;
  double c;
  double test_double;
  
  
  // metabolic rates
  NumericVector X = NumericVector(nb_s);
  // consumption rates
  NumericVector w = NumericVector(nb_s);
  // assimilation efficiencies
  NumericVector e = NumericVector(nb_s);
  // plants half saturation density
  NumericVector K1 = NumericVector(nb_b);
  // growth rates of plants
  NumericVector r = NumericVector(nb_b);
  // maximal nutrient level
  NumericVector S = NumericVector(nb_b);
  // body masses:
  NumericVector BM = NumericVector(nb_s);
  // vector of biomasses
  NumericVector bioms = NumericVector(nb_s + nb_n);
  // r*G for plants, as there is no need to compute it at each ODE call
    
  NumericVector test;
  
  // vector of derivatives
  NumericVector dB = NumericVector(nb_n + nb_s);
  // 
  

  NumericMatrix b = NumericMatrix(nb_s,nb_s);
  // handling times
  NumericMatrix h = NumericMatrix(nb_s,nb_s);
  // functional response
  NumericMatrix F = NumericMatrix(nb_s,nb_s - nb_b);
  
  // internal variables for optimisation
  // index of plants (optimisation purpose)
  IntegerVector plants = Range(0, nb_b - 1);
  IntegerVector animals = Range(nb_b, nb_s - 1);
  IntegerVector plants_bioms = Range(nb_n, nb_b - 1 + nb_n);
  IntegerVector animals_bioms = Range( nb_b + nb_n, nb_s - 1 + nb_n);
  IntegerVector non_nut = Range(0, nb_s - 1);
  NumericVector G = NumericVector(nb_b);
  IntegerVector::iterator cons;
  IntegerVector::iterator cons2;
  IntegerVector::iterator res;
  NumericVector uptake = NumericVector(nb_b);
  double out = 0;
  int i = 0;
  
  parameters(int s, int b, int n):
    nb_s(s), nb_b(b), nb_n(n) {}

  void print(){
    Rcout << "nb_s:"  << std::endl << nb_s << std::endl;
    Rcout << "nb_b:"  << std::endl << nb_b << std::endl;
    Rcout << "plants: " << plants << std::endl; 
    Rcout << "bioms: " << bioms << std::endl; 
    Rcout << "bioms plants: " << bioms[plants] << std::endl; 
    Rcout << "G: " << G << std::endl; 
    Rcout << "Gplant: " << G[plants] << std::endl;
    Rcout << "dbplant " << dB[plants] << std::endl;
    Rcout << "r[plants]" << r[plants] << std::endl;
    test_double = pow(0,1.5);
    test = F(_, 0);
    Rcout << "test " << test << std::endl;
  }
  
  // NumericVector F_rate_vec(int pred){
  //   return (w(pred)*b(_,pred)*pow(bioms[Range(nb_n, nb_s + nb_n)],1+q)) / 
  //   (1 + c*bioms(pred) + w(pred)*h(_,pred)*sum(b(_,pred)*pow(bioms[Range(nb_n, nb_s + nb_n -1)],1+q)));
  // }
  
  double F_rate(int prey, int pred, NumericVector bioms){
    double tot = 0;
    int i;
    // double res = 0;
    //     return ((w(pred)*b(prey,pred)*pow(bioms(prey + 1),1+q)) / 
    // ((1 + c*bioms(pred + 1) + w(pred)*h(prey,pred)*sum(b(_,pred)*pow(bioms[Range(1, nb_s + 1-1)],1+q)))*BM(pred)));

    for (i=0; i<nb_s; i++){
      tot += b(i,pred) * pow(bioms[i+1], 1+q);
    }
    return ((w(pred)*b(prey,pred)*pow(bioms(prey + 1),1+q)) / 
    ((1 + c*bioms(pred + 1) + w(pred)*h(prey,pred)*tot)*BM(pred)));
  }
  
  // NumericVector ODE(double t, NumericVector bioms){
  NumericVector ODE(NumericVector bioms, double t){

    bioms[bioms < 0.0000001] = 0.0;
    // Define the matrix of feeding interactions
    for (res = non_nut.begin(); res != non_nut.end(); res++){
      for (cons = animals.begin(); cons != animals.end(); cons++){
        F(*res, *cons - nb_b) = F_rate(*res, *cons, bioms);
      }
    }

    // G = pmin(bioms(0)/(K1 + bioms(0)), bioms(1)/(K2 + bioms(1)));
    // G = bioms(0)/(K1 + bioms(0));
    uptake =  r[plants]*bioms[plants_bioms]*(bioms(0)/(K1 + bioms(0))); 
    // uptake = uptake[plants]*G[plants];

    // derivates for plants
    // Rcout << "plants: ";
    for (res = plants.begin(); res != plants.end(); res++){
      out = 0;
      for (cons = animals.begin(); cons != animals.end(); cons++){
        out += bioms[*cons + 1] * F(*res, *cons - nb_b);
      }
      dB[*res + 1] = uptake[*res] - out - X[*res]*bioms[*res + 1];
    }

    // Rcout << "animals: ";
    // derivative for animals
    for (cons = animals.begin(); cons != animals.end(); cons++){
      out = 0;
      for (cons2 = animals.begin(); cons2 != animals.end(); cons2++){
        out += bioms[*cons2 + 1]*F(*cons,*cons2-nb_b);
      }
      // Rcout << " out done: " << *cons;
      dB[*cons + 1] = sum(e * F(_,*cons-nb_b)) * bioms[*cons + 1] - out - X[*cons]*bioms[*cons + 1];
      // Rcout << " db done";
    }

    // Rcout << "end animals: ";
    dB[bioms < 0.0000001] = 0;
    // derivate for nutrients
    // Rcout << "nuts: ";
    dB(0) = D * (S(0) - bioms(0)) - v1*sum(uptake);
    return dB;
  }

};



RCPP_MODULE(ParamModule){
using namespace Rcpp;
  class_<parameters>("parameters")
  
  .constructor<int, int, int>("constructor")
  
  .method("print", &parameters::print)
  .method("ODE", &parameters::ODE)
  
  .field("nb_s", &parameters::nb_s)
  .field("nb_b", &parameters::nb_b)
  .field("nb_n", &parameters::nb_n)
  .field("BM", &parameters::BM)
  .field("K1", &parameters::K1)
  .field("D", &parameters::D)
  .field("S", &parameters::S)
  .field("r", &parameters::r)
  .field("X", &parameters::X)
  .field("e", &parameters::e)
  .field("w", &parameters::w)
  .field("b", &parameters::b)
  .field("c", &parameters::c)
  .field("h", &parameters::h)
  .field("q", &parameters::q)
  .field("v1", &parameters::v1)
  .field("bioms", &parameters::bioms)
  .field("dB", &parameters::dB)
  .field("D", &parameters::D)
  .field("F", &parameters::F)
  .field("uptake", &parameters::uptake)
  ;  
}
// [[Rcpp::export]]
NumericVector mult_sliced_vecs(NumericVector aa, NumericVector bb, NumericVector cc, int begin, int end){
  // slicing different elements, but same NUMBER of elements
  IntegerVector selection = Range(begin, end);
  IntegerVector selection2 = Range(begin+2, end+2);
  NumericVector y;
  Rcout << "aa[selection]: " << aa[selection] << std::endl;
  Rcout << " bb[selection] " <<  bb[selection] << std::endl;
  y = aa[selection] * bb[selection];
  y = y[selection]*cc[selection2];
  //y[selection] = y[selection]*aa[selection];
  return y;
}
