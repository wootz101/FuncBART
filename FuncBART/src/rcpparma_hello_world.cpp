// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"
#include "rtnorm.h"




// [[Rcpp::export]]
Rcpp::List functionalBART(Rcpp::List& x,
                     Rcpp::NumericVector y,
                     size_t p,
                     size_t n,
                     size_t funcP,
                     size_t intvls,
                     Rcpp::List& x_test,
                     size_t n_test,
                     const Rcpp::List& iXinfo,
                     Rcpp::IntegerVector& nc,
                     size_t nkeeptrain,
                     double binaryOffset,
                     int iter,
                     int burn,      //BART Priors
                     double alpha,
                     double beta,
                     double tau,    //All the DART stuff after
                     bool dart,
                     double theta,
                     double omega,
                     double a,
                     double b,
                     double rho,
                     bool aug
){

  Rcpp::List ret;

  //tree input x[i]
  ret["var1"] = x.size();

  std::vector<double*> ix;
  double *t_x;

  std::vector<double*> ix_test;
  double *t_x_test;




#define TRDRAW(a, b) trdraw[a][b]

  Rcpp::NumericVector  yv(y);
  double *iy = &yv[0];

  double* iz = new double[n];

  Rprintf("ix %.4f %.4f" ,iy[0], " ", iy[5], " \n");

  // Priors
  //double alpha{0.95};
  //double beta{2};
  //double tau{0.5};

  //number of cuts
  int *numcut = &nc[0];


  //The total trees desired
  size_t m_total = x.size();
  //Vector of bart objects
  std::vector<bart>  m_bart(m_total);
  //random number generator
  arn gen;

  //init
  for(size_t k=0; k<n; k++) {
    if(iy[k]==0) iz[k]= -rtnorm(0., binaryOffset, 1., gen);
    else iz[k]=rtnorm(0., -binaryOffset, 1., gen);
    printf("%.4f", iz[k] , "\n");
  }

  std::vector<double*> trdraw(nkeeptrain);

  //Output cube of results
  Rcpp::NumericMatrix yhats0(m_total, n);

  Rcpp::NumericMatrix yhats(m_total, n);

  Rcpp::NumericMatrix yhats_test(m_total, n_test);

  //std::vector<double*> tedraw(nkeeptest);

  //for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
  //for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

  arma::cube tot_varprb = arma::zeros<arma::cube>(iter,p, m_total);
  arma::cube tot_varcnt = arma::zeros<arma::cube>(iter,p, m_total);
  //Rcpp::NumericMatrix tot_varprb(m_total, p);
  //Rcpp::NumericMatrix tot_varcnt(m_total, p);


  std::vector<double> varprb (p,0.);
  std::vector<size_t> varcnt (p,0);

  //Probs and Counts for functional predictors
  std::vector<double> varprb_funcP (funcP,0.);
  std::vector<size_t> varcnt_funcP (funcP,0);

  Rcpp::NumericMatrix tot_varprb_funcP(m_total, funcP);
  Rcpp::NumericMatrix tot_varcnt_funcP(m_total, funcP);

  //Probs and Counts for intervals
  std::vector<double> varprb_intvls (intvls,0.);
  std::vector<size_t> varcnt_intvls (intvls,0);

  Rcpp::NumericMatrix tot_varprb_intvls(m_total, intvls);
  Rcpp::NumericMatrix tot_varcnt_intvls(m_total, intvls);





  /////////////////////////////////////////////////////////////////////////////////////
  //MAIN LOOP TO LOOP THROUGH ALL TREES FOR SET UP!

  for(int i=0; i < m_total; i++){
    Rcpp::NumericVector xv(Rcpp::as<Rcpp::NumericVector>(x[i]));
    ix[i] = &xv[0];
    t_x = &xv[0];
    Rprintf("ix %.4f %.4f" ,ix[i][0], " ", ix[i][5], " \n");

    // Setting a signle tree for this bart structure. Ensemble of M bart structures
    m_bart[i].setm(1);


    if(iXinfo.size()>0){
      Rcpp::NumericMatrix Xinfo(Rcpp::as<Rcpp::NumericMatrix>(iXinfo[0]));
      xinfo _xi;
      _xi.resize(p);
      for(size_t i=0;i<p;i++) {
        _xi[i].resize(numcut[i]);
        printf("numcut %i", numcut[i], " ' \n'");
        printf("\n");
        //Rcpp::IntegerVector cutpts(Xinfo[i]);
        for(size_t j=0;j<numcut[i];j++)
          _xi[i][j]=Xinfo(i, j);
      }

      m_bart[i].setxinfo(_xi);
    }

    m_bart[i].setprior(alpha,beta,tau);
    m_bart[i].setdata(p, n, funcP, intvls, t_x, iz, numcut);
    m_bart[i].setdart(a,b,rho,aug,dart,theta,omega); //OPTIONAL RIGHT NOW


    /*
     for(size_t j=0;j<n;j++) yhats0(i,j) =m_bart[i].f(j);
     //tree result before the draw
     Rprintf("__________BEFORE_DRAW__________\n");
     Rprintf("%.4f", m_bart[i].f(1));
     Rprintf("____________________\n");
     //tree result before the draw
     Rprintf("__________BEFORE_DRAW__________\n");
     Rprintf("%.4f", m_bart[i].f(150));
     Rprintf("____________________\n");


     m_bart[i].draw(1., gen);
     //m_bart[i].draw(1., gen);
     //m_bart[i].draw(1., gen);
     //m_bart[i].draw(1., gen);

     for(size_t j=0;j<n;j++) yhats(i,j) =m_bart[i].f(j);

     //tree result AFTER the draw
     Rprintf("__________AFTER_DRAW__________\n");
     Rprintf("%.4f", m_bart[i].f(1));
     Rprintf("____________________\n");
     //tree result AFTER the draw
     Rprintf("__________AFTER_DRAW__________\n");
     Rprintf("%.4f", m_bart[i].f(150));
     Rprintf("____________________\n");

     */
    //Get initial values
    for(size_t j=0;j<n;j++) yhats0(i,j) = m_bart[i].f(j);

    Rprintf("\n");
    Rprintf("\n");

  }

  Rprintf("END INITIALIZATION: \n");


  /////////////////////////////////////////////////////////////
  // CREATE A LOOP FOR UPDATING EACH TREE, ONE AT A TIME

  //
  for(int it=0; it<(iter+burn); it++){


    //Rprintf("START ITER: \n");
    for(int i=0; i < m_total; i++){
      //Rprintf("START M_ITER: \n");


      // Metropolis hastings ratio


      if(dart){
        m_bart[i].startdart();
      }
      //m_bart[i].draw(funcP, intvls, 1., gen);


      if(it>(burn-1)){
        //get probabilites
        varcnt = m_bart[i].getnv();
        varprb = m_bart[i].getpv();

        for(int j =0; j<p; j++){
          tot_varcnt(it-burn,j, i) = varcnt[j];
          tot_varprb(it-burn,j, i) = varprb[j];
        }

      }

      // get Residuals




      //Get final predictions
      if(it ==(iter+burn-1)){
        for(size_t j=0;j<n;j++) yhats(i,j) =m_bart[i].f(j);
        Rprintf("TREE NUMBER: \n");
        Rprintf("%i", i+1);
        m_bart[i].pr();

        //MAKE PREDICTION ON TEST SET

        Rcpp::NumericVector xv_test(Rcpp::as<Rcpp::NumericVector>(x_test[i]));
        ix_test[i] = &xv_test[0];
        t_x_test = &xv_test[0];

        double *fhattest = new double[n_test];
        m_bart[i].predict(p, n_test, t_x_test, fhattest);
        Rprintf(  " \n");
        Rprintf(  " \n");

        Rprintf("TEST VALUES %.4f %.4f" , fhattest[0], " ", fhattest[1], " \n");
        Rprintf(  " \n");
        Rprintf(  " \n");

        for(size_t j=0;j<n_test;j++) yhats_test(i,j)= fhattest[j];

        delete[] fhattest;


        //get probabilites of functional predictors
        varcnt_funcP = m_bart[i].getnv_funcP();
        varprb_funcP = m_bart[i].getpv_funcP();
        for(int j =0; j<funcP; j++){
          tot_varcnt_funcP(i, j) = varcnt_funcP[j];
          tot_varprb_funcP(i, j) = varprb_funcP[j];
        }

        //get probabilites of intervals
        varcnt_intvls = m_bart[i].getnv_intvls();
        varprb_intvls = m_bart[i].getpv_intvls();
        for(int j =0; j<intvls; j++){
          tot_varcnt_intvls(i, j) = varcnt_intvls[j];
          tot_varprb_intvls(i, j) = varprb_intvls[j];
        }




      }
    }
  }


  ret["yhats"] = yhats;
  ret["yhats0"] = yhats0;
  ret["yhats_test"] = yhats_test;

  ret["varcount"] = tot_varcnt;
  ret["probcount"] = tot_varprb;

  //return values for functional predictors and intervals
  ret["vc_funcP"] =tot_varcnt_funcP;
  ret["vp_funcP"] =tot_varprb_funcP;

  ret["vc_intvls"] =tot_varcnt_intvls;
  ret["vp_intvls"] =tot_varprb_intvls;



  //return data structures (using C++)


  //random number generation
  //arn gen(n1, n2);



  return ret;

}






