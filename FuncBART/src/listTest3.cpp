
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


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//


// and we can use Rcpp::List to return both at the same time
//

//FEBRUARY 20th


// [[Rcpp::export]]
Rcpp::List functionalBART_res(Rcpp::List& x,
                          Rcpp::NumericVector y,
                          size_t p,
                          size_t n,
                          size_t funcP,
                          size_t intvls,
                          Rcpp::List& x_test,
                          size_t n_test,

                          double nu,      /// hyper paramaters to update sigma value
                          double lambda,  /// hyper paramaters to update sigma value

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

                          bool dart_fp,
                          double theta_fp,
                          double omega_fp,
                          double a_fp,
                          double b_fp,
                          double rho_fp,
                          bool aug_fp,

                          bool dart_int,
                          double theta_int,
                          double omega_int,
                          double a_int,
                          double b_int,
                          double rho_int,
                          bool aug_int


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

  //Rprintf("ix %.4f %.4f" ,iy[0], " ", iy[5], " \n");

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
  arn gen_fp;
  arn gen_int;

  //init
  for(size_t k=0; k<n; k++) {
    if(iy[k]==0) iz[k]= -rtnorm(0., binaryOffset, 1., gen);
    else iz[k]=rtnorm(0., -binaryOffset, 1., gen);
    //printf("%.4f", iz[k] , "\n");
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

  arma::cube tot_varprb_fP = arma::zeros<arma::cube>(iter,funcP, m_total);
  arma::cube tot_varcnt_fP = arma::zeros<arma::cube>(iter,funcP, m_total);

  arma::cube tot_varprb_iv = arma::zeros<arma::cube>(iter,intvls, m_total);
  arma::cube tot_varcnt_iv = arma::zeros<arma::cube>(iter,intvls, m_total);

  arma::cube yhatsCube = arma::zeros<arma::cube>(iter,n, m_total);
  arma::cube yhatsCube_test = arma::zeros<arma::cube>(iter,n_test, m_total);
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
  Rprintf("START INITIALIZATION: \n");

  for(int i=0; i < m_total; i++){
    Rcpp::NumericVector xv(Rcpp::as<Rcpp::NumericVector>(x[i]));
    ix[i] = &xv[0];
    t_x = &xv[0];
    //Rprintf("ix %.4f %.4f" ,ix[i][0], " ", ix[i][5], " \n");

    // Setting a signle tree for this bart structure. Ensemble of M bart structures
    m_bart[i].setm(1);


    if(iXinfo.size()>0){
      Rcpp::NumericMatrix Xinfo(Rcpp::as<Rcpp::NumericMatrix>(iXinfo[0]));
      xinfo _xi;
      _xi.resize(p);
      for(size_t k=0;k<p;k++) {
        _xi[k].resize(numcut[k]);
        //printf("numcut %i", numcut[i], " ' \n'");
        //printf("\n");
        //Rcpp::IntegerVector cutpts(Xinfo[i]);
        for(size_t j=0;j<numcut[k];j++)
          _xi[k][j]=Xinfo(k, j);
      }

      m_bart[i].setxinfo(_xi);
    }

    m_bart[i].setprior(alpha,beta,tau);
    m_bart[i].setdata(p, n, funcP, intvls, t_x, iz, numcut);

    m_bart[i].setdart_fp(a_fp,b_fp,rho_fp,aug_fp,dart_fp,theta_fp,omega_fp); //OPTIONAL RIGHT NOW
    m_bart[i].setdart_int(a_int,b_int,rho_int,aug_int,dart_int,theta_int,omega_int);

    //Get initial values
    for(size_t j=0;j<n;j++) yhats0(i,j) = m_bart[i].f(j);
  }

  Rprintf("END INITIALIZATION: \n");


  /////////////////////////////////////////////////////////////
  // CREATE A LOOP FOR UPDATING EACH TREE, ONE AT A TIME

  //Initialize Residuals
  double restemp{0.0};
  double sig_update{1.0};
  std::vector<double> sigma_values(iter+burn,0.);

  //std::vector<double> r_vector(n,0.);
  //Rcpp::NumericMatrix r_mat(m_total, n, 0.0);


// begin ITERATIVE LEARNING
Rprintf("START MCMC: \n");


  for(int it=0; it<(iter+burn); it++){



    //Friendly output status

    if(it%100==0) {
      printf("done %zu (out of %lu)\n",it,iter+burn);
      //Rprintf(" GLOBAL: %i", global_nv[0], " \n");
      }

    //n-length vector adding all tree values
    std::vector<double> all_fhat(n, 0.);
    std::vector<double> r_vector(n,0.);

    //std::vector<double> r_i(n,0.);
    // residual for sigma
    double rs_tot{0.0};

    std::vector<size_t> global_nv (p,0);

    for(int i=0; i < m_total; i++){
      // Get global NV count first
      std::vector<size_t> temp_nv (p,0);
      temp_nv = m_bart[i].getnv();
      for(size_t j=0; j<p; j++){
        global_nv[j] += temp_nv[j];
      }
    }

    if(it%100==0) {
      //printf("done %zu (out of %lu)\n",it,iter+burn);
      Rprintf(" GLOBAL: %i", global_nv[0], " \n");
      Rprintf(" GLOBAL: %i", global_nv[1], " \n");
      Rprintf(" GLOBAL: %i", global_nv[2], " \n");
      Rprintf(  "    \n");
    }


    for(int i=0; i < n; i++){
      // Get Residuals first
      for(size_t j=0; j<m_total; j++){
        all_fhat[i] += m_bart[j].f(i)/m_total;
      }
    }



    //START TREE LOOP
    for(int i=0; i < m_total; i++){
      // Get Residuals first
      for(size_t j=0; j<n; j++){
        restemp=(iz[j] - all_fhat[j]);
        r_vector[j]= restemp; }


      // Get Global Count from all trees first


      //Rprintf(" observation: %i", n-1, " \n");
      //Rprintf(" fitted values: %.4f" , all_fhat[n-1], " \n");

      //Rprintf(" residuals: %.4f" , restemp, " \n");
      //Rprintf(  " \n");

      //add up all the trees R_ti becumes R_i
      //for(int i=0; i < m_total; i++){

      //}


      // Metropolis hastings ratio
      if( (it==(burn/2-1)) &&dart){
        m_bart[i].startdart();
      }

      if( (it==(burn/2-1)) &&dart_fp){
        m_bart[i].startdart_fp();
      }

      if( (it==(burn/2-1))&&dart_int){
        m_bart[i].startdart_int();
      }


      m_bart[i].draw_global(funcP, intvls, 1.0, gen, gen_fp, gen_int, r_vector, global_nv);
      //m_bart[i].draw(funcP, intvls, sig_update, gen, gen_fp, gen_int, r_vector);


      if(it>(burn-1)){
        //get probabilites
        varcnt = m_bart[i].getnv();
        varprb = m_bart[i].getpv();

        varcnt_funcP = m_bart[i].getnv_funcP();
        varprb_funcP = m_bart[i].getpv_funcP();

        varcnt_intvls = m_bart[i].getnv_intvls();
        varprb_intvls = m_bart[i].getpv_intvls();



        for(int j =0; j<p; j++){
          tot_varcnt(it-burn,j, i) = varcnt[j];
          tot_varprb(it-burn,j, i) = varprb[j];

        }

        //get functional Predictor Probability
        for(int j =0; j<funcP; j++){
          tot_varcnt_fP(it-burn,j, i) = varcnt_funcP[j];
          tot_varprb_fP(it-burn,j, i) = varprb_funcP[j];
        }

        //get interval probabilites
        for(int j =0; j<intvls; j++){
          tot_varcnt_iv(it-burn,j, i) = varcnt_intvls[j];
          tot_varprb_iv(it-burn,j, i) = varprb_intvls[j];
        }

        //get the Yhat values
        for(size_t j=0;j<n;j++) yhatsCube(it-burn,j, i) =m_bart[i].f(j)/m_total;

        //get the Yhat_test values
        Rcpp::NumericVector xv_test(Rcpp::as<Rcpp::NumericVector>(x_test[i]));
        ix_test[i] = &xv_test[0];
        t_x_test = &xv_test[0];

        double *fhattest = new double[n_test];
        m_bart[i].predict(p, n_test, t_x_test, fhattest);

        for(size_t j=0;j<n_test;j++) yhatsCube_test(it-burn,j, i)= fhattest[j]/m_total;

        delete[] fhattest;


      }

      // get Residuals






      //Get final predictions
      if(it ==(iter+burn-1)){
        for(size_t j=0;j<n;j++) yhats(i,j) =m_bart[i].f(j)/m_total;
        Rprintf("TREE NUMBER: \n");
        Rprintf("%i", i+1);
        //m_bart[i].pr();

        //MAKE PREDICTION ON TEST SET

        Rcpp::NumericVector xv_test(Rcpp::as<Rcpp::NumericVector>(x_test[i]));
        ix_test[i] = &xv_test[0];
        t_x_test = &xv_test[0];

        double *fhattest = new double[n_test];
        m_bart[i].predict(p, n_test, t_x_test, fhattest);


       // Rprintf("TEST VALUES %.4f %.4f" , fhattest[0], " ", fhattest[1], " \n");


        for(size_t j=0;j<n_test;j++) yhats_test(i,j)= fhattest[j]/m_total;

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

    for(size_t j=0; j<n; j++){
      rs_tot += r_vector[j]*r_vector[j];
    }
    //Update sigma value.. draw from posterior
    sig_update = sqrt( (nu*lambda + rs_tot)/gen.chi_square(n+nu) ) ;
//Rprintf("SIGMAAAAA %.4f" , sig_update, " \n");
    //Rprintf(  " \n");
    sigma_values[it] = sig_update;


  }
  Rprintf("END MCMC: \n");

  ret["yhats"] = yhats;
  ret["yhats0"] = yhats0;
  ret["yhats_test"] = yhats_test;


  ret["yhatsCube_test"] = yhatsCube_test;
  ret["yhatsCube"] = yhatsCube;

  ret["vc_funcP_Cube"]=tot_varcnt_fP;
  ret["vp_funcP_Cube"]=tot_varprb_fP;

  ret["vc_intvls_Cube"]=tot_varcnt_iv;
  ret["vp_intvls_Cube"]=tot_varprb_iv;

  ret["varcount"] = tot_varcnt;
  ret["probcount"] = tot_varprb;

  //return values for functional predictors and intervals
  ret["vc_funcP"] =tot_varcnt_funcP;
  ret["vp_funcP"] =tot_varprb_funcP;

  ret["vc_intvls"] =tot_varcnt_intvls;
  ret["vp_intvls"] =tot_varprb_intvls;


  ret["sigma_values"] = sigma_values;




  //return data structures (using C++)


  //random number generation
  //arn gen(n1, n2);



  return ret;

}






