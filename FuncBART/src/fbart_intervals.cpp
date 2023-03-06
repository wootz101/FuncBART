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
Rcpp::List functionalBART_intervals(Rcpp::List& x,
                              Rcpp::NumericVector y,
                              size_t p,
                              size_t n,
                              size_t funcP,
                              size_t intvls,
                              Rcpp::IntegerVector& interval_sizes,
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
  //ret["var1"] = x.size();

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
  printf("M %i", m_total, " ' \n' ");
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



  Rcpp::List list_interval_prob(m_total);

  //arma::cube tot_varprb = arma::zeros<arma::cube>(iter,p, m_total);
  //arma::cube tot_varcnt = arma::zeros<arma::cube>(iter,p, m_total);

  arma::cube tot_varprb_fP = arma::zeros<arma::cube>(iter,funcP, m_total);
  arma::cube tot_varcnt_fP = arma::zeros<arma::cube>(iter,funcP, m_total);
  std::vector<double> varprb_funcP (funcP,0.);
  std::vector<size_t> varcnt_funcP (funcP,0);

  arma::cube yhatsCube = arma::zeros<arma::cube>(iter,n, m_total);
  arma::cube yhatsCube_test = arma::zeros<arma::cube>(iter,n_test, m_total);






  /////////////////////////////////////////////////////////////////////////////////////
  //MAIN LOOP TO LOOP THROUGH ALL TREES FOR SET UP!
  Rprintf("START INITIALIZATION: \n");

  for(int i=0; i < m_total; i++){
    //size_t p = funcP*interval_sizes[i];

    Rcpp::NumericVector xv(Rcpp::as<Rcpp::NumericVector>(x[i]));
    ix[i] = &xv[0];
    t_x = &xv[0];
    //Rprintf("ix %.4f %.4f" ,ix[i][0], " ", ix[i][5], " \n");

    // Setting a signle tree for this bart structure. Ensemble of M bart structures
    m_bart[i].setm(1);


    //Set up probability list

    list_interval_prob[i] = Rcpp::NumericMatrix (iter, interval_sizes[i]);

    if(iXinfo.size()>0){
      Rcpp::NumericMatrix Xinfo(Rcpp::as<Rcpp::NumericMatrix>(iXinfo[i]));
      xinfo _xi;

      _xi.resize((funcP*interval_sizes[i]));


      for(size_t k=0;k<(funcP*interval_sizes[i]);k++) {

        _xi[k].resize(numcut[k]);

        //printf("numcut %i", numcut[k], " ' \n'");
        //printf("\n");
        //Rcpp::IntegerVector cutpts(Xinfo[i]);
        for(size_t j=0;j<numcut[k];j++)
          _xi[k][j]=Xinfo(k, j);
      }

      m_bart[i].setxinfo(_xi);

    }

   m_bart[i].setprior(alpha,beta,tau);
   m_bart[i].setdata(funcP*interval_sizes[i], n, funcP, interval_sizes[i], t_x, iz, numcut);


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
  //double sig_update{1.0};
  //std::vector<double> sigma_values(iter+burn,0.);

  //std::vector<double> r_vector(n,0.);
  //Rcpp::NumericMatrix r_mat(m_total, n, 0.0);


  // begin ITERATIVE LEARNING
  Rprintf("START MCMC: \n");




  for(int it=0; it<(iter+burn); it++){



    //Friendly output status

    if(it%100==0) printf("done %zu (out of %lu)\n",it,iter+burn);

    //n-length vector adding all tree values
    std::vector<double> all_fhat(n, 0.);
    std::vector<double> r_vector(n,0.);

    //std::vector<double> r_i(n,0.);
    // residual for sigma
    //double rs_tot{0.0};

    for(int i=0; i < n; i++){
      // Get Residuals first
      for(size_t j=0; j<m_total; j++){
        all_fhat[i] += m_bart[j].f(i)/m_total;
      }
    }


    std::vector<size_t> global_fp_nv (p,0);

    for(int i=0; i < m_total; i++){
      // Get global NV count first
      std::vector<size_t> temp_nv (p,0);
      temp_nv = m_bart[i].getnv();
      for(size_t j=0;j<funcP;j++){
        for(size_t k=0; k<interval_sizes[i]; k++){
          global_fp_nv[j] += temp_nv[(j)*interval_sizes[i]+k];
        }
      }
    }



    if(it%100==0) {
      //printf("done %zu (out of %lu)\n",it,iter+burn);
      Rprintf(" GLOBAL: %i", global_fp_nv[0], " \n");
      Rprintf(" GLOBAL: %i", global_fp_nv[1], " \n");
      Rprintf(" GLOBAL: %i", global_fp_nv[2], " \n");
      Rprintf(  "    \n");
    }


    for(int i=0; i < m_total; i++){
      // Get Residuals first
      for(size_t j=0; j<n; j++){
        restemp=(iz[j] - all_fhat[j]);
        r_vector[j]= restemp;
      }


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



      // Metropolis hastings ratio


       // m_bart[i].draw(funcP, interval_sizes[i], 1.0, gen, gen_fp, gen_int, r_vector);
      m_bart[i].draw_global_int(funcP, interval_sizes[i], 1.0, gen, gen_fp, gen_int, r_vector, global_fp_nv);
      //m_bart[i].draw(funcP, intvls, sig_update, gen, gen_fp, gen_int, r_vector);

      if(it>(burn-1)){
        //get probabilites
        std::vector<double> varprb_intvls (interval_sizes[i],0.);
        varprb_intvls = m_bart[i].getpv_intvls();

        varcnt_funcP = m_bart[i].getnv_funcP();
        varprb_funcP = m_bart[i].getpv_funcP();

        //get interval probabilites
        for(int j =0; j<interval_sizes[i]; j++){

          Rcpp::NumericMatrix temp_intvls = list_interval_prob[i];
          temp_intvls(it-burn,j) = varprb_intvls[j];
          list_interval_prob[i] = temp_intvls;
        }



        //get functional Predictor Probability
        for(int j =0; j<funcP; j++){
          tot_varcnt_fP(it-burn,j, i) = varcnt_funcP[j];
          tot_varprb_fP(it-burn,j, i) = varprb_funcP[j];
        }


        //get the Yhat values
        for(size_t j=0;j<n;j++) yhatsCube(it-burn,j, i) =m_bart[i].f(j)/m_total;

        //get the Yhat_test values

        Rcpp::NumericVector xv_test(Rcpp::as<Rcpp::NumericVector>(x_test[i]));
        ix_test[i] = &xv_test[0];
        t_x_test = &xv_test[0];

        double *fhattest = new double[n_test];
        m_bart[i].predict(funcP*interval_sizes[i], n_test, t_x_test, fhattest);

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
        m_bart[i].predict(funcP*interval_sizes[i], n_test, t_x_test, fhattest);


        for(size_t j=0;j<n_test;j++) yhats_test(i,j)= fhattest[j]/m_total;

        delete[] fhattest;



      }
    }

   // for(size_t j=0; j<n; j++){
    //  rs_tot += r_vector[j]*r_vector[j];
    //}
    //Update sigma value.. draw from posterior
   // sig_update = sqrt( (nu*lambda + rs_tot)/gen.chi_square(n+nu) ) ;
    //Rprintf("SIGMAAAAA %.4f" , sig_update, " \n");
    //Rprintf(  " \n");
    //sigma_values[it] = sig_update;


  }



  Rprintf("END MCMC: \n");

  ret["interval_prob"] = list_interval_prob;

  ret["yhats"] = yhats;
  ret["yhats0"] = yhats0;
  ret["yhats_test"] = yhats_test;


  ret["yhatsCube_test"] = yhatsCube_test;
  ret["yhatsCube"] = yhatsCube;

  ret["vc_funcP_Cube"]=tot_varcnt_fP;
  ret["vp_funcP_Cube"]=tot_varprb_fP;

  //return values for functional predictors and intervals

 // ret["sigma_values"] = sigma_values;




  //return data structures (using C++)


  //random number generation
  //arn gen(n1, n2);



  return ret;

}
