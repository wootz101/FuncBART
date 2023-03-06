/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#include "bart.h"

//--------------------------------------------------
//constructor
bart::bart():m(200),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false) {}
bart::bart(size_t im):m(im),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false) {}
bart::bart(const bart& ib):m(ib.m),t(m),pi(ib.pi),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false)
{
   this->t = ib.t;
}
bart::~bart()
{
   if(allfit) delete[] allfit;
   if(r) delete[] r;
   if(ftemp) delete[] ftemp;
}

//--------------------------------------------------
//operators
bart& bart::operator=(const bart& rhs)
{
   if(&rhs != this) {

      this->t = rhs.t;
      this->m = t.size();

      this->pi = rhs.pi;

      p=0;n=0;x=0;y=0;
      xi.clear();

      if(allfit) {delete[] allfit; allfit=0;}
      if(r) {delete[] r; r=0;}
      if(ftemp) {delete[] ftemp; ftemp=0;}

   }
   return *this;
}
//--------------------------------------------------
//get,set
void bart::setm(size_t m)
{
   t.resize(m);
   this->m = t.size();

   if(allfit && (xi.size()==p)) predict(p,n,x,allfit);
}

//--------------------------------------------------
void bart::setxinfo(xinfo& _xi)
{
  //printf("setinfo 1 \n");
   size_t p=_xi.size();
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
     size_t nc=_xi[i].size();
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = _xi[i][j];
   }
}
//--------------------------------------------------
/*

void bart::setdataM(size_t p, size_t n, std::vector<double*> x_vector, double *y, size_t numcut, size_t m)
{

  for(int i; i< m; i++){
    double *_t_x = x_vector[i];
    setdata(p,n, _t_x, y, numcut);
  }
  printf("setDataM 1 \n");
}
*/

void bart::setdata(size_t p, size_t n, size_t funcP, size_t intvls, double *x, double *y, size_t numcut)
{
  //printf("setData 1 \n");
  int* nc = new int[p];
  for(size_t i=0; i<p; ++i) nc[i]=numcut;

  this->setdata(p, n, funcP, intvls, x, y, nc);
  delete [] nc;
}


void bart::setdata(size_t p, size_t n, size_t funcP, size_t intvls, double *x, double *y, int *nc)
{
 // printf("setData 2 \n");
   this->p=p; this->n=n; this->x=x; this->y=y;
   if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);

   if(allfit) delete[] allfit;
   allfit = new double[n];
   predict(p,n,x,allfit);

   if(r) delete[] r;
   r = new double[n];

   if(ftemp) delete[] ftemp;
   ftemp = new double[n];

   di.n=n; di.p=p; di.x = &x[0]; di.y=r;

   for(size_t j=0;j<p;j++){
     nv.push_back(0);
     pv.push_back(1/(double)p);
   }

   for(size_t j=0;j<funcP;j++){
     nv_funcP.push_back(0);
     pv_funcP.push_back(1/(double)funcP);
   }

   for(size_t j=0;j<intvls;j++){
     nv_intvls.push_back(0);
     pv_intvls.push_back(1/(double)intvls);
   }
}

//Set data for Fucntional DATA


/* ORIGNIAL CODE
void bart::setdata(size_t p, size_t n, double *x, double *y, size_t numcut)
{
  printf("setData 1 \n");
  int* nc = new int[p];
  for(size_t i=0; i<p; ++i) nc[i]=numcut;

  this->setdata(p, n, x, y, nc);
  delete [] nc;
}


void bart::setdata(size_t p, size_t n, double *x, double *y, int *nc)
{
  printf("setData 2 \n");
  this->p=p; this->n=n; this->x=x; this->y=y;
  if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);

  if(allfit) delete[] allfit;
  allfit = new double[n];
  predict(p,n,x,allfit);

  if(r) delete[] r;
  r = new double[n];

  if(ftemp) delete[] ftemp;
  ftemp = new double[n];

  di.n=n; di.p=p; di.x = &x[0]; di.y=r;

  for(size_t j=0;j<p;j++){
    nv.push_back(0);
    pv.push_back(1/(double)p);
  }
}


 */



/*
void bart::setdataM(size_t p, size_t n, std::vector<double*> x, double *y, int *nc)
{
  printf("setDataM 2 \n");
  this->p=p; this->n=n; this->x=x; this->y=y;
  if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);

  if(allfit) delete[] allfit;
  allfit = new double[n];
  predict(p,n,x,allfit);

  if(r) delete[] r;
  r = new double[n];

  if(ftemp) delete[] ftemp;
  ftemp = new double[n];

  di.n=n; di.p=p; di.x = &x[0]; di.y=r;
  for(size_t j=0;j<p;j++){
    nv.push_back(0);
    pv.push_back(1/(double)p);
  }
}
*/


//--------------------------------------------------
void bart::predict(size_t p, size_t n, double *x, double *fp)
//uses: m,t,xi
{
   double *fptemp = new double[n];

   for(size_t j=0;j<n;j++) fp[j]=0.0;
   for(size_t j=0;j<m;j++) {
      fit(t[j],xi,p,n,x,fptemp);
      for(size_t k=0;k<n;k++) fp[k] += fptemp[k];
   }

   delete[] fptemp;
}
//--------------------------------------------------
void bart::draw(size_t funcP, size_t intvls, double sigma, rn& gen, rn& gen_fp, rn& gen_int, std::vector<double> r_vector)
{

  //double temp_fit{0.0};

   for(size_t j=0;j<m;j++) {
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) {
         allfit[k] = allfit[k]-ftemp[k]; //allfit[k]-ftemp[k]; // allfit[k]-ftemp[k]
         r[k] =   r_vector[k] + ftemp[k] ;  //r_vector[k] + ftemp[k] Keeps the one tree inquestion to be updated  y[k]-allfit[k] //;
         //Rprintf(" DRAW Response: %.4f" , y[k], " \n");
         //Rprintf(" DRAW fitVal: %.4f" , allfit[k], " \n");
         //Rprintf(" DRAW residuals: %.4f" , r[k], " \n");
         //Rprintf( " \n");

      }
      bd(t[j],xi,di,pi,sigma,nv,pv,aug,gen);
      drmu(t[j],xi,di,pi,sigma,gen);
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
   }


   if(dartOn) {
     // ADD Fucntional Selection Portion 'draw_r'

     //DART for FunctionalPredictors FP
     //if(dartOn_fp){
       for(size_t j=0;j<funcP;j++){
         int temp_funcP{0};
         for(size_t k=0; k<intvls; k++){
           temp_funcP += nv[(j)*intvls+k];
         }
         nv_funcP[j] += temp_funcP;
       }

       draw_s(nv_funcP,lpv_funcP,theta_fp,gen_fp);
       draw_theta0(const_theta_fp,theta_fp,lpv_funcP,a_fp,b_fp,rho_fp,gen_fp);

       for(size_t j=0;j<funcP;j++) pv_funcP[j]=::exp(lpv_funcP[j]);

     //}

     //DART for intervals intvls
    // if(dartOn_int){
       for(size_t j=0;j<intvls;j++){
         int temp_intvls{0};

         for(size_t k=0; k<funcP; k++){
           temp_intvls += nv[(j)+k*funcP];
         }
         nv_intvls[j] += temp_intvls;
       }

       draw_s(nv_intvls,lpv_intvls,theta_int,gen_int);
       draw_theta0(const_theta_int,theta_int,lpv_intvls,a_int,b_int,rho_int,gen_int);


       for(size_t j=0;j<intvls;j++) pv_intvls[j]=::exp(lpv_intvls[j]);
// }

     // Combine the probabilities of funcP and intvls
     //draw_s(nv,lpv,theta,gen);
     //draw_theta0(const_theta,theta,lpv,a,b,rho,gen);
     //updates the overall probabilities of each predictor variable


     for(size_t j=0;j<funcP;j++){
       for(size_t k=0; k<intvls; k++){
         pv[(j)*intvls+k] =::exp( lpv_funcP[j] + lpv_intvls[k]); // ::exp( lpv_funcP[j])+::exp( lpv_intvls[k]);

       }
     }


     //for(size_t j=0;j<p;j++) pv[j]=::exp(lpv[j]);
   }

}

/*ORIGINAL VERSION

 void bart::draw(double sigma, rn& gen)
 {
 for(size_t j=0;j<m;j++) {
 fit(t[j],xi,p,n,x,ftemp);
 for(size_t k=0;k<n;k++) {
 allfit[k] = allfit[k]-ftemp[k];
 r[k] = y[k]-allfit[k];
 }
 bd(t[j],xi,di,pi,sigma,nv,pv,aug,gen);
 drmu(t[j],xi,di,pi,sigma,gen);
 fit(t[j],xi,p,n,x,ftemp);
 for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
 }
 if(dartOn) {

 draw_s(nv,lpv,theta,gen);
 draw_theta0(const_theta,theta,lpv,a,b,rho,gen);
 for(size_t j=0;j<p;j++) pv[j]=::exp(lpv[j]);
 }
 }*/



//--------------------------------------------------
//public functions
void bart::pr() //print to screen
{
   cout << "*****bart object:\n";
   cout << "m: " << m << std::endl;
   cout << "t[0]:\n " << t[0] << std::endl;
   //cout << "t[m-1]:\n " << t[m-1] << std::endl;
   cout << "prior and mcmc info:\n";
   pi.pr();
   if(dart){
     cout << "*****dart prior (On):\n";
     cout << "a: " << a << std::endl;
     cout << "b: " << b << std::endl;
     cout << "rho: " << rho << std::endl;
     cout << "augmentation: " << aug << std::endl;
   }
   else cout << "*****dart prior (Off):\n";
   if(p) cout << "data set: n,p: " << n << ", " << p << std::endl;
   else cout << "data not set\n";
}
