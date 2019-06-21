#include <stdlib.h>
#include <math.h>
#include <stdio.h>
/////////   c code for MB Lasso Regression via adaptive shooting   ////////////


// neiborhood selection via sparse Regression for partial correlation estimation;
// using the adaptive shooting method;

/// Variables:
/// n: Number of observations;
/// p: Number of variables;
/// lambda1: lasso panelty;
/// lambda2: Elastic Net panelty;
/// sigma_sr: vector of length p, the ith element is the square root of sigma^i
/// Y_m: matrix of n X p
// Remark: each column of Y_m will be standardized to mean zero and norm 1 in the code


void MBshootOne (int index, int n, int P, int n_iter, float lambda1, float lambda2, float * Y_m, float * Beta, float * Esig);

void MBshoot(int* NN, int * PP, float * L1, float * L2, float *Y_data, float * sigma, int * N_iter, float *Rho_output)
{
  int i,j,k;
  int n_iter;
  int n, p;

  float x;
  float lambda1;
  float lambda2;

  float *Y_m;
  float *Beta;
  float *Esig;
  float temp, temp1, temp2;
  float *meanx;
  float *normx;

  n_iter=*N_iter;
  n=*NN;
  p=*PP;
  lambda1=*L1;
  lambda2=*L2;
  temp=0;
  temp1=0;
  temp2=0;


  meanx=(float *) malloc(p*sizeof(float));
  normx=(float *) malloc(p*sizeof(float));
  Y_m=(float *) malloc(n*p*sizeof(float));

  Beta=(float *) malloc(p*p*sizeof(float));
  Esig=(float *) malloc(p*sizeof(float));

  for(i=0; i<n; i++)
    for(j=0; j<p; j++)
      Y_m[i*p+j]=Y_data[i*p+j];


 ////// normalize each column of Y_m into mean=0 and norm=1
    for(j=0;j<p;j++)  meanx[j] = 0 ;
    for(j=0;j<p;j++)  for(i=0;i<n;i++) meanx[j] = meanx[j]+ Y_m[i*p+j] ;
    for(j=0;j<p;j++)  meanx[j] = meanx[j]/n  ;
    for(j=0;j<p;j++)  for(i=0;i<n;i++) Y_m[i*p+j] = Y_m[i*p+j] - meanx[j];

    for(j=0;j<p;j++)  normx[j] = 0 ;
    for(j=0;j<p;j++)  for(i=0;i<n;i++) normx[j] = normx[j]+ Y_m[i*p+j]*Y_m[i*p+j] ;
    for(j=0;j<p;j++)  normx[j] = sqrt(normx[j]) ;
    for(j=0;j<p;j++)  for(i=0;i<n;i++) Y_m[i*p+j] = Y_m[i*p+j]/ normx[j];


//////////////////////////////     fit model for each column     //////////////////////////////

 for(i=0;i<p;i++)
   MBshootOne(i, n, p, n_iter, lambda1, lambda2, Y_m, Beta, Esig);

for(i=0; i<p; i++)
 for(j=0; j<=i; j++)
   {
     if(j==i)
       Rho_output[i*p+j]=1;
     else
      {
        temp1=Beta[i*p+j];
        temp2=Beta[j*p+i];
        temp=0;
        if(temp1>0 && temp2>0)
          temp=sqrt(temp1*temp2);
        if(temp1<0 && temp2<0)
          temp=-sqrt(temp1*temp2);
        Rho_output[i*p+j]=temp;
        Rho_output[j*p+i]=temp;
      }
   }

 /////// convert Esig into scale correspoding to sqrt(Y)<--- columns have sd=1, can then calculate inverse.
 for(i=0; i<p; i++)
  {
     sigma[i]=1/(Esig[i]*(n-1));
  }

 free(Beta);
 free(Y_m);
 free(Esig);
 free(meanx);
 free(normx);
} //end of MBmain


/////////////////////////////////////////
//////////// function for one regression
/////////////////////////////////////////

void MBshootOne (int index, int n, int P, int n_iter, float lambda1, float lambda2, float * Y_m, float * Beta, float * Esig)
{
   int p;
   int i, j, k;
   int iter;
   int iter_count;
   int *pick;
   int change_i;
   int cur_j, jrow, nrow_pick;

   float *X_Train;
   float *Y_Train;
   float temp, temp1, temp2;
   float *beta_tilda;
   float *beta_new;
   float eps1;       //   tolerant value;
   float *Yhat;
   float *E;
   float beta_next, beta_change;
   float *beta_old;
   float *beta_last;
   float tmp1, tmp2;
   float maxdif;

   iter=0;
   iter_count=0;
   p=P-1;
   eps1 = 1e-6;
   maxdif=-100.0;
   tmp1=0;
   tmp2=1;
   nrow_pick=0;

   pick=(int *) malloc(p*sizeof(int));
   beta_tilda=(float *) malloc(p*sizeof(float));
   beta_new=(float *) malloc(p*sizeof(float));
   beta_old=(float *) malloc(p*sizeof(float));
   beta_last=(float *) malloc(p*sizeof(float));

   Yhat=(float *) malloc(n*sizeof(float));
   E=(float *) malloc(n*sizeof(float));

   X_Train= (float *) malloc(n*p*sizeof(float));
   Y_Train= (float *) malloc(n*sizeof(float));

  /////////////////////////////////////// initial X_train, Y_train
   for(i=0; i<n; i++)
      Y_Train[i]=Y_m[i*P+index];

   for(i=0; i<n; i++)
     for(j=0; j<P; j++)
       {
          if(j<index)
            {
              X_Train[i*p+j]=Y_m[i*P+j];
            }
           if(j>index)
            {
              X_Train[i*p+j-1]=Y_m[i*P+j];
            }
       }

 ////////////////////////////////////////////// begin regression

 for(j=0;j<p;j++)
  {
   beta_tilda[j]=0;
   for(i = 0;i<n;i++)
   {
   beta_tilda[j] = beta_tilda[j] + X_Train[i*p+j] * Y_Train[i] ;
   }
  }

 // Get beta_new and shrink;
   for(j=0;j<p;j++)
   {
     if ( beta_tilda[j] > 0 )
          temp = beta_tilda[j] - lambda1;
     else
          temp = - beta_tilda[j] - lambda1;
     if(temp <=0 )
          beta_new[j] = 0;
     else
      {
          beta_new[j] =  temp /(1+lambda2);
       if ( beta_tilda[j] < 0 )
          beta_new[j] = (-1) * beta_new[j];
       }
  }

///////////////////////// End of Step 0 ///////////////////////////////////
///////////////////////// Get initial E

   for(k=0; k<n; k++)
       {
         Yhat[k]=0;
         for(i=0; i<p; i++)
           Yhat[k]=Yhat[k]+X_Train[k*p+i]*beta_new[i];
         E[k]=Y_Train[k]-Yhat[k];
       }

///////////////////////// update one beta

 for(i=0;i<p;i++)
     beta_old[i]=beta_new[i];

   k = 0;
   for(j=0;j<p;j++)
    {
     if( beta_old[j]>  eps1  || beta_old[j] < -eps1 )
       {
        pick[k]= j;
        k = k + 1;
        if(k>0)
          break;
        }
      if(k>0)
         break;
     }

if(k>0) //otherwise, converge at the first step.
 {
   cur_j=pick[0];

     beta_next = 0;
     for(i = 0;i<n;i++)
         beta_next = beta_next + X_Train[i*p+cur_j] * E[i] ;

     beta_next=beta_next+beta_old[cur_j];

     // Get beta_new and shrink;
     if ( beta_next > 0 )
        temp = beta_next - lambda1;
     else
        temp = - beta_next - lambda1;

     if(temp < 0 )
         beta_new[cur_j] = 0;
     else
      {
         beta_new[cur_j] =  temp /(1+lambda2);
       if(beta_next < 0 )
         beta_new[cur_j] = (-1) * beta_new[cur_j];
       }
     beta_change=beta_old[cur_j]-beta_new[cur_j];
     change_i=cur_j;


 ////////////////////////     Step 1, 2, .....  //////////////////////////////

    for ( iter = 0; iter < n_iter; iter++ ) {


   for(j=0;j<p;j++)
      beta_last[j] = beta_new[j];

   for(j=0;j<p;j++)
      pick[j] = 0 ;

   // Get active Set;
   k = 0;
   for(j=0;j<p;j++)
     if( beta_last[j]>  eps1  || beta_last[j] < -eps1 )
       {
        pick[k]= j;
        k = k + 1;
        }

   nrow_pick = k;



  if(nrow_pick>0) // if all the beta equals 0 here, stop the program
  {
   /////////////////////// loop for active set ///////////////////
   for(jrow =0; jrow< nrow_pick; jrow++) {
     cur_j = pick[jrow];
     ////////////////////    Update   ////////////////////
     beta_old[change_i] = beta_new[change_i];
     for(k=0;k<n;k++)
         E[k] = E[k]+X_Train[k*p+change_i]*beta_change;
     beta_next = 0;
     for(i = 0;i<n;i++)
         beta_next = beta_next + X_Train[i*p+cur_j] * E[i] ;
     beta_next=beta_next+beta_old[cur_j];

     /////////   Get new beta for (cur_j)  //////////
     // Get beta_new and shrink;
     if ( beta_next > 0 )
        temp = beta_next - lambda1;
     else
        temp = - beta_next - lambda1;
     if(temp < 0 )
        beta_new[cur_j] = 0;
     else
      {
        beta_new[cur_j] = temp /(1+lambda2);
       if ( beta_next < 0 )
        beta_new[cur_j] = (-1) * beta_new[cur_j];
       }
     beta_change=beta_old[cur_j]-beta_new[cur_j];
     change_i=cur_j;
     iter_count=iter_count+1 ;

     }   ///////////////////////  End of loop for active set ///////////////////



     //////////////////  converge on the active set ?  //////////////////////
   // Get the max difference between beta_last and beta_new :  maxdif

     maxdif = -100.0; tmp1 = 0; tmp2 = 0;
     for(j=0;j<p;j++)
       { tmp1 = beta_last[j]- beta_new[j] ;
         tmp2 = -tmp1 ;
         temp=tmp1;
         if(tmp2>tmp1)
            temp = tmp2;
         if( temp > maxdif ) maxdif = temp;
       }
   } // end of if(nrow_pick>0)


   if ( maxdif < 1e-6 || nrow_pick<1)
     {

      for(j=0;j<p;j++)
         beta_last[j] = beta_new[j];

     /////////////////////// loop for all set //////////////////////
     for ( cur_j = 0; cur_j < p ; cur_j++ )
      {
        beta_old[change_i] = beta_new[change_i];

        ////////////////////    Update   ////////////////////

        for(k=0;k<n;k++)
            E[k] = E[k]+X_Train[k*p+change_i]*beta_change;

        beta_next = 0;
        for(i = 0;i<n;i++)
            beta_next = beta_next + X_Train[i*p+cur_j] * E[i] ;

        beta_next=beta_next+beta_old[cur_j];

        /////////   Get new beta for (cur_j)  //////////

        // Get beta_new and shrink;
        if ( beta_next > 0 )
           temp = beta_next - lambda1;
        else
           temp = - beta_next - lambda1;
        if(temp < 0 )
           beta_new[cur_j] = 0;
        else
         {
           beta_new[cur_j] = temp /(1+lambda2);
          if ( beta_next < 0 )
           beta_new[cur_j] = (-1) * beta_new[cur_j];
          }
        beta_change=beta_old[cur_j]-beta_new[cur_j];
        change_i=cur_j;

        iter_count=iter_count+1 ;

     }///////////////////////  End of loop for all set ///////////////////

     maxdif = -100.0; tmp1 = 0; tmp2 = 0;

     // Get the max difference between beta_last and beta_new :  maxdif
     for(j=0;j<p;j++)   { tmp1 = beta_last[j]- beta_new[j] ;
     tmp2 = -tmp1 ;
     temp = tmp1;
     if(tmp2>tmp1)
       temp=tmp2;
     if( temp > maxdif ) maxdif = temp;
     }

     if( maxdif <1e-6)   break ; // jump out the iteration of iter

   } // end of if(maxdif<1e-6)


  }// end of iteration of iter

  for(k=0;k<n;k++)
               E[k] = E[k]+X_Train[k*p+change_i]*beta_change;

}// end of if(k>0)
 /////////////////////   End of Step 1,2,...........///////////////////////////

 /////////////////     Output Esig ///////////////////
 Esig[index]=0;
 for(i=0; i<n; i++)
   {
     Esig[index]=Esig[index]+E[i]*E[i];
   }
 Esig[index]=Esig[index]/n;

 ////////////////      Output  Beta      /////////////////
 for(i=0; i<P; i++)
   {
      if(i<index)
         Beta[index*P+i]=beta_new[i];
      if(i>index)
         Beta[index*P+i]=beta_new[i-1];
      if(i==index)
         Beta[index*P+i]=0;
   }


 free(X_Train);
 free(Y_Train);
 free(beta_new);
 free(beta_tilda);
 free(beta_old);
 free(beta_last);
 free(E);

 }  // End of MBshooting();

