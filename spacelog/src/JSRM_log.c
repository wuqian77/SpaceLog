#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
//static double *vector(int n)
//{
//  double *v;
//  v = Calloc(n, double);
//  return v;
//}

//static void free_vector(double *v)
//{
//  Free(v);
//}

/////////   c code of Joint Lasso Regression for Rpackage space.   ////////////

// Joint Sparse Regression model; using the adaptive shooting method;

/// Variables:
/// n: Number of observations;
/// p: Number of variables;
/// lambda1: lasso panelty;
/// lambda2: Elastic Net panelty if lambda2 !=0
/// lambda3: log panelty if lambda3!=0
/// sigma_sr: vector of length p, the ith element is the square root of sigma^i
/// Y_data: vector of length n*p (converted from input data matrix nXp by row);

/// Remark: the columns of Y_data do not need to be standardized beforehand, such that we can play with weight from outside.
/// In the code, each column Y_data will be centered to mean 0  (but not re-scaled)
///   .C("JSRM",
///   as.integer(n),
///   as.integer(p),
///   as.single(lambda1),
///   as.single(lambda2),
///   as.single(lambda3),
///   as.single(Y_data),
///   as.single(sigma_sr),
///   as.integer(n_iter),
///   iter.count=as.integer(iter_count),
///   beta.estimate=as.single(Beta_output)
///   )


void JSRM_log(int* NN, int * PP, float * L1, float * L2, float * L3, float *Y_data, float *sigma_sr, int * N_iter,
              int *Iter_count, float *Beta_output, float *Beta_1, float *Beta_2, float *Beta_3, float *Beta_4)
{
 /////////////////////////  In put data Y.m /////////////////////////////

 int n, p, n_iter;
 float lambda1, lambda2, lambda3;
 float *Y_m;

 float *meanx;
 float *normx;
 float *B_s;
 float *B;
 float *A;
 float *beta_new, *beta_old, *beta_last;
 float *Yhat_m;
 float *E_m;
 //float *A2;
 float *AA;
 int *pick;


 float temp, temp1, threshold;
 float  Aij, Aji;
 float  beta_next;
 //float temp2;
 float beta_change;
 float tmp1, tmp2;
 float maxdif;

 int i,j,k;
 int iter_count;
 int iter, nbeta;
 int nrow_pick, jrow;
 int change_i, change_j;
 int cur_i,cur_j;

 float eps1;       //   tolerant value;


 iter=0;
 eps1 = 1e-6;
 maxdif = -100.0;
 iter_count=0;

 n=*NN;
 p=*PP;
 n_iter=*N_iter;
 nbeta = p*(p-1)/2;

 lambda1=*L1;
 lambda2=*L2;
 lambda3=*L3;

 meanx=(float *) malloc(p*sizeof(float));
 normx=(float *) malloc(p*sizeof(float));
 Y_m=(float *) malloc(n*p*sizeof(float));

 for(i=0; i<n; i++)
   for(j=0; j<p; j++)
     Y_m[i*p+j]=Y_data[i*p+j];


 ////// normalize each column of Y_m into mean=0 and norm=1

    for(j=0;j<p;j++)
        meanx[j] = 0 ;
    for(j=0;j<p;j++)
        for(i=0;i<n;i++)
           meanx[j] = meanx[j]+ Y_m[i*p+j] ;
    for(j=0;j<p;j++)
        meanx[j] = meanx[j]/n  ;
    for(j=0;j<p;j++)
        for(i=0;i<n;i++)
            Y_m[i*p+j] = Y_m[i*p+j] - meanx[j];

    for(j=0;j<p;j++)  normx[j] = 0 ;
    for(j=0;j<p;j++)  for(i=0;i<n;i++) normx[j] = normx[j]+ Y_m[i*p+j]*Y_m[i*p+j] ;

    
 //////////////////////////////     Step 0      //////////////////////////////
 ////////////////////////////// Get initial value (Equation 3)

 B_s=(float *) malloc(p*p*sizeof(float));
 B=(float *) malloc(p*p*sizeof(float));
 A =(float *) malloc(p*p*sizeof(float));
 beta_new=(float *) malloc(p*p*sizeof(float));

  for(i = 0;i<p;i++)
    {
      for(j=0;j<p;j++)
        B[i*p+j] = (sigma_sr[i]/sigma_sr[j]);
   }

  

 for(i = 0;i<p;i++)
     for(j=0;j<(i+1);j++)
      {
      B_s[i*p+j]=B[i*p+j]*B[i*p+j]*normx[i] +  B[j*p+i]*B[j*p+i]*normx[j];
      B_s[j*p+i]=B_s[i*p+j];
      }

  


 for(i = 0;i<p;i++)
     for(j=0;j<(i+1);j++)
      {
       temp=0;
       for(k = 0; k < n ; k++ )
        {
	 temp = temp + Y_m[k*p+i]*Y_m[k*p+j];
	}
       A[i*p+j]=temp*B[j*p+i];
       A[j*p+i]=temp*B[i*p+j];
      }

    

 for(i = 0;i<p;i++) {
   for(j=0;j<(i+1);j++) {
     temp1 =   A[i*p+j] + A[j*p+i]  ;

    // double s=0;
    // double v=1;

     if ( temp1 > 0 )
         temp = temp1 - lambda1;
     else
         temp = - temp1 - lambda1;
     if(temp < 0 )
         beta_new[i*p+j] = 0;
     else
       {
         beta_new[i*p+j] =  temp /(B_s[i*p+j]*(1+lambda2));
         if ( temp1 < 0)
         beta_new[i*p+j] = (-1) * beta_new[i*p+j];
        }

     if (lambda3 <= 0)
         beta_new[j*p+i]=beta_new[i*p+j] ;
     else {

       threshold = lambda1/lambda3 ;
       if (temp1 > 0) temp = temp1 - threshold;
       else if (temp1 <= 0) temp = - temp1 - threshold;


       /// threshold = lambda1/(fabs(beta_old)+lambda3) ; set the initial value beta_old=0

       if (fabs(temp1) <= threshold) beta_new[i*p+j]=0;
       else
       {
         beta_new[i*p+j] = temp/(B_s[i*p+j]);
          if (temp1 < 0)
            beta_new[i*p+j] = (-1) * beta_new[i*p+j];
       }
       beta_new[j*p+i]=beta_new[i*p+j] ;
     }
    }
   }


 for(i = 0;i<p;i++)
     { beta_new[i*p+i] = 0 ; }

 for(i = 0;i<p ;i++)
 {
   for(j=0;j<p;j++)
     Beta_1[i*p+j]=beta_new[i*p+j];
 }

 free(A);
 ///////////////////////// End of Step 0 ///////////////////////////////////

 E_m=(float *) malloc(n*p*sizeof(float));
 beta_old=(float *) malloc(p*p*sizeof(float));
 beta_last=(float *) malloc(p*p*sizeof(float));
 pick=(int *) malloc(nbeta*2*sizeof(int));
 AA =(float *) malloc(p*p*sizeof(float));
/////////////////////// Step 1:   Get Initial E
/////////////////////// Equation (4)-(6)

 Yhat_m=(float *) malloc(n*p*sizeof(float));

 for(k=0; k<n; k++)
   for(j=0; j<p; j++)
     {
       Yhat_m[k*p+j]=0;
       for(i=0; i<p; i++)
         Yhat_m[k*p+j]=Yhat_m[k*p+j]+Y_m[k*p+i]*beta_new[i*p+j]*B[i*p+j];
         E_m[k*p+j]=Y_m[k*p+j]-Yhat_m[k*p+j];
     }

   ///printf("%d",printMatrice(beta_new));
   
 free(Yhat_m);
///////////////////////// Step 2: update one beta  //////////////////////////////


  for(i = 0;i<nbeta;i++)
    {
     for(j=0;j<2;j++)
        pick[i*2+j] = 0 ;
    }

  ///////// Get active Set;
     k = 0;
     for(j = p-1; j>=1; j--)
      {
       for(i = j-1; i>=0; i--)
       {
         if( beta_new[i*p+j]>  eps1  || beta_new[i*p+j] < -eps1 )
  	      {pick[k*2+0] =i; pick[k*2+1] = j;
  	       k = k + 1;
      	  break;
       	 }
       }
       if (k>0)
          break;
      }

    

 if(k>0) //otherwise, converge to 0 at the initial step.
 {

    for(i = 0;i<p;i++)
     for(j=0;j<p;j++)
      {
           beta_old[i*p+j] = beta_new[i*p+j];
       }

  ///////// change one beta
     cur_i = pick[0];
     cur_j = pick[1];

  //////// Equation (7), (8)
     Aij=0;
     Aji=0;
          for(k=0; k<n; k++)
             {
               Aij=Aij+E_m[k*p+cur_j]*Y_m[k*p+cur_i];
               Aji=Aji+E_m[k*p+cur_i]*Y_m[k*p+cur_j];
             }
     Aij=Aij*B[cur_i*p+cur_j];
     Aji=Aji*B[cur_j*p+cur_i];

    
  //////// Equation (10)
     beta_next=(Aij+Aji)/B_s[cur_i*p+cur_j]+beta_old[cur_i*p+cur_j];

    

      ///shrink beta_next
      temp1 = beta_next;
     // double s=0;
     // double v=1;

      if ( temp1 > 0 )
            temp = temp1 - lambda1/B_s[cur_i*p+cur_j];
      else
            temp = - temp1 - lambda1/B_s[cur_i*p+cur_j];
      if(temp < 0 )
            temp = 0;
      else
        {
            temp =  temp /(1+lambda2);
         if(temp1 < 0 )
            temp = (-1) * temp;
        }

          if (lambda3 <= 0) {
          beta_change=beta_old[cur_i*p+cur_j]-temp;

          beta_new[cur_i*p+cur_j] = temp;
          beta_new[cur_j*p+cur_i] = temp;

          change_i=cur_i;
          change_j=cur_j;
          }
          else{
            threshold = lambda1/((fabs(beta_old[cur_i*p+cur_j])+lambda3)*B_s[cur_i*p+cur_j]) ;

            if (temp1 > 0)
              temp = temp1 - threshold;
            else
              temp = - temp1 - threshold;

            if (fabs(temp1) < threshold){
              beta_new[cur_i*p+cur_j]=0;
              beta_new[cur_j*p+cur_i]=0;
              beta_change=beta_old[cur_i*p+cur_j];
              change_i=cur_i;
              change_j=cur_j;
            }
           else {
            beta_new[cur_i*p+cur_j] = temp;
            beta_new[cur_j*p+cur_i] = temp;
            beta_change=beta_old[cur_i*p+cur_j] - temp;
            change_i=cur_i;
            change_j=cur_j;
            if (temp1 < 0) {
              beta_new[cur_i*p+cur_j] = (-1)*temp ;
              beta_new[cur_j*p+cur_i] = (-1)*temp ;
              beta_change=beta_old[cur_i*p+cur_j] + temp;
              change_i=cur_i;
              change_j=cur_j;
            }
          }
        }


        

          for(i = 0;i<p ;i++)
          {
            for(j=0;j<p;j++)
              Beta_2[i*p+j]=beta_new[i*p+j];
          }

          for(i = 0;i<p ;i++)
          {
            for(j=0;j<p;j++)
              Beta_3[i*p+j]=B_s[i*p+j];
          }


 //////////////////////////////////////////////////////////////////////////////////
 ///////////////////////// Step 3: begin to iterate  //////////////////////////////

 for ( iter = 0; iter < n_iter; iter++ ){// iteration for all loop

   printf("\n\n start iteraction : \n\n") ;
   printf("%d", iter);

     for(i = 0;i<nbeta;i++)
     {
      for(j=0;j<2;j++)
        pick[i*2+j] = 0 ;
      }

    /////////// update beta_last
    for(i = 0;i<p;i++)
    {
       for(j=0;j<p;j++)
         beta_last[i*p+j] = beta_new[i*p+j];
    }

   // Get active Set;
   k = 0;
   for(j = p-1; j>=1; j--)
     for(i = j-1; i>=0; i--)
     {
       if( beta_new[i*p+j]>  eps1  || beta_new[i*p+j] < -eps1 )
	       {pick[k*2] =i; pick[k*2+1] = j; k = k + 1; }
     }


   nrow_pick = k;

   printf("\n\n number of non-zero beta : \n\n") ;
   printf("%d", k);


 if(nrow_pick>0)  // otherwise, go to all loop directly.
  {
   /////////////////////// loop for active set ///////////////////
   for(jrow = 0; jrow<=nrow_pick - 1; jrow++)
   {

     ///printf("\n\n loop for active set : \n\n") ;
     ///printf("%f", jrow);

      cur_i = pick[jrow*2];
      cur_j = pick[jrow*2+1];

     // cout << "cur_i:" <<  cur_i << " cur_j: " << cur_j << '\n' ;

     beta_old[change_i*p+change_j]=beta_new[change_i*p+change_j];
     beta_old[change_j*p+change_i]=beta_new[change_j*p+change_i];


     /////////// Update Residue   /////////////////
     /////////// Equation (11)

    for(k=0; k<n; k++)
     {
       E_m[k*p+change_i]=E_m[k*p+change_i]+ Y_m[k*p+change_j]*beta_change*B[change_j*p+change_i];
       E_m[k*p+change_j]=E_m[k*p+change_j]+ Y_m[k*p+change_i]*beta_change*B[change_i*p+change_j];
     }

 ///Yhat_m[k*p+j]=Yhat_m[k*p+j]+Y_m[k*p+i]*beta_new[i*p+j]*B[i*p+j];



    ///////////// update beta
    //////////// Equation (12)
     Aij=0;
     Aji=0;
     for(k=0; k<n; k++)
        {
          Aij=Aij+E_m[k*p+cur_j]*Y_m[k*p+cur_i];
          Aji=Aji+E_m[k*p+cur_i]*Y_m[k*p+cur_j];
        }
     Aij=Aij*B[cur_i*p+cur_j];
     Aji=Aji*B[cur_j*p+cur_i];



     beta_next=(Aij+Aji)/B_s[cur_i*p+cur_j]+beta_old[cur_i*p+cur_j];


    ///shrink beta_next
           temp1 = beta_next;
           if ( beta_next > 0 )
                 temp = temp1 - lambda1/B_s[cur_i*p+cur_j];
           else
                 temp = - temp1 - lambda1/B_s[cur_i*p+cur_j];
           if(temp < 0 )
                 temp = 0;
           else
             {
                 temp =  temp /(1+lambda2);
              if(temp1 < 0 )
                 temp = (-1) * temp;
              }

          if (lambda3 <= 0){
            beta_change=beta_old[cur_i*p+cur_j]-temp;

          beta_new[cur_i*p+cur_j] = temp;
          beta_new[cur_j*p+cur_i] = temp;

          change_i=cur_i;
          change_j=cur_j;
          }
          else {
            threshold = lambda1/((fabs(beta_old[cur_i*p+cur_j])+lambda3)*B_s[cur_i*p+cur_j]) ;

            if (temp1 > 0)
              temp = temp1 - threshold;
            else
              temp = - temp1 - threshold;

            if (fabs(temp1) < threshold){
            beta_new[cur_i*p+cur_j]=0;
            beta_new[cur_j*p+cur_i]=0;
            beta_change=beta_old[cur_i*p+cur_j];
            change_i=cur_i;
            change_j=cur_j;
            }
            else
            {
              beta_new[cur_i*p+cur_j] = temp;
              beta_new[cur_j*p+cur_i] = temp;
              beta_change=beta_old[cur_i*p+cur_j] - temp;
              change_i=cur_i;
              change_j=cur_j;
              if (temp1 < 0){
              beta_new[cur_i*p+cur_j] = (-1)*temp ;
              beta_new[cur_j*p+cur_i] = (-1)*temp ;
              beta_change=beta_old[cur_i*p+cur_j] + temp;
              change_i=cur_i;
              change_j=cur_j;
            }
          }
        }
     iter_count=iter_count+1 ;



 }
   ///////////////////////  End of loop for active set ///////////////////

   //////////////////  converge on the active set ?  //////////////////////
   // Get the max difference between beta_last and beta_new :  maxdif
   maxdif=-100;
   for(i = 0;i<p;i++)
     for(j=0;j<p;j++)
      {
      tmp1 = beta_last[i*p+j]- beta_new[i*p+j] ;
      tmp2 = -tmp1 ;
      temp=tmp1;
      if(tmp2>tmp1)
        temp=tmp2;
      //temp = max(tmp1,tmp2);
      if( temp > maxdif )
         maxdif = temp;
     }


   }// end of if(nrow_pick>0)



   // if convergent on active set, move to full set
   if ( maxdif < 1e-6 || nrow_pick<1)
     {

       for(i = 0;i<p;i++)
	     for(j=0;j<p;j++)
	       beta_last[i*p+j] = beta_new[i*p+j];

       /////////////////////// loop for all set //////////////////////

       for( cur_i = 0; cur_i < p-1 ; cur_i++ )
	    for ( cur_j = cur_i + 1; cur_j < p ; cur_j++ )
	    {
           beta_old[change_i*p+change_j]=beta_new[change_i*p+change_j];
	   beta_old[change_j*p+change_i]=beta_new[change_j*p+change_i];



	   if(beta_change< -eps1 || beta_change> eps1)
	   {

	     ///////////    Update E_m   /////////////////
	     //////////// Equation (11)
	    for(k=0; k<n; k++)
	         {
	           E_m[k*p+change_i]=E_m[k*p+change_i]+ Y_m[k*p+change_j]*beta_change*B[change_j*p+change_i];
	           E_m[k*p+change_j]=E_m[k*p+change_j]+ Y_m[k*p+change_i]*beta_change*B[change_i*p+change_j];
            }

	   }

	 ///////////// update beta
     //////////// Equation (12)
        Aij=0;
	      Aji=0;
	      for(k=0; k<n; k++)
	         {
	           Aij=Aij+E_m[k*p+cur_j]*Y_m[k*p+cur_i];
	           Aji=Aji+E_m[k*p+cur_i]*Y_m[k*p+cur_j];
	         }
	      Aij=Aij*B[cur_i*p+cur_j];
	      Aji=Aji*B[cur_j*p+cur_i];

	      beta_next=(Aij+Aji)/B_s[cur_i*p+cur_j]+beta_old[cur_i*p+cur_j];

	    ///shrink beta_next
	           temp1 = beta_next;
	           if ( beta_next > 0 )
	                 temp = temp1 - lambda1/B_s[cur_i*p+cur_j];
	           else
	                 temp = - temp1 - lambda1/B_s[cur_i*p+cur_j];
	           if(temp < 0 )
	                 temp = 0;
	           else
	             {
	                 temp =  temp /(1+lambda2);
	              if(temp1 < 0 )
	                 temp = (-1) * temp;
	             }

            if (lambda3 <= 0){
             beta_change=beta_old[cur_i*p+cur_j]-temp;

             beta_new[cur_i*p+cur_j] = temp;
             beta_new[cur_j*p+cur_i] = temp;

             change_i=cur_i;
             change_j=cur_j;
           }
           else {
             threshold = lambda1/((fabs(beta_old[cur_i*p+cur_j])+lambda3)*B_s[cur_i*p+cur_j]) ;

             if (temp1 > 0)
               temp = temp1 - threshold;
             else
               temp = - temp1 - threshold;

             if (fabs(temp1) < threshold){
               beta_new[cur_i*p+cur_j]=0;
               beta_new[cur_j*p+cur_i]=0;
               beta_change=beta_old[cur_i*p+cur_j];
               change_i=cur_i;
               change_j=cur_j;
             }
             else
             {
               beta_new[cur_i*p+cur_j] = temp;
               beta_new[cur_j*p+cur_i] = temp;
               beta_change=beta_old[cur_i*p+cur_j] - temp;
               change_i=cur_i;
               change_j=cur_j;
               if (temp1 < 0){
               beta_new[cur_i*p+cur_j] = (-1)*temp ;
               beta_new[cur_j*p+cur_i] = (-1)*temp ;
               beta_change=beta_old[cur_i*p+cur_j] + temp;
               change_i=cur_i;
               change_j=cur_j;
             }
           }
	       }

           iter_count=iter_count+1 ;


	    }
   ///////////////////////  End of loop for all set ///////////////////
       // Get the max difference between beta_last and beta_new :  maxdif
       maxdif = -100.0;
       for(i = 0;i<p;i++)
	     for(j=0;j<p;j++)
	       { tmp1 = beta_last[i*p+j]- beta_new[i*p+j] ;
	 		tmp2 = -tmp1 ;
	 		temp = tmp1;
	 		if(tmp2>tmp1)
	 		  temp=tmp2;

			 if( temp > maxdif)
			    maxdif = temp;
			 }

       if( maxdif <1e-06)
          break ;
     }   // end of if(maxdif<1e-6)

  }// end of iter
}// end of if(k>0)
 /////////////////////   End of Step 1,2,...........///////////////////////////


 ////////////////      Output  Beta to file2        /////////////////

 for(i = 0;i<p ;i++)
 {
   for(j=0;j<p;j++)
     Beta_4[i*p+j]=beta_last[i*p+j];
 }

 for(i = 0;i<p ;i++)
 {
   for(j=0;j<p;j++)
     Beta_output[i*p+j]=beta_new[i*p+j];
 }

 *Iter_count=iter_count;

 free(beta_new);
 free(beta_old);
 free(beta_last);
 free(E_m);
 free(B_s);
 free(B);
 free(Y_m);
 free(pick);
 free(meanx);
 free(normx);
 }  // End of jsrm.shoot();



