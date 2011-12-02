/*
  proj02.c: project 2 for ARC
  autor: xjorda01@stud.fit.vutbr.cz
   
  Je lepsi projizdeni po radcich, je to rychlejsi kvuli rozlozeni poli v pamenti.
  A jeste se inkrementuje jenom jeden register coz je rychlejssi nez aby jsme zadavali nova hodnota do registru.

 P.S.  Nestihl jsem otestovat 1024x1024 na Merlinu :( Nevim pokud pojede pod minutu :(

*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
/* definition of matrix functions  - uses MATRIX_SIZE */
#include "matrix.h"
#include <omp.h>
#define RADIUS 3 
#define C_TEP_KAPACITA 1000
  matrix_t *even_matrix = NULL;
  matrix_t *odd_matrix  = NULL ;
  matrix_t *temp_matrix = NULL;
  double precision = 0.0;

  /*promene */ 
  int i = 0,j = 0;
  int xx = 0, yy = 0;
  int xx1 = 0, xx2 = 0;
  int yy1 = 0, yy2 = 0;
  double temp_temp = 0.0;
  double temp_diff = 0.0;
  double distance = 0.0;
  double maxDiff = 100;
  int m = 0;
#ifdef _OPENMP
  double * p_max_Diff = NULL;
#endif
  /*********/

void usage(char *prog)
{
  fprintf(stderr, "usage : %s filename precision\n", prog);
  exit(1);
}
int main(int argc, char *argv[]){

  FILE *fp;

  even_matrix = alloc_matrix();
  odd_matrix = alloc_matrix();

  if(argc != 3)
  	usage(argv[0]);
    
  if(!(fp=fopen(argv[1], "r"))){
    perror(argv[1]);
  	exit(1);
  }

  precision = atof(argv[2]);

  printf("precision: %f\n", precision);
  srand(time(NULL));
  

#ifdef _OPENMP
  omp_set_num_threads(omp_get_num_procs());
  p_max_Diff = malloc(sizeof(double)* omp_get_num_procs());
#endif

  read_matrix(fp,even_matrix);
  read_matrix(fp,odd_matrix);

   
  printf("input matrix:\n");
/*  print_matrix(odd_matrix);*/
  

  /* solve start */
  /* zde umistete svuj kod, ostatni casti 
    (mimo hlavicku souboru, definici promenych
     a pripadne hlavickove soubory nemente! */
  /* prvek matice a_ji je (*matrix)[j][i] */
 


while( maxDiff > precision ){

   maxDiff = 0.0;
#ifdef _OPENMP
   for(m=0;m<omp_get_num_procs();m++)
        *(p_max_Diff+m) = maxDiff;
      
#endif
   /*   */
#pragma omp parallel
{
  #pragma omp  for schedule(static) private(i,j,m,xx,xx1,xx2,yy,yy1,yy2,temp_temp,temp_diff,distance)
   for(i= 1; i<MATRIX_SIZE - 1; i++)
    for(j= 1; j<MATRIX_SIZE - 1; j++){
      /* Miss the cell that has value 100 */ 
      if((*even_matrix)[j][i] != 100.0){ 
        
         /* Q */
         temp_temp = 0.0;
         
         /* xx */
         if( j > RADIUS){
           xx1= j  - RADIUS;
         }
         else{
          xx1 = 0;
         }

         if( (MATRIX_SIZE - j) < RADIUS ){
          xx2 = MATRIX_SIZE - 1;
         }
         else{
          xx2 = j + RADIUS;
         }

         /* yy */
         if( i > RADIUS){
           yy1= i - RADIUS;
         }
         else{
          yy1 = 0;
         }

         if( (MATRIX_SIZE - i) < RADIUS ){
           yy2 = MATRIX_SIZE - 1;
         }
         else{
          yy2 = i + RADIUS;
         }
        /* Count Q for this cell*/ 
        for(yy = yy1; yy <= yy2; yy++) 
          for(xx = xx1; xx <= xx2; xx++) {
               /* Count distance */
               distance = sqrt( pow((xx - j),2) + pow((yy - i),2) );

               if( ( distance != 0 ) && ( distance < RADIUS ) )
                   /* Q = k * dT1/dx */
                   temp_temp = temp_temp + (double) ((*even_matrix)[xx][yy] - (*even_matrix)[j][i] ) / distance;
      
         }/* for */

        /* Q = c * sT2/dt */
       temp_diff =  temp_temp / C_TEP_KAPACITA;
       /**/
       (*odd_matrix)[j][i]= (temp_diff ) + (*even_matrix)[j][i];
       /* find  max dT2 */   
//    #pragma omp critical
#ifdef _OPENMP
       if ( fabs( temp_diff ) > *(p_max_Diff + omp_get_thread_num() ))
           *( p_max_Diff + omp_get_thread_num() ) = fabs( temp_diff );
#else 
       if ( fabs( temp_diff ) > maxDiff )
           maxDiff = fabs( temp_diff );
#endif
     }/* if 100 - na misto continue */
    }/* for */  
 }/* Pragma Parallel*/
#ifdef _OPENMP
   for(m=0;m<omp_get_num_procs();m++){
       if( *(p_max_Diff+m)>maxDiff)
          maxDiff = *(p_max_Diff+m);
      }
#endif
    /* Swap the matrix */
    temp_matrix = odd_matrix;
    odd_matrix  = even_matrix;
    even_matrix = temp_matrix;
  }/*while*/ 
  /* solve end */
  
  printf("\noutput matrix:\n");
  print_matrix(even_matrix);

  fclose(fp);

  free_matrix(even_matrix);
  free_matrix(odd_matrix);

#ifdef _OPENMP
  free(p_max_Diff);
#endif   

  return 0;
}






