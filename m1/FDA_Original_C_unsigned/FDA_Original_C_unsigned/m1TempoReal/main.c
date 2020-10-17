#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pgmfiles.h"
#include "diff2d.h"
#include <time.h>
#define minhaFreq 2700000000
//gcc -o fda pgmtolist.c pgmfiles.c diff2d.c main.c -lm
void calculaVetor(float ht, float lambda, long alfa, long beta, float vet[]){
    vet[0] = 0;
    float qN = 0.0;
    for (int i=1; i<=alfa; i++){
        vet[i] = (1.0 - exp(-8.0 * ht * dco(i, 0, lambda))) / 8.0;
    }
}


void main (int argc, char **argv) {
  char   row[80];
  float  **matrix;
  int i, j;
  FILE   *inimage, *outimage;
  long   imax;
  float  lambda;
  int result;
  eightBitPGMImage *PGMImage;

  double primeiroTempo;
  double segundoTempo;
  double terceiroTempo;
  double quartoTempo;
  /* ---- read image name  ---- */
  PGMImage = (eightBitPGMImage *) malloc(sizeof(eightBitPGMImage));

  if (!argv[1])
  {
	  printf("name of input PGM image file (with extender): ");
    scanf("%s", PGMImage->fileName);
  }
  else
  {
    strcpy(PGMImage->fileName, argv[1]);
  }

  result = read8bitPGM(PGMImage);

  if(result < 0)
    {
      printPGMFileError(result);
      exit(result);
    }
  float vet[PGMImage->x];
  // calculando os valores na funcao
  /* ---- allocate storage for matrix ---- */
  unsigned __int64 i_anterior, i_posterior;
  unsigned __int64 i_primeiro, i_segundo;
  unsigned __int64 i_terceiro;


  //i_anterior = __rdtsc();
  double tempoAnterior = clock();
  matrix = (float **) malloc (PGMImage->x * sizeof(float *));
  if (matrix == NULL)
    {
      printf("not enough storage available\n");
      exit(1);
    }
  for (i=0; i<PGMImage->x; i++)
    {
      matrix[i] = (float *) malloc (PGMImage->y * sizeof(float));
      if (matrix[i] == NULL)
        {
	  printf("not enough storage available\n");
	  exit(1);
        }
    }
  double tempoPosterior = clock();
  //i_posterior = __rdtsc();
  primeiroTempo = (tempoPosterior - tempoAnterior) / minhaFreq;
  //i_primeiro = i_posterior - i_anterior;
  //printf("\nCycle in function allocate storage for matrix: %I64d \n", i_primeiro);
  printf("\nTime in function allocate storage for matrix: %f \n", primeiroTempo);

  /* ---- read image data into matrix ---- */
  //i_anterior = __rdtsc();
  tempoAnterior = clock();

  for (i=0; i<PGMImage->x; i++)
    for (j=0; j<PGMImage->y; j++)
      matrix[i][j] = (float) *(PGMImage->imageData + (i*PGMImage->y) + j);

  tempoPosterior = (double)clock();
  //i_posterior = __rdtsc();
  segundoTempo = (tempoPosterior - tempoAnterior) / minhaFreq;
  printf("\nTime in function read image data into matrix: %f \n", segundoTempo);
  //i_segundo = i_posterior - i_anterior;
  //printf("\nCycle in function read image data into matrix: %I64d \n", i_segundo);


  /* ---- process image ---- */
  printf("contrast paramter lambda (>0) : ");
  //~ gets(row);  sscanf(row, "%f", &lambda);
  scanf("%f", &lambda);
  printf("number of iterations: ");
  //~ gets(row);  sscanf(row, "%ld", &imax);
  scanf("%ld", &imax);


  // calculando a formula
  calculaVetor(0.5, lambda, PGMImage->x, PGMImage->y, vet);


  //i_anterior = __rdtsc();
  tempoAnterior = (double)clock();

  for (i=1; i<=imax; i++)
    {
      printf("iteration number: %3ld \n", i);
      diff2d (0.5, lambda, PGMImage->x, PGMImage->y, matrix);
    }

  tempoPosterior = (double)clock();
  //i_posterior = __rdtsc();
  terceiroTempo = (tempoPosterior - tempoAnterior) / minhaFreq;
  printf("\nTime in for the diffLUT: %.30lf \n", terceiroTempo);
  //i_terceiro = i_posterior - i_anterior;
  //printf("\nCycle in for the diffLUT: %I64d \n", i_terceiro);


  /* copy the Result Image to PGM Image/File structure */
  for (i=0; i<PGMImage->x; i++)
    for (j=0; j<PGMImage->y; j++)
      *(PGMImage->imageData + i*PGMImage->y + j) = (char) matrix[i][j];

  /* ---- write image ---- */
  if (!argv[2])
  {
    printf("name of output PGM image file (with extender): ");
    scanf("%s", PGMImage->fileName);
  }
  else
  {
    strcpy(PGMImage->fileName, argv[2]);
  }

  write8bitPGM(PGMImage);
  /* ---- disallocate storage ---- */

  for (i=0; i<PGMImage->x; i++)
    free(matrix[i]);
  free(matrix);

  free(PGMImage->imageData);
  free(PGMImage);
}


