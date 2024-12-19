/***************************************************************************
Programa que realiza o agrupamento de partículas de acordo com o raio de
percolação fornecido, usando o algoritmo Friends of Friends.
Ultima atualização 12/01/2009
Autor: Renata S. Rocha Ruiz
******************************************************************************
Paralelizado com OpenACC
Por: Ana Solórzano
******************************************************************************/

//---------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <sys/time.h>
//---------------------------------------------------------------------------

int  *igru, *iden, N;
float  *x, *y, *z, *v1, *v2, *v3;

//---------------------------------------------------------------------------
long getTime(void)
{
	struct timeval time;
	gettimeofday(&time, (struct timezone *) NULL);
	return time.tv_sec*1000000 + time.tv_usec;
}

//---------------------------------------------------------------------------
/********************* Lendo o arquivo de dados de entrada *****************/
//---------------------------------------------------------------------------
bool LeDados(char *fn)
  {
  int m;
  float vx, vy, vz;
  FILE  *fp;
  fp = fopen(fn,"rt");

  fscanf (fp, "%d", &N);

//********* alocando memória************//
  iden  = new int [N+1];  // indentificador da partícula
  igru  = new int [N+1];  // índice do grupo
  x  = new float [N+1];  // posição x da partícula;
  y  = new float [N+1];  // posição y da partícula;
  z  = new float [N+1];  // posição z da partícula
  v1 = new float [N+1];  // coordenada x da velocidade
  v2 = new float [N+1];  // coordenada y da velocidade
  v3 = new float [N+1];  // coordenada z da velocidade

  for (int i = 0 ; i < N ; i++)
    {
    //fscanf (fp, "%f  %f %f %f  %f %f %f ",&x[i], &y[i],&z[i],&v1[i],&v2[i],&v3[i],&iden[i]);
    fscanf (fp, "%d %f %f %f %f %f %f %d ", &m , &x[i], &y[i],&z[i],&v1[i],&v2[i],&v3[i],&iden[i]);
//		printf("%d %f %f %f %f %f %f %d", m , x[i], y[i],z[i],v1[i],v2[i],v3[i],iden[i]);  
    igru[i] = 0;
    }

  fclose (fp);
  return false;
  }
 //---------------------------------------------------------------------------
/**************************** Realiza o agrupamento ******************************/
//---------------------------------------------------------------------------

//void Friends(float *x, float *y, float *z, int *iden)
void Friends(float rperc)
  {

  float dist;
  int k = 0;
  int i = 0;
  int j, l;

  #pragma acc data copyin(x[0:N+1], y[0:N+1], z[0:N+1]) create(igru[0:N+1])
  {
  #pragma acc loop reduction(+:k) private(igru)
  for (i = 0; i < N; i++){
      k++;
      while(igru[i] != 0)i++;
    	igru[i] = k;
    	for (j = i ; j < N ; j++){
      	    if(igru[j] == k){
	    for (l = (i + 1) ; l < N ; l++){
                if (igru[l] == 0){
                    dist = (x[j] - x[l])*(x[j] - x[l]) + (y[j] - y[l])*(y[j] - y[l]) + (z[j] - z[l])*(z[j] - z[l]) ;
	            if (dist <= rperc){
                    igru[l] = k;
                    }
               }
            }

           }
       }
     }
  }


/********************escrevendo arquivo de saida ************************/

  char str1[10],str2[10];

  FILE *fp;

  int num;
  char resultado;
  num = sprintf(str2, "%.2lf", rperc);

  strcpy(str1, "Grupos_RP");
  strcat(str1,str2);

  fp = fopen(str1,"w");

  fprintf(fp, "%d %d \n", N, k);

  for (i = 0 ; i < N ; i++)
  fprintf(fp,"%4d % 10d %4d % 10.6e % 10.6e % 10.6e % 10.6e % 10.6e % 10.6e \n", i,iden[i],igru[i],x[i], y[i],z[i],
  v1[i], v2[i],v3[i]);

  fclose (fp);

  printf("Numero total de grupos: ");
  printf("%d \n", k);

  /********* Calculando a quantidade de particulas por grupo ************/
  int nn, si, mult;
  int *Ngr;

  Ngr  = new int [k+1];
    for ( nn = 1; nn <= k; nn++ )
      {
      mult = 0;
      for (si = 0; si < N; si++)
        if(igru[si] == nn) mult =  mult + 1;
        Ngr[nn] = mult;
       }


  int cont1 = 0;
  for (nn = 1 ; nn <= k ; nn++)
    if (Ngr[nn] > 1) cont1++;

    printf("Grupos com massa maior que 1: %d \n", cont1);

  delete Ngr;

  }

//---------------------------------------------------------------------------
/**************************** Limpa a memória ******************************/
//---------------------------------------------------------------------------
void LimpaMemoria(void)
  {
  delete iden;
  delete igru;
  delete x;
  delete y;
  delete z;
  delete v1;
  delete v2;
  delete v3;
  }
//---------------------------------------------------------------------------
/***************************************************************************/
//---------------------------------------------------------------------------

main(int argc, char **argv)
{
  float  local_v[100], rperc ;
  char *Arg1;

  if(argc != 3 )
    {
    puts( "Por favor entre com o nome do arquivo de dados e o raio de percolação! ");
    getchar();
    return(0);
    }

  puts ("Iniciando...");
  Arg1 = argv[1];
  sscanf(argv[2], "%f", &rperc);

  long start = getTime();
  LeDados(Arg1);
	//Raio de percolação ao quadrado, para verificar distância entre particulas
  Friends(rperc*rperc);
  LimpaMemoria();
  long stop = getTime();

	//Tempo em milisec
  printf("Tempo: %ld\n", (long)(stop-start));

  return 0;
}
