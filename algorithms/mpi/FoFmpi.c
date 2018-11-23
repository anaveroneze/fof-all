#include <stdio.h>
#include "mpi.h"
#include <string.h>
#include <stdlib.h>

#define MAX_LOCAL_ORDER 8388760
#define MAX_HOST 120
#define NMax 16777500
#define BUFFSIZE 67208864

//#define MAX_LOCAL_ORDER 200000
//#define NMax 400000
//#define BUFFSIZE 900000

main(int argc, char* argv[])
  {
  float *local_x, *local_y, *local_z, *x1, *y1, *z1, *v1, *v2, *v3;
  double startwtime = 0.0, endwtime;
  int *iden, *local_id, *igru;

  local_x = calloc(MAX_LOCAL_ORDER,sizeof(float));
  local_y = calloc(MAX_LOCAL_ORDER,sizeof(float));
  local_z = calloc(MAX_LOCAL_ORDER,sizeof(float));

  igru = calloc(NMax,sizeof(int));

  if(local_x == NULL || local_y == NULL || local_z == NULL)
    {
    printf("Nao foi possivel alocar memoria.\n");
    exit(1);
    }

  x1 = calloc(NMax,sizeof(float));
  y1 = calloc(NMax,sizeof(float));
  z1 = calloc(NMax,sizeof(float));
  v1 = calloc(NMax,sizeof(float));
  v2 = calloc(NMax,sizeof(float));
  v3 = calloc(NMax,sizeof(float));
  iden = calloc(NMax,sizeof(int));

  if(x1 == NULL || y1 == NULL || z1 == NULL || v1 == NULL || v2 == NULL || v3 == NULL ||iden == NULL){
    printf("Nao foi possivel alocar memoria.\n");
    exit(1);
  }
  float rperc;
  int    N; // ordem do vetor
  int    elem;  /* = n/p */
  int    P; // numero de processadores
  int    my_rank; // identificação do processador

  int aLinhas[MAX_HOST]; // Vetor com o numero de linhas que serão processadas por cada processador
  int aLinhasDisp[MAX_HOST];  // Vetor com a linha inicial que será processada por cada processador
  int linhasprocessador ;// Número de linhas que serão processadas pelo processador
  int numlinhas, linhainicial, resto, processador;
  MPI_Status status;
  void Friends(float rperc, float x[], float y[], float z[], int n1, int my_rank);
  void TrInterface(float rperc, float x[], float y[], float z[], float v1[], float v2[], float v3[], int iden[], int igru[],int qtp[], int LF[], int N, int P, int Tgr);

  if(argc != 2 ){
    puts( "Por favor entre com o nome do arquivo de dados! ");
    getchar();
    exit(1);
  }

  /*******Inicializando a MPI ********/
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int size ;
  char *buff;
  MPI_Buffer_attach( malloc(BUFFSIZE), BUFFSIZE);
  /* a buffer of 10000 bytes can now be used by MPI_Bsend */
  /* Buffer size reduced to zero */
  //MPI_Buffer_attach( buff, size);
  /* Buffer of 10000 bytes available again */

  if (my_rank == 0){
    // 1 - SE É O PROCESSADOR MESTRE LÊ O INPUT DO USUÁRIO. ABRE E LÊ O ARQUIVO
    printf("Raio de Percolação (em Mpc):\n");
    scanf("%f", &rperc);
    startwtime = MPI_Wtime();
    int i;
    float m, tp, p1;
    FILE  *fp;
    fp = fopen(argv[1],"r");
    fscanf (fp, "%d", &N);
    int auxvariable = 1;

    for( i = 0 ; i < N ; i++){
      fscanf (fp, " %f %f %f %f %f %f %f %d", &tp, &x1[i], &y1[i],&z1[i],&v1[i],&v2[i],&v3[i], &auxvariable);
      //printf("%f %f %f %f %f %f %f %d\n", tp, x1[i], y1[i], z1[i], v1[i], v2[i], v3[i], auxvariable);
      iden[i] = (int)(tp);
    }
    fclose(fp);

    // 2 - FAZ O BALANCEAMENTO DE CARGA ENTRE OS PROCESSADORES E ACHA O RESTO
    elem = N/P;
    printf("\nTotal partículas = %d\nPartículas por processador = %d \n",N, elem);
    resto = N%P;
    printf("Resto = %d \n\n", resto);

    linhainicial = 0; // linhainicial = 1;
    for(processador = 0; processador < P; processador++){
      if(processador < resto) // SE AINDA TEM RESTO O PROCESSADOR IRÁ LER UMA LINHA A MAIS
        numlinhas = elem + 1;
      else
        numlinhas = elem; // número de linhas que cada processador vai calcular

      aLinhasDisp[processador] = linhainicial; // Vetor com a linha inicial para cada processador.
      aLinhas[processador] = numlinhas; // Vetor com o número de linhas correspondente a cada processador.
      //ENVIA AO PROCESSADOR O NUMERO DE LINHAS QUE SE DEVE CALCULAR
      MPI_Send(&numlinhas, 1, MPI_INT, processador, 0, MPI_COMM_WORLD);
      //ENVIA AO PROCESSADOR O DESLOCAMENTO NO VETOR DE PARTICULAS
      MPI_Send(&linhainicial, 1, MPI_INT, processador, 0, MPI_COMM_WORLD);

      printf("Linha inicial = %d Numero Linhas = %d \n", linhainicial, numlinhas);
      printf("P = %d Resto = %d \n", processador, resto );
      linhainicial = linhainicial + numlinhas;
    }
  }

  MPI_Bcast(&rperc, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

  linhasprocessador = 0;
  linhainicial = 0;
  MPI_Recv(&linhasprocessador, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
  MPI_Recv(&linhainicial, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

  //**** Divide os vetores x, y e z entre todos os processadores ****/
  MPI_Scatterv(x1, aLinhas, aLinhasDisp, MPI_FLOAT, local_x, linhasprocessador, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(y1, aLinhas, aLinhasDisp, MPI_FLOAT, local_y, linhasprocessador, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(z1, aLinhas, aLinhasDisp, MPI_FLOAT, local_z, linhasprocessador, MPI_FLOAT, 0, MPI_COMM_WORLD) ;
  /*******Cada processador realiza o agrupamento usando o FOF *********/
  Friends(rperc, local_x, local_y,local_z, linhasprocessador, my_rank);

  /*******O processador mestre escreve os resultados *********/
  if (my_rank == 0){
    int t, tt,i;
    int *temp, *temp2, *qtp;
    int ind_proc, incremento = 0;
    temp = calloc(linhasprocessador,sizeof(int));
    temp2 = calloc(N,sizeof(int));
    qtp = calloc(N,sizeof(int));
    if(temp == NULL || temp2 == NULL || qtp == NULL){
      printf("Não foi possível alocar memória.\n");
      exit(1);
    }

    int  temp1 = 0, max[P], LF[P], maxL = 0, Linhas = 0;
    for (tt = 0 ; tt < P; tt++)
      max[tt] = 0;

    int jj = 0;
    int ii = 0;
    int Tgr = 0;
    for( ind_proc = 0 ; ind_proc < P; ind_proc++ ){
      MPI_Recv(&Linhas, 1,MPI_INT, ind_proc, 0, MPI_COMM_WORLD, &status);// recebendo numero de linhas de cada proc.
      MPI_Recv(temp, Linhas, MPI_INT, ind_proc, 0, MPI_COMM_WORLD, &status);// recebendo indice de grupo de cada proc.

      MPI_Recv(&temp1, 1, MPI_INT, ind_proc, 0, MPI_COMM_WORLD, &status);// recebendo numero de grupos encontrados.
      MPI_Recv(temp2, temp1, MPI_INT, ind_proc, 0, MPI_COMM_WORLD, &status);// recebendo qtde de particulas do grupo.

      printf("Numero de grupos = %d linha final = %d my_rank = %d \n", temp1, maxL, ind_proc);
      printf("Partículas por grupo:\n");
      for (i = 1 ; i < temp1; i++){
        ii++;
        qtp[ii] = temp2[i];
        //printf("grupo i = %d part = %d \n", i, temp2[i]);
      }

      temp1 = temp1 -1;
      max[ind_proc] = temp1;
      Tgr += temp1;// soma o total de grupos para fazer a atualização
      maxL += Linhas ;//soma as linhas para saber a linha final de cada proc.
      LF[ind_proc] = maxL;//armazena a linha final de cada proc. para fazer o trat. da interface

      if(ind_proc != 0){
        incremento += max[ind_proc-1] ;
        for (t = 0 ; t < Linhas; t++){
          igru[jj] = temp[t] + incremento; //temp[P] num de particulas por P + incremento em relação ao processador
          jj++;
        }
        printf("TOTAL: igr = %d jj = %d \n\n", igru[jj-1], jj-1);
      }
      else{
         for (t = 0 ; t < Linhas; t++){
           igru[jj] = temp[t]; //se é o primeiro processador ele vai recever o numero de particulas do indice 0
           jj++;
          }
         printf("TOTAL: igr = %d jj = %d \n\n", igru[jj-1], jj-1);
      }
    }
    printf("Quantidade Total de grupos: %d \n", Tgr-1);

    endwtime = MPI_Wtime();
    double tfof;
    tfof = (endwtime - startwtime);
    printf("Tempo de processamento FoF: %lfs \n",tfof );

    double T_init = 0, T_final;
    T_init = MPI_Wtime();
    printf("Fazendo o pós processamento...\n");
    TrInterface(rperc, x1, y1, z1, v1, v2, v3, iden, igru, qtp, LF, N, P, Tgr);
    //IDEN: vetor identificador das partículas
    //IGRU: vetor com número de grupos
    //QTP: vetor com quantidade de particulas por grupo ?
    //LF: vetor com a linha final de cada P
    //N: numero total de particulas
    //P: num processadores
    //TGR: quantidade total de grupos
    T_final = MPI_Wtime();

    double tpos;
    tpos = (T_final - T_init);
    printf ("Tempo de pós-processamento: %lfs \n",tpos );

    free (temp);
    free (temp2);
    free (qtp);
    }

    MPI_Buffer_detach( &buff, &size);
    MPI_Finalize();

/**************************** Libera a memória ******************************/
  free (iden);
  free (x1);
  free (y1);
  free (z1);
  free (v1);
  free (v2);
  free (v3);
  free (local_x);
  free (local_y);
  free (local_z);
  free (igru);

  
//---------------------------------------------------------------------------
}  /* main */

//PÓS PROCESSAMENTO
/*****************************************************************/
void TrInterface(float rperc, float x[], float y[], float z[], float v1[], float v2[], float v3[], int iden[], int igru[],int qtp[], int LF[], int N, int P, int Tgr){
  int *aux, *t1;
  aux = calloc(P,sizeof(int)); //aux tamanho do número de processadores P
  t1 = calloc(N,sizeof(int)); //aux tamanho do número de partículas N
  if(aux == NULL || t1 == NULL){
      printf("Nao foi possivel alocar memoria.\n");
      exit(1);
  }
  float lim_i, lim_s, dist;
  int t, l, kk, i, j, k, jj, aa;
  aux[0] = 0;

  for(i = 1 ; i <= P ; i++){ // definindo a linha final de cada processador.
    aux[i] = LF[i-1]-1;
    printf("P%d aux = %d \n", i, aux[i]);
  }

  printf("Número de processadores: %d\n", P);
  for (kk = 1; kk < P; kk++){//para cada processador
    lim_i = x[aux[kk]] - rperc; // definindo o limite inferior da interface
    lim_s = x[aux[kk]] + rperc; // definindo o limite superior da interface
    t = 0;
    printf("Limite inferior: %f\nLimite superior %f EM RELACAO A X\n\n", lim_i, lim_s);
    printf("Inicio for: %d  Fim for: %d\n", aux[kk-1], aux[kk+1]);//LINHA INICIAL E FINAL DO PROCESSADOR NUM KK
    for (i = aux[kk-1]; i < aux[kk+1]; i++)
      if (lim_i <= x[i] && x[i] <= lim_s){ //ENCONTRA AS PARTÍCULAS PERTENCENTES A CADA CHUNK
        t1[t++] = i;
      //  printf("I: %d\n", i);
      }
      printf("Particulas dentro do limite EM RELACAO A X = %d \nNo total de %d particulas\n\n", t, aux[kk+1] - aux[kk-1]);
      for (j = 0; j < t; j++){
        for (k = (j+1); k < t; k++){
          if(igru[t1[j]] != igru[t1[k]]){
            dist = (x[t1[j]] - x[t1[k]])*(x[t1[j]] - x[t1[k]]) + (y[t1[j]] - y[t1[k]])*(y[t1[j]] - y[t1[k]]) + (z[t1[j]] - z[t1[k]])*(z[t1[j]] - z[t1[k]]);
            if (dist*dist <= 0.1){//PQ 0.1
              for (l = aux[kk-1]; l < aux[kk+1] ; l++){
                if (igru[l] == igru[t1[k]] && l != t1[k] && igru[l] != igru[t1[j]]){
                  printf("MUDOU: igru[%d] para igru[%d]\n", igru[l], igru[t1[j]]);
                  igru[l] = igru[t1[j]];
                }
              }
              printf("igru[%d] = %d igru[%d] = %d \n", t1[k], igru[t1[k]], t1[j], igru[t1[j]]);
              qtp[igru[t1[k]]] = 0; //o grupo deixa de existir
              igru[t1[k]] = igru[t1[j]];// o grupo antigo recebe a nova denominação
              qtp[igru[t1[j]]] = qtp[igru[t1[j]]] + 1;//adiciona 1 part.
            }
          }
        }
      }
     printf("Processador = %d num particulas = %d \n", kk, t);
   }
//  free (t1);
  free (aux);

/********************************************************************************/
//------Escrevendo o arquivo de dados com o índice do grupo de cada partícula-----//
/********************************************************************************/
  char str1[10],str2[10];
  FILE *fp;
  int num = sprintf(str2, "%d", P);
  char resultado;

  strcpy(str1, "GrupoNp");
  strcat(str1,str2);
  fp = fopen(str1,"w");
  //fprintf(fp, "%d  %d \n", N, (Tgr));

  int *t5;
  t5 = calloc((Tgr),sizeof(int));
  int s, t2 = 0;
  printf("Calculando numero de grupos com massa maior q um \n");
  for (s = 1 ; s <= Tgr; s++)
    if(qtp[s] > 1){
      t5[t2] = s;
      t2++;
      //printf ("igru = %d  %d \n",s, qtp[s]);
    }
  //printf ("igru = %d  %d \n",s, qtp[s]);
  printf("Numero de grupos massa maior q um = %d\n", t2);

/********************escrevendo arquivo de saida ************************/
FILE *fp1;
fp1 = fopen("GruposCSel.txt","w");
fprintf(fp1, "%d  %d \n", N, t2);
fprintf(fp, "%d  %d \n", N, Tgr);
int ii;
  for (ii = 0 ; ii < t2 ; ii++)
    fprintf(fp1,"%d \n",t5[ii] );

  for (ii = 0 ; ii < N ; ii++)
  fprintf(fp,"%d %d %d % 10.6e % 10.6e % 10.6e % 10.6e % 10.6e % 10.6e \n", ii,iden[ii],igru[ii],x[ii], y[ii],z[ii], v1[ii], v2[ii],v3[ii]);

  fclose (fp);
  fclose (fp1);
}

//****************************************************************************/
void Friends(float rperc, float x[], float y[], float z[], int n1, int my_rank){
  int *igru;
  int *Ngr;

  igru = calloc((n1),sizeof(int));
  Ngr = calloc((n1),sizeof(int));
    if(igru == NULL || Ngr == NULL)
      {
      printf("Nao foi possivel alocar memoria.\n");
      exit(1);
      }

  float dist;
  int i = 0;
  int k = 0;
  int j, l, mult;
    //printf("Raio de percolação: %f \n", rperc);
  for (i = 0; i < n1; i++)
    {
     //if (my_rank == 1) printf ("k = %d \n", k);
    k++;

    while (igru[i] != 0 && i < (n1-1) )i++;

    igru[i] = k;
    mult = 1;

    for (j = i ; j < n1 ; j++)
      if(igru[j] == k)
        for (l = (i + 1) ; l < n1 ; l++)
          if (igru[l] == 0)
            {
            dist = (x[j] - x[l])*(x[j] - x[l]) + (y[j] - y[l])*(y[j] - y[l]) + (z[j] - z[l])*(z[j] - z[l]);
            if (dist*dist <= rperc)
              {
              igru[l] = k;
              mult =  mult + 1;
              }
            }
        Ngr[k] = mult;
       // printf("igru = %d mult = %d \n", igru[i], Ngr[k]);
        }
   //printf("k = %d Linhas = %d \n", k, n1);
  /********** Enviando os dados para o processador mestre **************/
    //printf(" Enviando o indice \n");
    k = k+1;
    MPI_Send(&n1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Bsend(igru, n1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&k, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Bsend(Ngr, k, MPI_INT, 0, 0, MPI_COMM_WORLD);

    //printf(" Envio concluido \n");

  free (igru);
  //free (Ngr);
  }
