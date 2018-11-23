/***************************************************************************
Programa que realiza o agrupamento de partículas de acordo com o raio de
percolação fornecido, usando o algoritmo Friends of Friends.
Ultima atualização 12/01/2009
Autor: Renata S. Rocha Ruiz
 ******************************************************************************/

using namespace std;
//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#include "No.h"
#include "grupo.h"
#include "Corpo.h"
#include "common.h"
#include "Segmento.h"
#include "Tupla.h"

#include <vector>
#include <string>


int *igru, *iden, N;
float *x, *y, *z, *v1, *v2, *v3;

float max_x, min_x, max_y, min_y, max_z, min_z;

timeval tStart, tEnd;
timeval tempo_inicio, tempo_fim;

double ti = 0, tf = 0, tempo = 0;

double tempo_relabel = 0;
double tempo_verifica_grupo = 0;
double tempo_d = 0;

float raio = 1.7; //1,7 p testes controlados

int No::n_grupos = 1;

vector<int> DivEntrada[8];

vector<Tupla*> t1, t2;
vector<int> r1, r2;

vector<grupo*> grupos_total;



//---------------------------------------------------------------------------
/********************* Lendo o arquivo de dados de entrada *****************/
//---------------------------------------------------------------------------

bool LeDados(char *fn) {
    int m;
    float vx, vy, vz;
    FILE *fp;
    fp = fopen(fn, "rt");

    fscanf(fp, "%d", &N);

    //********* alocando memória************//
    iden = new int [N + 1]; // indentificador da partícula
    igru = new int [N + 1]; // índice do grupo
    x = new float [N + 1]; // posição x da partícula;
    y = new float [N + 1]; // posição y da partícula;
    z = new float [N + 1]; // posição z da partícula
    v1 = new float [N + 1]; // coordenada x da velocidade
    v2 = new float [N + 1]; // coordenada y da velocidade
    v3 = new float [N + 1]; // coordenada z da velocidade
    int valueaux = 0;
    fscanf(fp, "%d %f %f %f %f %f %f %d ", &m, &x[0], &y[0], &z[0], &v1[0], &v2[0], &v3[0], &valueaux); //O: Valores de velocidade nao utilizados.
    igru[0] = 0;
 iden[0] = 0;

    //set valores minimos e máximos
    max_x = min_x = x[0];
    max_y = min_y = y[0];
    max_z = min_z = z[0];

    for (int i = 1; i < N; i++) {
        fscanf(fp, "%d %f %f %f %f %f %f %d", &m, &x[i], &y[i], &z[i], &v1[i], &v2[i], &v3[i], &valueaux);
	iden[i] = i;
        igru[i] = 0;
        if (x[i] > max_x)
            max_x = x[i];
        if (x[i] < min_x)
            min_x = x[i];
        if (y[i] > max_y)
            max_y = y[i];
        if (y[i] < min_y)
            min_y = y[i];
        if (z[i] > max_z)
            max_z = z[i];
        if (z[i] < min_z)
            min_z = z[i];
    }
    printf("\n Maximos e Minimos:\n\t x_max : %f \t y_max: %f \t z_max: %f \n\t x_min : %f \t y_min: %f \t z_min: %f \n", max_x, max_y, max_z, min_x, min_y, min_z);
    fclose(fp);
    return false;
}
//---------------------------------------------------------------------------
/**************************** Realiza o agrupamento ******************************/
//---------------------------------------------------------------------------

//void Friends(float *x, float *y, float *z, int *iden)

void Friends() {
    printf("\n Tempo n2 original\n");
    float dist, rperc;
    int i = 0;
    int k = 0;
    int j, l;

    //        printf("Raio de Percolação (em Mpc):");
    //        scanf("%f", &rperc);
    rperc = raio;
    printf("\nRaio = %f\n", rperc);

    //ORIGINAL! dont touch it!
    for (i = 0; i < N; i++) {
        k++; //set k pra 1, 0 é usado para quando a particula nao possui grupo
        while (igru[i] != 0)i++; // já tem designação, então passa para o próximo!
        igru[i] = k; //set identificar do grupo como indice k

        for (j = i; j < N; j++)//verifica todos os indices do vetor seguintes ao i.
        {
            if (igru[j] == k) //verifica se o vetor de identificação por grupos ja não poussui um valor previamente calculado
            {
                for (l = (i + 1); l < N; l++)//a partir do indice i+1 percorre o restante do vetor com o indice l
                {
                    if (igru[l] == 0)//caso o igrup ainda não tenha sido definido E a distancia entre as particulas j e l seja menor que rperc
                    {
                        dist = sqrt((x[j] - x[l])*(x[j] - x[l]) + (y[j] - y[l])*(y[j] - y[l]) + (z[j] - z[l])*(z[j] - z[l]));

                        //                        printf("\ndistancia : %d e %d é %f",iden[j],iden[l],dist);
                        if (dist <= rperc) {
                            igru[l] = k; //set o identificador do grupo como k.
                        }
                    }
                }
            }
        }
        //O:i=l;?
    }


    //    for (i = 0; i < N; i++) {
    //        k++; //set k pra 1, 0 é usado para quando a particula nao possui grupo
    //        while (igru[i] != 0)i++; // já tem designação, então passa para o próximo!
    //        igru[i] = k; //set identificar do grupo como indice k
    //
    //        for (j = i; j < N; j++)//verifica todos os indices do vetor seguintes ao i.
    //        {
    //            if (igru[j] == k) //verifica se o vetor de identificação por grupos ja não poussui um valor previamente calculado
    //            {
    //                for (l = (i + 1); l < N; l++)//a partir do indice i+1 percorre o restante do vetor com o indice l
    //                {
    //                    if (igru[l] == 0)//caso o igrup ainda não tenha sido definido E a distancia entre as particulas j e l seja menor que rperc
    //                    {
    //                        dist = sqrt((x[j] - x[l])*(x[j] - x[l]) + (y[j] - y[l])*(y[j] - y[l]) + (z[j] - z[l])*(z[j] - z[l]));
    //
    //                        //                        printf("\ndistancia : %d e %d é %f",iden[j],iden[l],dist);
    //                        if (dist <= rperc) {
    //                            igru[l] = k; //set o identificador do grupo como k.
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //        //O:i=l;?
    //    }


    //        for (i = 0; i < N; i++) {
    //            printf("\n no %d é grupo %d", iden[i], igru[i]);
    //        }

    //PRINT DE GRUPOS E SEUS ELEMENTOS (QUANDO O GRUPO TEM MAIS DE 1 CORPO)
    //    int contadornovo;
    //    for (int i = 0; i < N; i++) {
    //        //            printf("\nID: %d grupo %d", iden[i], igru[i]);
    //        int ngrupo = igru[i];
    //        contadornovo = 0;
    //        for (int j = i + 1; j < N; j++) {
    //            if (ngrupo == igru[j]) {
    //                contadornovo++;
    //                printf("\tcorpo %d pertence ao grupo %d\n", iden[j], ngrupo);
    //            }
    //        }
    //        if (contadornovo > 0) {
    //            printf("\tcorpo %d pertence ao grupo %d\n", iden[i], ngrupo);
    //            printf("grupo %d tem %d corpos\n", ngrupo, contadornovo + 1);
    //        }
    //    }

    /********************escrevendo arquivo de saida ************************/

    //    char str1[10], str2[10];
    //
    //    FILE *fp;
    //
    //    int num;
    //    char resultado;
    //    num = sprintf(str2, "%.2lf", rperc);
    //
    //    strcpy(str1, "Grupos_RP");
    //    strcat(str1, str2);
    //
    //    fp = fopen(str1, "w");
    //
    //    fprintf(fp, "%d %d \n", N, k);
    //
    //    for (i = 0; i < N; i++)
    //        fprintf(fp, "%4d % 10d %4d % 10.6e % 10.6e % 10.6e % 10.6e % 10.6e % 10.6e \n", i, iden[i], igru[i], x[i], y[i], z[i],
    //            v1[i], v2[i], v3[i]);
    //
    //    fclose(fp);
    //
    //    printf("Numero total de grupos: ");
    //    printf("%d \n", k);
    //
    //    printf("Numero de nos: %d\n", N);

    /********* Calculando a quantidade de particulas por grupo ************/
    int nn, si, mult;
    int *Ngr;

    Ngr = new int [k + 1];
    for (nn = 1; nn <= k; nn++) {
        mult = 0;
        for (si = 0; si < N; si++)
            if (igru[si] == nn) mult = mult + 1;
        Ngr[nn] = mult;
    }


    int cont1 = 0;
    for (nn = 1; nn <= k; nn++)
        if (Ngr[nn] > 1) cont1++;

    printf("Grupos com massa maior que 1: %d \n", cont1);

    //    printf("Tempo de execução: %ld segundos e %ld microsegundos\n", tEnd.tv_sec - tStart.tv_sec, usec);

    delete Ngr;

}

//---------------------------------------------------------------------------
/**************************** Limpa a memória ******************************/
//---------------------------------------------------------------------------

void LimpaMemoria(void) {
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

//-------------------------CRIA ARVORE------------------------------------

void CriaArvore(Segmento *s, int quadro) {
    gettimeofday(&tempo_inicio, NULL);
    No::n_grupos = 0;
    printf("\n----------------------------------");
    printf("\nnlogn:");
    printf("\nRaio utilizado: %.2f", raio);



    //lista dos grupos
    s->l_grupos = new vector<grupo*>;

    //primeiro grupo    
    grupo *g;
    g = new grupo(1);
    //primeiro corpo
    printf("\ttamanho da entrada : %d", DivEntrada[quadro].size());
    if (DivEntrada[quadro].size() <= 0) {
        printf("\n\n ENTRADA ZERADA\n\n");
        return;
    }

    max_x = x[DivEntrada[quadro].at(0)];
    min_x = x[DivEntrada[quadro].at(0)];
    max_y = y[DivEntrada[quadro].at(0)];
    min_y = y[DivEntrada[quadro].at(0)];
    max_z = z[DivEntrada[quadro].at(0)];
    min_z = z[DivEntrada[quadro].at(0)];
    for (int i = 1; i < DivEntrada[quadro].size(); i++) {
        int atual = DivEntrada[quadro].at(i);
        if (x[i] > max_x)
            max_x = x[atual];
        if (x[i] < min_x)
            min_x = x[atual];
        if (y[i] > max_y)
            max_y = y[atual];
        if (y[i] < min_y)
            min_y = y[atual];
        if (z[i] > max_z)
            max_z = z[atual];
        if (z[i] < min_z)
            min_z = z[atual];
    }
    printf("\nLIMITES Segmento: ");
    printf("\n\t MAX \t\t MIN \nX:\t %f\t %f\nY:\t %f\t %f\nZ:\t %f\t %f\n", max_x, min_x, max_y, min_y, max_z, min_z);


    int frst = DivEntrada[quadro].at(0);

    Corpo *primeiro;
    primeiro = new Corpo(iden[frst], x[frst], y[frst], z[frst]);
    primeiro->igru = 1;

    s->l_grupos->push_back(g);
    s->l_grupos->at(0)->add_no(primeiro);
    No::n_grupos += 1;
    //inicia arvore
    s->arvore = new No(primeiro, max_x, min_x, max_y, min_y, max_z, min_z);
    //printf("\nprimeiro no: id %d\n", s->arvore->c->id);
    s->arvore->populado = true;
    if (primeiro->x < s->arvore->f_xmin || primeiro->x > s->arvore->f_xmax ||
            primeiro->y < s->arvore->f_ymin || primeiro->y > s->arvore->f_ymax ||
            primeiro->z < s->arvore->f_zmin || primeiro->z > s->arvore->f_zmax) {
        s->arvore->fronteira.push_back(primeiro);
    }

    float a, b, c;
    //adiciona demais corpos
    unsigned long int i;

    for (i = 1; i < DivEntrada[quadro].size(); i++) {
        frst = DivEntrada[quadro].at(i);
        //printf("\ni:%d , index: %d e  id no : %d", i, frst, iden[frst]);

        Corpo *seguinte = new Corpo(iden[frst], x[frst], y[frst], z[frst]);
        s->arvore->add(s->l_grupos, seguinte, max_x, min_x, max_y, min_y, max_z, min_z);
    }

    /*
        if (s->arvore->fronteira.size() != 0) {
            printf("\n size fronteira: %d", s->arvore->fronteira.size());
            for (int i = 0; i < s->arvore->fronteira.size(); i++) {
                printf("\n\t id:%d", s->arvore->fronteira.at(i)->id);
            }
        }
     */
    gettimeofday(&tempo_fim, NULL);


    tf = (double) tempo_fim.tv_usec + ((double) tempo_fim.tv_sec * (1000000.0));
    ti = (double) tempo_inicio.tv_usec + ((double) tempo_inicio.tv_sec * (1000000.0));
    tempo = (tf - ti) / 1000;

    //printf("\n Numero de grupos : %d\n", No::n_grupos);

    //    printf("\ntestando acesso a nó %d\n", s->l_grupos->at(0)->lista->at(0)->id);
    //    int gp = 0;
    //    for (int x = 0; x < s->l_grupos->size(); x++) {
    //        //printf("\nG:%d size : %d", x, s->l_grupos->at(x)->lista->size());
    //        if (s->l_grupos->at(x)->lista->size() > 1) {
    //            gp++;
    //        }
    //    }
    //
    //    printf("\n Grupos com massa maior que 1 : %d", gp);
}

void DivideEntrada() {
    printf("Dividindo entrada:");
    float meioX = min_x + ((max_x - min_x) / 2);
    float meioY = min_y + ((max_y - min_y) / 2);
    float meioZ = min_z + ((max_z - min_z) / 2);
    //    printf("\nMX: %f\nMY: %f\nMZ: %f", meioX, meioY, meioZ);
    int quadro;
    for (int contador = 0; contador < N; contador++) {
        quadro = 0;
        if (x[contador] > meioX) {
            quadro += 1;
        }
        if (y[contador] > meioY) {
            quadro += 2;
        }
        if (z[contador] > meioZ) {
            quadro += 4;
        }
        DivEntrada[quadro].push_back(contador);
    }

    /*printf("\n\nEntrada dividida:\n");
    for (int c = 0; c < 8; c++) {
        printf("\nQuadro %d", c);
        int i;
        for (i = 0; i < DivEntrada[c].size(); i++) {
            printf("\n\tcorpo index : %d   \t id: %d", DivEntrada[c].at(i), iden[DivEntrada[c].at(i)]);
        }
        printf("\n\tnumero total corpos : %d", i);
    }*/
}

float verif(Corpo* n1, Corpo* n2) {
    return sqrt(pow(n1->x - n2->x, 2) + pow(n1->y - n2->y, 2) + pow(n1->z - n2->z, 2));
}

main(int argc, char **argv) {
    float local_v[100];
    char *Arg1;
    if (argc != 2) {
        puts("Por favor entre com o nome do arquivo de dados! ");
        getchar();
        exit(1);
    }

    //puts("Iniciando...");
    Arg1 = argv[1];

    LeDados(Arg1);
    //printf("\n Entrada : %d dados", N);

    gettimeofday(&tempo_inicio, NULL);
    //printf("\nraio = %f", raio);


    for (int i = 0; i < N; i++) {
        DivEntrada[0].push_back(i);
    }
    Segmento s;
    CriaArvore(&s, 0);

    /*for(int i=0;i<s.l_grupos->size();i++){
         printf("\n Grupo %d",i);
         for(int j=0;j<s.l_grupos->at(i)->lista->size();j++){
             printf("\n\t ID: %d",s.l_grupos->at(i)->lista->at(j)->id);
         }
     }
     */
    int cont = 0;
    for (int x = 0; x < s.l_grupos->size(); x++) {
        //printf("\nG %d size: %d", x, s.l_grupos->at(x)->lista->size());
        if (s.l_grupos->at(x)->lista->size() == 0) {
            s.l_grupos->erase(s.l_grupos->begin() + x);
            x--;
            //                 printf("\napagou");
        } else {
            if (s.l_grupos->at(x)->lista->size() > 1) {
                cont++;
            }
        }
    }
    //      IMPRIME TEMPO 

    gettimeofday(&tempo_fim, NULL);
    tf = (double) tempo_fim.tv_usec + ((double) tempo_fim.tv_sec * (1000000.0));
    ti = (double) tempo_inicio.tv_usec + ((double) tempo_inicio.tv_sec * (1000000.0));
    tempo = (tf - ti) / 1000;
    printf("\nTempo total :    %.5f (ms)\n", tempo);


    //        printf("\n-------------------IMPRIME GRUPOS-----------------");
    //        int cont2 = 0;
    //        int cont3 = 0;
    //        for (int i = 0; i < s.l_grupos->size(); i++) {
    //            if (s.l_grupos->at(i)->lista->size() == 1) {
    //                cont2++;
    //                cont3++;
    //            } else {
    //                if (s.l_grupos->at(i)->lista->size() > 1) {
    //                    printf("\n G %d size: %d", i, s.l_grupos->at(i)->lista->size());
    //                    cont3 += s.l_grupos->at(i)->lista->size();
    //                    for (int j = 0; j < s.l_grupos->at(i)->lista->size(); j++) {
    //                        printf("\n\t id : %d", s.l_grupos->at(i)->lista->at(j)->id);
    //                    }
    //                }
    //            }
    //        }
    printf("\nNumero de grupos: %d", s.l_grupos->size());
    printf("\nNumero de grupos + 1 elemento: %d", cont);
    //        printf("\nNumero de grupos exatamente 1 elemento: %d", cont2);
    //        printf("\nNumero de corpos nos grupos: %d", cont3);

    gettimeofday(&tempo_fim, NULL);
    tf = (double) tempo_fim.tv_usec + ((double) tempo_fim.tv_sec * (1000000.0));
    ti = (double) tempo_inicio.tv_usec + ((double) tempo_inicio.tv_sec * (1000000.0));
    tempo = (tf - ti) / 1000;
    printf("\nTempo total :    %.5f (ms)\n", tempo);

    printf("\nFIM EXECUCAO\n\n");
    return 0;
}
