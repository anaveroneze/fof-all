/*
 * File:   No.cpp
 * Author: Otavio Migliavacca Madalosso
 *
 * Created on 2 de Julho de 2013, 14:26
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <vector>

#include "No.h"
#include "grupo.h"
#include "Corpo.h"
#include "common.h"

#define debug false

int overhead = 0;

timeval ts, te, ttotal;
timeval t_relabel_s, t_relabel_e;

No::No(Corpo *corp, float _xmax, float _xmin, float _ymax, float _ymin,
        float _zmax, float _zmin) {

    menor_raio = false;
    c = corp;


    xmax = _xmax;
    xmin = _xmin;
    ymax = _ymax;
    ymin = _ymin;
    zmax = _zmax;
    zmin = _zmin;

    if (xmax - xmin < raio && ymax - ymin < raio && zmax - zmin < raio) {
        identificador = corp->igru;
    }


    f_xmax = xmax - raio;
    f_xmin = xmin + raio;
    f_ymax = ymax - raio;
    f_ymin = ymin + raio;
    f_zmax = zmax - raio;
    f_zmin = zmin + raio;

    populado = true;


    for (int i = 0; i < 8; i++) {
        filhos[i] = NULL;
    }
}

float ver(float pos_x, float pos_y, float pos_z, No* a) {
    float dist = sqrt(pow(pos_x - a->c->x, 2) + pow(pos_y - a->c->y, 2) + pow(pos_z - a->c->z, 2));
    return dist;
}

float ver(Corpo* n1, Corpo* n2) {
    return sqrt(pow(n1->x - n2->x, 2) + pow(n1->y - n2->y, 2) + pow(n1->z - n2->z, 2));
}

void No::relabel(vector<grupo*> *l, int g1, int g2) {
    grupo *aux1, *aux2;
    if (g1 > g2) {
        int a = g1;
        g1 = g2;
        g2 = a;
    }
    aux1 = l->at(g1 - 1);
    aux2 = l->at(g2 - 1);
#if debug==true
    printf("\nRemover grupo %d", aux2->getI_gru());
#endif
    for (int i = 0; i < aux2->lista->size(); i++) {
        Corpo *c;
        c = aux2->lista->at(i);
        c->igru = g1;
        aux1->add_no(c);
    }
    l->at(g2 - 1)->lista->clear();
#if debug == true
    printf("Numero de grupos %d\n", l->size());
#endif

    ////    printf("\nTempo 1 : %d segundos e %d microssegundos", te.tv_sec - ts.tv_sec, abs(te.tv_usec-ts.tv_usec));
    //    ttotal.tv_sec+=te.tv_sec - ts.tv_sec;
    //    ttotal.tv_usec+=abs(te.tv_usec-ts.tv_usec);
    //    printf("\nTempo: %d seg e %d usec",ttotal.tv_sec,ttotal.tv_usec);
}

void No::print_fronteira() {
#if debug==true
    printf("\nFronteira com %d corpos", fronteira.size());
    if (fronteira.size() != 0) {
        for (int i = 0; i < fronteira.size(); i++) {

            printf("\n\tid: %d", fronteira.at(i)->id);
        }
    }
#endif
}

void No::add(vector<grupo*> *l, Corpo *novo, float xmax, float xmin, float ymax, float ymin, float zmax,
        float zmin) {
    int qx, qy, qz, quadro;
    float dx, dy, dz, dist;
    //limites das areas de fronteiras, qualquer objeto maior que o maximo ou
    //menor que o minimo pertencerá a fronteira
    float lim_max_x, lim_min_x;
    float lim_max_y, lim_min_y;
    float lim_max_z, lim_min_z;
    bool novo_fronteira;
#if debug==true
    printf("Numero de grupos %d\n", l->size());
    printf("\n DADOS DO CORPO:");
    printf("\n\t ID :%d", novo->id);
    printf("\n\t POS X:%f", novo->x);
    printf("\n\t POS Y:%f", novo->y);
    printf("\n\t POS Z:%f", novo->z);
    if (novo->igru != -1) {
        printf("\n\t IGRU:%d", novo->igru);
    }
#endif
    //    centro das coordenadas

    dx = ((xmax + xmin) / 2);
    dy = ((ymax + ymin) / 2);
    dz = ((zmax + zmin) / 2);

#if debug==true
    //        printf("\n AREA DO CUBO");
    //        printf("\n\t xmax : %f ,xmin : %f ", xmax, xmin);
    //        printf("\n\t ymax : %f ,ymin : %f ", ymax, ymin);
    //        printf("\n\t zmax : %f ,zmin : %f ", zmax, zmin);
//    printf("\nAREA FRONTEIRA:");
//    printf("\n\t fxmax : %f ,fxmin : %f ", f_xmax, f_xmin);
//    printf("\n\t fymax : %f ,fymin : %f ", f_ymax, f_ymin);
//    printf("\n\t fzmax : %f ,fzmin : %f ", f_zmax, f_zmin);

#endif

    //verifica se novo nó encontra-se em uma área de fronteira
    if (novo->x <= f_xmin || novo->x >= f_xmax || novo->y <= f_ymin || novo->y >= f_ymax
            || novo->z <= f_zmin || novo->z >= f_zmax) {
#if debug==true
        printf("\nO novo está na área de fronteira");
#endif
        fronteira.push_back(novo);

        novo_fronteira = true;
    } else {
        novo_fronteira = false;
#if debug==true
        printf("\nO novo nó não está na área de fronteira");
#endif
    }


    //caso o no atual seja populado
    if (populado) {
#if debug==true
        printf("\n-v\tNo OCUPADO, realocando antigo:");
#endif
        //calcula subdivisão do novo sub no
        if (c->x <= dx) {
            qx = 0;
            lim_max_x = dx;
            lim_min_x = xmin;
        } else {
            qx = 1;
            lim_max_x = xmax;
            lim_min_x = dx;
        }
        if (c->y <= dy) {
            qy = 0;
            lim_max_y = dy;
            lim_min_y = ymin;
        } else {
            qy = 1;
            lim_max_y = ymax;
            lim_min_y = dy;
        }
        if (c->z <= dz) {
            qz = 0;
            lim_max_z = dz;
            lim_min_z = zmin;
        } else {
            qz = 1;
            lim_max_z = zmax;
            lim_min_z = dz;
        }

        //index da subdivisão ao qual o nó será alocado
        quadro = 1 * qx + 2 * qy + 4 * qz;
#if debug==true
        printf("\n-v\tAdd em filhos[%d]", quadro);
#endif

        //caso não exista a subdivisão, cria um novo objeto.
        No *no_velho = new No(c, lim_max_x, lim_min_x,
                lim_max_y, lim_min_y, lim_max_z, lim_min_z);
        filhos[quadro] = no_velho;
        //verifica se está na rea de fronteira do novo no criado
        if (c->x < no_velho->f_xmin || c->x > no_velho->f_xmax || c->y < no_velho->f_ymin || c->y > no_velho->f_ymax
                || c->z < no_velho->f_zmin || c->z > no_velho->f_zmax) {
#if debug==true
            printf("\n-v\tO  no (id : %d )está na área de fronteira", c->id);
#endif
            filhos[quadro]->fronteira.push_back(c);
        }
        //        filhos[quadro]->populado=true;
        c = NULL;
        populado = false;
    } else {
#if debug==true
        printf("\nNó pai, adicionando filho");
#endif
    }

    int verF = false;
    /****************************NAO POPULADO**********************************/
    if (novo->x < dx) {
        qx = 0;
        lim_max_x = dx;
        lim_min_x = xmin;
        if (novo->x > dx - raio) {
            verF = true;
        }
    } else {
        qx = 1;
        lim_max_x = xmax;
        lim_min_x = dx;
        if (novo->x < dx + raio) {
            verF = true;
        }
    }
    if (novo->y < dy) {
        qy = 0;
        lim_max_y = dy;
        lim_min_y = ymin;
        if (novo->y > dy - raio) {
            verF = true;
        }
    } else {
        qy = 1;
        lim_max_y = ymax;
        lim_min_y = dy;
        if (novo->y < dy + raio) {
            verF = true;
        }
    }
    if (novo->z < dz) {
        qz = 0;
        lim_max_z = dz;
        lim_min_z = zmin;
        if (novo->z > dz - raio) {
            verF = true;
        }
    } else {
        qz = 1;
        lim_max_z = zmax;
        lim_min_z = dz;
        if (novo->z < dz + raio) {
            verF = true;
        }
    }
    quadro = 1 * qx + 2 * qy + 4 * qz;
#if debug==true
    printf("\n-f\tAdd em filhos[%d]", quadro);
#endif
    //verifica se o novo corpo está na fronteira do filho ao qual será adicionado.
    if (filhos[quadro] != NULL) {
        if (novo->x < filhos[quadro]->f_xmin || novo->x > filhos[quadro]->f_xmax ||
                novo->y < filhos[quadro]->f_ymin || novo->y > filhos[quadro]->f_ymax ||
                novo->z < filhos[quadro]->f_zmin || novo->z > filhos[quadro]->f_zmax) {
            verF = true;
#if debug==true
            printf("\n Sim, esta\n");
#endif
        }
#if debug==true
        else {
            printf("\n nao\n");
        }
#endif
    }


    //verifica agrupamento com vizinhos (incluindo fronteiras)
    bool achou = false;
    for (int i = 0; i < 8; i++) { //verifica todos corpos vizinhos ao cubo que irá add
        if (i != quadro && filhos[i] != NULL) {
#if debug==true
            printf("\nFilhos[%d] Populado: %s", i, (filhos[i]->populado) ? "True" : "False");
#endif
            if (filhos[i]->populado) {
#if debug==true
                printf("\n-f\t\tFilho Vizinho");
#endif
                //verificar corpo vizinho
                dist = ver(novo, filhos[i]->c);
#if debug==true
                printf("\n-f\tDistancia entre nó e filhos[%d](id = %d) = %f", i, filhos[i]->c->id, dist);
#endif
                if (dist <= raio) {
#if debug==true
                    printf("\n-f\t O nó pertence ao grupo %d", filhos[i]->c->igru);
#endif
                    if (novo->igru == -1) {
#if debug==true
                        printf("\nSEM GRUPO");
#endif
                        novo->igru = filhos[i]->c->igru;
                    } else {
                        if (novo->igru != filhos[i]->c->igru) {
#if debug==true
                            printf("\n\n\n-f\t **************pertence a 2 grupos"); //se o nó ja está em outro grupo, e se estiver
                            printf("\n-f\tPertence aos grupos : %d e %d", novo->igru, filhos[i]->c->igru);
#endif
                            relabel(l, novo->igru, filhos[i]->c->igru);
                            if (novo->igru > filhos[i]->c->igru) {
                                novo->igru = filhos[i]->c->igru;
                            }
                        }
                    }
                }
            } else {
                if (novo_fronteira || verF) {

                    //NOVO TEMPO START
                    gettimeofday(&ts, NULL);
#if debug==true
                    printf("\n-f\t\tArea fronteira do no vizinho");
                    //verificar area de fronteira
                    printf("\n-f\t\tTamanho da lista: %d", filhos[i]->fronteira.size());
#endif
                    for (int j = 0; j < filhos[i]->fronteira.size(); j++) {
                        Corpo *f = filhos[i]->fronteira.at(j);
#if debug==true
                        printf("\n-f\t\tVerifica distancia com o no %d", f->id);
#endif
                        dist = ver(f, novo);
#if debug==true
                        printf(": %f", dist);
#endif
                        if (dist <= raio) {
#if debug==true
                            printf("\n-f\t O nó pertence ao grupo %d", f->igru);
                            if (f->igru == -1) {
                                printf("\n\n\n\n------------------------__FAIL__-------------------\n\n\n");
                            }
                            printf("\n-f\t novo tem igru: %d", novo->igru);
#endif
                            if (novo->igru == -1) {
                                novo->igru = f->igru;
                                //                                l->at(novo->igru - 1)->add_no(novo); //LINHA ADD //add2 x
                            } else {
                                if (novo->igru != f->igru) {
                                    //printf("\n\n\n-f\t Pertence a 2 grupos"); //se o nó ja está em outro grupo, e se estiver
                                    //printf("\n-f\tO nó %d pertence aos grupos : %d e %d", novo->id, novo->igru, f->igru);
#if debug==true
                                    printf("\n\n\n-f\t Pertence a 2 grupos"); //se o nó ja está em outro grupo, e se estiver
                                    printf("\n-f\tPertence aos grupos : %d e %d", novo->igru, f->igru);
#endif
                                    relabel(l, novo->igru, f->igru);
                                    if (novo->igru > f->igru) {
                                        novo->igru = f->igru;
                                    }
                                }
#if debug==true
                                printf("\n-f\t novo tem igru: %d", novo->igru);
#endif
                                //                                l->at(novo->igru - 1)->add_no(novo);
                            }
                        }
                    }

                }
            }
        }
    }

    if (filhos[quadro] == NULL) {
        if (novo->igru == -1) { //só cria novo grupo caso quando for criar o no, ele ainda não tiver


#if debug==true
            printf("\nADICIONA NOVO GRUPO"); //seu IGRU.
#endif
            No::n_grupos++; //incrementa variável estatica
            novo->igru = No::n_grupos; //define igru como o novo grupo
            grupo *novo_grupo = new grupo(No::n_grupos); //adiciona nova lista de nos(grupo novo)
            l->push_back(novo_grupo);


#if debug==true
            printf("\n Adicionado novo grupo, igru %d", No::n_grupos);
#endif
        }
        No *novo_no = new No(novo, lim_max_x, lim_min_x,
                lim_max_y, lim_min_y, lim_max_z, lim_min_z);
        if (novo->x < novo_no->f_xmin || novo->x > novo_no->f_xmax ||
                novo->y < novo_no->f_ymin || novo->y > novo_no->f_ymax ||
                novo->z < novo_no->f_zmin || novo->z > novo_no->f_zmax) {


#if debug==true
            printf("\n-f\t O novo corpo está na FRONTEIRA de seu novo nó");
#endif
            novo_no->fronteira.push_back(novo);
        } else {

#if debug==true
            printf("\n-f\t O novo corpo NAO está na FRONTEIRA de seu novo nó");
#endif
        }
        filhos[quadro] = novo_no;

        gettimeofday(&t_relabel_s, NULL);

        l->at(novo->igru - 1)->add_no(novo);

        gettimeofday(&t_relabel_e, NULL);
        double tfim = (double) t_relabel_e.tv_usec + ((double) t_relabel_e.tv_sec * (1000000.0));
        double tinicio = (double) t_relabel_s.tv_usec + ((double) t_relabel_s.tv_sec * (1000000.0));
        tempo_verifica_grupo += (tfim - tinicio) / 1000;
    } else {
        //passa para proximo sub no
        filhos[quadro]->add(l, novo, lim_max_x, lim_min_x, lim_max_y, lim_min_y, lim_max_z,
                lim_min_z);
    }
}

No::No(const No & orig) {

}

No::~No() {

}

