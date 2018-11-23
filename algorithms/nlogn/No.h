/* 
 * File:   No.h
 * Author: otavio
 *
 * Created on 2 de Julho de 2013, 14:26
 */

#ifndef NO_H
#define	NO_H
#include <vector>
#include "grupo.h"
#include "Corpo.h"
using namespace std;

class No {
public:
    float xmax, xmin, ymax, ymin, zmax, zmin;
    float f_xmax, f_xmin, f_ymax, f_ymin, f_zmax, f_zmin;
    bool populado,menor_raio;
    Corpo *c;
    vector<Corpo*> fronteira;
    vector<Corpo*> conteudo;
    No* filhos[8];
    int identificador;
    static int n_grupos;

    No(Corpo *c, float xmax, float xmin, float ymax, float ymin, float zmax,
            float zmin);
    No(const No& orig);

    void print_fronteira();

    void print_no(No a);

    void relabel(vector<struct grupo*> *l, int g1, int g2);

    void add(vector<struct grupo*> *l, Corpo *n, float xmax, float xmin, float ymax, float ymin, float zmax,
            float zmin);


    virtual ~No();




private:

};

#endif	/* NO_H */

