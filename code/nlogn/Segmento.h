/* 
 * File:   Segmento.h
 * Author: otavio
 *
 * Created on 10 de Março de 2014, 11:09
 */

#ifndef SEGMENTO_H
#define	SEGMENTO_H
#include "grupo.h"
#include "No.h"

class Segmento {
public:
    vector<grupo*> *l_grupos;
    No* arvore;
    Segmento();
    Segmento(const Segmento& orig);
    virtual ~Segmento();
private:

};

#endif	/* SEGMENTO_H */

