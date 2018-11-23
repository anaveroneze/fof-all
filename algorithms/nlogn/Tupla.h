/* 
 * File:   Tupla.h
 * Author: otavio
 *
 * Created on 25 de Março de 2014, 14:04
 */

#ifndef TUPLA_H
#define	TUPLA_H

class Tupla {
public:
    int ngrupos,segmento;
    Tupla(int n,int s);
    Tupla(const Tupla& orig);
    virtual ~Tupla();
private:

};

#endif	/* TUPLA_H */

