/* 
 * File:   grupo.h
 * Author: otavio
 *
 * Created on 11 de Setembro de 2013, 12:41
 */

#ifndef GRUPO_H
#define	GRUPO_H
#include "No.h"
#include "Corpo.h"
#include <vector>

using namespace std;

class grupo {
    int i_gru;
public:
    vector<struct Corpo*> *lista;
    int getI_gru() {
        return i_gru;
    }
    void list_ids();
    void lista_grupos();
    grupo(int i);
    void add_no(struct Corpo *a);
    grupo(const grupo& orig);
    virtual ~grupo();
private:
};

#endif	/* GRUPO_H */

