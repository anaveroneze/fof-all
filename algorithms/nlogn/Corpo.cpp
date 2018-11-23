/* 
 * File:   Corpo.cpp
 * Author: otavio
 * 
 * Created on 1 de Novembro de 2013, 10:23
 */

#include "Corpo.h"

Corpo::Corpo(int _id, float _x, float _y, float _z) {
    id=_id;
    x=_x;
    y=_y;
    z=_z;
    igru=-1;
}

Corpo::Corpo(const Corpo& orig) {
}

Corpo::~Corpo() {
}

