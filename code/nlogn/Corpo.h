/* 
 * File:   Corpo.h
 * Author: otavio
 *
 * Created on 1 de Novembro de 2013, 10:23
 */

#ifndef CORPO_H
#define	CORPO_H

class Corpo {
public:
    int id, igru;
    float x, y, z;
    Corpo(int id, float x, float y, float z);
    Corpo(const Corpo& orig);
    virtual ~Corpo();
private:

};

#endif	/* CORPO_H */

