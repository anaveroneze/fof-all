/* 
 * File:   grupo.cpp
 * Author: otavio
 * 
 * Created on 11 de Setembro de 2013, 12:41
 */

#include "grupo.h"
#include "No.h"
#include <vector>
#include <stdio.h>
#include "Corpo.h"

using namespace std;

grupo::grupo(int i) {
    i_gru=i;
    lista = new vector<struct Corpo*>();
}

void grupo::add_no(Corpo *a){
    lista->push_back(a);
}

void grupo::list_ids(){
//    printf("\n IGRU = %d",i_gru);
//    if(lista.size()==0){
//        printf("\nGRUPO VAZIO");
//    }
//    else{
//        int i = 0;
//        while(lista.size()>i){
//                 printf("\n no id %d",lista[i]->id);
//                 i++;
//        }
//    }
}
grupo::grupo(const grupo& orig) {
}

grupo::~grupo() {
}

