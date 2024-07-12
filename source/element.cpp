#include "element.h"

element::element():number(0), mass(0){}

element::~element(){
    delete [] atoms;
}