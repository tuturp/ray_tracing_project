#include "floor.h"

Floor::Floor(){}

Floor::~Floor(){}

Floor::Floor(int Zf, QImage solImage):
    Zf(Zf),
    solImage(solImage){}

int Floor:: getR(){return R;}

int Floor:: getG(){return G;}

int Floor:: getB(){return B;}

int Floor::getZf(){return Zf;}
