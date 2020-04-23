#include "background.h"

Background::Background(){}

Background::Background( QImage backImage):
    backImage(backImage){}

int Background:: getR(){return R;}

int Background:: getG(){return G;}

int Background:: getB(){return B;}
