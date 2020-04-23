#include "plan.h"

plan::plan(){}

plan::plan(double a,double b,double c,double d,int R,int G,int B):
    rgb{R,G,B},normal{a,b,c},offset(d)
{}

int plan::getR(){return rgb[0];};
int plan::getG(){return rgb[1];};
int plan::getB(){return rgb[2];};
std::vector<double> plan::getNrm(){return normal;};
double plan::getOffset(){return offset;};

std::vector<double> plan::intersectPlanVect(std::vector<double> origin,std::vector<double> vect)
{
    double k=-1*(normal[0]*origin[0]+normal[1]*origin[1]+normal[2]*origin[2]+offset)/(vect[0]+vect[1]+vect[2]);
    std::vector<double> intersection{origin[0]+k*vect[0],origin[1]+k*vect[1],origin[2]+k*vect[2]};
    return intersection;
};
