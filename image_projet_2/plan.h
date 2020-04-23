#ifndef PLAN_H
#define PLAN_H
#include <vector>
#include <string>

class plan
{
public:
    plan();
    plan(double a,double b,double c,double d,int R,int G,int B);
    ~plan();
    int getR();
    int getG();
    int getB();
    std::vector<double> getNrm();
    double getOffset();
    std::vector<double> intersectPlanVect(std::vector<double> origin,std::vector<double> vect);

private:
    int rgb[3];
    std::string typeObj="plan";
    std::vector<double> normal;
    double offset;
};

#endif // PLAN_H
