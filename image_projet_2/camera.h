#ifndef CAMERA_H
#define CAMERA_H

#include <vector>
#include <cmath>


using namespace std;

class camera
{
private:

    int resW,resH;
    std::vector<double> vectC; //camera origin
    std::vector<double> vectG; //grid center
    //int fov,alpha,beta;
    std::vector<double> normal;//grid normal
    std::vector<double> vectW;//base grid plan : (vectW,vectH)
    std::vector<double> vectH;



public:
    camera();
    camera(int h,int w,double xC,double yC,double zC, int fov, int alpha,int beta);
    ~camera();
    int getWidth();
    int getHeight();
    std::vector<double> getC();
    std::vector<double> getG();
    std::vector<double> getNrm();
    std::vector<double> getVectW();
    std::vector<double> getVectH();

    std::vector<double> getPixel(int i,int j);

};

#endif // CAMERA_H
