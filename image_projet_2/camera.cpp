#include "camera.h"
#include <cmath>

camera::camera(): //default camera 1920*1080 fov=pi/2 nrm=x
    resW(1920),resH(1080),
    vectC{0,0,0},
    vectG{resW/2/tan(M_PI/2/2),0,0},
    normal{1,0,0},
    vectW{0, -1, 0},
    vectH{0, 0, -1}
    {}

camera::camera(int h,int w,double xC,double yC,double zC, int fov, int alpha,int beta):
    resW(w),resH(h),
    vectC{xC,yC,zC},
    vectG{xC+w/2/tan(fov/2)*cos(alpha)*(1+sin(beta)), yC+w/2/tan(fov/2)*sin(alpha)*(1+sin(beta)), zC+w/2/tan(fov/2)*cos(beta)},
    normal{cos(alpha)*(1+sin(beta)), sin(alpha)*(1+sin(beta)), cos(beta)},
    vectW{sin(alpha), -cos(alpha), 0},
    vectH{-cos(alpha)*cos(beta), sin(alpha)*(1+sin(beta)), sin(beta)}
    {}

camera::~camera(){};

int camera::getWidth(){return resW;};
int camera::getHeight(){return resH;};
std::vector<double> camera::getC(){return vectC;}
std::vector<double> camera::getG(){return vectG;}
std::vector<double> camera::getNrm(){return normal;};
std::vector<double> camera::getVectW(){return vectW;};
std::vector<double> camera::getVectH(){return vectH;};

std::vector<double> camera::getPixel(int i,int j)//returns location of (i,j) pixel of grid
{
    std::vector<double> pixel=vectG;
    for (int k=0;k<3;k++)
    {
        pixel[k]+=(i-resW/2)*vectW[k]+(j-resH/2)*vectH[k];
    }
    return pixel;
};






