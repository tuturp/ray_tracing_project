#ifndef OBJET_H
#define OBJET_H

#include <string>




class Object{
    private:
    int centreX;
    int centreY;
    int centreZ;
    int r;
    int RGB[3];
    float indice=1.0;
    std::string typeObj="sphere";



    public:
    Object();
    ~Object();
    Object(int centreX,int centreY, int centreZ , int r , int R , int G, int B);
    Object(int centreX,int centreY, int centreZ , int r , int R , int G, int B, float indice);

    void objectTrace();
    int getCentreX();
    int getCentreY();
    int getCentreZ();
    int getRayon();
    int getR();
    int getG();
    int getB();
    float getIndice();
    void setRGB(int R , int G , int B);

};

#endif // OBJET_H
