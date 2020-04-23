#include "calcul.h"
#include<QtCore>
#include<iostream>
#include <QImage>
#include<mainwindow.h>






Calcul::Calcul(QMutex *mu, QImage *imageI, vector<Object*>* listObject, vector<Source*>* listSource, camera* eye, Background* back, Floor * sol,int nbthread):
    mutex(mu),
    image(imageI),
    listObject(listObject),
    listSource(listSource),
    eye(eye),
    back(back),
    sol(sol),
    nbthread(nbthread){}

void Calcul::run(){
    for (int i = (this->id)*eye->getWidth()/nbthread; i<(this->id+1)*eye->getWidth()/nbthread;i++){
        for (int j= 0 ; j<eye->getHeight();j++){
            std::vector<double> pixel = eye->getPixel(i,j);
            Rayon ray(eye, pixel[0],pixel[1], pixel[2], listObject, listSource, sol, back,i,j);
            std::vector<int> rgb=RayInterRec(&ray,listObject,listSource,eye, sol, back);
            ray.updateRGB(rgb[0],rgb[1],rgb[2]);
            mutex->lock();
            image->setPixel(i,j,qRgba(ray.getR(),ray.getG(),ray.getB(),255));
            mutex->unlock();
            ray.~Rayon();
            //delete ray;
        }

    }

}
