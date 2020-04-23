#ifndef CALCUL_H
#define CALCUL_H
#include<QtCore>
#include<mainwindow.h>


class Calcul:public QThread
{
private:
    QMutex * mutex;
public:
    Calcul(QMutex *mu, QImage *image, vector<Object*>* listObject, vector<Source*>* listSource, camera* eye, Background* back, Floor * sol,int nbthread);
    void run();
    QImage getImage();
    int id;
    int nbthread;
    QImage * image;
    vector<Object*>* listObject;
    vector<Source*>* listSource;
    camera* eye;
    Background* back;
    Floor* sol;

};

#endif // CALCUL_H
