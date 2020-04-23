#include "floor.h"
#include "background.h"
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include<vector>
#include"objet.h"
#include"rayon.h"
#include"view.h"
#include<iostream>
#include<camera.h>
#include"calcul.h"
#include <fstream>

using namespace std;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    //Floor sol(75,200,200,200);
    //Background back(0, 0, 120);
    QImage backImage("ciel1.png");
    QImage  solImage("sol1.jpg");
    Background back(backImage);

    class camera eye;

    Floor sol(eye.getHeight()/20, solImage);
    int sizeX = 2502;
    int sizeY = 1500;


    vector<Object*> listObject;
    vector<Source*> listSource;

    /*listObject.push_back(new Object(200,-20,-13,50,161, 6, 132,0.1));
    listObject.push_back(new Object(200,150,-13,20,4, 160, 111,0.5));
    listObject.push_back(new Object(200,160,-60,23,64, 160, 70));


    listSource.push_back(new Source(-1250,100,750,50));*/
    descripTxt(&listObject, &listSource);
    QImage imageI(sizeX, sizeY, QImage::QImage::Format_RGBA8888); //création de l'image

    QMutex mu;//création d'un mutex pour éviter les conflits entre les différents thread la modifiant

    vector<Calcul*>  listThread;
    int nbthread = 3;//nombre de threads à creer (3 optimal)
    for(int i=0; i<nbthread; i++){ //création des threads
        listThread.push_back(new Calcul(&mu , &imageI, &listObject,&listSource, &eye, &back, &sol, nbthread));
        listThread[listThread.size()-1]->id = i;
        listThread[listThread.size()-1]->start();
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();//début chrono les threads travaillent

    for(int i=0; i<nbthread; i++){//attente de la conception totale de l'image par les différent threads
        listThread[i]->wait();
    }

    end = std::chrono::system_clock::now();//fin chrono
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Temps assemblage image : " << elapsed_seconds.count()<< std::endl;


    imageI.save("image.png");
    QGraphicsScene *graphic = new QGraphicsScene(this);
    graphic->addPixmap(QPixmap::fromImage(imageI));

    ui->graphicsView->setScene(graphic);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::descripTxt(vector<Object*>* listObject, vector<Source*>* listSource){

    // Opens 'objet.txt'.
    std::ifstream obj("objet.txt");
    std::ifstream src("source.txt");
    std::string lineS;
    std::string line;


    while(obj){
            getline(obj, line);
            if(line=="Objet"){
                std::string typeObject;
                int centreX, centreY,centreZ ,r , R ,G, B,trans1, trans2;
                while(obj>>typeObject>>centreX>>centreY>>centreZ>>r>>R>>G>>B>>trans1>>trans2){
                    if (typeObject=="sphere"){
                        listObject->push_back(new Object(centreX,centreY,centreZ,r,R,G,B,trans1+0.1*trans2));
                    }
                }
            }

    }

    obj.close();


    while(src){
            getline(src, lineS);
            if(lineS=="Source"){
                std::string typeObject;
                int centreX, centreY,centreZ ,r ;
                while(src>>typeObject>>centreX>>centreY>>centreZ>>r){
                    if (typeObject=="Source"){
                        listSource->push_back(new Source(centreX,centreY,centreZ,r));
                    }
                }
            }

    }
    src.close();


}


