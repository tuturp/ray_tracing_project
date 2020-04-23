#ifndef BACKGROUND_H
#define BACKGROUND_H

#include<QImage>

class Background
{

private:
    int R,G,B;
public:
    Background();
    Background(QImage backImage);
    int  getR();
    int getG();
    int getB();
    QImage backImage;
};

#endif // BACKGROUND_H
