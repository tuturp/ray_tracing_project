#ifndef FLOOR_H
#define FLOOR_H

#include<QImage>
class Floor
{
private:
    int Zf, R, G, B;
public:
    Floor();
    ~Floor();
    Floor(int Zf, QImage solImage);
    int getR();
    int getG();
    int getB();
    int getZf();
    QImage solImage;
};

#endif // FLOOR_H
