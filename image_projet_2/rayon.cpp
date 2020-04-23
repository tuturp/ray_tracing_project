#include "rayon.h"

#include<iostream>
#include <cmath>




Rayon::Rayon(){}
Rayon::~Rayon(){}

Rayon::Rayon(camera* eye, double VectX, double VectY, double VectZ, vector<Object*>* listObject, vector<Source*>* listSource, Floor* sol, Background* back, int i ,int j) :
    origineRayX(eye->getC()[0]),
    origineRayY(eye->getC()[1]),
    origineRayZ(eye->getC()[2]),
    i(i),
    j(j)
{
    VectDir[0]=VectX-origineRayX;
    VectDir[1]=VectY-origineRayY;
    VectDir[2]=VectZ-origineRayZ;
    //RayInter( listObject , listSource, eye, sol, back);

}



int Rayon::getOrigineRayX(){return origineRayX;}
int Rayon::getOrigineRayY(){return origineRayY;}
int Rayon::getOrigineRayZ(){return origineRayZ;}
int Rayon::getInter(){return inter;}
int Rayon::getR(){return R;}
int Rayon::getG(){return G;}
int Rayon::getB(){return B;}
void Rayon::updateRGB(int r,int g, int b){R=r;G=g;B=b;}


void Rayon::RayInter(vector<Object*>* listObject, vector<Source*>* listSource, camera* eye, Floor* sol, Background* back){
    Source * testSource = (*listSource)[(listSource->size())-1];
    int Vx = VectDir[0];
    int Vy = VectDir[1];
    int Vz = VectDir[2];
    int Xs = testSource->getCentreX();
    int Ys = testSource->getCentreY();
    int Zs = testSource->getCentreZ();
    double solTestN = -(eye->getC()[2]+sol->getZf());
    double T = 0;
    if (Vz!=0){T = solTestN/Vz;}
    if (T>0){
        R= sol->getR();
        G = sol->getG();
        B = sol->getB();
        for (unsigned int s = 0 ; s< listObject->size();s++){
                Object* testObjetShadow= (*listObject)[s];
                double Xss = Vx*T + eye->getC()[0];
                double Yss = Vy*T + eye->getC()[1];
                double Zss = Vz*T + eye->getC()[2];
                int Xo = testObjetShadow->getCentreX();
                int Yo = testObjetShadow->getCentreY();
                int Zo = testObjetShadow->getCentreZ();
                int r = testObjetShadow->getRayon();
                double ass = (Xs-Xss)*(Xs-Xss)+(Ys-Yss)*(Ys-Yss)+(Zs-Zss)*(Zs-Zss);
                double bss = 2*((Xs-Xss)*(Xss-Xo)+(Ys-Yss)*(Yss-Yo)+(Zs-Zss)*(Zss-Zo));
                double css = (Xss-Xo)*(Xss-Xo)+(Yss-Yo)*(Yss-Yo)+(Zss-Zo)*(Zss-Zo)-r*r;
                double detSs = bss*bss-4*ass*css;
                if(detSs>0){
                    double solTS1 = (-bss-sqrt(detSs))/(2*ass);
                    double solTS2 = (-bss+sqrt(detSs))/(2*ass);
                    if( solTS1>0 && solTS2>0 && (R!=0 && G!=0 && B!=0)){
                        R=floor(0.5+R/2);
                        G= floor(0.5+G/2);
                        B= floor(0.5+B/2);
                    }
                }
        }


    }else {
        R= back->getR();
        G = back->getG();
        B = back->getB();
    }
    double solT1=0;
    double solT2 = 0;
    for(unsigned int i= 0; i< listObject->size(); i++){
        Object * testObjet= (*listObject)[i];
        int Xr = origineRayX;
        int Yr = origineRayY;
        int Zr = origineRayZ;
        int Xo = testObjet->getCentreX();
        int Yo = testObjet->getCentreY();
        int Zo = testObjet->getCentreZ();
        int r = testObjet->getRayon();
        double a = Vx*Vx+Vy*Vy+Vz*Vz;
        double b = 2*(Vx*(Xr-Xo)+Vy*(Yr-Yo)+Vz*(Zr-Zo));
        double c = (Xr-Xo)*(Xr-Xo)+(Yr-Yo)*(Yr-Yo)+(Zr-Zo)*(Zr-Zo)-r*r;
        double detX = b*b-4*a*c;
        if (detX>0){
            double currentT1 = (-b - sqrt(detX))/(2*a);
            double currentT2 = (-b + sqrt(detX))/(2*a);
            if (solT1==0){solT1 = currentT1;solT2 = currentT2;}
            else{
                if (currentT1<solT1 || currentT2<solT2){
                    solT1 = currentT1;
                    solT2 = currentT2;
                }
            }
            double solX;
            double solY;
            double solZ;
            if (((T>0 && (solT1< T || solT2<T))|| T<=0)&&solT1==currentT1){
                if(solT1<solT2){
                     solX =Vx*solT1 + Xr;
                     solY = Vy*solT1 + Yr;
                     solZ = Vz*solT1 + Zr;
                }else{
                    solX =Vx*solT2 + Xr;
                    solY = Vy*solT2 + Yr;
                    solZ = Vz*solT2 + Zr;
                }
                double N[3] = {(solX-Xo)/r,(solY-Yo)/r,(solZ-Zo)/r};
                double NormeN = sqrt((testSource->getCentreX()-solX)*(testSource->getCentreX()-solX) + (testSource->getCentreY()-solY)*(testSource->getCentreY()-solY)+(testSource->getCentreZ()-solZ)*(testSource->getCentreZ()-solZ));
                double L[3] = {(testSource->getCentreX()-solX)/NormeN , (testSource->getCentreY()-solY)/NormeN ,(testSource->getCentreZ()-solZ)/NormeN };
                double NormeV= sqrt((eye->getC()[0]-solX)*(eye->getC()[0]-solX) + (eye->getC()[1]-solY)*(eye->getC()[1]-solY)+(eye->getC()[2]-solZ)*(eye->getC()[2]-solZ));
                double V[3] = {(eye->getC()[0]-solX)/NormeV , (eye->getC()[1]-solY)/NormeV ,(eye->getC()[1]-solZ)/NormeV };
                double NormeH = sqrt((L[0]+V[0])*(L[0]+V[0])+(L[1]+V[1])*(L[1]+V[1])+(L[2]+V[2])*(L[2]+V[2]));
                double H[3]= {(L[0]+V[0])/NormeH,(L[1]+V[1])/NormeH,(L[2]+V[2])/NormeH};
                double hf = N[0]*H[0]+N[1]*H[1]+N[2]*H[2];
                double fctr = N[0]*L[0]+N[1]*L[1]+N[2]*L[2];
                double ka = 0.2;
                double kd = 0.8;
                double ks = 0.5;
                int n = 100;
                double specAdd = ks*pow(hf,n);
                if (ka*testObjet->getR()+kd*fctr*testObjet->getR()+specAdd*testObjet->getR()>0){
                    if (ka*testObjet->getR()+kd*fctr*testObjet->getR()+specAdd*testObjet->getR()<255){
                        R= floor(0.5+ka*testObjet->getR()+kd*fctr*testObjet->getR()+specAdd*testObjet->getR());
                    }else {R=255;}
                }else{R=0;}

                if(ka*testObjet->getG()+kd*fctr*testObjet->getG()+specAdd*testObjet->getR()>0){
                    if (ka*testObjet->getG()+kd*fctr*testObjet->getG()+specAdd*testObjet->getR()<255){
                        G= floor(0.5+ka*testObjet->getG()+kd*fctr*testObjet->getG()+specAdd*testObjet->getR());
                    }else {G=255;}
                }else{G=0;}

                if(ka*testObjet->getB()+kd*fctr*testObjet->getB()+specAdd*testObjet->getR()>0){
                    if (ka*testObjet->getB()+kd*fctr*testObjet->getB()+specAdd*testObjet->getR()<255){
                        B= floor(0.5+ka*testObjet->getB()+kd*fctr*testObjet->getB()+specAdd*testObjet->getR());
                    }else {B=255;}
                }else {B=0;}

                for (unsigned int s = 0 ; s< listObject->size();s++){
                    if (s!=i){
                        Object* testObjetShadow= (*listObject)[s];
                        Xo = testObjetShadow->getCentreX();
                        Yo = testObjetShadow->getCentreY();
                        Zo = testObjetShadow->getCentreZ();
                        r = testObjetShadow->getRayon();
                        double as = (Xs-solX)*(Xs-solX)+(Ys-solY)*(Ys-solY)+(Zs-solZ)*(Zs-solZ);
                        double bs = 2*((Xs-solX)*(solX-Xo)+(Ys-solY)*(solY-Yo)+(Zs-solZ)*(solZ-Zo));
                        double cs = (solX-Xo)*(solX-Xo)+(solY-Yo)*(solY-Yo)+(solZ-Zo)*(solZ-Zo)-r*r;
                        double detS = bs*bs-4*as*cs;
                        if(detS>0){
                            double solTS1 = (-bs-sqrt(detS))/(2*as);
                            double solTS2 = (-bs+sqrt(detS))/(2*as);
                            if( solTS1>0 && solTS2>0 && (R!=0 && G!=0 && B!=0)){
                                R=floor(0.5+ka*testObjet->getR());
                                G= floor(0.5+ka*testObjet->getG());
                                B= floor(0.5+ka*testObjet->getB());
                            }
                        }
                    }
                }

            }
        }
    }
}


std::vector<int> RayInterRec(Rayon* ray,vector<Object*>* listObject, vector<Source*>* listSource, camera* eye, Floor* sol, Background* back){
    int R,G,B;
    Source * testSource = (*listSource)[(listSource->size())-1];
    int Vx = ray->VectDir[0];
    int Vy = ray->VectDir[1];
    int Vz = ray->VectDir[2];
    int Xs = testSource->getCentreX();
    int Ys = testSource->getCentreY();
    int Zs = testSource->getCentreZ();
    double solTestN = -(eye->getC()[2]+sol->getZf());
    double T = 0;
    if (Vz!=0){T = solTestN/Vz;}
    if (Vz<0 && T>0){
        ray->pix =  sol->solImage.pixel(ray->i,ray->j);
        R = qRed(ray->pix);
        G = qGreen(ray->pix);
        B = qBlue(ray->pix);
        for (unsigned int s = 0 ; s< listObject->size();s++){
                Object* testObjetShadow= (*listObject)[s];
                double Xss = Vx*T + eye->getC()[0];
                double Yss = Vy*T + eye->getC()[1];
                double Zss = Vz*T + eye->getC()[2];
                int Xo = testObjetShadow->getCentreX();
                int Yo = testObjetShadow->getCentreY();
                int Zo = testObjetShadow->getCentreZ();
                int r = testObjetShadow->getRayon();
                double ass = (Xs-Xss)*(Xs-Xss)+(Ys-Yss)*(Ys-Yss)+(Zs-Zss)*(Zs-Zss);
                double bss = 2*((Xs-Xss)*(Xss-Xo)+(Ys-Yss)*(Yss-Yo)+(Zs-Zss)*(Zss-Zo));
                double css = (Xss-Xo)*(Xss-Xo)+(Yss-Yo)*(Yss-Yo)+(Zss-Zo)*(Zss-Zo)-r*r;
                double detSs = bss*bss-4*ass*css;
                if(detSs>0){
                    double solTS1 = (-bss-sqrt(detSs))/(2*ass);
                    double solTS2 = (-bss+sqrt(detSs))/(2*ass);
                    if( solTS1>0 && solTS2>0 && (R!=0 && G!=0 && B!=0)){
                        R=floor(0.5+R/2);
                        G= floor(0.5+G/2);
                        B= floor(0.5+B/2);
                    }
                }
        }


    }if(Vz>=0){
        ray->pix=  back->backImage.pixel(ray->i,ray->j);

        R = 0;
                G = 0;
                B = 0;
        /*R = qRed(ray->pix);
        G = qGreen(ray->pix);
        B = qBlue(ray->pix);*/
    }
    double solT1=0;
    double solT2 = 0;
    for(unsigned int i= 0; i< listObject->size(); i++){
        Object * testObjet= (*listObject)[i];
        int Xr = ray->getOrigineRayX();
        int Yr = ray->getOrigineRayY();
        int Zr = ray->getOrigineRayZ();
        int Xo = testObjet->getCentreX();
        int Yo = testObjet->getCentreY();
        int Zo = testObjet->getCentreZ();
        int r = testObjet->getRayon();
        double a = Vx*Vx+Vy*Vy+Vz*Vz;
        double b = 2*(Vx*(Xr-Xo)+Vy*(Yr-Yo)+Vz*(Zr-Zo));
        double c = (Xr-Xo)*(Xr-Xo)+(Yr-Yo)*(Yr-Yo)+(Zr-Zo)*(Zr-Zo)-r*r;
        double detX = b*b-4*a*c;
        if (detX>0){
            double currentT1 = (-b - sqrt(detX))/(2*a);
            double currentT2 = (-b + sqrt(detX))/(2*a);
            if (solT1==0){solT1 = currentT1;solT2 = currentT2;}
            else{
                if (currentT1<solT1 || currentT2<solT2){
                    solT1 = currentT1;
                    solT2 = currentT2;
                }
            }
            double solX;
            double solY;
            double solZ;
            if (((Vz>0 && (solT1< Vz || solT2<Vz))|| Vz<=0)&&solT1==currentT1){//((Vz>0 && (solT1< Vz || solT2<Vz))|| Vz<=0)&&solT1==currentT1  ((T>0 && (solT1< T || solT2<T))|| T<=0)&&solT1==currentT1
                if(solT1<solT2){
                     solX =Vx*solT1 + Xr;
                     solY = Vy*solT1 + Yr;
                     solZ = Vz*solT1 + Zr;
                }else{
                    solX =Vx*solT2 + Xr;
                    solY = Vy*solT2 + Yr;
                    solZ = Vz*solT2 + Zr;
                }
                double N[3] = {(solX-Xo)/r,(solY-Yo)/r,(solZ-Zo)/r};
                double NormeN = sqrt((testSource->getCentreX()-solX)*(testSource->getCentreX()-solX) + (testSource->getCentreY()-solY)*(testSource->getCentreY()-solY)+(testSource->getCentreZ()-solZ)*(testSource->getCentreZ()-solZ));
                double L[3] = {(testSource->getCentreX()-solX)/NormeN , (testSource->getCentreY()-solY)/NormeN ,(testSource->getCentreZ()-solZ)/NormeN };
                double NormeV= sqrt((eye->getC()[0]-solX)*(eye->getC()[0]-solX) + (eye->getC()[1]-solY)*(eye->getC()[1]-solY)+(eye->getC()[2]-solZ)*(eye->getC()[2]-solZ));
                double V[3] = {(eye->getC()[0]-solX)/NormeV , (eye->getC()[1]-solY)/NormeV ,(eye->getC()[1]-solZ)/NormeV };
                double NormeH = sqrt((L[0]+V[0])*(L[0]+V[0])+(L[1]+V[1])*(L[1]+V[1])+(L[2]+V[2])*(L[2]+V[2]));
                double H[3]= {(L[0]+V[0])/NormeH,(L[1]+V[1])/NormeH,(L[2]+V[2])/NormeH};
                double hf = N[0]*H[0]+N[1]*H[1]+N[2]*H[2];
                double fctr = N[0]*L[0]+N[1]*L[1]+N[2]*L[2];
                double ka = 0.2;
                double kd = 0.8;
                double ks = 0.5;
                int n = 100;
                double specAdd = ks*pow(hf,n);
                if (ka*testObjet->getR()+kd*fctr*testObjet->getR()+specAdd*testObjet->getR()>0){
                    if (ka*testObjet->getR()+kd*fctr*testObjet->getR()+specAdd*testObjet->getR()<255){
                        R= floor(0.5+ka*testObjet->getR()+kd*fctr*testObjet->getR()+specAdd*testObjet->getR());
                    }else {R=255;}
                }else{R=0;}

                if(ka*testObjet->getG()+kd*fctr*testObjet->getG()+specAdd*testObjet->getR()>0){
                    if (ka*testObjet->getG()+kd*fctr*testObjet->getG()+specAdd*testObjet->getR()<255){
                        G= floor(0.5+ka*testObjet->getG()+kd*fctr*testObjet->getG()+specAdd*testObjet->getR());
                    }else {G=255;}
                }else{G=0;}

                if(ka*testObjet->getB()+kd*fctr*testObjet->getB()+specAdd*testObjet->getR()>0){
                    if (ka*testObjet->getB()+kd*fctr*testObjet->getB()+specAdd*testObjet->getR()<255){
                        B= floor(0.5+ka*testObjet->getB()+kd*fctr*testObjet->getB()+specAdd*testObjet->getR());
                    }else {B=255;}
                }else {B=0;}
                //shadow
                for (unsigned int s = 0 ; s< listObject->size();s++){
                    if (s!=i){
                        Object* testObjetShadow= (*listObject)[s];
                        Xo = testObjetShadow->getCentreX();
                        Yo = testObjetShadow->getCentreY();
                        Zo = testObjetShadow->getCentreZ();
                        r = testObjetShadow->getRayon();
                        double as = (Xs-solX)*(Xs-solX)+(Ys-solY)*(Ys-solY)+(Zs-solZ)*(Zs-solZ);
                        double bs = 2*((Xs-solX)*(solX-Xo)+(Ys-solY)*(solY-Yo)+(Zs-solZ)*(solZ-Zo));
                        double cs = (solX-Xo)*(solX-Xo)+(solY-Yo)*(solY-Yo)+(solZ-Zo)*(solZ-Zo)-r*r;
                        double detS = bs*bs-4*as*cs;
                        if(detS>0){
                            double solTS1 = (-bs-sqrt(detS))/(2*as);
                            double solTS2 = (-bs+sqrt(detS))/(2*as);
                            if( solTS1>0 && solTS2>0 && (R!=0 && G!=0 && B!=0)){
                                R=floor(0.5+ka*testObjet->getR());
                                G= floor(0.5+ka*testObjet->getG());
                                B= floor(0.5+ka*testObjet->getB());
                            }
                        }
                    }
                }
                //
                if (testObjet->getIndice()!=1.0)
                {
                    // 1e deviation
                    double angle1=asin(testObjet->getIndice()*asin(acos(hf)));
                    double sinAngle1=sin(M_PI-angle1);
                    double cosAngle1=cos(M_PI-angle1);
                    double axe1[3]={N[1]*Vz-N[2]*Vy,N[2]*Vx-N[0]*Vz,N[0]*Vy-N[1]*Vx};
                    double nV1[3]={N[0]*(axe1[0]*axe1[0]*(1-cosAngle1)+cosAngle1)+N[1]*(axe1[0]*axe1[1]*(1-cosAngle1)-axe1[2]*sinAngle1)+N[2]*(axe1[0]*axe1[2]*(1-cosAngle1)+axe1[1]*sinAngle1),
                                   N[1]*(axe1[1]*axe1[1]*(1-cosAngle1)+cosAngle1)+N[2]*(axe1[2]*axe1[1]*(1-cosAngle1)-axe1[0]*sinAngle1)+N[0]*(axe1[0]*axe1[1]*(1-cosAngle1)+axe1[2]*sinAngle1),
                                   N[2]*(axe1[2]*axe1[2]*(1-cosAngle1)+cosAngle1)+N[0]*(axe1[0]*axe1[2]*(1-cosAngle1)-axe1[1]*sinAngle1)+N[1]*(axe1[1]*axe1[2]*(1-cosAngle1)+axe1[0]*sinAngle1)};
                    // 2e deviation


                    double coeff1 = nV1[0]*nV1[0]+nV1[1]*nV1[1]+nV1[2]*nV1[2];
                    double coeff2 = 2*(nV1[0]*(solX-Xo)+nV1[1]*(solY-Yo)+nV1[2]*(solZ-Zo));
                    double coeff3 = (solX-Xo)*(solX-Xo)+(solY-Yo)*(solY-Yo)+(solZ-Zo)*(solY-Zo)-r*r;
                    double delta = coeff2*coeff2-4*coeff1*coeff3;
                    double nX=0,nY=0,nZ=0;
                    double nN[3];
                    if (delta>0)
                    {
                        double sol1 = (-coeff2 - sqrt(delta))/(2*coeff1);
                        double sol2 = (-coeff2 + sqrt(delta))/(2*coeff1);

                            if(sol1>sol2){
                                 nX =nV1[0]*sol1 + solX;
                                 nY = nV1[1]*sol1 + solY;
                                 nZ = nV1[2]*sol1 + solZ;
                            }
                            else{
                                nX =nV1[0]*sol2 + solX;
                                nY = nV1[1]*sol2 + solY;
                                nZ = nV1[2]*sol2 + solY;
                            }
                            //double nN[3] = {(nX-Xo)/r,(nY-Yo)/r,(nZ-Zo)/r};
                            nN[0]=(nX-Xo)/r;
                            nN[1]=(nY-Yo)/r;
                            nN[2]=(nZ-Zo)/r;

                    }


                    double angle2=asin(testObjet->getIndice()*asin(acos(nV1[0]*nN[0]+nV1[1]*nN[1]+nV1[2]*nN[2])));
                    double sinAngle2=sin(M_PI-angle2);
                    double cosAngle2=cos(M_PI-angle2);
                    double axe2[3]={nN[1]*nV1[2]-nN[2]*nV1[1],nN[2]*nV1[0]-nN[0]*nV1[2],nN[0]*nV1[1]-nN[1]*nV1[0]};
                    double nV2[3]={nN[0]*(axe2[0]*axe2[0]*(1-cosAngle2)+cosAngle2)+nN[1]*(axe2[0]*axe2[1]*(1-cosAngle2)-axe2[2]*sinAngle2)+nN[2]*(axe2[0]*axe2[2]*(1-cosAngle2)+axe2[1]*sinAngle2),
                                   nN[1]*(axe2[1]*axe2[1]*(1-cosAngle2)+cosAngle2)+nN[2]*(axe2[2]*axe2[1]*(1-cosAngle2)-axe2[0]*sinAngle2)+nN[0]*(axe2[0]*axe2[1]*(1-cosAngle2)+axe2[2]*sinAngle2),
                                   nN[2]*(axe2[2]*axe2[2]*(1-cosAngle2)+cosAngle2)+nN[0]*(axe2[0]*axe2[2]*(1-cosAngle2)-axe2[1]*sinAngle2)+nN[1]*(axe2[1]*axe2[2]*(1-cosAngle2)+axe2[0]*sinAngle2)};




                    camera nEye(0,0,nX,nY,nZ,0,0,0);

                    Rayon nRay(&nEye,nV2[0]+nX,nV2[1]+nY,nV2[2]+nZ,listObject,listSource,sol,back,ray->i,ray->j);

                    std::vector<int> nRGB=RayInterRec(&nRay,listObject,listSource,&nEye,sol,back);
                    if(3*R/4+nRGB[0]/4<255){R=3*R/4+nRGB[0]/4;}else{R=255;}
                    if(3*G/4+nRGB[1]/4<255){G=3*G/4+nRGB[1]/4;}else{G=255;}
                    if(3*B/4+nRGB[2]/4<255){B=3*B/4+nRGB[2]/4;}else{B=255;}
                }

            }
        }
    }
    std::vector<int> rgb{R,G,B};
    return rgb;
}



