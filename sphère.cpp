#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>     // std::cout
#include <algorithm>    
#include <math.h>    
 
class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        //constructeur. Par défaut, le point est à l'origine. 
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};
 
Vector operator+(const Vector& a, const Vector& b) {
    //addition de deux vecteurs
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {  //multiplication terme à terme, utilisée pour les sphères colorées 
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {   //produit scalaire 
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {   //produit vectoriel 
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
 
class Ray {
public:
	Ray(const Vector& O, const Vector& C) : origine(O), direction(C){};
	Vector origine;
	Vector direction;
};

class Source { 
	public:
		Source(const double intensite, const Vector& origine) : I(intensite), O(origine){};
		double I;
		Vector O;
};

class Sphere {  //une origine et un rayon
    public:
        Sphere(const Vector& origine, const double rayon, const Vector couleur) : O(origine), R(rayon) , C(couleur) {};
    bool intersect(const Ray& r, Vector& P, Vector& N, double& t) {
    //P la position, N la normale 
    double a = 1;
    double b = dot(r.direction,r.origine-O);
    double c = (r.origine - O).norm2()-R*R;

    double delta = b*b - a*c;
    if (delta<0) {return false;}
    else{
        double x1 = (-b -sqrt(delta))/(a);
        double x2 = (-b +sqrt(delta))/(a);
        if (x2 < 0) {return false;}  //pas d'intersection
        else { 
            if (x1>0){
                t=x1;
            }
            else{
                t=x2;
            }
            P = r.origine + t*r.direction;
            N = (P-O);
            N.normalize();
        };
        return true;
    };
};
    Vector O;
    double R;
    Vector C;
};

class Scene {  // ensemble de sphères 
public: 
    Scene() {};
    void ajoutersphere(const Sphere& s) {spheres.push_back(s);}
    std::vector<Sphere> spheres;

    bool intersection(const Ray& r, Vector& P, Vector& N, int& sphere_inter_id) {   //pour toutes les sphères
        //P la position, N la normale 
        bool inter = false;
        double min_t = 10E10;
        for (int i=0; i<spheres.size();i++){
            Vector Pprime, Nprime;
            double t;
            bool interprime = spheres[i].intersect(r,Pprime,Nprime,t);
            if (interprime){  //il y a bien une intersection avec cette sphère
                inter = true;
                if (t < min_t){   //la sphère est plus proche 
                    min_t = t; 
                    P = Pprime;
                    N = Nprime;
                    sphere_inter_id = i;
                }
            }
        }
        return inter;
    };
};
 
int main() {
    int W = 512;
    int H = 512;
 
    Vector source(-10, 20, 40);
    double intensite = 100000000;

	Sphere s(Vector(0,0,0),10, Vector(0,0,1));   //spère bleue 
	double fov=60*M_PI/180.;   //champ visuel 
	double tanfov2 = tan(fov/2);

    Sphere s1(Vector(0,-2000-10,0), 2000, Vector(0,1,1));  //sol blanc, de telle sorte à ce que la sphre soit posée sur le sol
    Sphere s2(Vector(0,2000+25,0), 2000, Vector(0,0,1));  //plafond
    Sphere s3(Vector(-2000-25,0,0), 2000, Vector(1,1,1));  //mur à gauche 
    Sphere s4(Vector(2000+25,0,0), 2000, Vector(1,1,1));  //mur à droite 
    Sphere s5(Vector(0,0,-2000-50), 2000, Vector(1,1,1));  //mur au fond 

    Scene scene;
    scene.ajoutersphere(s);
    scene.ajoutersphere(s1);
    scene.ajoutersphere(s2);
    scene.ajoutersphere(s3);
    scene.ajoutersphere(s4);
    scene.ajoutersphere(s5);

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {   // on parcourt l'image 
        for (int j = 0; j < W; j++) {
            Vector C(0,0,55);  //origine du vecteur vision
			Vector u(j-W/2+0.5, i-H/2+0.5, -W/(2*tan(fov/2.)));   // direction du vecteur vision
			u.normalize();
			Ray r(C,u);  //rayon de la vision

            Vector P,N; //vecteurs position et normal
            int sphere_inter_id = -1;
            Vector eclairage = Vector(0,0,0);
            double value0;
            double value1;
            double value2;
			if (scene.intersection(r,P,N, sphere_inter_id) == true){  //une intersection est trouvée 
                Vector A = source - P;
                double B = (source - P).norm2();  //vecteur unitaire qui part de P et se dirige vers la lumière 
                A.normalize();
                eclairage = scene.spheres[sphere_inter_id].C*fmax(0,dot(A,N))*intensite/B;
                for (int k=0; k<3; k++){
                    eclairage[k]=pow(eclairage[k], 1.0/2.2);
                };
				value0 = eclairage[0];  //coordonnée rouge
            	value1 = eclairage[1];  //coordonnée verte
            	value2 = eclairage[2];  //coordonnée bleue 

                //vecteur secondaire d'origine P et se dirigeant vers la lumière 
                Vector A2 = source - P;
                A2.normalize();
                Ray r2(P+0.005*A2,A2);
                Vector P2,N2;
                sphere_inter_id = -1;
                if (scene.intersection(r2,P2,N2, sphere_inter_id) == true){  //une intersection est trouvée 
                    double distance = (P-P2).norm();
                    double distance2 = (P-source).norm();
                    if (distance < distance2){
                        value0 = 0;  //coordonnée rouge
                        value1 = 0;  //coordonnée verte
                        value2 = 0;  //coordonnée bleue 
                    }
                };
				image[((H-i-1) * W + j) * 3 + 0] = fmin(255, (fmax(0,value0)));  //coordonnée rouge
            	image[((H-i-1) * W + j) * 3 + 1] = fmin(255,fmax(0,value1));  //coordonnée verte
            	image[((H-i-1) * W + j) * 3 + 2] = fmin(255,fmax(0,value2));  //coordonnée bleue 
			}
			else {
				image[((H-i-1) * W + j) * 3 + 0] = 0;
            	image[((H-i-1) * W + j) * 3 + 1] = 0;
            	image[((H-i-1) * W + j) * 3 + 2] = 0;
			};
        }
    };
    stbi_write_png("sphère.png", W, H, 3, &image[0], 0);
 
    return 0;
}