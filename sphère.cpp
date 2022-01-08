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
        Vector O;
        double R;
        Vector C;
};

bool intersect(const Ray& r, Sphere& s, Vector& P, Vector& N) {
    //P la position, N la normale 
    double a = 1;
    double b = 2*dot(r.direction,r.origine-s.O);
    double c = (r.origine - s.O).norm2()-s.R*s.R;

    double delta = b*b - 4*a*c;
    if (delta<0) {return false;}
        else{
            double x1 = (-b -sqrt(delta))/(2*a);
            double x2 = (-b +sqrt(delta))/(2*a);
            if (x2 < 0) {return false;}  //pas d'intersection
            else { 
                double t;
                if (x1>0){
                    t=x1;
                }
                else{
                    t=x2;
                }
                P = r.origine + t*r.direction;
                N = (P-s.O);
                N.normalize();
            };
            return true;
    };
};
 
int main() {
    int W = 512;
    int H = 512;
 
    Vector source(-15, 20, -30);
    double intensite = 500000;

	Sphere s(Vector(0,0,-55),20, Vector(0,0,1));   //spère bleue 
	double fov=60*M_PI/180.;   //champ visuel 
	double tanfov2 = tan(fov/2);

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {   // on parcourt l'image 
        for (int j = 0; j < W; j++) {
            Vector C(0,0,0);  //origine du vecteur vision
			Vector u(j-W/2+0.5, i-H/2+0.5, -W/(2*tan(fov/2.)));   // direction du vecteur vision
			u.normalize();
			Ray r(C,u);  //rayon de la vision

            Vector P,N; //vecteurs position et normal
            Vector eclairage = Vector(0,0,0);
			if (intersect(r,s,P,N) == true){
                Vector A = source - P;
                double B = (source - P).norm2();
                A.normalize();
                eclairage = s.C*fmax(0,dot(A,N))*intensite/B;
				image[((H-i-1) * W + j) * 3 + 0] = fmin(255, fmax(0,eclairage[0]));  //coordonnée rouge
            	image[((H-i-1) * W + j) * 3 + 1] = fmin(255,fmax(0,eclairage[1]));  //coordonnée verte
            	image[((H-i-1) * W + j) * 3 + 2] = fmin(255,fmax(0,eclairage[2]));  //coordonnée bleue 
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