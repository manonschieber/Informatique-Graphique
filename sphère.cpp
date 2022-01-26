#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>     // std::cout
#include <algorithm>    
#include <math.h>    
#include <random>
#include <iomanip>

using std::cout;
using std::endl;
using std::setprecision;

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

Vector operator-(const Vector& b) {
    return Vector(- b[0], - b[1], - b[2]);
}
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
    Sphere(const Vector& origine, const double rayon, const Vector couleur, bool isMirror = false, bool isTransparent = false) : O(origine), R(rayon) , C(couleur), isMirror(isMirror), isTransparent(isTransparent){};
    bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const {
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
    bool isMirror;
    bool isTransparent;
};

class Scene {  // ensemble de sphères 
public: 
    Scene() {};
    void ajoutersphere(const Sphere& s) {spheres.push_back(s);}
    std::vector<Sphere> spheres;
    Vector source;
    double intensite;
    bool intersection(const Ray& r, Vector& P, Vector& N, int& sphere_inter_id) const {   //pour toutes les sphères
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

Vector getColor(const Ray &r, const Scene &scene, int nombre_rebond){   //renvoie la couleur du pixel 
    if (nombre_rebond==0) return Vector(0,0,0);

    Vector P,N; //vecteurs position et normal
    int sphere_inter_id;
    Vector intensite_pixel(0,0,0);

	if (scene.intersection(r,P,N, sphere_inter_id) == true){  //une intersection est trouvée 
        if (scene.spheres[sphere_inter_id].isMirror){  //sphère miroir
            Vector direction_miroir = r.direction - 2*dot(N,r.direction)*N;
            Ray rayon_miroir(P+0.001*N, direction_miroir);
            intensite_pixel = getColor(rayon_miroir, scene, nombre_rebond -1);  
        };
        if (scene.spheres[sphere_inter_id].isTransparent){  //sphère transparente
            double indice1=1;   //air
            double indice2 = 1.3;   //verre
            Vector normale_transparence(N);
            if (dot(r.direction, N) >0 ){   //on est en train de sortir de la sphère 
                indice1 = 1.3;
                indice2=1;
                normale_transparence = -N;
            }
            double nombre_racine = 1 - indice1/indice2*indice1/indice2*(1-dot(normale_transparence,r.direction)*dot(normale_transparence,r.direction));
            if (nombre_racine>0){
                Vector direction_refracte = (indice1/indice2)*(r.direction - dot(r.direction, normale_transparence)*normale_transparence) - normale_transparence * sqrt(nombre_racine);
                Ray rayon_refracte(P-0.001*normale_transparence, direction_refracte);
                intensite_pixel = getColor(rayon_refracte, scene, nombre_rebond -1);  //un rebond de moins    
            ;}
        };
        //CONTRIBUTION DIRECTE : pas de sphère transparente, ni miroir
        Vector A2 = scene.source - P;
        A2.normalize();
        Ray r2(P+0.005*A2,A2);
        Vector P2,N2;
        int sphere_inter_id2;

        if (scene.intersection(r2,P2,N2, sphere_inter_id2) == true){  //une intersection est trouvée 
            double distance = (P-P2).norm();
            double distance2 = (scene.source-P).norm();
            if (distance < distance2){
                intensite_pixel = Vector(0,0,0);
            }
            else{
                Vector A = scene.source - P;
                double B = (scene.source - P).norm2();  //vecteur unitaire qui part de P et se dirige vers la lumière 
                A.normalize();
                intensite_pixel = scene.spheres[sphere_inter_id].C/M_PI*scene.intensite*fmax(0,dot(A,N))/B;
            };
        };

        //CONTRIBUTION INDIRECTE
        std::random_device rd;
        std::default_random_engine eng(rd());
        std::uniform_real_distribution<float> distr(0, 1);
        setprecision(6);
        double t1=distr(eng);
        double t2=distr(eng);

        double x = cos(2*M_PI*t1)*sqrt(1-t2);
        double y = sin(2*M_PI*t1)*sqrt(1-t2);
        double z = sqrt(t2);
        Vector direction_aleatoire = Vector(x,y,z);

        Vector T1;
        if (abs(N[0]) < abs(N[1]) && abs(N[0]) < abs(N[2])) {
            T1 = Vector(0, N[2], -N[1]);
        } else{
            if (abs(N[1]) < abs(N[0]) && abs(N[1]) < abs(N[2])){
                T1 = Vector(N[2], 0, -N[0]);
            } else{
                T1 = Vector(-N[2], N[0], 0);
            }
        };
        T1.normalize();
        Vector T2 = cross(N,T1);

        direction_aleatoire = direction_aleatoire[0]*T1+direction_aleatoire[1]*T2+direction_aleatoire[2]*N;
        Ray rayon_aleatoire(P+0.001*N, direction_aleatoire);
        intensite_pixel= intensite_pixel + getColor(rayon_aleatoire, scene, nombre_rebond -1)*scene.spheres[sphere_inter_id].C;  

    };
    return intensite_pixel;
};

int main() {
    int W = 512;
    int H = 512;

	Sphere s(Vector(-12,0,0),10, Vector(1,1,1), false, false);   //spère bleue 
	Sphere sbis(Vector(12,0,0),10, Vector(0,0,1), false, false);
	double fov=60*M_PI/180.;   //champ visuel 
	double tanfov2 = tan(fov/2);


    Sphere s1(Vector(0,-2000-10,0), 2000, Vector(0,1,1));  //sol blanc, de telle sorte à ce que la sphre soit posée sur le sol
    Sphere s2(Vector(0,2000+25,0), 2000, Vector(0,0,1));  //plafond
    Sphere s3(Vector(-2000-25,0,0), 2000, Vector(1,1,1));  //mur à gauche 
    Sphere s4(Vector(2000+25,0,0), 2000, Vector(1,1,1));  //mur à droite 
    Sphere s5(Vector(0,0,-2000-25), 2000, Vector(1,1,1));  //mur au fond 

    Scene scene;
    scene.ajoutersphere(s);
    scene.ajoutersphere(sbis);
    scene.ajoutersphere(s1);
    scene.ajoutersphere(s2);
    scene.ajoutersphere(s3);
    scene.ajoutersphere(s4);
    scene.ajoutersphere(s5);
 
    scene.source = Vector(-10, 20, 40);
    scene.intensite = 100000000;

    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for 
    for (int i = 0; i < H; i++) {   // on parcourt l'image 
        for (int j = 0; j < W; j++) {
            Vector C(0,0,55);  //origine du vecteur vision
			Vector u(j-W/2+0.5, i-H/2+0.5, -W/(2*tan(fov/2.)));   // direction du vecteur vision
			u.normalize();
			Ray r(C,u);  //rayon de la vision

            Vector color(0,0,0);
            for (int k=0; k<10; k++){  //on envoie 5 rayons au lieu d'un seul 
                color = color + getColor(r, scene, 5)/5;
            };

			image[((H-i-1) * W + j) * 3 + 0] = fmin(255, (fmax(0,pow(color[0], 1.0/2.2))));  //coordonnée rouge
            image[((H-i-1) * W + j) * 3 + 1] = fmin(255,fmax(0,pow(color[1], 1.0/2.2)));  //coordonnée verte
            image[((H-i-1) * W + j) * 3 + 2] = fmin(255,fmax(0,pow(color[2],1.0/2.2)));  //coordonnée bleue 
			}
        }
    stbi_write_png("sphère.png", W, H, 3, &image[0], 0);
 
    return 0;
}