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
#define _USE_MATH_DEFINES
#include <cmath>

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

std::random_device rd;
std::default_random_engine eng(rd());
std::uniform_real_distribution<float> distr(0, 1);

Vector randomcos(const Vector&N){
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
    return direction_aleatoire;
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

class Object { 
public:
	Object() {};
    virtual bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const = 0;
    Vector albedo;
    bool isMirror;
    bool isTransparent;
};

class Sphere : public Object{  //une origine et un rayon
public:
    Sphere(const Vector& origine, const double rayon, const Vector couleur, bool isMirror = false, bool isTransparent = false) : O(origine), R(rayon){
        albedo = couleur; 
        isMirror = isMirror; 
        isTransparent = isTransparent;
    };
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
};

class Triangle : public Object { 
public:
    Triangle(const Vector& A, const Vector& B, const Vector& C, const Vector couleur, bool isMirror = false, bool isTransparent = false) : A(A), B(B), C(C){
        albedo = couleur;
        isMirror = isMirror;
        isTransparent = isTransparent;
    };
        bool intersect(const Ray& r, Vector& P, Vector& Normale, double& t) const {
            Normale = cross(B-A, C-A);
            Normale.normalize();
            t = dot(C - r.origine, Normale)/dot(r.direction, Normale);
            P = r.origine + t*r.direction;
            if (t<0){
                return false;
            } else {
                // formule de Cramer pour calculer les coordonnées barycentriques
                Vector v0 = B-A;
                Vector v1 = C-A;
                Vector v2 = P-A;

                double dot00 = dot(v0, v0);
                double dot01 = dot(v0, v1);
                double dot02 = dot(v0, v2);
                double dot11 = dot(v1, v1);
                double dot12 = dot(v1, v2);

                double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
                double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
                double v = (dot00 * dot12 - dot01 * dot02) * invDenom;
                double w = 1 - u - v;
                // Check if point is in triangle
                if (u<0 || u>1){ return false; };
                if (v<0 || v>1){ return false; };
                if (w<0 || w>1){ return false; };
                return true;             
        }

    };
    Vector A;
    Vector B;
    Vector C;
};


class Scene {  // ensemble de sphères 
public: 
    Scene() {};
    void ajoutersphere(const Sphere& s) {objects.push_back((const Object*)&s);}
    void ajoutertriangle(const Triangle& t) {objects.push_back((const Object*)&t);}
    std::vector<const Object*> objects;
    Sphere *lumiere;
    double intensite;
    bool intersection(const Ray& r, Vector& P, Vector& N, int& sphere_inter_id) const {   //pour toutes les sphères
        //P la position, N la normale
        bool inter = false;
        double min_t = 10E10;
        for (int i=0; i<objects.size();i++){
            Vector Pprime, Nprime;
            double t;
            bool interprime = objects[i]->intersect(r,Pprime,Nprime,t);
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
    Vector P,N; //vecteurs position et normal
    int sphere_inter_id;
    Vector intensite_pixel(0,0,0);

	if (scene.intersection(r,P,N, sphere_inter_id) == true){  //une intersection est trouvée 
        if (scene.objects[sphere_inter_id]->isMirror && nombre_rebond>0){  //sphère miroir
            Vector direction_miroir = r.direction - 2*dot(N,r.direction)*N;
            Ray rayon_miroir(P+0.001*N, direction_miroir);
            intensite_pixel = getColor(rayon_miroir, scene, nombre_rebond -1);  
        };
        if (scene.objects[sphere_inter_id]->isTransparent && nombre_rebond>0){  //sphère transparente
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
        Vector axe = P-scene.lumiere->O;
        axe.normalize();
        Vector direction_aleatoire = randomcos(axe);
        Vector point_aleatoire = direction_aleatoire * scene.lumiere->R + scene.lumiere->O;
        Vector wi = (point_aleatoire - P);
        wi.normalize();
        double distance = (point_aleatoire - P).norm();

        Ray r2(P+0.005*N,wi);
        Vector P2,N2;
        int sphere_inter_id2;

        bool isInter = scene.intersection(r2,P2,N2, sphere_inter_id2);
        double distance2 = (scene.lumiere->O-P).norm();
        if (isInter && 0.99*distance > distance2){  //une intersection est trouvée 
            intensite_pixel = Vector(0,0,0);
        } else {
            intensite_pixel = (scene.intensite)/(4*M_PI*distance)*fmax(0,dot(N,wi))*dot(direction_aleatoire, -wi)*dot(axe, direction_aleatoire)*scene.objects[sphere_inter_id]->albedo;
        };

        // //CONTRIBUTION INDIRECTE
        Vector direction_aleatoire2 = randomcos(N);
        Ray rayon_aleatoire(P+0.001*N, direction_aleatoire2);
        intensite_pixel= intensite_pixel + getColor(rayon_aleatoire, scene, nombre_rebond -1)*scene.objects[sphere_inter_id]->albedo;  
    };
    return intensite_pixel;
};

int main() {
int W = 512;
    int H = 512;

    Sphere lumiere(Vector(-10, 40, 40),15, Vector(1,1,1));
	// Sphere s(Vector(-12,0,0),10, Vector(1,0,1), false, false);   //spère bleue 
	// Sphere sbis(Vector(12,0,0),10, Vector(0,0,1), false, false);
	double fov=60*M_PI/180.;   //champ visuel 
	double tanfov2 = tan(fov/2);

    Sphere s1(Vector(0,-2000-10,0), 2000, Vector(0,1,1));  //sol blanc, de telle sorte à ce que la sphre soit posée sur le sol
    Sphere s2(Vector(0,2000+25,0), 2000, Vector(0,0,1));  //plafond
    Sphere s3(Vector(-2000-25,0,0), 2000, Vector(1,1,1));  //mur à gauche 
    Sphere s4(Vector(2000+25,0,0), 2000, Vector(1,1,1));  //mur à droite 
    Sphere s5(Vector(0,0,-2000-25), 2000, Vector(1,1,1));  //mur au fond 
    
    Triangle triangle(Vector(-12, 0, 0), Vector(12, 0, 0),Vector(0, 7, 0), Vector(1,0,0));

    Scene scene;
    scene.ajoutersphere(lumiere);
    //scene.ajoutersphere(s);
    //scene.ajoutersphere(sbis);
    scene.ajoutersphere(s1);
    scene.ajoutersphere(s2);
    scene.ajoutersphere(s3);
    scene.ajoutersphere(s4);
    scene.ajoutersphere(s5);
    scene.ajoutertriangle(triangle);
 
    scene.lumiere = &lumiere;
    scene.intensite = 100000000;
    Vector position_camera(0,0,55);  //origine du vecteur vision
    double focus = 55;  // tout ce qui est avant ou après cette distance là sera plus floue
    int nbrayons = 2; 

    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for 
    for (int i = 0; i < H; i++) {   // on parcourt l'image 
        for (int j = 0; j < W; j++) {
            Vector color(0,0,0);
            for (int k=0; k<nbrayons ; k++){  //on envoie 5 rayons au lieu d'un seul 
                double depth = H/(2*tan(fov*0.5));
                double r1 = drand48();
                double r2 = drand48();
                double x = sqrt(-2*log(r1))*cos(2*M_PI*r2)*0.5;
                double y = sqrt(-2*log(r1))*sin(2*M_PI*r2)*0.5;

                double dx_profondeur = (drand48()-0.5)*5;
                double dy_profondeur = (drand48()-0.5)*5;

                Vector u(j-W/2+0.5+x, i-H/2+0.5+y, -W/(2*tan(fov/2.)));   // direction du vecteur vision
                u.normalize();

                Vector mise_au_point = position_camera + focus*u;
                Vector nouvelle_origine = position_camera + Vector(dx_profondeur,dy_profondeur,0);

                Vector a = mise_au_point - nouvelle_origine;
                a.normalize();
                Ray r(nouvelle_origine,a);  //rayon de la vision
                color = color + getColor(r, scene, 5)/nbrayons;
            };

			image[((H-i-1) * W + j) * 3 + 0] = fmin(255, (fmax(0,pow(color[0], 1.0/2.2))));  //coordonnée rouge
            image[((H-i-1) * W + j) * 3 + 1] = fmin(255,fmax(0,pow(color[1], 1.0/2.2)));  //coordonnée verte
            image[((H-i-1) * W + j) * 3 + 2] = fmin(255,fmax(0,pow(color[2],1.0/2.2)));  //coordonnée bleue 
			}
        }
    stbi_write_png("sphère.png", W, H, 3, &image[0], 0);
 
    return 0;
};