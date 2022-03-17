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
#include <string>
#include <map>

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
            double x = dot(r.direction, Normale);
            if(x==0){return false;}
            t = dot(C - r.origine, Normale)/x;
            P = r.origine + t*r.direction;
            if (t<0){ return false;}
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

    };
    const Vector A;
    const Vector B;
    const Vector C;
};

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
    int uvi, uvj, uvk;  // indices within the uv coordinates array
    int ni, nj, nk;  // indices within the normals array
    int group;       // face group
};

class Box {
public: 
    Box(){};
    Box(const Vector& bmin, const Vector& bmax): bmin(bmin), bmax(bmax) {};
    Vector bmin, bmax;
    bool intersect(const Ray& r) const {
        // on va tester les trois plans possibles de cette box pour trouver les intersections
        double t1x = (bmin[0] - r.origine[0])/r.direction[0];
        double t2x = (bmax[0] - r.origine[0])/r.direction[0];
        double tminx = std::min(t1x,t2x);
        double tmaxx = std::max(t1x,t2x);

        double t1y = (bmin[1] - r.origine[1])/r.direction[1];
        double t2y = (bmax[1] - r.origine[1])/r.direction[1];
        double tminy = std::min(t1y,t2y);
        double tmaxy = std::max(t1y,t2y);

        double t1z = (bmin[2] - r.origine[2])/r.direction[2];
        double t2z = (bmax[2] - r.origine[2])/r.direction[2];
        double tminz = std::min(t1z,t2z);
        double tmaxz = std::max(t1z,t2z);

        if(std::min(std::min(tmaxx, tmaxy), tmaxz)- std::max(std::max(tminx,tminy),tminz)>0){
            return true;
        } else{
            return false;
        }
    };
};

class TriangleMesh : public Object {
public:
  ~TriangleMesh() {}
    TriangleMesh() {};
    TriangleMesh(const char* obj, double scaling, const Vector& offset, const Vector& couleur, bool isMirror = false, bool isTransparent = false) {
        albedo = couleur;
        isMirror = isMirror;
        isTransparent = isTransparent;
        readOBJ(obj);
        for (int i = 0; i < vertices.size(); i++) {
            vertices[i] = vertices[i] * scaling + offset;
        }
    };
    
    void readOBJ(const char* obj) {
 
        char matfile[255];
        char grp[255];
 
        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;
 
            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());
 
            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }
 
            if (line[0] == 'v' && line[1] == ' ') {
                Vector vec;
 
                Vector col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));
 
                    vertices.push_back(vec);
                    vertexcolors.push_back(col);
 
                } else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;
 
                char* consumedline = line + 1;
                int offset;
 
                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        } else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }
 
                consumedline = consumedline + offset;
 
                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    } else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;                             
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            } else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                } else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }
 
            }
 
        }
        fclose(f);

        bb.bmax = vertices[0];
        bb.bmin = vertices[0];
        for (int i=1; i<vertices.size();i++){
            for (int j=0; j<3;j++){
                bb.bmin[j]  = std::min(bb.bmin[j], vertices[i][j]);
                bb.bmax[j]  = std::max(bb.bmax[j], vertices[i][j]);
            }
        }
    };
 
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    
    bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const{
        if (!bb.intersect(r)) return false;

        t = 1e99;
        bool has_inter = false;
        for(int i=0; i<indices.size(); i++){
            int i0 = indices[i].vtxi;
            int i1 = indices[i].vtxj;
            int i2 = indices[i].vtxk;
            Vector x = vertices[i0];
            Vector y = vertices[i1];
            Vector z = vertices[i2];

            Triangle triangle(vertices[i0],vertices[i1],vertices[i2], albedo, isMirror, isTransparent);
            Vector localP, localN;
            double localt;
            if (triangle.intersect(r,localP, localN, localt)){
                has_inter = true;
                if (localt < t){
                    t = localt;
                    P = localP;
                    N = localN;
                }
            }
        }
        return has_inter;
    };

private:
    Box bb;
};

class Sphere : public Object{  //une origine et un rayon
public:
    Sphere(const Vector& origine, const double rayon, const Vector couleur, bool isMirror = false, bool isTransparent = false) : O(origine), R(rayon){
        albedo = couleur; 
        this -> isMirror = isMirror; 
        this -> isTransparent = isTransparent;
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

class Scene {  // ensemble de sphères 
public: 
    Scene() {};
    void ajoutersphere(const Sphere& s) {objects.push_back((const Object*)&s);}
    void ajoutertriangle(const Triangle& t) {objects.push_back((const Object*)&t);}
    void ajoutergeometry(const TriangleMesh& g) {objects.push_back((const Object*)&g);}

    std::vector<const Object*> objects;
    Sphere *lumiere;
    double intensite;

    bool intersection(const Ray& r, Vector& P, Vector& N, int& sphere_inter_id) const {   //gestion des intersections dans le cas où l'on a plusieurs sphères
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

Vector getColor(const Ray &r, const Scene &scene, int nombre_rebond){   //fonction qui renvoie la couleur du pixel 
    Vector P,N; //vecteurs position et normal
    int sphere_inter_id;
    Vector intensite_pixel(0,0,0);
    if (nombre_rebond == 0){return Vector(0,0,0);} //condition d'arrêt pour éviter une récursion infinie
    if (nombre_rebond < 0){return Vector(0,0,0);}
    
    bool a = scene.intersection(r,P,N, sphere_inter_id);
	if (a == true){  //une intersection est trouvée 
        if (scene.objects[sphere_inter_id]->isMirror){  //surface miroir ou spéculaire
            Vector direction_miroir = r.direction - 2*dot(N,r.direction)*N;
            Ray rayon_miroir(P+0.001*N, direction_miroir);
            intensite_pixel = getColor(rayon_miroir, scene, nombre_rebond -1);  
        } else {
            if (scene.objects[sphere_inter_id]->isTransparent){  //surface transparente
                double indice1=1;   //air
                double indice2 = 1.3;   //verre
                Vector normale_transparence(N);
                if (dot(r.direction, N) >0 ){   //on est en train de sortir de la sphère 
                    indice1 = 1.3;
                    indice2=1;
                    normale_transparence = -N;
                } else{
                    indice1 = 1;
                    indice2=1.3;
                    normale_transparence = N;
                }
                double nombre_racine = 1 - std::pow(indice1/indice2,2)*(1-std::pow(dot(r.direction,normale_transparence),2));
                if (nombre_racine>0){
                    Vector direction_refracte = (indice1/indice2)*(r.direction -dot(r.direction, normale_transparence)*normale_transparence)-normale_transparence*sqrt(nombre_racine);
                    Ray rayon_refracte(P-0.01*normale_transparence, direction_refracte);
                    intensite_pixel = getColor(rayon_refracte, scene, nombre_rebond -1);  //un rebond de moins    
                ;}
            } else{
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
                    intensite_pixel = std::pow(scene.intensite,2.2)/(4*M_PI*distance)*fmax(0,dot(N,wi))*dot(direction_aleatoire, -wi)*dot(axe, direction_aleatoire)*scene.objects[sphere_inter_id]->albedo;
                };
            }
        }

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

    double fov=60*M_PI/180.;   //champ visuel 
	double tanfov2 = tan(fov/2);
    Vector position_camera(0,0,55);  //origine du vecteur vision
    double focus = 55;  // tout ce qui est avant ou après cette distance là sera plus floue
    int nbrayons = 10; 
    int nbrebonds = 4;

    Sphere lumiere(Vector(-10, 40, 40),15, Vector(1,1,1));
	Sphere s(Vector(-12,0,0),10, Vector(0,0,1),false,false);   //spère rouge
	Sphere sbis(Vector(12,0,0),10, Vector(0,1,0),false,false);  //sphère verte

    Sphere s1(Vector(0,-2000-10,0), 2000, Vector(1,0,0));  //sol rouge, la sphère est posée dessus
    Sphere s2(Vector(0,2000+25,0), 2000, Vector(0,0,1));  //plafond bleu
    Sphere s3(Vector(-2000-25,0,0), 2000, Vector(1,0,1));  //mur à gauche rose
    Sphere s4(Vector(2000+25,0,0), 2000, Vector(1,0,1));  //mur à droite rose 
    Sphere s5(Vector(0,0,-2000-25), 2000, Vector(1,0,1));  //mur au fond rose
     

    //TriangleMesh chien("dog/13463_Australian_Cattle_Dog_v3.obj", 1, Vector(0,0,0), Vector(0,1,0));

    Scene scene;
    scene.ajoutersphere(lumiere);
    //scene.ajoutergeometry(chien);
    scene.ajoutersphere(s);
    scene.ajoutersphere(sbis);
    scene.ajoutersphere(s1);
    scene.ajoutersphere(s2);
    scene.ajoutersphere(s3);
    scene.ajoutersphere(s4);
    scene.ajoutersphere(s5);
 
    scene.lumiere = &lumiere;
    scene.intensite = 5000;

    std::vector<unsigned char> image(W * H * 3, 0);
    #pragma omp parallel for 
    for (int i = 0; i < H; i++) {   // on parcourt l'image 
        for (int j = 0; j < W; j++) {
            Vector color(0,0,0);
            for (int k=0; k<nbrayons ; k++){  //on envoie nbrayons 
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
                color = color + getColor(r, scene, nbrebonds)/nbrayons;
            };

			image[((H-i-1) * W + j) * 3 + 0] = fmin(255, (fmax(0,pow(color[0], 1.0/2.2))));  //coordonnée rouge
            image[((H-i-1) * W + j) * 3 + 1] = fmin(255,fmax(0,pow(color[1], 1.0/2.2)));  //coordonnée verte
            image[((H-i-1) * W + j) * 3 + 2] = fmin(255,fmax(0,pow(color[2],1.0/2.2)));  //coordonnée bleue 
			}
        }
    stbi_write_png("sphère.png", W, H, 3, &image[0], 0);
 
    return 0;
};