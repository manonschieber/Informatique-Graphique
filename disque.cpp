#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
 
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
 
class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
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
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {   //produit scalaire 
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {   //produit vectoriel 
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
 
class Ray {  //une origine et une direction 
	public:
		Ray(const Vector& origine, const Vector& direction) : C(origine), u(direction){};
		Vector C;
		Vector u;
    	// ...
};

class Sphere {  //une origine et un rayon
public:
	Sphere(const Vector& O, const double R) : O(O), R(R) {};
	bool intersect(const Ray& r) const {
		double a = 1;
		double b = 2*dot(r.u,r.C-O);
		double c = (r.C - O).norm2()-R*R;
		double delta = b*b - 4*a*c;
		return delta>=0; 
	};
	Vector O;
	double R;
    // ...
};
 
int main() {
    int W = 512;
    int H = 512;
 
	Sphere s(Vector(0,0,10),1);
	double fov=60*M_PI/180.;   //champ visuel 
	double tanfov2 = tan(fov/2);

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
			Vector u(j-W/2+0.5, i-H/2+0.5, -W/(2*tan(fov/2.)));
			u.normalize();
			Vector C(0,0,0);
			Ray r(C,u);
			if (s.intersect(r) == true){
				image[(i * W + j) * 3 + 0] = 255;
            	image[(i * W + j) * 3 + 1] = 255;
            	image[(i * W + j) * 3 + 2] = 255;
			}
			else {
				image[(i * W + j) * 3 + 0] = 0;
            	image[(i * W + j) * 3 + 1] = 0;
            	image[(i * W + j) * 3 + 2] = 0;
			}
        }
    }
    stbi_write_png("disque.png", W, H, 3, &image[0], 0);
 
    return 0;
}