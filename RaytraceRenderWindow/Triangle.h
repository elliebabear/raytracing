#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Cartesian3.h"
#include "Homogeneous4.h"
#include "Material.h"
#include "Ray.h"
#include <cmath>


class Triangle
{
public:
    Homogeneous4 verts[3];
    Homogeneous4 normals[3];
    Homogeneous4 colors[3];
    Cartesian3 uvs[3];
    Material *shared_material;
    Triangle();
    float intersect(Ray r);
    Cartesian3 baricentric(Cartesian3 o);
    Homogeneous4 BlinnPhong(Homogeneous4 lpos, Homogeneous4 lc, Cartesian3 bc, Ray r);
    Homogeneous4 directLight(Homogeneous4 lpos, Homogeneous4 lc, Cartesian3 bc, Ray r);
};

#endif // TRIANGLE_H
