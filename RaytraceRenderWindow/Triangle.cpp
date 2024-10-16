#include "Triangle.h"

Triangle::Triangle()
{
    shared_material= nullptr;
}

float Triangle::intersect(Ray r)
{
    float result = -100.0f;
    const float epsilon = 0.00001;
    Cartesian3 e01 = verts[1].Point() - verts[0].Point();
    Cartesian3 e02 = verts[2].Point() - verts[0].Point();
    Cartesian3 s = r.origin - verts[0].Point();
    auto s1 = r.direction.cross(e02);
    auto s2 = s.cross(e01);
    auto se = 1.0f / s1.dot(e01);
    float t = s2.dot(e02) * se;
    float b1 = s1.dot(s) * se;
    float b2 = s2.dot(r.direction) * se;
    if(t+epsilon>0 && b1+epsilon>0 && b2+epsilon>0 && 1-b1-b2+epsilon>0)
        return t;
    return result;
}

Cartesian3 Triangle::baricentric(Cartesian3 o)
{//input intersection o = og + dir * t
    Cartesian3 bc;
    Cartesian3 e01 = verts[1].Point() - verts[0].Point();
    Cartesian3 e02 = verts[2].Point() - verts[0].Point();
    Cartesian3 st = e01.cross(e02);
    float S = 0.5f * st.length();

    Cartesian3 e0o = o - verts[0].Point();
    Cartesian3 st0 = e01.cross(e0o);
    float S0 = 0.5f * st0.length();

    Cartesian3 e12 = verts[2].Point() - verts[1].Point();
    Cartesian3 e1o = o - verts[1].Point();
    Cartesian3 st1 = e12.cross(e1o);
    float S1 = 0.5f * st1.length();

    bc.x = S1 / S;//*vo
    bc.z = S0 / S;//*v2
    bc.y = 1.0f - S1/S - S0/S;//*v1

    return bc;
}

Homogeneous4 Triangle::BlinnPhong(Homogeneous4 lpos, Homogeneous4 lc, Cartesian3 bc, Ray ray)
{
    Cartesian3 hitPoint = (verts[0].Point()*bc.x) + (verts[1].Point()*bc.y) + (verts[2].Point()*bc.z);
    Cartesian3 n = (normals[0].Vector()*bc.x) + (normals[1].Vector()*bc.y) + (normals[2].Vector()*bc.z);
    n = n.unit();
    Cartesian3 l = (lpos.Point() - hitPoint).unit();
    Cartesian3 v = (ray.origin - hitPoint).unit();

    //diffuse = k*I*max(0,n*l)
    Cartesian3 kd(lc.x*shared_material->diffuse.x , lc.y*shared_material->diffuse.y , lc.z*shared_material->diffuse.z);
    Cartesian3 diffuseRef = kd * std::max(0.0f,n.dot(l));

    //specular = k*I*max(0,n*h)^shinness    h - half vector
    Cartesian3 h = (v + l).unit();
    Cartesian3 ks(lc.x*shared_material->specular.x , lc.y*shared_material->specular.y , lc.z*shared_material->specular.z);
    Cartesian3 specularRef = 0.4f * ks * pow(std::max(0.0f,n.dot(h)),shared_material->shininess);

    //ambient = k*I
    Cartesian3 ambientRef(lc.x*shared_material->ambient.x , lc.y*shared_material->ambient.y , lc.z*shared_material->ambient.z);
    Cartesian3 Lightall = diffuseRef + specularRef + ambientRef + shared_material->emissive;
    //Cartesian3 Lightall = diffuseRef + specularRef + shared_material->emissive;
    Homogeneous4 color(Lightall);

    return color;
}

Homogeneous4 Triangle::directLight(Homogeneous4 lpos, Homogeneous4 lc, Cartesian3 bc, Ray ray)
{
    Cartesian3 hitPoint = (verts[0].Point()*bc.x) + (verts[1].Point()*bc.y) + (verts[2].Point()*bc.z);
    Cartesian3 n = (normals[0].Vector()*bc.x) + (normals[1].Vector()*bc.y) + (normals[2].Vector()*bc.z);
    n = n.unit();
    Cartesian3 l = (lpos.Point() - hitPoint).unit();
    Cartesian3 v = (ray.origin - hitPoint).unit();

    //diffuse = k*I*max(0,n*l)
    Cartesian3 kd(lc.x*shared_material->diffuse.x , lc.y*shared_material->diffuse.y , lc.z*shared_material->diffuse.z);
    Cartesian3 diffuseRef = kd * std::max(0.0f,n.dot(l));

    //specular = k*I*max(0,n*h)^shinness    h - half vector
    Cartesian3 h = (v + l).unit();
    Cartesian3 ks(lc.x*shared_material->specular.x , lc.y*shared_material->specular.y , lc.z*shared_material->specular.z);
    Cartesian3 specularRef = 0.4f * ks * pow(std::max(0.0f,n.dot(h)),shared_material->shininess);

    //ambient = k*I
    //Cartesian3 ambientRef(lc.x*shared_material->ambient.x , lc.y*shared_material->ambient.y , lc.z*shared_material->ambient.z);

    Cartesian3 Lightall = diffuseRef + specularRef + shared_material->emissive;
    Homogeneous4 color(Lightall);

    return color;
}


