#ifndef SCENE_H
#define SCENE_H
#include "ThreeDModel.h"
#include "Triangle.h"


class Scene
{
public:
    std::vector<ThreeDModel>* objects;
    std::vector<Triangle> triangles;
    RenderParameters* rp;
    Material *default_mat;
    struct CollisionInfo
    {
        Triangle trl;
        float t;
    };

    Scene(std::vector<ThreeDModel> *texobjs,RenderParameters *renderp);
    void updateScene();
    Matrix4 getModelview();
    CollisionInfo closestTriangle(Ray r);
};

#endif // SCENE_H
