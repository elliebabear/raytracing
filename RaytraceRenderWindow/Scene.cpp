#include "Scene.h"

void Scene::updateScene()
{
    triangles.clear(); //Clear the list so it can be populated again
    for (int i = 0;i< int(objects->size());i++)
    {
        typedef unsigned int uint;
        ThreeDModel obj = objects->at(uint(i));
        for (uint face = 0; face < obj.faceVertices.size(); face++)
        {
            for (uint triangle = 0; triangle < obj.faceVertices[face].size()-2; triangle++)
            {
                Triangle t;
                for (uint vertex = 0; vertex < 3; vertex++)
                {
                    uint faceVertex = 0;
                    if (vertex != 0)
                        faceVertex = triangle + vertex;
                    //this is our vertex before any transformations. (world space)
                    Homogeneous4 v = Homogeneous4(obj.vertices[obj.faceVertices [face][faceVertex]].x,
                                                  obj.vertices[obj.faceVertices [face][faceVertex]].y,
                                                  obj.vertices[obj.faceVertices [face][faceVertex]].z);
                    //apply modelview
                    v = getModelview() * v;
                    t.verts[vertex] = v;
                    Homogeneous4 n = Homogeneous4(obj.normals[obj.faceNormals [face][faceVertex]].x,
                                                  obj.normals[obj.faceNormals [face][faceVertex]].y,
                                                  obj.normals[obj.faceNormals [face][faceVertex]].z,0.0f);
                    n = getModelview() * n;
                    t.normals[vertex] = n;
                    Cartesian3 tex = Cartesian3(obj.textureCoords[obj.faceTexCoords[face][faceVertex]].x,
                                                obj.textureCoords[obj.faceTexCoords[face][faceVertex]].y,0.0f);
                    t.uvs[vertex] = tex;
                    t.colors[vertex] = Cartesian3( 0.7f, 0.7f, 0.7f);
                }
                if(obj.material== nullptr)
                {
                    t.shared_material = default_mat;
                }else{
                    t.shared_material = obj.material;
                }
                triangles.push_back(t);
            }
        }
    }
}

Scene::Scene(std::vector<ThreeDModel> *texobjs,RenderParameters *renderp)
{
    objects = texobjs;
    rp = renderp;
    Cartesian3 ambient = Cartesian3(0.5f,0.5f,0.5f);
    Cartesian3 diffuse = Cartesian3(0.5f,0.5f,0.5f);
    Cartesian3 specular = Cartesian3(0.5f,0.5f,0.5f);
    Cartesian3 emissive = Cartesian3(0,0,0);
    float shininess = 1.0f;
    default_mat = new Material(ambient,diffuse,specular,emissive,shininess);
}

Matrix4 Scene::getModelview()
{
    Matrix4 result;
    result.SetIdentity();//fill in the diagonal with 1
    result.SetTranslation({rp->xTranslate, rp->yTranslate, rp->zTranslate});

    result =  result * rp->rotationMatrix;
    return result;
}

Scene::CollisionInfo Scene::closestTriangle(Ray r)
{
    Scene::CollisionInfo ci;
    ci.t = -1000;
    for(Triangle tri:triangles)
    {
        const float epsilon = 0.00001;
        float temp = tri.intersect(r);
        if(temp + epsilon> 0)
        {
            if(ci.t <  0)
            {
                ci.trl = tri;
                ci.t = temp;
            }
            if(ci.t > 0 && temp < ci.t)
            {
                ci.trl = tri;
                ci.t = temp;
            }
        }
    }
    return ci;
}








