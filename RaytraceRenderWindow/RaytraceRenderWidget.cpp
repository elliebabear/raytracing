//////////////////////////////////////////////////////////////////////
//
//  University of Leeds
//  COMP 5812M Foundations of Modelling & Rendering
//  User Interface for Coursework
////////////////////////////////////////////////////////////////////////


#include <math.h>
#include <random>
// include the header file
#include "RaytraceRenderWidget.h"

#define N_THREADS 16
#define N_LOOPS 100
#define N_BOUNCES 5
#define TERMINATION_FACTOR 0.35f

// constructor
RaytraceRenderWidget::RaytraceRenderWidget
        (   
        // the geometric object to show
        std::vector<ThreeDModel>      *newTexturedObject,
        // the render parameters to use
        RenderParameters    *newRenderParameters,
        // parent widget in visual hierarchy
        QWidget             *parent
        )
    // the : indicates variable instantiation rather than arbitrary code
    // it is considered good style to use it where possible
    : 
    // start by calling inherited constructor with parent widget's pointer
    QOpenGLWidget(parent),
    // then store the pointers that were passed in
    texturedObjects(newTexturedObject),
    renderParameters(newRenderParameters)
    { // constructor
        scene = new Scene(texturedObjects,renderParameters);
        QTimer *timer = new QTimer(this);
        connect(timer, &QTimer::timeout, this, &RaytraceRenderWidget::forceRepaint);
        timer->start(30);
    } // constructor


// destructor
RaytraceRenderWidget::~RaytraceRenderWidget()
    { // destructor
    // empty (for now)
    // all of our pointers are to data owned by another class
    // so we have no responsibility for destruction
    // and OpenGL cleanup is taken care of by Qt
    } // destructor                                                                 

// mouse-handling
void RaytraceRenderWidget::mousePressEvent(QMouseEvent *event)
    { // RaytraceRenderWidget::mousePressEvent()
    // store the button for future reference
    int whichButton = int(event->button());
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0f * event->x() - size) / size;
    float y = (size - 2.0f * event->y() ) / size;

    
    // and we want to force mouse buttons to allow shift-click to be the same as right-click
    unsigned int modifiers = event->modifiers();
    
    // shift-click (any) counts as right click
    if (modifiers & Qt::ShiftModifier)
        whichButton = Qt::RightButton;
    
    // send signal to the controller for detailed processing
    emit BeginScaledDrag(whichButton, x,y);
    } // RaytraceRenderWidget::mousePressEvent()
    
void RaytraceRenderWidget::mouseMoveEvent(QMouseEvent *event)
    { // RaytraceRenderWidget::mouseMoveEvent()
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0f * event->x() - size) / size;
    float y = (size - 2.0f * event->y() ) / size;
    
    // send signal to the controller for detailed processing
    emit ContinueScaledDrag(x,y);
    } // RaytraceRenderWidget::mouseMoveEvent()
    
void RaytraceRenderWidget::mouseReleaseEvent(QMouseEvent *event)
    { // RaytraceRenderWidget::mouseReleaseEvent()
    // scale the event to the nominal unit sphere in the widget:
    // find the minimum of height & width   
    float size = (width() > height()) ? height() : width();
    // scale both coordinates from that
    float x = (2.0f * event->x() - size) / size;
    float y = (size - 2.0f * event->y() ) / size;
    
    // send signal to the controller for detailed processing
    emit EndScaledDrag(x,y);
    } // RaytraceRenderWidget::mouseReleaseEvent()

// called when OpenGL context is set up
void RaytraceRenderWidget::initializeGL()
    { // RaytraceRenderWidget::initializeGL()
	// this should remain empty
    } // RaytraceRenderWidget::initializeGL()

// called every time the widget is resized
void RaytraceRenderWidget::resizeGL(int w, int h)
    { // RaytraceRenderWidget::resizeGL()
    // resize the render image
    frameBuffer.Resize(w, h);
    } // RaytraceRenderWidget::resizeGL()
    
// called every time the widget needs painting
void RaytraceRenderWidget::paintGL()
    { // RaytraceRenderWidget::paintGL()
    // set background colour to white
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    // and display the image
    glDrawPixels(frameBuffer.width, frameBuffer.height, GL_RGBA, GL_UNSIGNED_BYTE, frameBuffer.block);
    } // RaytraceRenderWidget::paintGL()

//1.1
void RaytraceRenderWidget::Raytrace()
{
    scene->updateScene();
    raytracingThread = std::thread(&RaytraceRenderWidget::RaytraceThread,this);
    raytracingThread.detach();
}

void RaytraceRenderWidget::RaytraceThread()
{//1.1
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
    int sampleNum = 1;
    if(renderParameters->monteCarloEnabled) sampleNum = 5;
    frameBuffer.clear(RGBAValue(0.0f,0.0f,0.0f,1.0f));
    #pragma omp parallel for schedule(dynamic)
    for(int j = 0; j < frameBuffer.height; j++)
    {
        for(int i = 0; i < frameBuffer.width; i++)
        {
            Homogeneous4 color(0,0,0,0);
            for (int k = 0; k < sampleNum; ++k) {
                float x = i+distribution(generator);
                float y = j+distribution(generator);
                Ray camRay = calculateRay(x,y,renderParameters->orthoProjection);
                Homogeneous4 sampleColor = calculateLight(camRay,7,1.0f);
                color = color + sampleColor;
            }
            color = color / sampleNum;
            float gamma = 2.2f;
            color.x = pow(color.x,1/gamma);
            color.y = pow(color.y,1/gamma);
            color.z = pow(color.z,1/gamma);
            frameBuffer[j][i]=RGBAValue(color.x*255.0f,color.y*255.0f,color.z*255.0f,255.0f);
        }
    }
    std::cout<< "Done!" <<std::endl;
}

Homogeneous4 RaytraceRenderWidget::calculateLight(Ray ray, int depth, float nowIOR)
{
    Homogeneous4 reflectcolor(0,0,0,0);
    Homogeneous4 refractcolor(0,0,0,0);
    Homogeneous4 result(0,0,0,0);
    if(depth==0) return result;
    Scene::CollisionInfo hit = scene->closestTriangle(ray);
    const float epsilon = 0.00001;
    if(hit.t >= 0)
    {
        Cartesian3 o = ray.origin + ray.direction * hit.t;
        Cartesian3 BariC = hit.trl.baricentric(o);
        Cartesian3 hitNor = (hit.trl.normals[0].Vector() * BariC.x) +
                            (hit.trl.normals[1].Vector() * BariC.y) +
                            (hit.trl.normals[2].Vector() * BariC.z);
        hitNor = hitNor.unit();
        Cartesian3 hitPoint = (hit.trl.verts[0].Point()* BariC.x) +
                              (hit.trl.verts[1].Point()* BariC.y) +
                              (hit.trl.verts[2].Point()* BariC.z);
        hitPoint = hitPoint + hitNor * 0.00001f;
        if(renderParameters->phongEnabled)
        {
            for(Light *l:renderParameters->lights)
            {
                Homogeneous4 lpos = scene->getModelview() * l->GetPositionCenter();
                Homogeneous4 lc = l->GetColor();
                if(renderParameters->shadowsEnabled)
                {
                    Cartesian3 tritolight = (lpos.Point() - hitPoint).unit();
                    Ray shadowray(hitPoint,tritolight);
                    Scene::CollisionInfo shadowhit = scene->closestTriangle(shadowray);
                    float tl = (lpos.Point().x - shadowray.origin.x) / shadowray.direction.x;
                    float d = -(hitNor.x*hitPoint.x + hitNor.y*hitPoint.y + hitNor.z*hitPoint.z);
                    float distance = (hitNor.x*lpos.Point().x + hitNor.y*lpos.Point().y
                                      + hitNor.z*lpos.Point().z + d)/ sqrt(hitNor.x*hitNor.x + hitNor.y*hitNor.y+ hitNor.z*hitNor.z);
                    if( shadowhit.t>=epsilon && shadowhit.t+2.0f*epsilon<tl && distance>2.0f*epsilon){
                        Cartesian3 ambientRef(lc.x*hit.trl.shared_material->ambient.x ,
                                              lc.y*hit.trl.shared_material->ambient.y ,
                                              lc.z*hit.trl.shared_material->ambient.z);
                        result = result+hit.trl.shared_material->emissive + ambientRef;
//                        if(renderParameters->monteCarloEnabled){
//                            result = result+hit.trl.shared_material->emissive;
//                        }else{
//                            result = result+hit.trl.shared_material->emissive + ambientRef;
//                        }
                    }
                    else{
                        if(renderParameters->monteCarloEnabled){
                            result = result+hit.trl.directLight(lpos, lc, BariC, ray);
                        }else{
                            result = result+hit.trl.BlinnPhong(lpos, lc, BariC, ray);
                        }
                        //result = result+hit.trl.BlinnPhong(lpos, lc, BariC, ray);
                    }
                }else{
                    if(renderParameters->monteCarloEnabled){
                        result = result+hit.trl.directLight(lpos, lc, BariC, ray);
                    }else{
                        result = result+hit.trl.BlinnPhong(lpos, lc, BariC, ray);
                    }
                    //result = result+hit.trl.BlinnPhong(lpos, lc, BariC, ray);
                }
            }

            if(renderParameters->monteCarloEnabled && !hit.trl.shared_material->isLight()){
                Ray sampleRay = randomRay(hitNor, hitPoint);
                int newDepth = depth - 1;
                Homogeneous4 sampleColor = calculateLight(sampleRay,newDepth,nowIOR);
                Cartesian3 indirectColor(sampleColor.x*hit.trl.shared_material->ambient.x ,
                                         sampleColor.y*hit.trl.shared_material->ambient.y ,
                                         sampleColor.z*hit.trl.shared_material->ambient.z);
                result = result+indirectColor;
            }

            if(renderParameters->refractionEnabled && renderParameters->fresnelRendering
                && renderParameters->reflectionEnabled && hit.trl.shared_material->transparency>0.0f){
                Ray mirrorRay = reflectRay(ray,hitNor,hitPoint);
                Ray refractRay = refract(ray,hitNor,hitPoint,nowIOR,hit.trl.shared_material->indexOfRefraction);
                int newDepth = depth - 1;
                reflectcolor =calculateLight(mirrorRay, newDepth, nowIOR);
                refractcolor = calculateLight(refractRay, newDepth, hit.trl.shared_material->indexOfRefraction);
                float kr = fresnel(ray,hitNor,nowIOR,hit.trl.shared_material->indexOfRefraction);
                result = kr * hit.trl.shared_material->reflectivity*reflectcolor + (1.0f-hit.trl.shared_material->reflectivity)*result
                         + (1-kr)*hit.trl.shared_material->transparency*refractcolor + (1-kr)*(1.0f-hit.trl.shared_material->transparency)*result;
            }
            if(renderParameters->refractionEnabled && !renderParameters->fresnelRendering && hit.trl.shared_material->transparency>0.01f){
                Ray refractRay = refract(ray,hitNor,hitPoint,nowIOR,hit.trl.shared_material->indexOfRefraction);
                int newDepth = depth - 1;
                refractcolor = calculateLight(refractRay, newDepth, hit.trl.shared_material->indexOfRefraction);
                result =hit.trl.shared_material->transparency*refractcolor + (1.0f-hit.trl.shared_material->transparency)*result;
                if(renderParameters->fresnelRendering){
                    float kr = fresnel(ray,hitNor,nowIOR,hit.trl.shared_material->indexOfRefraction);
                    result = (1.0f-kr)*hit.trl.shared_material->transparency*refractcolor + (1.0f-hit.trl.shared_material->transparency)*result;
                }
            }
            if(renderParameters->reflectionEnabled && !renderParameters->fresnelRendering && hit.trl.shared_material->reflectivity>0.0f){
                Ray mirrorRay = reflectRay(ray,hitNor,hitPoint);
                int newDepth = depth - 1;
                reflectcolor =calculateLight(mirrorRay, newDepth, nowIOR);
                result = hit.trl.shared_material->reflectivity*reflectcolor + (1.0f-hit.trl.shared_material->reflectivity)*result;
                if(renderParameters->fresnelRendering){
                    float kr = fresnel(ray,hitNor,nowIOR,hit.trl.shared_material->indexOfRefraction);
                    result = kr * hit.trl.shared_material->reflectivity*reflectcolor + (1.0f-hit.trl.shared_material->reflectivity)*result;
                }
            }
        }
        if(renderParameters->interpolationRendering)
        {
            result = Homogeneous4(abs(hitNor.x),abs(hitNor.y),abs(hitNor.z),1);
        }
    }
    return result;
}

Ray RaytraceRenderWidget::reflectRay(Ray ray, Cartesian3 normal, Cartesian3 hitpoint)
{
    Cartesian3 dir = (ray.direction.unit() - 2.0f*ray.direction.dot(normal)*normal).unit();
    hitpoint = hitpoint + normal * 0.00001f;
    Ray result(hitpoint,dir);
    return result;
}

Ray RaytraceRenderWidget::refract(Ray ray, Cartesian3 normal, Cartesian3 hitpoint, float nowIOR, float newIOR)
{
    float angle = ray.direction.unit().dot(normal.unit());
    if(angle>1.0f)   angle=1.0f;
    if(angle<-1.0f)  angle=-1.0f;
    float anglei = nowIOR;
    float anglet = newIOR;
    Cartesian3 n = normal;
    if(angle<0){
        angle = -angle;
    }else{
        n= -1.0f*normal;
    }
    float a = anglei / anglet;
    float k = 1.0f - a * a * (1.0f - angle * angle);
    if(k>=0.0f) {
        Cartesian3 dir=(a*ray.direction.unit() + (a*angle-sqrtf(k))* n.unit()).unit();
        hitpoint = hitpoint - n * 0.00005f;
        Ray result(hitpoint,dir);
        return result;
    }else{
        Ray result(hitpoint,Cartesian3(0.0f,0.0f,0.0f));
        return result;
    }
}

Ray RaytraceRenderWidget::randomRay(Cartesian3 normal, Cartesian3 hitpoint)
{
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
    float r1 = distribution(generator);
    float r2 = distribution(generator);
    Cartesian3 a, b;
    if (std::fabs(normal.x) > std::fabs(normal.y))
        a = Cartesian3(normal.z, 0, -normal.x) / sqrtf(normal.x * normal.x + normal.z * normal.z);
    else
        a = Cartesian3(0, -normal.z, normal.y) / sqrtf(normal.y * normal.y + normal.z * normal.z);
    b = normal.cross(a);
    float sin = sqrtf(1 - r1 * r1);
    float p = 2 * M_PI * r2;
    float x = sin * cosf(p);
    float z = sin * sinf(p);
    Cartesian3 sample = Cartesian3(x, r1, z);
    Cartesian3 dir = Cartesian3(
        sample.x * b.x + sample.y * normal.x + sample.z * a.x,
        sample.x * b.y + sample.y * normal.y + sample.z * a.y,
        sample.x * b.z + sample.y * normal.z + sample.z * a.z).unit();
    Ray result(hitpoint,dir);
    return result;
}

float RaytraceRenderWidget::fresnel(Ray ray, Cartesian3 normal, float nowIOR, float newIOR)
{
    float angle = ray.direction.unit().dot(normal.unit());
    if(angle<0){angle = -angle;}
    if(angle>1.0f)   angle=1.0f;
    if(angle<-1.0f)  angle=-1.0f;
    float anglei = nowIOR;
    float anglet = newIOR;
    float a = (anglei-anglet)/(anglei+anglet);
    float kr = a*a+(1.0f-a*a)*(1.0f-angle)*(1.0f-angle)*(1.0f-angle)*(1.0f-angle)*(1.0f-angle);
    return kr;
}

Ray RaytraceRenderWidget::calculateRay(int pixelx, int pixely, bool orthoP)
{
    float aspectRatio = float(frameBuffer.width) / float(frameBuffer.height);
    float x = (float(pixelx) / float(frameBuffer.width) - 0.5f) * 2;
    float y = (float(pixely) / float(frameBuffer.height) - 0.5f) * 2;
    Cartesian3 og;
    Cartesian3 dir;

    if(aspectRatio > 1.0)
    {
        x = x * aspectRatio;
    }
    else
    {
        y = y / aspectRatio;
    }

    if(orthoP)
    {
        og.x = x; og.y = y; og.z = 1;
        dir.x = 0.0; dir.y = 0.0; dir.z = -1;
    }
    else
    {
        og.x = 0; og.y = 0; og.z = 1;
        dir.x = x; dir.y = y; dir.z = -1;
    }
    dir = dir.unit();
    Ray r(og,dir);
    return r;
}

void RaytraceRenderWidget::forceRepaint()
{//1.1
    update();
}
