#ifndef SQUARE_H
#define SQUARE_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySquareIntersection{
    bool intersectionExists;
    float t;
    float u,v;
    Vec3 intersection;
    Vec3 normal;
    RaySquareIntersection() : intersectionExists(false), t(FLT_MAX) {}
};


class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    Square() : Mesh() {}
    Square(Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
           float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void setQuad( Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
                  float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector , upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector*width;
        m_up_vector = m_up_vector*height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;                                      vertices[0].u = uMin; vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;                     vertices[1].u = uMax; vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;       vertices[2].u = uMax; vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;                        vertices[3].u = uMin; vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;


    }

    RaySquareIntersection intersect(const Ray &ray) const {
        RaySquareIntersection intersection;

        Vec3 a = vertices[0].position;
        //2 points quelconques V1 et V2
        Vec3 v1 = vertices[1].position - vertices[0].position, v2 = vertices[3].position - vertices[0].position;

        Vec3 ori = ray.origin();
        Vec3 dir = ray.direction();
        Vec3 normale = Vec3::cross(v1, v2);
        normale.normalize();
        float Distance = Vec3::dot(a,normale);

        //calcule de la solution t
        float t = (Distance - Vec3::dot(ori, normale))/ Vec3::dot(dir, normale);
        Vec3 p = ori + t * dir;

        if(t > 0){
            Vec3 p0 = Vec3::cross(vertices[1].position - vertices[0].position, p - vertices[0].position);
            Vec3 p1 = Vec3::cross(vertices[2].position - vertices[1].position, p - vertices[1].position);
            Vec3 p2 = Vec3::cross(vertices[3].position - vertices[2].position, p - vertices[2].position);
            Vec3 p3 = Vec3::cross(vertices[0].position - vertices[3].position, p - vertices[3].position);
            if((Vec3::dot(p0, normale) > 0) == (Vec3::dot(p1, normale) > 0) && (Vec3::dot(p0, normale) > 0) == (Vec3::dot(p2, normale) > 0) && (Vec3::dot(p0, normale) > 0) == (Vec3::dot(p3, normale) > 0)){
                intersection.t = t;
                intersection.intersectionExists = true;
                intersection.intersection = p; 
                intersection.normal = normale;
            }
        }

        return intersection;
    }
};
#endif // SQUARE_H
