#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"

struct RayTriangleIntersection{
    bool intersectionExists;
    float t;
    float w0,w1,w2;
    unsigned int tIndex;
    Vec3 intersection;
    Vec3 normal;
};

class Triangle {
private:
    Vec3 m_c[3] , m_normal;
    float area;
public:
    Triangle() {}
    Triangle( Vec3 const & c0 , Vec3 const & c1 , Vec3 const & c2 ) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
    }
    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross( m_c[1] - m_c[0] , m_c[2] - m_c[0] );
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }
    void setC0( Vec3 const & c0 ) { m_c[0] = c0; } // remember to update the area and normal afterwards!
    void setC1( Vec3 const & c1 ) { m_c[1] = c1; } // remember to update the area and normal afterwards!
    void setC2( Vec3 const & c2 ) { m_c[2] = c2; } // remember to update the area and normal afterwards!
    Vec3 const & normal() const { return m_normal; }
    Vec3 projectOnSupportPlane( Vec3 const & p ) const {
        Vec3 result;
        //TODO completer
        return result;
    }
    float squareDistanceToSupportPlane( Vec3 const & p ) const {
        float result;
        //TODO completer
        return result;
    }
    float distanceToSupportPlane( Vec3 const & p ) const { return sqrt( squareDistanceToSupportPlane(p) ); }
    bool isParallelTo( Line const & L ) const {
        return Vec3::dot(m_normal, L.direction()) == 0;
    }
    Vec3 getIntersectionPointWithSupportPlane( Line const & L ) const {
        Vec3 result;

        Plane plane = Plane(m_c[0], m_normal);

        if(!isParallelTo(L)) 
            result = plane.getIntersectionPoint(L);
        
        return result;
    }
    void computeBarycentricCoordinates( Vec3 const & p , float & u0 , float & u1 , float & u2 ) const {
        u0 = Vec3::dot(m_normal, Vec3::cross(m_c[1] - m_c[0], p - m_c[0]));
        u1 = Vec3::dot(m_normal, Vec3::cross(m_c[2] - m_c[1], p - m_c[1]));
        u2 = Vec3::dot(m_normal, Vec3::cross(m_c[0] - m_c[2], p - m_c[2]));
    }

    RayTriangleIntersection getIntersection( Ray const & ray ) const {
        RayTriangleIntersection result;
        result.intersectionExists = false;

        if(!isParallelTo(ray)) {

            Vec3 p = getIntersectionPointWithSupportPlane(ray);
            float w0, w1, w2;
            computeBarycentricCoordinates(p, w0, w1, w2);

            if(w0 <= 1 && w0 >= 0 && w1 <= 1 && w1 >= 0 && w2 <= 1 && w2 >= 0) {
                
                Plane plane = Plane(m_c[0], m_normal);

                Vec3 center = plane.center();
                Vec3 normal = plane.normal();

                Vec3 ori = ray.origin();
                Vec3 dir = ray.direction();

                float t = (Vec3::dot(center, normal) - Vec3::dot(ori, normal)) / Vec3::dot(dir, normal);

                result.intersectionExists = true;
                result.t = t;
                result.w0 = w0;
                result.w1 = w1;
                result.w2 = w2;
                result.intersection = p;
                result.normal = m_normal;

            }

        return result;
        }
    }
};
#endif
