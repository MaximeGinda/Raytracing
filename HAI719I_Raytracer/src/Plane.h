#ifndef PLANE_H
#define PLANE_H
#include "Vec3.h"
#include "Line.h"
class Plane {
private:
    Vec3 m_center , m_normal;
public:
    Plane() {}
    Plane( Vec3 const & c , Vec3 const & n ) {
        m_center = c;
        m_normal = n; m_normal.normalize();
    }
    void setCenter( Vec3 const & c ) { m_center = c; }
    void setNormal( Vec3 const & n ) { m_normal = n; m_normal.normalize(); }
    Vec3 const & center() const { return m_center; }
    Vec3 const & normal() const { return m_normal; }
    Vec3 project( Vec3 const & p ) const {
        Vec3 result;
        //TODO completer
        return result;
    }
    float squareDistance( Vec3 const & p ) const { return (project(p) - p).squareLength(); }
    float distance( Vec3 const & p ) const { return sqrt( squareDistance(p) ); }
    bool isParallelTo( Line const & L ) const {
        return Vec3::dot(m_normal, L.direction()) == 0;
    }
    Vec3 getIntersectionPoint( Line const & L ) const {
        // you should check first that the line is not parallel to the plane!
        Vec3 result;
        if(!isParallelTo(L)) {
            Vec3 a = m_center;
            Vec3 n = m_normal;

            float D = Vec3::dot(a, n);
            Vec3 o = L.origin();
            Vec3 d = L.direction();

            float t = (D - Vec3::dot(o, n)) / Vec3::dot(d, n);

            result = o + t * d;
        }
        return result;
    }
};
#endif
