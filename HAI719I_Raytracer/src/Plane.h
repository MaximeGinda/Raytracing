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
        Vec3 result;

        if(!isParallelTo(L)) {

            Vec3 center = m_center;
            Vec3 normal = m_normal;

            Vec3 ori = L.origin();
            Vec3 dir = L.direction();

            float t = (Vec3::dot(center, normal) - Vec3::dot(ori, normal)) / Vec3::dot(dir, normal);

            result = ori + t * dir;
        }
        return result;
    }
};
#endif
