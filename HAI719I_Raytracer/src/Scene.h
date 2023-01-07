#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"


#include <GL/glut.h>


enum LightType {
    LightType_Spherical,
    LightType_Quad
};


struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Mesh quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

};

// struct Node {
//     std::vector<double> point;
//     KDNode *left;
//     KDNode *right;

//     KDNode(const std::vector<double> &pt) : point(pt), left(nullptr), right(nullptr) {}
    
// };

// struct Tree {   
//     Node *root;
//     Tree() : root(nullptr) {}

//     void insert(const std::vector<double> &point) {
//         root = insert(root, point, 0);
//     }

//     Node* insert((Node *node, const std::vector<double> &point, int depth){
//         if (node == nullptr) return new KDNode(point);

//         int dim = depth % point.size();

//         if (point[dim] < node->point[dim]) {
//             node->left = insert(node->left, point, depth + 1);
//         } else {
//             node->right = insert(node->right, point, depth + 1);
//         }

//         return node;
//     }
// };

struct BoundingBox {
    std::array<double, 3> min;
    std::array<double, 3> max;

    BoundingBox() : min({0, 0, 0}), max({0, 0, 0}) {}

    BoundingBox(const std::array<double, 3> &min, const std::array<double, 3> &max) : min(min), max(max) {}

    bool intersects(const Ray &ray) const {
        // Code de l'intersection de la bounding box et du rayon
        // ...
    }

    void expand(const BoundingBox &other) {
        for (int i = 0; i < 3; i++) {
            min[i] = std::min(min[i], other.min[i]);
            max[i] = std::max(max[i], other.max[i]);
        }
    }

    void draw(const BoundingBox &box) {
        std::array<double, 24> vertices = {
            box.min[0], box.min[1], box.min[2],
            box.max[0], box.min[1], box.min[2],
            box.max[0], box.max[1], box.min[2],
            box.min[0], box.max[1], box.min[2],
            box.min[0], box.min[1], box.max[2],
            box.max[0], box.min[1], box.max[2],
            box.max[0], box.max[1], box.max[2],
            box.min[0], box.max[1], box.max[2]
        };

        std::array<unsigned int, 36> indices = {
            0, 1, 2, 2, 3, 0,
            4, 5, 6, 6, 7, 4,
            4, 5, 1, 1, 0, 4,
            6, 7, 3, 3, 2, 6,
            1, 5, 6, 6, 2, 1,
            4, 0, 3, 3, 7, 4
        };

        glBegin(GL_LINE_STRIP);
        for (unsigned int i : indices) {
            glVertex3d(vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]);
        }
        glEnd();
    }


    std::pair<std::array<float, 3>, std::array<float, 3>> getBounds(const Mesh &mesh) {
        std::array<float, 3> min = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
        std::array<float, 3> max = {std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};

        for (unsigned int i : mesh.triangles_array) {
            const MeshVertex &vertex = mesh.vertices[i];
            for (int j = 0; j < 3; j++) {
                min[j] = std::min(min[j], vertex.position[j]);
                max[j] = std::max(max[j], vertex.position[j]);
            }
        }

        return {min, max};
    }
};

struct RaySceneIntersection{
    bool intersectionExists;
    unsigned int typeOfIntersectedObject;
    unsigned int objectIndex;
    float t;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;
    RaySceneIntersection() : intersectionExists(false) , t(FLT_MAX) {}
};



class Scene {
    std::vector< Mesh > meshes;
    std::vector< Sphere > spheres;
    std::vector< Square > squares;
    std::vector< Light > lights;

public:


    Scene() {
    }

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for( unsigned int It = 0 ; It < meshes.size() ; ++It ) {
            Mesh const & mesh = meshes[It];
            mesh.draw();

            
            BoundingBox box;
            std::pair<std::array<float, 3>, std::array<float, 3>> bounds = box.getBounds(mesh);
            std::array<float, 3> min = bounds.first;
            std::array<float, 3> max = bounds.second;

            BoundingBox box1(min, max);
            box.expand(box1);
            box.draw(box);
        }
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            sphere.draw();
        }
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) {
            Square const & square = squares[It];
            square.draw();
        }
    }

    // Tree buildKDTree(std::vector<std::vector<double>> points) {
    //     Tree tree;

    //     for (const auto &point : points) {
    //         tree.insert(point);
    //     }

    //     return tree;
    // }

    // renvoie un float de l'intersection la plus proche
    float searchFirstIntersection(Ray const &ray)
    {

        size_t meshesSize = meshes.size();
        for (size_t i = 0; i < meshesSize; i++)
        {
            // RayTriangleIntersection rayMesh = meshes[i].intersect(ray);
            // if (rayMesh.intersectionExists) return rayMesh.t;
        }

        size_t spheresSize = spheres.size();
        for (size_t i = 0; i < spheresSize; i++)
        {
            RaySphereIntersection raySphere = spheres[i].intersect(ray);
            if (raySphere.intersectionExists && spheres[i].material.type != Material_Glass) return raySphere.t;
            // if (raySphere.intersectionExists) return raySphere.t;
        }

        size_t squaresSize = squares.size();
        for (size_t i = 0; i < squaresSize; i++)
        {
            RaySquareIntersection raySquare = squares[i].intersect(ray);
            if (raySquare.intersectionExists && squares[i].material.type != Material_Glass) return raySquare.t;
        }

        return FLT_MAX;
    }

    float calculateCoef(int l_num, int echant, Vec3 intersect)
    {
        int nb_ombre = 0;
        float x0 = lights[l_num].pos[0] - lights[l_num].radius/2,
              y0 = lights[l_num].pos[1] - lights[l_num].radius/2,
              z0 = lights[l_num].pos[2] - lights[l_num].radius/2;

        Vec3 Lvec;
        Ray omb;

        float pasX, pasZ;
        for (int x = 0; x < echant; x++)
        {
            pasX = (float)(rand() / (float)(RAND_MAX / (lights[l_num].radius)));
            pasZ = (float)(rand() / (float)(RAND_MAX / (lights[l_num].radius)));
            Lvec = Vec3(x0 + pasX, y0, z0 + pasZ) - intersect;
            Lvec.normalize();

            omb = Ray(intersect, Lvec);

            // on cherche l'intersection la plus proche
            float ombre = searchFirstIntersection(omb);
            if (ombre < 1 && ombre > 0.00001) nb_ombre++;

        }

        return (float)nb_ombre / echant; 
    }

    Vec3 deapthOfField(RaySceneIntersection result, Vec3 color, Vec3 intersect){

        // Calcul de la distance de mise au point et du rayon de confusion
        float focus_distance = 4; // distance de mise au point en mètres
        float aperture_size = 1; // taille de l'ouverture en millimètres
        float blur_radius = (1.0 / aperture_size) * focus_distance; // rayon de confusion en mètres

        Vec3 blur_color = Vec3(0.3,0.3,0.3) ;

        // Calcul de la distance de l'objet à la distance de mise au point
        float distance_to_focus = abs(result.t - focus_distance);

        // Ajout de flou au pixel si nécessaire
        if (distance_to_focus < blur_radius) {

            float blur_amount = (blur_radius - distance_to_focus) / blur_radius;

            for(int i = 0; i < 3; i++)
                color[i] = (color[i] * (1.0 - blur_amount)) + blur_amount * blur_color[i];;
        }

        return color;
    }

    RaySceneIntersection computeIntersection(Ray const & ray, float znear) {
        RaySceneIntersection result;

         //On regarde toutes les meshes
        size_t meshesSize = meshes.size();
        for (size_t i = 0; i < meshesSize; i++)
        {
            RayTriangleIntersection rmi = this->meshes[i].intersect(ray);
            if (rmi.intersectionExists){
                // Est-ce que c'est le plus proche ?
                if (rmi.t > znear && rmi.t < result.t) {
                    result.intersectionExists = rmi.intersectionExists;
                    result.typeOfIntersectedObject = 0;
                    result.objectIndex = i;
                    result.t = rmi.t;
                    result.rayMeshIntersection = rmi; 
                }
            }
        }

        // On regarde toutes les spheres
        size_t spheresSize = spheres.size();
        for(size_t i = 0; i < spheresSize; i++){
            RaySphereIntersection rsi = spheres[i].intersect(ray);
            if (rsi.intersectionExists){
                // Est-ce que c'est le plus proche ?
                if(rsi.t < result.t && rsi.t > znear){
                    result.intersectionExists = rsi.intersectionExists;
                    result.typeOfIntersectedObject = 1;
                    result.objectIndex = i;
                    result.t = rsi.t;    
                    result.raySphereIntersection = rsi;     
                }
            }
        }  

        //On regarde tous les carrés
        size_t squaresSize = squares.size();
        for(size_t i = 0; i < squaresSize; i++){
            RaySquareIntersection rsi = squares[i].intersect(ray);
            if (rsi.intersectionExists){
                // Est-ce que c'est le plus proche ?
                if(rsi.t < result.t && rsi.t > znear){
                    result.intersectionExists = rsi.intersectionExists;
                    result.typeOfIntersectedObject = 2;
                    result.objectIndex = i;
                    result.t = rsi.t;    
                    result.raySquareIntersection = rsi;     
                }
            }
        }                
        return result;
    }

    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces, float znear ) {

        RaySceneIntersection raySceneIntersection = computeIntersection(ray, znear);
        Ray rayShadow;
        Vec3 color;
        Vec3 inter;

        size_t lightsSize = lights.size();

         if(raySceneIntersection.intersectionExists){
            switch(raySceneIntersection.typeOfIntersectedObject){
                case 1: { // SPHERE
                    for(size_t lnum = 0; lnum < lightsSize; lnum++){
                        Vec3 newOrigin = raySceneIntersection.raySphereIntersection.intersection;
                        Vec3 newNorm = raySceneIntersection.raySphereIntersection.normal;


                        if(spheres[raySceneIntersection.objectIndex].material.type == Material_Diffuse_Blinn_Phong){

                            Vec3 n = raySceneIntersection.raySphereIntersection.normal;
                            n.normalize();

                            Vec3 l = lights[lnum].pos - raySceneIntersection.raySphereIntersection.intersection;
                            l.normalize();

                            Vec3 r  = 2 * (Vec3::dot(l,n)) * n - l;
                            r.normalize();

                            Vec3 v = ray.origin() - raySceneIntersection.raySphereIntersection.intersection;
                            v.normalize();

                            float theta = std::max(Vec3::dot(n,l), 0.f);
                            float alpha = std::max(Vec3::dot(r,v), 0.f);

                            float compAmbianteR = spheres[raySceneIntersection.objectIndex].material.ambient_material[0] * lights[lnum].material[0], 
                                  compAmbianteG = spheres[raySceneIntersection.objectIndex].material.ambient_material[1] * lights[lnum].material[1], 
                                  compAmbianteB = spheres[raySceneIntersection.objectIndex].material.ambient_material[2] * lights[lnum].material[2];
                            float compDiffuseR = spheres[raySceneIntersection.objectIndex].material.diffuse_material[0] * lights[lnum].material[0] * theta,
                                  compDiffuseG = spheres[raySceneIntersection.objectIndex].material.diffuse_material[1] * lights[lnum].material[1] * theta, 
                                  compDiffuseB = spheres[raySceneIntersection.objectIndex].material.diffuse_material[2] * lights[lnum].material[2] * theta;
                            float compSpeculaireR = spheres[raySceneIntersection.objectIndex].material.specular_material[0] * lights[lnum].material[0] * pow(alpha, this->spheres[raySceneIntersection.objectIndex].material.shininess),
                                  compSpeculaireG = spheres[raySceneIntersection.objectIndex].material.specular_material[1] * lights[lnum].material[1] * pow(alpha, this->spheres[raySceneIntersection.objectIndex].material.shininess),
                                  compSpeculaireB = spheres[raySceneIntersection.objectIndex].material.specular_material[2] * lights[lnum].material[2] * pow(alpha, this->spheres[raySceneIntersection.objectIndex].material.shininess);

                            
                            color[0] = compAmbianteR + compDiffuseR + compSpeculaireR;

                            // pour éviter de corrompre l'image
                            if(color[0] <= 0.00001){
                                color[0] = 0;
                            }

                            color[1] = compAmbianteG + compDiffuseG + compSpeculaireG;

                            if(color[1] <= 0.00001){
                                color[1] = 0;
                            }

                            color[2] = compAmbianteB + compDiffuseB + compSpeculaireB;

                            if(color[2] <= 0.00001){
                                color[2] = 0;
                            }
                        }
                        else if(spheres[raySceneIntersection.objectIndex].material.type == Material_Glass){

                            if(NRemainingBounces > 0){

                                
                                float ind_med = spheres[raySceneIntersection.objectIndex].material.index_medium;
                                float dot = Vec3::dot(newNorm, ray.direction());
                                float k = (ind_med * ind_med) * (dot * dot);

                                Vec3 newDir;
                                if(k >= 0.0){
                                    newDir = ind_med * ray.direction() + ((ind_med * (dot + sqrt(k))) * newNorm);
                                    newDir.normalize();
                                }
                                
                                Ray newRay = Ray(newOrigin,newDir);
                                color = rayTraceRecursive(newRay,NRemainingBounces-1,0.f);
                            }
                        }
                        else if(spheres[raySceneIntersection.objectIndex].material.type == Material_Mirror){

                            if(NRemainingBounces > 0){

                                Vec3 L = ray.origin() - newOrigin;
                                L.normalize();

                                Vec3 newDir = 2 * Vec3::dot(newNorm, L) * newNorm - L;
                                newDir.normalize();

                                Ray newRay = Ray(newOrigin,newDir);
                                color = rayTraceRecursive(newRay, NRemainingBounces-1,0.000001f);
                            }
                        }
                        else if(spheres[raySceneIntersection.objectIndex].material.type == texture){
                            
                            
                        }
                        
                        inter = raySceneIntersection.raySphereIntersection.intersection;
                    }
                } 
                    break;
                case 2: { // SQUARE 
                    for(size_t lnum = 0; lnum < lightsSize; lnum++){

                        Vec3 newOrigin = raySceneIntersection.raySquareIntersection.intersection;
                        Vec3 newNorm = raySceneIntersection.raySquareIntersection.normal;

                        if(squares[raySceneIntersection.objectIndex].material.type == Material_Diffuse_Blinn_Phong){
                            
                            Vec3 n = raySceneIntersection.raySquareIntersection.normal;
                            n.normalize();

                            Vec3 l = lights[0].pos - raySceneIntersection.raySquareIntersection.intersection;
                            l.normalize();

                            Vec3 r  = 2 * (Vec3::dot(l,n)) * n - l;
                            r.normalize();

                            Vec3 v = ray.origin() - raySceneIntersection.raySquareIntersection.intersection;
                            v.normalize();

                            float theta = std::max(Vec3::dot(n,l), 0.f);
                            float alpha = std::max(Vec3::dot(r,v), 0.f);

                            float compAmbianteR = squares[raySceneIntersection.objectIndex].material.ambient_material[0] * lights[lnum].material[0], 
                                  compAmbianteG = squares[raySceneIntersection.objectIndex].material.ambient_material[1] * lights[lnum].material[1], 
                                  compAmbianteB = squares[raySceneIntersection.objectIndex].material.ambient_material[2] * lights[lnum].material[2];
                            float compDiffuseR = squares[raySceneIntersection.objectIndex].material.diffuse_material[0] * lights[lnum].material[0] * theta,
                                  compDiffuseG = squares[raySceneIntersection.objectIndex].material.diffuse_material[1] * lights[lnum].material[1] * theta, 
                                  compDiffuseB = squares[raySceneIntersection.objectIndex].material.diffuse_material[2] * lights[lnum].material[2] * theta;
                            float compSpeculaireR = squares[raySceneIntersection.objectIndex].material.specular_material[0] * lights[lnum].material[0] * pow(alpha, this->squares[raySceneIntersection.objectIndex].material.shininess),
                                  compSpeculaireG = squares[raySceneIntersection.objectIndex].material.specular_material[1] * lights[lnum].material[1] * pow(alpha, this->squares[raySceneIntersection.objectIndex].material.shininess),
                                  compSpeculaireB = squares[raySceneIntersection.objectIndex].material.specular_material[2] * lights[lnum].material[2] * pow(alpha, this->squares[raySceneIntersection.objectIndex].material.shininess);

                            
                            color[0] = compAmbianteR + compDiffuseR + compSpeculaireR;

                            // pour éviter de corrompre l'image
                            if(color[0] <= 0.00001){
                                color[0] = 0;
                            }

                            color[1] = compAmbianteG + compDiffuseG + compSpeculaireG;

                            if(color[1] <= 0.00001){
                                color[1] = 0;
                            }

                            color[2] = compAmbianteB + compDiffuseB + compSpeculaireB;
                            
                            if(color[2] <= 0.00001){
                                color[2] = 0;
                            }
                        }
                        else if(squares[raySceneIntersection.objectIndex].material.type == Material_Glass){
     
                            if(NRemainingBounces > 0){

                                float ind_med = squares[raySceneIntersection.objectIndex].material.index_medium;
                                float dot = Vec3::dot(newNorm, ray.direction());
                                float k = (ind_med * ind_med) * (dot * dot);

                                Vec3 newDir;
                                if(k >= 0.0){
                                    newDir = ind_med * ray.direction() + ((ind_med * (dot + sqrt(k))) * newNorm);
                                    newDir.normalize();
                                }
                                
                                Ray newRay = Ray(newOrigin,newDir);
                                color = rayTraceRecursive(newRay,NRemainingBounces-1, 0.f);
                            }
                        }
                        else if(squares[raySceneIntersection.objectIndex].material.type == Material_Mirror){

                            if(NRemainingBounces > 0){

                                Vec3 L = ray.origin() - newOrigin;
                                L.normalize();

                                Vec3 newDir = 2 * Vec3::dot(newNorm, L) * newNorm - L;
                                newDir.normalize();

                                Ray newRay = Ray(newOrigin,newDir);
                                color = rayTraceRecursive(newRay, NRemainingBounces-1,0.000001f);

                            }
                        }
                    }

                    inter = raySceneIntersection.raySquareIntersection.intersection;
                    
                }
                    break;
                case 0 : { // Mesh
                    for(size_t lnum = 0; lnum < lightsSize; lnum++){
                        if(meshes[raySceneIntersection.objectIndex].material.type == Material_Diffuse_Blinn_Phong){
                            
                            Vec3 n = raySceneIntersection.rayMeshIntersection.normal;
                            n.normalize();

                            Vec3 l = lights[0].pos - raySceneIntersection.rayMeshIntersection.intersection;
                            l.normalize();

                            Vec3 r  = 2 * (Vec3::dot(l,n)) * n - l;
                            r.normalize();

                            Vec3 v = ray.origin() - raySceneIntersection.rayMeshIntersection.intersection;
                            v.normalize();

                            float theta = std::max(Vec3::dot(n,l), 0.f);
                            float alpha = std::max(Vec3::dot(r,v), 0.f);

                            float compAmbianteR = meshes[raySceneIntersection.objectIndex].material.ambient_material[0] * lights[lnum].material[0], 
                                  compAmbianteG = meshes[raySceneIntersection.objectIndex].material.ambient_material[1] * lights[lnum].material[1], 
                                  compAmbianteB = meshes[raySceneIntersection.objectIndex].material.ambient_material[2] * lights[lnum].material[2];
                            float compDiffuseR = meshes[raySceneIntersection.objectIndex].material.diffuse_material[0] * lights[lnum].material[0] * theta,
                                  compDiffuseG = meshes[raySceneIntersection.objectIndex].material.diffuse_material[1] * lights[lnum].material[1] * theta, 
                                  compDiffuseB = meshes[raySceneIntersection.objectIndex].material.diffuse_material[2] * lights[lnum].material[2] * theta;
                            float compSpeculaireR = meshes[raySceneIntersection.objectIndex].material.specular_material[0] * lights[lnum].material[0] * pow(alpha, this->meshes[raySceneIntersection.objectIndex].material.shininess),
                                  compSpeculaireG = meshes[raySceneIntersection.objectIndex].material.specular_material[1] * lights[lnum].material[1] * pow(alpha, this->meshes[raySceneIntersection.objectIndex].material.shininess),
                                  compSpeculaireB = meshes[raySceneIntersection.objectIndex].material.specular_material[2] * lights[lnum].material[2] * pow(alpha, this->meshes[raySceneIntersection.objectIndex].material.shininess);

                            
                            color[0] = compAmbianteR + compDiffuseR + compSpeculaireR;

                            // pour éviter de corrompre l'image
                            if(color[0] <= 0.00001){
                                color[0] = 0;
                            }

                            color[1] = compAmbianteG + compDiffuseG + compSpeculaireG;

                            if(color[1] <= 0.00001){
                                color[1] = 0;
                            }

                            color[2] = compAmbianteB + compDiffuseB + compSpeculaireB;
                            
                            if(color[2] <= 0.00001){
                                color[2] = 0;
                            }
                        }
                    }
                }   
                    break;
                default:
                    break;
            }
        }

        // Shadow
        for(size_t i = 0; i < lightsSize; i++){
            float coeff = calculateCoef(i, 10, inter);
            color *= 1 - coeff;
        }

        // Deapth of field
        color = deapthOfField(raySceneIntersection, color, inter);

        return color;
    }

    Vec3 rayTrace( Ray const & rayStart ) {
        return rayTraceRecursive(rayStart, 5, 4.9f);
    }

    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0. , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Diffuse_Blinn_Phong;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
    }

    void setup_single_sphere_text() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        ppmLoader::ImageRGB imageRGB;

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0. , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = texture;
            s.material.imageRGB = imageRGB;
        }

        // ppmLoader::load_ppm(imageRGB,"s2.ppm");
    }

    void setup_single_mesh()
    {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

         {
            lights.resize(lights.size() + 1);
            Light &light = lights[lights.size() - 1];
            light.pos = Vec3(0.0, 1.5, 0.0);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;
        }   
        {
            meshes.resize(meshes.size() + 1);
            Mesh &m = meshes[meshes.size() - 1];
            m.loadOFF("data/suzanne.off");
            m.centerAndScaleToUnit();
            m.material.type = Material_Diffuse_Blinn_Phong;
            m.material.diffuse_material = Vec3(1., 0., 0.);
            m.material.specular_material = Vec3(1., 1., 1.);
            m.material.shininess = 16;
            m.build_arrays();
        }
        
    }


    // void setup_single_mesh()
    // {
    //     meshes.clear();
    //     spheres.clear();
    //     squares.clear();
    //     lights.clear();

    //      {
    //         lights.resize(lights.size() + 1);
    //         Light &light = lights[lights.size() - 1];
    //         light.pos = Vec3(0.0, 1.5, 0.0);
    //         light.radius = 2.5f;
    //         light.powerCorrection = 2.f;
    //         light.type = LightType_Spherical;
    //         light.material = Vec3(1, 1, 1);
    //         light.isInCamSpace = false;
    //     }   

    //     { //Left Wall

    //         squares.resize(squares.size() + 1);
    //         Square &s = squares[squares.size() - 1];
    //         s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
    //         s.scale(Vec3(2., 2., 1.));
    //         s.translate(Vec3(0., 0., -2.));
    //         s.rotate_y(90);
    //         s.build_arrays();
    //         s.material.diffuse_material = Vec3(1., 0., 0.);
    //         s.material.specular_material = Vec3(1., 0., 0.);
    //         s.material.shininess = 16;
    //     }

    //     { //Back Wall
    //         squares.resize(squares.size() + 1);
    //         Square &s = squares[squares.size() - 1];
    //         s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
    //         s.scale(Vec3(2., 2., 1.));
    //         s.translate(Vec3(0., 0., -2.));
    //         s.build_arrays();
    //         s.material.diffuse_material = Vec3(1., 1., 1.);
    //         s.material.specular_material = Vec3(1., 1., 1.);
    //         s.material.shininess = 16;
    //     }

    //     { //Right Wall
    //         squares.resize(squares.size() + 1);
    //         Square &s = squares[squares.size() - 1];
    //         s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
    //         s.translate(Vec3(0., 0., -2.));
    //         s.scale(Vec3(2., 2., 1.));
    //         s.rotate_y(-90);
    //         s.build_arrays();
    //         s.material.diffuse_material = Vec3(0.0, 1.0, 0.0);
    //         s.material.specular_material = Vec3(0.0, 1.0, 0.0);
    //         s.material.shininess = 16;
    //     }

    //     { //Front Wall
    //     	squares.resize( squares.size() + 1 );
    //         Square & s = squares[squares.size() - 1];
    //         s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
    //         s.translate(Vec3(0., 0., -2.));
    //         s.scale(Vec3(2., 2., 1.));
    //         s.rotate_y(180);
    //         s.build_arrays();
    //         s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
    //         s.material.specular_material = Vec3( 1.0,1.0,1.0 );
    //         s.material.shininess = 16;
    //     }

    //     { //Floor
    //         squares.resize(squares.size() + 1);
    //         Square &s = squares[squares.size() - 1];
    //         s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
    //         s.translate(Vec3(0., 0., -2.));
    //         s.scale(Vec3(2., 2., 1.));
    //         s.rotate_x(-90);
    //         s.build_arrays();
    //         s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
    //         s.material.specular_material = Vec3(1.0, 1.0, 1.0);
    //         s.material.shininess = 16;
    //     }

    //     { //Ceiling
    //         squares.resize(squares.size() + 1);
    //         Square &s = squares[squares.size() - 1];
    //         s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
    //         s.translate(Vec3(0., 0., -2.));
    //         s.scale(Vec3(2., 2., 1.));
    //         s.rotate_x(90);
    //         s.build_arrays();
    //         s.material.diffuse_material = Vec3(1.0, 1.0, 1.0);
    //         s.material.specular_material = Vec3(1.0, 1.0, 1.0);
    //         s.material.shininess = 16;
    //     }

    //     {
    //         meshes.resize(meshes.size() + 1);
    //         Mesh &m = meshes[meshes.size() - 1];
    //         m.loadOFF("data/suzanne.off");
    //         m.centerAndScaleToUnit();
    //         m.material.type = Material_Diffuse_Blinn_Phong;
    //         m.material.diffuse_material = Vec3(1., 0., 0.);
    //         m.material.specular_material = Vec3(1., 1., 1.);
    //         m.material.shininess = 16;
    //         m.build_arrays();
    //     }
    // }

    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.8,0.8,0.8 );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
    }

    void setup_cornell_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();  
            // s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.,0.,1. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            // s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3(1.,1.,.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Front Wall
        	squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }


        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            // s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            // s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

    }

};



#endif
