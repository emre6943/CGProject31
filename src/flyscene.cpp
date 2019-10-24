#include "flyscene.hpp"
#include <GLFW/glfw3.h>


void Flyscene::initialize(int width, int height) {
    // initiliaze the Phong Shading effect for the Opengl Previewer
    phong.initialize();


    // set the camera's projection matrix
    flycamera.setPerspectiveMatrix(60.0, width / (float) height, 0.1f, 100.0f);
    flycamera.setViewport(Eigen::Vector2f((float) width, (float) height));

    // load the OBJ file and materials
    Tucano::MeshImporter::loadObjFile(mesh, materials,
                                      "resources/models/dodgeColorTest.obj");


    // normalize the model (scale to unit cube and center at origin)
    mesh.normalizeModelMatrix();

    // pass all the materials to the Phong Shader
    for (int i = 0; i < materials.size(); ++i)
        phong.addMaterial(materials[i]);



    // set the color and size of the sphere to represent the light sources
    // same sphere is used for all sources
    lightrep.setColor(Eigen::Vector4f(1.0, 1.0, 0.0, 1.0));
    lightrep.setSize(0.15);

    // create a first ray-tracing light source at some random position
    lights.push_back(Eigen::Vector3f(-1.0, 1.0, 1.0));

    //set the color and the size of our debugOrbRep
    debugOrbRep.setColor(Eigen::Vector4f(1.0, 0.0, 0.0, 1.0));
    debugOrbRep.setSize(.08);

    // scale the camera representation (frustum) for the ray debug
    camerarep.shapeMatrix()->scale(0.2);

    // the debug ray is a cylinder, set the radius and length of the cylinder
    ray.setSize(0.005, 10.0);

    // craete a first debug ray pointing at the center of the screen
    createDebugRay(Eigen::Vector2f(width / 2.0, height / 2.0));

    glEnable(GL_DEPTH_TEST);

    // for (int i = 0; i<mesh.getNumberOfFaces(); ++i){
    //   Tucano::Face face = mesh.getFace(i);
    //   for (int j =0; j<face.vertex_ids.size(); ++j){
    //     std::cout<<"vid "<<j<<" "<<face.vertex_ids[j]<<std::endl;
    //     std::cout<<"vertex "<<mesh.getVertex(face.vertex_ids[j]).transpose()<<std::endl;
    //     std::cout<<"normal "<<mesh.getNormal(face.vertex_ids[j]).transpose()<<std::endl;
    //   }
    //   std::cout<<"mat id "<<face.material_id<<std::endl<<std::endl;
    //   std::cout<<"face   normal "<<face.normal.transpose() << std::endl << std::endl;
    // }



}

void Flyscene::paintGL(void) {

    // update the camera view matrix with the last mouse interactions
    flycamera.updateViewMatrix();
    Eigen::Vector4f viewport = flycamera.getViewport();

    // clear the screen and set background color
    glClearColor(0.9, 0.9, 0.9, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // position the scene light at the last ray-tracing light source
    scene_light.resetViewMatrix();
    scene_light.viewMatrix()->translate(-lights.back());

    // render the scene using OpenGL and one light source
    phong.render(mesh, flycamera, scene_light);


    // render the ray and camera representation for ray debug
    ray.render(flycamera, scene_light);
    camerarep.render(flycamera, scene_light);

    // render ray tracing light sources as yellow spheres
    for (int i = 0; i < lights.size(); ++i) {
        lightrep.resetModelMatrix();
        lightrep.modelMatrix()->translate(lights[i]);
        lightrep.render(flycamera, scene_light);
    }

    //render the debugOrbRep's as red spheres
    for (int i = 0; i < debugOrbs.size(); i++) {
        debugOrbRep.resetModelMatrix();
        debugOrbRep.modelMatrix()->translate(debugOrbs[i]);
        debugOrbRep.render(flycamera, scene_light);
    }

    // render coordinate system at lower right corner
    flycamera.renderAtCorner();
}

void Flyscene::simulate(GLFWwindow *window) {
    // Update the camera.
    // NOTE(mickvangelderen): GLFW 3.2 has a problem on ubuntu where some key
    // events are repeated: https://github.com/glfw/glfw/issues/747. Sucks.
    float dx = (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS ? 1.0 : 0.0) -
               (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS ? 1.0 : 0.0);
    float dy = (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS ||
                glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS
                ? 1.0
                : 0.0) -
               (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS ||
                glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS
                ? 1.0
                : 0.0);
    float dz = (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS ? 1.0 : 0.0) -
               (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS ? 1.0 : 0.0);
    flycamera.translate(dx, dy, dz);
}


void Flyscene::raytraceScene(int width, int height) {
    std::cout << "ray tracing ..." << std::endl;

    // if no width or height passed, use dimensions of current viewport
    Eigen::Vector2i image_size(width, height);
    if (width == 0 || height == 0) {
        image_size = flycamera.getViewportSize();
    }

    // create 2d vector to hold pixel colors and resize to match image size
    vector<vector<Eigen::Vector3f>> pixel_data;
    pixel_data.resize(image_size[1]);
    for (int i = 0; i < image_size[1]; ++i)
        pixel_data[i].resize(image_size[0]);

    // origin of the ray is always the camera center
    Eigen::Vector3f origin = flycamera.getCenter();
    Eigen::Vector3f screen_coords;

    // for every pixel shoot a ray from the origin through the pixel coords
    for (int j = 0; j < image_size[1]; ++j) {
        for (int i = 0; i < image_size[0]; ++i) {
            // create a ray from the camera passing through the pixel (i,j)
            screen_coords = flycamera.screenToWorld(Eigen::Vector2f(i, j));
            // launch raytracing for the given ray and write result to pixel data
            pixel_data[i][j] = traceRay(origin, screen_coords);
        }
    }

    // write the ray tracing result to a PPM image
    Tucano::ImageImporter::writePPMImage("result.ppm", pixel_data);
    std::cout << "ray tracing done! " << std::endl;
}


//Intersection with plane
auto intersectPlane(Eigen::Vector3f start, Eigen::Vector3f to, Eigen::Vector3f normal, Eigen::Vector4f p_onPlane) {
    Eigen::Vector3f p = Eigen::Vector3f(p_onPlane.x(), p_onPlane.y(), p_onPlane.z());
    float distance = p.dot(normal);
    Eigen::Vector4f point;
    struct result {
        bool inter;
        float t;
        Eigen::Vector4f point;
    };

    if (to.dot(normal) == 0) {
        return result{false, 0, p_onPlane};
    }

    float t = (distance - (start.dot(normal))) / (to.dot(normal));
    point = Eigen::Vector4f((start + t * to), p_onPlane.w());
    return result{true, t, point};
}

// a method used in triangle intersection
float sign(Eigen::Vector3f p1, Eigen::Vector4f p2, Eigen::Vector4f p3) {
    return ((p1.x() - p3.x()) * (p2.y() - p3.y())) - ((p2.x() - p3.x()) * (p1.y() - p3.y()));
}

//Intersection with triangle
//give normalized vectors
auto intersectTriange(Eigen::Vector3f start, Eigen::Vector3f to, Tucano::Face face, Tucano::Mesh mesh) {
    Eigen::Vector3f color;
    Eigen::Vector3f facenormal = face.normal;
    std::vector<Eigen::Vector4f> vecs;

    struct result {
        bool inter;
        Tucano::Face face;
        Tucano::Mesh mesh;
        Eigen::Vector4f hit;
    };

    for (int i = 0; i < 3; i++) {
        int vertexid = face.vertex_ids[i];
        vecs.push_back(mesh.getVertex(face.vertex_ids[i]));
    }


    auto intersectionPlane = intersectPlane(start, to, facenormal, vecs[0]);
    if (intersectionPlane.inter) {
        Eigen::Vector4f p_onPlane = intersectionPlane.point;

        float a = sign(p_onPlane, vecs[0], vecs[1]);
        float b = sign(p_onPlane, vecs[1], vecs[2]);
        float c = sign(p_onPlane, vecs[2], vecs[0]);


        bool has_neg = (a < 0) || (b < 0) || (c < 0);
        bool has_pos = (a > 0) || (b > 0) || (c > 0);

        return result{!(has_neg && has_pos), face, mesh, p_onPlane};
    }

    return result{false, face, mesh, vecs[0]};
}

//gives the min and max of the box
std::vector<Eigen::Vector3f> getBoxLimits(std::vector<Tucano::Face> box, Tucano::Mesh mesh) {

    std::vector<Eigen::Vector3f> vecs;
    Eigen::Vector4f vec;

    if (box.size() == 0) {
        return vecs;
    }

    std::vector<GLuint> vecsofface;

    for (int i = 0; i < box.size(); i++) {
        for (int a = 0; a < box[i].vertex_ids.size(); i++) {
            vecsofface.push_back(box[i].vertex_ids[a]);
        }
    }

    float mean_x = 0;
    float mean_y = 0;
    float mean_z = 0;

    vec = mesh.getVertex(vecsofface[0]);

    float min_x = vec.x();
    float max_x = vec.x();

    float min_y = vec.y();
    float max_y = vec.y();

    float min_z = vec.z();
    float max_z = vec.z();

    for (int i = 0; i < vecsofface.size(); i++) {
        Eigen::Vector4f v = mesh.getVertex(vecsofface[i]);

        mean_x = mean_x + v.x();
        mean_y = mean_y + v.y();
        mean_z = mean_z + v.z();

        if (min_x > v.x()) {
            min_x = v.x();
        }
        if (max_x < v.x()) {
            max_x = v.x();
        }
        if (min_y > v.y()) {
            min_y = v.y();
        }
        if (max_y < v.y()) {
            max_y = v.y();
        }
        if (min_z > v.z()) {
            min_z = v.z();
        }
        if (max_z < v.z()) {
            max_z = v.z();
        }
    }
    mean_x = mean_x / vecsofface.size();
    mean_y = mean_y / vecsofface.size();
    mean_z = mean_z / vecsofface.size();

    vecs.push_back(Eigen::Vector3f(min_x, min_y, min_z));
    vecs.push_back(Eigen::Vector3f(max_x, max_y, max_z));
    vecs.push_back(Eigen::Vector3f(mean_x, mean_y, mean_z));
    return vecs;
}

//recursuve box
auto createboxes(std::vector<Tucano::Face> box, Tucano::Mesh mesh, std::vector<std::vector<Tucano::Face>> boxes) {
    std::vector<Tucano::Face> box1;
    std::vector<Tucano::Face> box2;

    std::vector<Eigen::Vector3f> boxlim = getBoxLimits(box, mesh);
    float xdiff = boxlim[1].x() - boxlim[0].x();
    float ydiff = boxlim[1].y() - boxlim[0].y();
    float zdiff = boxlim[1].z() - boxlim[0].z();

    float cut;

    if (box.size() < 200) {
        boxes.push_back(box);
        return boxes;
    }

    if (xdiff > ydiff && xdiff > zdiff) {
        cut = boxlim[2].x();
        for (int i = 0; i < box.size(); i++) {
            std::vector<GLuint> vecsofface = box[i].vertex_ids;
            for (int a = 0; a < vecsofface.size(); a++) {
                if (mesh.getVertex(vecsofface[a]).x() < cut) {
                    // i changed this code so it compares the vertecis of the face, since we dont have an operated == operator for the face
                    if (box1.vector::back().vertex_ids[0] != mesh.getFace(i).vertex_ids[0] &&
                        box1.vector::back().vertex_ids[1] != mesh.getFace(i).vertex_ids[1] &&
                        box1.vector::back().vertex_ids[2] != mesh.getFace(i).vertex_ids[2]) {
                        box1.push_back(mesh.getFace(i));
                    }
                } else {
                    //changed this the same as the code above
                    if (box2.vector::back().vertex_ids[0] != mesh.getFace(i).vertex_ids[0] &&
                        box2.vector::back().vertex_ids[1] != mesh.getFace(i).vertex_ids[1] &&
                        box2.vector::back().vertex_ids[2] != mesh.getFace(i).vertex_ids[2]) {
                        box2.push_back(mesh.getFace(i));
                    }
                }
            }
        }
    } else if (ydiff > xdiff && ydiff > zdiff) {
        cut = boxlim[2].y();
        for (int i = 0; i < box.size(); i++) {
            std::vector<GLuint> vecsofface = box[i].vertex_ids;
            for (int a = 0; a < vecsofface.size(); a++) {
                if (mesh.getVertex(vecsofface[a]).y() < cut) {
                    if (box1.vector::back().vertex_ids[0] != mesh.getFace(i).vertex_ids[0] &&
                        box1.vector::back().vertex_ids[1] != mesh.getFace(i).vertex_ids[1] &&
                        box1.vector::back().vertex_ids[2] != mesh.getFace(i).vertex_ids[2]) {
                        box1.push_back(mesh.getFace(i));
                    }
                } else {
                    if (box2.vector::back().vertex_ids[0] != mesh.getFace(i).vertex_ids[0] &&
                        box2.vector::back().vertex_ids[1] != mesh.getFace(i).vertex_ids[1] &&
                        box2.vector::back().vertex_ids[2] != mesh.getFace(i).vertex_ids[2]) {
                        box2.push_back(mesh.getFace(i));
                    }
                }
            }
        }
    } else {
        cut = boxlim[2].z();
        for (int i = 0; i < box.size(); i++) {
            std::vector<GLuint> vecsofface = box[i].vertex_ids;
            for (int a = 0; a < vecsofface.size(); a++) {
                if (mesh.getVertex(vecsofface[a]).z() < cut) {
                    if (box1.vector::back().vertex_ids[0] != mesh.getFace(i).vertex_ids[0] &&
                        box1.vector::back().vertex_ids[1] != mesh.getFace(i).vertex_ids[1] &&
                        box1.vector::back().vertex_ids[2] != mesh.getFace(i).vertex_ids[2]) {
                        box1.push_back(mesh.getFace(i));
                    }
                } else {
                    if (box2.vector::back().vertex_ids[0] != mesh.getFace(i).vertex_ids[0] &&
                        box2.vector::back().vertex_ids[1] != mesh.getFace(i).vertex_ids[1] &&
                        box2.vector::back().vertex_ids[2] != mesh.getFace(i).vertex_ids[2]) {
                        box2.push_back(mesh.getFace(i));
                    }
                }
            }
        }
    }

    createboxes(box1, mesh, boxes);
    createboxes(box2, mesh, boxes);

}

//first box
std::vector<std::vector<Tucano::Face>> firstBox(Tucano::Mesh mesh) {
    std::vector<Tucano::Face> box;
    std::vector<std::vector<Tucano::Face>> boxes;
    for (int i = 0; i < mesh.getNumberOfFaces(); i++) {
        box.push_back(mesh.getFace(i));
    }
    createboxes(box, mesh, boxes);
    return boxes;
}

// checks intersection with bounding box
auto intersectBox(Eigen::Vector3f start, Eigen::Vector3f to, std::vector<std::vector<Tucano::Face>> boxes,
                  std::vector<std::vector<Eigen::Vector3f>> boxbounds, int inter_node, Tucano::Mesh mesh) {
    struct result {
        bool inter;
        Tucano::Face face;
        Eigen::Vector4f hit;
    };


    float tmin_x = (boxbounds[inter_node][0].x() - start.x()) / to.x();
    float tmax_x = (boxbounds[inter_node][1].x() - start.x()) / to.x();

    float tmin_y = (boxbounds[inter_node][0].y() - start.y()) / to.y();
    float tmax_y = (boxbounds[inter_node][1].y() - start.y()) / to.y();

    float tmin_z = (boxbounds[inter_node][0].z() - start.z()) / to.z();
    float tmax_z = (boxbounds[inter_node][1].z() - start.z()) / to.z();

    float tin_x = std::min(tmin_x, tmax_x);
    float tout_x = std::max(tmin_x, tmax_x);

    float tin_y = std::min(tmin_y, tmax_y);
    float tout_y = std::max(tmin_y, tmax_y);

    float tin_z = std::min(tmin_z, tmax_z);
    float tout_z = std::max(tmin_z, tmax_z);

    float tin = 0;
    float tout = 1;

    if ((tin > tout) || (tout < 0)) {
        if (boxes.size() <= inter_node) {
            return result{false, boxes[0][0], Eigen::Vector4f(0, 0, 0, 0)};
        }
        return intersectBox(start, to, boxes, boxbounds, inter_node + 1, mesh);
    }

    for (int i = 0; i < boxes[inter_node].size(); i++) {
        auto ans = intersectTriange(start, to, boxes[inter_node][i], mesh);
        if (ans.inter) {
            return result{true, ans.face, ans.hit};
        }
    }
    return result{false, boxes[0][0], Eigen::Vector4f(0, 0, 0, 0)};

}

//INTERSECT RETURNS A STRUCTURE 
// intersect(start, to, mesh).inter -> returns true if there is a intersection
// intersect(start, to, mesh).face -> returns the face that it hit if it didnt hit returns a random face
// intersect(start, to, mesh).hit -> returns the vec3f that it hit if it didnt hit returns a random hit point
//intersect of one vector to the universe
// FOR THE ACCELERATION YOU NEED TO FEED INTERSECT METHOD WITH THOSE CREATE THEM BEFORE LOOPING THEY ARE EXPENSIVE
//std::vector<std::vector<Tucano::Face>> boxes = firstBox(mesh);
//std::vector<std::vector<Eigen::Vector3f>> boxbounds;
//for (int i = 0; i < boxes.size(); i++) {
//boxbounds.push_back(getBoxLimits(boxes[i], mesh));
//}
auto
intersect(Eigen::Vector3f start, Eigen::Vector3f to, Tucano::Mesh mesh, std::vector<std::vector<Tucano::Face>> boxes,
          std::vector<std::vector<Eigen::Vector3f>> boxbounds) {
    struct result {
        bool inter;
        Tucano::Face face;
        Eigen::Vector4f hit;
    };

    auto ans = intersectBox(start, to, boxes, boxbounds, 0, mesh);
    if (ans.inter) {
        return result{true, ans.face, Eigen::Vector3f(ans.hit.x(), ans.hit.y(), ans.hit.z())};
    }

    return result{false, boxes[0][0], Eigen::Vector3f(0, 0, 0, 0)};
}


Eigen::Vector3f reflect(const Eigen::Vector3f &I, const Eigen::Vector3f &N) {
    return I - 2 * I.dot(N) * N;
}

/* NOTE: whoever will implement hard shadow from a point, you need to check in the 2nd for loop, whether there is a face
 * between the light source and your current object. if so, you'll be in a shadow, otherwise, you're directly illuminated
Eigen::Vector3f directIllumination(const Eigen::Vector3f &I, const Eigen::Vector3f &N, std::vector<Eigen::Vector3f> lights,
                   Tucano::Mesh mesh) {
    std::vector<Eigen::Vector3f> lightvecs;
    for (int i = 0; i < lights.size(); i++) {
        Eigen::Vector3f current = (lights[i] - I);
        current = current.normalized();
        lightvecs.push_back(current);
    }
    for (int i = 0; i < lightvecs.size(); i++) {
        if (intersect(I, lightvecs[i], mesh).size() == 0) {

            float attenuation = 1 / (*dist);
            return lightColor * dot(N, L) * attenuation;

        }
    }
    if (!IsVisibile(I, L, dist)) return BLACK;

}
*/
Eigen::Vector3f shade(int level, Eigen::Vector3f start, Eigen::Vector3f to, Tucano::Mesh mesh,
                      Tucano::Effects::PhongMaterial phong, std::vector<Eigen::Vector3f> lights,
                      std::vector<std::vector<Tucano::Face>> boxes,
                      std::vector<std::vector<Eigen::Vector3f>> boxbounds) {
    //return empty vector which is just supposed to be black
    if (!intersect(start, to, mesh, boxes, boxbounds).inter) {
        return Eigen::Vector3f(0, 0, 0);
    }
    if (level == 0) {
        return start; //WE NEED to return color of current vertex we are at but we are just returning it's coords now just so it compiles
    }

}

Eigen::Vector3f recursiveraytracing(int level, Eigen::Vector3f start, Eigen::Vector3f to, Tucano::Mesh mesh,
                                    Tucano::Effects::PhongMaterial phong, std::vector<Eigen::Vector3f> lights,
                                    std::vector<std::vector<Tucano::Face>> boxes,
                                    std::vector<std::vector<Eigen::Vector3f>> boxbounds) {
    //return empty vector which is just supposed to be black
	auto intersection = intersect(start, to, mesh, boxes, boxbounds);
    if (!intersection.inter) {
        return Eigen::Vector3f(0, 0, 0);
    }
    if (level == 0) {
        return start; //WE NEED to return color of current vertex we are at but we are just returning it's coords now just so it compiles
    }
    return shade(level, intersection.hit,
                 reflect(intersection.hit - start,intersection.face.normal),
                 mesh, phong, lights, boxes, boxbounds); //either return just the color or after the shading
}


Eigen::Vector3f Flyscene::traceRay(Eigen::Vector3f &origin,
                                   Eigen::Vector3f &dest) {
    // just some fake random color per pixel until you implement your ray tracing
    // remember to return your RGB values as floats in the range [0, 1]!!!
    std::vector<std::vector<Tucano::Face>> boxes = firstBox(mesh);
    std::vector<std::vector<Eigen::Vector3f>> boxbounds;
    for (int i = 0; i < boxes.size(); i++) {
        boxbounds.push_back(getBoxLimits(boxes[i], mesh));
    }
    Eigen::Vector3f result = recursiveraytracing(1, origin, dest, mesh, phong, lights, boxes,
                                                 boxbounds);   //bounce one more ray

    return result;
}


void Flyscene::createDebugRay(const Eigen::Vector2f &mouse_pos) {
    ray.resetModelMatrix();
    // from pixel position to world coordinates
    Eigen::Vector3f screen_pos = flycamera.screenToWorld(mouse_pos);

    // direction from camera center to click position
    Eigen::Vector3f dir = (screen_pos - flycamera.getCenter()).normalized();

    // position and orient the cylinder representing the ray
    ray.setOriginOrientation(flycamera.getCenter(), dir);


    // place the camera representation (frustum) on current camera location,
    camerarep.resetModelMatrix();
    camerarep.setModelMatrix(flycamera.getViewMatrix().inverse());
}