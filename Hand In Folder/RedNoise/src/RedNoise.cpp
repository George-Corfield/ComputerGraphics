#include <CanvasTriangle.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <TextureMap.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <sstream>
#include <unordered_map>
#include <fstream>
#include <vector>
#include <typeinfo>
#include <glm/glm.hpp>

#define WIDTH 640
#define HEIGHT 480
#define SPECULAR_EXPONENT 256
// define colours for testing
#define WHITE Colour("White", 255, 255, 255)
#define BLACK Colour("Black", 0, 0, 0)
#define BLUE Colour("Blue", 0, 0, 255)
#define GREEN Colour("Green", 0, 255, 0)

// structure to return triangles and vertex normals in one
struct VertexProperties
{
    std::vector<ModelTriangle> modelPoints;
    std::vector<glm::vec3> Normals;
};

// structure to return colours and texture maps
struct Palettes
{
    std::unordered_map<std::string, Colour> colourPalette;
    std::vector<TextureMap> map;
};

// structure to hold vertices to calculate normals
struct ObjectVertex
{
    glm::vec3 vertex;
    int count;

    ObjectVertex(const glm::vec3 &v, int c) : vertex(v), count(c) {}
};

// creates 32 bit integer of colour
uint32_t returnColour(Colour c)
{
    return (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue;
}

// resolves colour into constituent components
Colour retrieveColour(uint32_t n)
{
    int r = (n & 0x00FF0000) >> 16;
    int g = (n & 0x0000FF00) >> 8;
    int b = (n & 0x000000FF);
    return Colour(r, g, b);
}

// Return vec3 of the colour components based on intensity between 0 and 255
glm::vec3 returnLightingColour(Colour c, float angleOfIncidenceIntensity, float proximityIntensity, float specularIntensity, float ambientIntensity)
{
    glm::vec3 objectColour(c.red, c.green, c.blue);
    glm::vec3 lightColour(WHITE.red, WHITE.green, WHITE.blue);
    glm::vec3 ambientComponent = ambientIntensity * objectColour;
    glm::vec3 diffuseComponent = angleOfIncidenceIntensity * proximityIntensity * objectColour;
    glm::vec3 specularComponent = specularIntensity * angleOfIncidenceIntensity * proximityIntensity * lightColour;
    return glm::clamp(ambientComponent + diffuseComponent + specularComponent, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(255.0f, 255.0f, 255.0f));
}

// sorts the triangle points vertically in place
void sortTrianglePointsVertical(CanvasTriangle &triangle)
{
    if (triangle.v0().y > triangle.v2().y)
    {
        std::swap(triangle.v0(), triangle.v2());
    }
    if (triangle.v0().y > triangle.v1().y)
    {
        std::swap(triangle.v0(), triangle.v1());
    }
    if (triangle.v1().y > triangle.v2().y)
    {
        std::swap(triangle.v1(), triangle.v2());
    }
}

// sorts triangle points 1 and 2 in place horizontally
void sortTrianglePointsHorizontal(CanvasTriangle &t)
{
    if (t.v1().x > t.v2().x)
    {
        std::swap(t.v1(), t.v2());
    }
}

// interpolates a point p on the line between a and c
CanvasPoint interpolatePoint(CanvasPoint a, CanvasPoint b, CanvasPoint c)
{
    float y1 = b.y - a.y;
    float y2 = c.y - a.y;
    float x2 = c.x - a.x;
    float d2 = c.depth - a.depth;
    float x1 = (y1 / y2) * x2;
    float d1 = (y1 / y2) * d2;
    return CanvasPoint(a.x + x1, b.y, a.depth + d1);
}

// interpolates a set of values between 2 numbers
std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues)
{
    float diff = to - from;
    float division = diff / (numberOfValues - 1);
    std::vector<float> out = {from};
    float current = from;
    for (size_t y = 0; y < numberOfValues - 1; y++)
    {
        out.push_back(current + division);
        current += division;
    }
    return out;
}

// interpolates a set of vec2's between 2 vec2's
std::vector<glm::vec2> interpolate2DVector(glm::vec2 from, glm::vec2 to, int n)
{
    glm::vec2 diff = to - from;
    glm::vec2 division = diff / float(n - 1);
    std::vector<glm::vec2> out = {from};
    glm::vec2 current = from;
    for (size_t y = 0; y < n - 1; y++)
    {
        out.push_back(current + division);
        current += division;
    }
    return out;
}

// calculates the normal to the triangle
glm::vec3 calculateNormal(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2)
{
    glm::vec3 e0 = v1 - v0;
    glm::vec3 e1 = v2 - v0;
    return glm::normalize(glm::cross(e0, e1));
}

// generates a random float between -1 and 1
float randomFloat()
{
    return ((float)rand() / (float)RAND_MAX) * 2.0f - 1.0f;
}

// generates a sample of lights in a sphere for soft shadowing
std::vector<glm::vec3> sampleLightSources(glm::vec4 light, int sampleSize)
{
    float radius = light.w;
    glm::vec3 center(light.x, light.y, light.z);
    std::vector<glm::vec3> out = {center};
    for (int i = 0; i < sampleSize - 1; i++)
    {
        glm::vec3 ratio(randomFloat() * radius, randomFloat() * radius, randomFloat() * radius);
        glm::vec3 newLight = center + ratio;
        out.push_back(newLight);
    }
    return out;
}

// draws a line between 2 canvas points
void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour c)
{
    float diffX = to.x - from.x;
    float diffY = to.y - from.y;
    float numberOfSteps = std::max(abs(diffX), abs(diffY));
    float xStepSize = diffX / numberOfSteps;
    float yStepSize = diffY / numberOfSteps;
    for (float i = 0.0; i <= numberOfSteps; i++)
    {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);
        window.setPixelColour(std::ceil(x), std::ceil(y), returnColour(c));
    }
}

// draws a raster between 2 points taking into account occlusion
void rasterise(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour c, float depthArr[HEIGHT][WIDTH])
{
    int y = std::round(from.y);
    std::vector<float> depth = interpolateSingleFloats(from.depth, to.depth, (std::ceil(to.x) - std::floor(from.x)) + 2);
    int start_x = std::floor(from.x);
    int end_x = std::ceil(to.x);
    for (int i = 0; i <= end_x - start_x; i++)
    {
        int x = start_x + i;
        if (x >= 0 && x < WIDTH)
        {
            if (depth[i] < depthArr[y][x])
            {
                window.setPixelColour(x, y, returnColour(c));
                depthArr[y][x] = depth[i];
            }
        }
    }
}

// draws an unfilled triangle
void drawStrokedTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour)
{
    drawLine(window, triangle.v0(), triangle.v1(), colour);
    drawLine(window, triangle.v1(), triangle.v2(), colour);
    drawLine(window, triangle.v2(), triangle.v0(), colour);
}

// fills a triangle with a colour using rasterise
void fillTriangle(DrawingWindow &window, CanvasTriangle t, Colour c, float depthArr[HEIGHT][WIDTH], int factor)
{
    sortTrianglePointsHorizontal(t);
    float initial_y = t.v0().y;
    float ydiff = std::round(std::abs(t.v0().y - t.v1().y)) + 1;
    std::vector<glm::vec2> left = interpolate2DVector(glm::vec2(t.v0().x, t.v0().depth), glm::vec2(t.v1().x, t.v1().depth), ydiff);
    std::vector<glm::vec2> right = interpolate2DVector(glm::vec2(t.v0().x, t.v0().depth), glm::vec2(t.v2().x, t.v2().depth), ydiff);
    for (int i = 0; i < ydiff; i++)
    {
        int y = initial_y + (i * factor);
        if (y >= 0 && y < HEIGHT - 1)
        {
            rasterise(window, CanvasPoint(left[i].x, y, left[i].y), CanvasPoint(right[i].x, y, right[i].y), c, depthArr);
        }
    }
}

// draws a filled triangle using fill triangle on 2 halves
void drawFilledTriangle(DrawingWindow &window, CanvasTriangle t, Colour c, float depthArr[HEIGHT][WIDTH])
{
    sortTrianglePointsVertical(t);
    CanvasPoint p = interpolatePoint(t.v0(), t.v1(), t.v2());
    CanvasTriangle t1(t.v0(), p, t.v1());
    CanvasTriangle t2(t.v2(), p, t.v1());
    fillTriangle(window, t1, c, depthArr, 1);
    fillTriangle(window, t2, c, depthArr, -1);
}

// projects a 3d point onto a 2d surface
CanvasPoint projectVertexOntoCanvasPoint(glm::vec3 camera, float focalLength, glm::vec3 vertexPosition, glm::mat3 cameraOrientation)
{
    float u, v;
    glm::vec3 pos = (vertexPosition - camera) * cameraOrientation;
    pos.x *= (float)320;
    pos.y *= (float)320;
    if (pos.z == 0)
        pos.z = -0.01;
    u = (focalLength * -1 * (pos.x / pos.z)) + (float)WIDTH / 2;
    v = (focalLength * (pos.y / pos.z)) + (float)HEIGHT / 2;

    return CanvasPoint(u, v, 1 / pos.z);
}

// detects if a point p is not on the screen based on width, height and depth
bool notOnTheScreen(CanvasPoint p)
{
    if ((p.x >= WIDTH || p.x < 0) || (p.y >= HEIGHT || p.y < 0) || (p.depth >= 0))
        return true;
    else
        return false;
}

// draws a wireframe scene
void drawUsingWireframes(DrawingWindow &window, std::vector<ModelTriangle> modelPoints, glm::vec3 &camera, float focalLength, glm::mat3 &cameraOrientation)
{
    for (size_t i = 0; i < modelPoints.size(); i++)
    {
        ModelTriangle t = modelPoints[i];
        CanvasTriangle cT(
            projectVertexOntoCanvasPoint(camera, focalLength, t.vertices[0], cameraOrientation),
            projectVertexOntoCanvasPoint(camera, focalLength, t.vertices[1], cameraOrientation),
            projectVertexOntoCanvasPoint(camera, focalLength, t.vertices[2], cameraOrientation));
        if (!(notOnTheScreen(cT.v0()) && notOnTheScreen(cT.v1()) && notOnTheScreen(cT.v2())))
        {
            if (t.colour.texture != -1)
                t.colour = WHITE;
            drawStrokedTriangle(window, cT, t.colour);
        }
    }
}

// draws a rasterised screen
void drawUsingRasterisation(DrawingWindow &window, std::vector<ModelTriangle> modelPoints, glm::vec3 &camera, float focalLength, glm::mat3 &cameraOrientation)
{

    float depthArr[HEIGHT][WIDTH] = {};
    for (size_t i = 0; i < modelPoints.size(); i++)
    {
        ModelTriangle t = modelPoints[i];
        CanvasTriangle cT(
            projectVertexOntoCanvasPoint(camera, focalLength, t.vertices[0], cameraOrientation),
            projectVertexOntoCanvasPoint(camera, focalLength, t.vertices[1], cameraOrientation),
            projectVertexOntoCanvasPoint(camera, focalLength, t.vertices[2], cameraOrientation));
        if (!(notOnTheScreen(cT.v0()) && notOnTheScreen(cT.v1()) && notOnTheScreen(cT.v2())))
        {
            if (t.colour.texture != -1)
                t.colour = WHITE;
            drawFilledTriangle(window, cT, t.colour, depthArr);
        }
    }
}

// calculates the barycentric coordinates of a point in a triangle
glm::vec3 calculateBarycentric(ModelTriangle t, glm::vec3 point)
{
    glm::vec3 e0 = t.vertices[1] - t.vertices[0];
    glm::vec3 e1 = t.vertices[2] - t.vertices[0];

    float area = 0.5f * glm::length(glm::cross(e0, e1));
    float u = glm::length(0.5f * glm::cross((point - t.vertices[0]), (point - t.vertices[1]))) / area;
    float v = glm::length(0.5f * glm::cross((point - t.vertices[0]), (point - t.vertices[2]))) / area;
    float w = 1.0f - (v + u);
    return glm::vec3(u, v, w);
}

// Gets the closest point in a direction from a source
RayTriangleIntersection getClosestValidIntersection(glm::vec3 source, glm::vec3 direction, std::vector<ModelTriangle> modelPoints, int index, float margin)
{
    RayTriangleIntersection solution;
    solution.distanceFromCamera = std::numeric_limits<float>::max();
    solution.triangleIndex = modelPoints.size() + 1;
    for (size_t i = 0; i < modelPoints.size(); i++)
    {
        ModelTriangle t = modelPoints[i];
        glm::vec3 e0 = t.vertices[1] - t.vertices[0];
        glm::vec3 e1 = t.vertices[2] - t.vertices[0];
        glm::vec3 SPVector = source - t.vertices[0];
        glm::mat3 DEMatrix(-direction, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
        if ((possibleSolution.y >= 0 && possibleSolution.y <= 1.0) && (possibleSolution.z >= 0 && possibleSolution.z <= 1.0) && ((possibleSolution.y + possibleSolution.z) <= 1.0) && (possibleSolution.x >= 0))
        {
            RayTriangleIntersection r(possibleSolution, possibleSolution.x, modelPoints[i], i);
            if (r.distanceFromCamera < solution.distanceFromCamera && r.triangleIndex != index && r.distanceFromCamera > margin)
            {
                solution = r;
            }
        }
    }
    return solution;
}

// uses closest interaction to detect if there is an object between a point and the light source
int hasShadow(std::vector<ModelTriangle> modelPoints, glm::vec3 light, glm::vec3 surfacePoint, int index)
{
    glm::vec3 d = glm::normalize(light - surfacePoint);
    float distance = glm::distance(surfacePoint, light);
    float margin = 0.01;
    RayTriangleIntersection t = getClosestValidIntersection(surfacePoint, d, modelPoints, index, margin);
    if (t.triangleIndex < modelPoints.size() && t.distanceFromCamera < distance)
        return 0;
    else
        return 1;
}

// calculates how close a point is to a light source based on 1/4pir^2
float calculateProximity(glm::vec3 light, glm::vec3 point)
{
    float distance = glm::distance(light, point);
    float diffuse = 10 / (4 * 3.14 * distance * distance);
    if (diffuse >= 1.0)
        return 1.0f;
    else
        return diffuse;
}

// calculates angle between point normal and light source
float calculateAngleOfIncidence(glm::vec3 point, glm::vec3 light, glm::vec3 normal)
{
    glm::vec3 direction = glm::normalize(light - point);
    float angle = glm::dot(direction, normal);
    if (angle <= 0.0)
        return 0.0f;
    else
        return angle;
}

// calculates specular degree beteen point, view and light source
float calculateSpecular(glm::vec3 point, glm::vec3 light, glm::vec3 camera, glm::vec3 normal)
{
    glm::vec3 rayOfIncidence = glm::normalize(point - light);
    glm::vec3 view = glm::normalize(point - camera);
    glm::vec3 rayOfReflection = glm::normalize(rayOfIncidence - (float)2.0 * normal * (glm::dot(rayOfIncidence, normal)));
    float brightness = std::abs(std::pow(glm::dot(view, rayOfReflection), SPECULAR_EXPONENT));
    if (brightness >= 1.0)
        return 1.0;
    else
        return brightness;
}

// calculates intensities and then uses returnLightingColour to generate a Colour for a point
glm::vec3 calculateLightingColour(glm::vec3 point, std::vector<glm::vec3> lights, glm::vec3 camera, glm::vec3 normal, RayTriangleIntersection r, std::vector<ModelTriangle> modelPoints, bool flatShading)
{
    float angleOfIncidenceIntensity = calculateAngleOfIncidence(point, lights[0], normal);
    float proximityIntensity = calculateProximity(lights[0], point);
    float specularIntensity = calculateSpecular(point, lights[0], camera, normal);
    float ambientIntensity = 0.3f;
    float shadowFactor = 0.0f;
    if (flatShading)
    {
        for (glm::vec3 light : lights)
        {
            shadowFactor += hasShadow(modelPoints, light, point, r.triangleIndex);
        }
        shadowFactor /= lights.size();
        return returnLightingColour(r.intersectedTriangle.colour, shadowFactor * angleOfIncidenceIntensity, proximityIntensity, shadowFactor * specularIntensity, ambientIntensity);
    }
    if (hasShadow(modelPoints, lights[0], point, r.triangleIndex) == 0)
        return returnLightingColour(r.intersectedTriangle.colour, angleOfIncidenceIntensity, proximityIntensity, 0.0f, ambientIntensity);
    return returnLightingColour(r.intersectedTriangle.colour, angleOfIncidenceIntensity, proximityIntensity, specularIntensity, ambientIntensity);
}

// flat shades the point, if it has a texture the texture is mapped onto the point
Colour flatShading(glm::vec3 point, std::vector<glm::vec3> lights, glm::vec3 camera, RayTriangleIntersection r, std::vector<ModelTriangle> modelPoints, std::vector<TextureMap> map)
{
    ModelTriangle t = r.intersectedTriangle;
    if (t.colour.texture != -1)
    {
        TextureMap m = map[t.colour.texture];
        r.intersectedTriangle.colour = WHITE;
        glm::vec3 barycentric = calculateBarycentric(t, point);
        glm::vec2 p = barycentric.x * t.texturePoints[2] + barycentric.y * t.texturePoints[1] + barycentric.z * t.texturePoints[0];
        glm::vec2 texture_point(std::floor(p.x * (float)m.width), std::floor(p.y * (float)m.height));
        Colour c(retrieveColour(m.pixels[texture_point.y * m.width + texture_point.x]));
        r.intersectedTriangle.colour = c;
    }
    glm::vec3 c = calculateLightingColour(point, lights, camera, r.intersectedTriangle.normal, r, modelPoints, true);
    return Colour(c.r, c.g, c.b);
}

// uses the phong shading model based on barycentric coordinates
Colour phongShading(glm::vec3 point, std::vector<glm::vec3> lights, glm::vec3 camera, glm::vec3 n0, glm::vec3 n1, glm::vec3 n2, RayTriangleIntersection r, std::vector<ModelTriangle> modelPoints)
{
    ModelTriangle t = r.intersectedTriangle;
    glm::vec3 barycentric = calculateBarycentric(t, point);

    glm::vec3 normalAtPoint = glm::normalize(barycentric.z * n0 + barycentric.y * n1 + barycentric.x * n2);

    glm::vec3 c = calculateLightingColour(point, lights, camera, normalAtPoint, r, modelPoints, false);
    return Colour(c.r, c.g, c.b);
}

// uses gouraud shading model to calculate lighting at the point
Colour gouraudShading(glm::vec3 point, std::vector<glm::vec3> lights, glm::vec3 camera, glm::vec3 n0, glm::vec3 n1, glm::vec3 n2, RayTriangleIntersection r, std::vector<ModelTriangle> modelPoints)
{
    ModelTriangle t = r.intersectedTriangle;
    glm::vec3 barycentric = calculateBarycentric(t, point);
    glm::vec3 v0 = calculateLightingColour(r.intersectedTriangle.vertices[0], lights, camera, n0, r, modelPoints, false);
    glm::vec3 v1 = calculateLightingColour(r.intersectedTriangle.vertices[1], lights, camera, n1, r, modelPoints, false);
    glm::vec3 v2 = calculateLightingColour(r.intersectedTriangle.vertices[2], lights, camera, n2, r, modelPoints, false);

    glm::vec3 c = barycentric.z * v0 + barycentric.y * v1 + barycentric.x * v2;
    return Colour(c.r, c.g, c.b);
}

// calculates what the reflection of a point would be to create a mirror effect on a surface
Colour mirrorShading(glm::vec3 point, std::vector<glm::vec3> lights, glm::vec3 camera, glm::vec3 n0, glm::vec3 n1, glm::vec3 n2, RayTriangleIntersection r, std::vector<ModelTriangle> modelPoints, std::vector<glm::vec3> vertexNormals, std::vector<TextureMap> map)
{
    glm::vec3 rayOfIncidence = glm::normalize(point - camera);
    glm::vec3 barycentric = calculateBarycentric(r.intersectedTriangle, point);
    glm::vec3 normalAtPoint = glm::normalize(barycentric.z * n0 + barycentric.y * n1 + barycentric.x * n2);
    glm::vec3 rayOfReflection = glm::normalize(rayOfIncidence - (float)2.0 * normalAtPoint * (glm::dot(rayOfIncidence, normalAtPoint)));
    RayTriangleIntersection reflectedTriangle = getClosestValidIntersection(point, rayOfReflection, modelPoints, r.triangleIndex, 0.05);
    if (reflectedTriangle.triangleIndex < modelPoints.size())
    {
        glm::vec3 reflectedPoint = point + reflectedTriangle.distanceFromCamera * rayOfReflection;
        if (reflectedTriangle.intersectedTriangle.shadingType == "phong")
            return phongShading(reflectedPoint,
                                lights,
                                camera,
                                vertexNormals[reflectedTriangle.intersectedTriangle.vertexIdx[0]],
                                vertexNormals[reflectedTriangle.intersectedTriangle.vertexIdx[1]],
                                vertexNormals[reflectedTriangle.intersectedTriangle.vertexIdx[2]],
                                reflectedTriangle,
                                modelPoints);
        else if (reflectedTriangle.intersectedTriangle.shadingType == "gouraud")
            return gouraudShading(reflectedPoint,
                                  lights,
                                  camera,
                                  vertexNormals[reflectedTriangle.intersectedTriangle.vertexIdx[0]],
                                  vertexNormals[reflectedTriangle.intersectedTriangle.vertexIdx[1]],
                                  vertexNormals[reflectedTriangle.intersectedTriangle.vertexIdx[2]],
                                  reflectedTriangle,
                                  modelPoints);
        else
            return flatShading(reflectedPoint, lights, camera, reflectedTriangle, modelPoints, map);
    }
    else
        return BLACK;
}

// projects a 2d point onto a 3d vector
glm::vec3 projectCanvasPointOntoVertex(CanvasPoint p, float focalLength, glm::vec3 camera, glm::mat3 cameraOrientation)
{
    float z = focalLength;
    float x = ((p.x - WIDTH / 2) / focalLength) * -1 * z;
    float y = ((p.y - HEIGHT / 2) / focalLength) * z;
    x /= 320;
    y /= 320;
    glm::vec3 out = camera - (cameraOrientation * glm::vec3(x, y, z));
    return out;
}

// uses raytracing to draw a 3d scene
void drawUsingRayTracing(DrawingWindow &window, std::vector<ModelTriangle> modelPoints, std::vector<glm::vec3> vertexNormals, std::vector<TextureMap> map, glm::vec3 &camera, float focalLength, std::vector<glm::vec3> lights, glm::mat3 cameraOrientation)
{
    float margin = 0.0;
    for (int y = 0; y < HEIGHT; y++)
    {
        for (int x = 0; x < WIDTH; x++)
        {
            glm::vec3 p = projectCanvasPointOntoVertex(CanvasPoint(x, y), focalLength, camera, cameraOrientation);
            glm::vec3 d = glm::normalize(p - camera);
            RayTriangleIntersection r = getClosestValidIntersection(camera, d, modelPoints, -1, margin);
            if (r.triangleIndex < modelPoints.size() && r.triangleIndex >= 0)
            {
                glm::vec3 point = camera + r.distanceFromCamera * d;
                Colour c;
                if (r.intersectedTriangle.shadingType == "phong")
                    c = phongShading(point, lights, camera, vertexNormals[r.intersectedTriangle.vertexIdx[0]], vertexNormals[r.intersectedTriangle.vertexIdx[1]], vertexNormals[r.intersectedTriangle.vertexIdx[2]], r, modelPoints);
                else if (r.intersectedTriangle.shadingType == "gouraud")
                    c = gouraudShading(point, lights, camera, vertexNormals[r.intersectedTriangle.vertexIdx[0]], vertexNormals[r.intersectedTriangle.vertexIdx[1]], vertexNormals[r.intersectedTriangle.vertexIdx[2]], r, modelPoints);
                else if (r.intersectedTriangle.shadingType == "mirror")
                    c = mirrorShading(point, lights, camera, vertexNormals[r.intersectedTriangle.vertexIdx[0]], vertexNormals[r.intersectedTriangle.vertexIdx[1]], vertexNormals[r.intersectedTriangle.vertexIdx[2]], r, modelPoints, vertexNormals, map);
                else
                    c = flatShading(point, lights, camera, r, modelPoints, map);
                window.setPixelColour(x, y, returnColour(c));
            }
            else
                window.setPixelColour(x, y, returnColour(BLACK));
        }
    }
}

// draws a scene based on the rendertype, default is Rasterisation
void draw(DrawingWindow &window, std::vector<ModelTriangle> modelPoints, std::vector<glm::vec3> vertexNormals, std::vector<TextureMap> map, glm::vec3 &camera, float focalLength, glm::mat3 &cameraOrientation, std::vector<glm::vec3> lights, int renderType)
{
    window.clearPixels();
    if (renderType == 1)
        drawUsingWireframes(window, modelPoints, camera, focalLength, cameraOrientation);
    else if (renderType == 2)
        drawUsingRasterisation(window, modelPoints, camera, focalLength, cameraOrientation);
    else if (renderType == 3)
        drawUsingRayTracing(window, modelPoints, vertexNormals, map, camera, focalLength, lights, cameraOrientation);
}

// Reads in an object file and returns the model triangles and vertex normals
VertexProperties readObjFile(std::string filePath, std::unordered_map<std::string, Colour> palette, float scale)
{
    std::vector<ModelTriangle> out;
    std::vector<ObjectVertex> vertices;
    std::vector<glm::vec3> vertexNormals;
    std::vector<glm::vec2> vertexTexture;
    Colour c = WHITE;
    std::string shadingType = "flat";
    bool hasVertexNormals = false;

    std::ifstream file(filePath);
    std::string line;
    std::string operation;
    while (getline(file, line))
    {
        std::vector<std::string> lineSplit = split(line, ' ');
        operation = lineSplit[0];
        if (operation == "v")
        {
            ObjectVertex vertex(
                glm::vec3(
                    std::stof(lineSplit[1]),
                    std::stof(lineSplit[2]),
                    std::stof(lineSplit[3])),
                0);

            vertex.vertex /= scale;
            vertices.push_back(vertex);
            if (!hasVertexNormals)
                vertexNormals.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
        }
        else if (operation == "f")
        {
            ModelTriangle t;
            for (int i = 1; i < 4; i++)
            {
                std::vector<std::string> points = split(lineSplit[i], '/');
                int v = std::stoi(points[0]) - 1;
                t.vertices[i - 1] = vertices[v].vertex;
                t.vertexIdx[i - 1] = v;
                if (points[1] != "")
                {
                    t.texturePoints[i - 1] = vertexTexture[std::stoi(points[1]) - 1];
                }
                vertices[v].count += 1;
            }
            glm::vec3 normal = calculateNormal(t.vertices[0], t.vertices[1], t.vertices[2]);
            t.normal = normal;
            t.colour = c;
            t.shadingType = shadingType;
            out.push_back(t);
            if (!hasVertexNormals)
            {
                for (int i = 0; i < 3; i++)
                {
                    int n = vertices[t.vertexIdx[i]].count;
                    if (n == 0)
                        vertexNormals[t.vertexIdx[i]] = normal;
                    else
                    {
                        vertexNormals[t.vertexIdx[i]] = glm::normalize((((float)n / (float)(n + 1)) * vertexNormals[t.vertexIdx[i]]) + ((float)1 / (float)(n + 1)) * normal);
                    }
                    vertices[t.vertexIdx[i]].count += 1;
                }
            }
        }
        else if (operation == "vn")
        {
            vertexNormals.push_back(
                glm::vec3(
                    std::stof(lineSplit[1]),
                    std::stof(lineSplit[2]),
                    std::stof(lineSplit[3])));
            hasVertexNormals = true;
        }
        else if (operation == "vt")
        {
            vertexTexture.push_back(
                glm::vec2(
                    std::stof(lineSplit[1]),
                    std::stof(lineSplit[2])));
        }
        else if (operation == "s")
        {
            shadingType = lineSplit[1];
        }
        else if (operation == "usemtl")
        {
            std::string mtlName = lineSplit[1];
            if (palette.find(mtlName) == palette.end())
            {
                c = WHITE;
            }
            else
            {
                c = palette[mtlName];
            }
        }
    }
    file.close();
    return VertexProperties{out, vertexNormals};
}

// Reads an mtl file and returns the colour dictionary and the texture map
Palettes readMtlFile(std::string filePath)
{
    std::unordered_map<std::string, Colour> colourPalette;
    std::vector<TextureMap> maps;
    std::ifstream file(filePath);
    std::string line;
    std::string operation;
    std::string mtlName;

    float r, g, b;

    while (getline(file, line))
    {
        std::vector<std::string> lineSplit = split(line, ' ');
        operation = lineSplit[0];

        if (operation == "newmtl")
        {
            mtlName = lineSplit[1];

            if (getline(file, line))
            {
                std::vector<std::string> lineSplit = split(line, ' ');
                r = std::stof(lineSplit[1]);
                g = std::stof(lineSplit[2]);
                b = std::stof(lineSplit[3]);
                Colour c(mtlName, round(r * 255), round(g * 255), round(b * 255));
                c.texture = -1;
                colourPalette[mtlName] = c;
            }
        }
        else if (operation == "map_Kd")
        {
            std::string fileLocation = "./models/" + lineSplit[1];
            TextureMap t(fileLocation);
            maps.push_back(t);
            colourPalette[mtlName].texture = maps.size() - 1;
        }
    }
    file.close();
    return Palettes{colourPalette, maps};
}

// calculates camera orientation for the camera to look at a point
glm::mat3 lookAt(glm::vec3 v, glm::vec3 c)
{
    glm::vec3 forward = glm::normalize((c - v));
    glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0, 1, 0), forward));
    glm::vec3 up = glm::normalize(glm::cross(forward, right));
    return glm::mat3(right, up, forward);
}

// generates an x rotation for a particular angle
glm::mat3 x_rotation(float angle)
{
    glm::mat3 rotation(1.0, 0.0, 0.0,
                       0.0, std::cos(angle), std::sin(angle),
                       0.0, -1.0 * std::sin(angle), std::cos(angle));
    return rotation;
}

// generates a y rotation for a particular angle
glm::mat3 y_rotation(float angle)
{
    glm::mat3 rotation(std::cos(angle), 0.0, -1.0 * std::sin(angle),
                       0.0, 1.0, 0.0,
                       std::sin(angle), 0.0, std::cos(angle));
    return rotation;
}

// generates a z rotation for a particular angle
glm::mat3 z_rotation(float angle)
{
    glm::mat3 rotation(std::cos(angle), std::sin(angle), 0.0,
                       -1.0 * std::sin(angle), std::cos(angle), 0.0,
                       0.0, 0.0, 1.0);
    return rotation;
}

// handles all events in a render
void handleEvent(SDL_Event event,
                 DrawingWindow &window,
                 glm::vec3 &camera,
                 glm::mat3 &cameraOrientation,
                 float &focalLength,
                 bool &orbit,
                 int &renderType,
                 glm::vec4 &light,
                 std::vector<glm::vec3> &sampleLights,
                 bool &lookAtEnabled)
{
    glm::vec3 x_translate(0.2, 0.0, 0.0);
    glm::vec3 y_translate(0.0, 0.2, 0.0);
    glm::vec3 z_translate(0.0, 0.0, 0.1);

    float angle = 0.087;

    if (event.type == SDL_KEYDOWN)
    {
        if (event.key.keysym.sym == SDLK_LEFT)
        {
            std::cout << "LEFT" << std::endl;
            camera += x_translate;
        }
        else if (event.key.keysym.sym == SDLK_RIGHT)
        {
            std::cout << "RIGHT" << std::endl;
            camera -= x_translate;
        }
        else if (event.key.keysym.sym == SDLK_UP)
        {
            std::cout << "UP" << std::endl;
            camera += y_translate;
        }
        else if (event.key.keysym.sym == SDLK_DOWN)
        {
            std::cout << "DOWN" << std::endl;
            camera -= y_translate;
        }
        else if (event.key.keysym.sym == SDLK_i)
        {
            std::cout << "FORWARD" << std::endl;
            camera -= z_translate;
        }
        else if (event.key.keysym.sym == SDLK_k)
        {
            std::cout << "BACKWARD" << std::endl;
            camera += z_translate;
        }
        else if (event.key.keysym.sym == SDLK_a)
        {
            std::cout << "-5 deg y-axis" << std::endl;
            camera = y_rotation(-1 * angle) * camera;
            if (lookAtEnabled)
                cameraOrientation = lookAt(glm::vec3(0, 0, 0), camera);
        }
        else if (event.key.keysym.sym == SDLK_d)
        {
            std::cout << "+5 deg y-axis" << std::endl;
            camera = y_rotation(angle) * camera;
            if (lookAtEnabled)
                cameraOrientation = lookAt(glm::vec3(0, 0, 0), camera);
        }
        else if (event.key.keysym.sym == SDLK_s)
        {
            std::cout << "-5 deg x-axis" << std::endl;
            camera = x_rotation(angle) * camera;
            if (lookAtEnabled)
                cameraOrientation = lookAt(glm::vec3(0, 0, 0), camera);
        }
        else if (event.key.keysym.sym == SDLK_w)
        {
            std::cout << "+5 deg x-axis" << std::endl;
            camera = x_rotation(-1 * angle) * camera;
            if (lookAtEnabled)
                cameraOrientation = lookAt(glm::vec3(0, 0, 0), camera);
        }
        else if (event.key.keysym.sym == SDLK_u)
        {
            std::cout << "Tilt up" << std::endl;
            cameraOrientation *= x_rotation(angle * 0.1f);
        }
        else if (event.key.keysym.sym == SDLK_j)
        {
            std::cout << "Tilt down" << std::endl;
            cameraOrientation *= x_rotation(-0.1f * angle);
        }
        else if (event.key.keysym.sym == SDLK_n)
        {
            std::cout << "Pan left" << std::endl;
            cameraOrientation *= y_rotation(0.1f * angle);
        }
        else if (event.key.keysym.sym == SDLK_m)
        {
            std::cout << "Pan right" << std::endl;
            cameraOrientation *= y_rotation(-0.1f * angle);
        }
        else if (event.key.keysym.sym == SDLK_o)
        {
            if (orbit)
            {
                std::cout << "Orbit Disabled" << std::endl;
                orbit = false;
            }
            else
            {
                std::cout << "Orbit Enabled" << std::endl;
                orbit = true;
            }
        }
        else if (event.key.keysym.sym == SDLK_0)
        {
            if (lookAtEnabled)
            {
                std::cout << "Look At Disabled" << std::endl;
                lookAtEnabled = false;
            }
            else
            {
                std::cout << "Look At Enabled" << std::endl;
                lookAtEnabled = true;
            }
        }
        else if (event.key.keysym.sym == SDLK_1)
        {
            std::cout << "WireFrame Render" << std::endl;
            renderType = 1;
        }
        else if (event.key.keysym.sym == SDLK_2)
        {
            std::cout << "Rasterising Render" << std::endl;
            renderType = 2;
        }
        else if (event.key.keysym.sym == SDLK_3)
        {
            std::cout << "Raytrace Render" << std::endl;
            renderType = 3;
        }
        else if (event.key.keysym.sym == SDLK_t)
        {
            std::cout << "Move Light up" << std::endl;
            light += glm::vec4(0.0, 0.1, 0.0, 0.0);
            sampleLights = sampleLightSources(light, sampleLights.size());
        }
        else if (event.key.keysym.sym == SDLK_g)
        {
            std::cout << "Move Light down" << std::endl;
            light += glm::vec4(0.0, -0.1, 0.0, 0.0);
            sampleLights = sampleLightSources(light, sampleLights.size());
        }
        else if (event.key.keysym.sym == SDLK_f)
        {
            std::cout << "Move Light left" << std::endl;
            light += glm::vec4(-0.1, 0.0, 0.0, 0.0);
            sampleLights = sampleLightSources(light, sampleLights.size());
        }
        else if (event.key.keysym.sym == SDLK_h)
        {
            std::cout << "Move Light right" << std::endl;
            light += glm::vec4(0.1, 0.0, 0.0, 0.0);
            sampleLights = sampleLightSources(light, sampleLights.size());
        }
        else if (event.key.keysym.sym == SDLK_v)
        {
            std::cout << "Move Light out" << std::endl;
            light += glm::vec4(0.0, 0.0, 0.1, 0.0);
            sampleLights = sampleLightSources(light, sampleLights.size());
        }
        else if (event.key.keysym.sym == SDLK_b)
        {
            std::cout << "Move Light in" << std::endl;
            light += glm::vec4(0.0, 0.0, -0.1, 0.0);
            sampleLights = sampleLightSources(light, sampleLights.size());
        }
        else if (event.key.keysym.sym == SDLK_8)
        {
            std::cout << "Hard Shadows" << std::endl;
            sampleLights = sampleLightSources(light, 1);
        }
        else if (event.key.keysym.sym == SDLK_9)
        {
            std::cout << "Soft Shadows" << std::endl;
            sampleLights = sampleLightSources(light, 10);
        }
        else if (event.key.keysym.sym == SDLK_p)
        {
            std::cout << "-------------" << std::endl;
            std::cout << "light Pos: " << light.x << ";" << light.y << ";" << light.z << std::endl;
            std::cout << "Camera Pos: " << camera.x << ";" << camera.y << ";" << camera.z << std::endl;
            std::cout << "Camera Orientation: " << std::endl;
            std::cout << "Camera Right: " << cameraOrientation[0].x << ";" << cameraOrientation[0].y << ";" << cameraOrientation[0].z << std::endl;
            std::cout << "Camera Up: " << cameraOrientation[1].x << ";" << cameraOrientation[1].y << ";" << cameraOrientation[1].z << std::endl;
            std::cout << "Camera Forward: " << cameraOrientation[2].x << ";" << cameraOrientation[2].y << ";" << cameraOrientation[2].z << std::endl;
        }
    }
    else if (event.type == SDL_MOUSEBUTTONDOWN)
    {
        int mouseX, mouseY;
        SDL_GetMouseState(&mouseX, &mouseY);
        std::cout << mouseX << ";" << mouseY << std::endl;
    }
}

int main(int argc, char *argv[])
{
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
    SDL_Event event;
    Palettes palettes = readMtlFile("./models/material.mtl");
    std::unordered_map<std::string, Colour> colourPalette = palettes.colourPalette;
    std::vector<TextureMap> map = palettes.map;
    VertexProperties v = readObjFile("./models/coursework-box.obj", colourPalette, 3);
    std::vector<ModelTriangle> modelPoints = v.modelPoints;
    std::vector<glm::vec3> vertexNormals = v.Normals;
    glm::vec3 camera(0.0, 0.0, 4.0);
    glm::vec4 light(0.0, 0.5, 0.5, 0.05); // describes (x,y,z) of light and the radius for soft shadows
    std::vector<glm::vec3> sampledLights = sampleLightSources(light, 1);
    glm::mat3 cameraOrientation = glm::mat3(1.0, 0.0, 0.0,
                                            0.0, 1.0, 0.0,
                                            0.0, 0.0, 1.0);
    float focalLength = 2.0;
    bool orbit = false;
    int renderType = 2;
    bool lookAtEnabled = true;

    while (true)
    {
        // We MUST poll for events - otherwise the window will freeze !
        if (window.pollForInputEvents(event))
            handleEvent(event, window, camera, cameraOrientation, focalLength, orbit, renderType, light, sampledLights, lookAtEnabled);
        if (orbit)
        {
            camera = y_rotation(0.005) * camera;
            cameraOrientation = lookAt(glm::vec3(0, 0, 0), camera);
        }
        draw(window, modelPoints, vertexNormals, map, camera, focalLength, cameraOrientation, sampledLights, renderType);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }
}