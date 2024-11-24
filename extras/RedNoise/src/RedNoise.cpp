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
#include <glm/glm.hpp>

#define WIDTH 640
#define HEIGHT 480
#define SPECULAR_EXPONENT 256
#define WHITE Colour("White", 255, 255, 255)
#define BLACK Colour("Black", 0, 0, 0)
#define BLUE Colour("Blue", 0, 0, 255)
#define GREEN Colour("Green", 0, 255, 0)

struct vertexNormalPair
{
    std::vector<ModelTriangle> modelPoints;
    std::vector<glm::vec3> vertexNormals;
};

struct objectVertex
{
    glm::vec3 vertex;
    int count;

    objectVertex(const glm::vec3 &v, int c) : vertex(v), count(c) {}
};

int normaliseColour(float channel)
{
    return int((channel) * (255.0f / (714.0f)));
}

uint32_t returnColour(Colour c)
{
    return (255 << 24) + (c.red << 16) + (c.green << 8) + c.blue;
}

glm::vec3 returnLightingColour(Colour c, float diffuseBrightness, float specularBrightness, float ambientBrightness)
{
    int red = normaliseColour(c.red * (diffuseBrightness + ambientBrightness) + WHITE.red * specularBrightness);
    int green = normaliseColour(c.green * (diffuseBrightness + ambientBrightness) + WHITE.green * specularBrightness);
    int blue = normaliseColour(c.blue * (diffuseBrightness + ambientBrightness) + WHITE.blue * specularBrightness);
    return glm::vec3(red, green, blue);
}

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

void sortTrianglePointsHorizontal(CanvasTriangle &t)
{
    if (t.v1().x > t.v2().x)
    {
        std::swap(t.v1(), t.v2());
    }
}

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

void drawStrokedTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour)
{
    drawLine(window, triangle.v0(), triangle.v1(), colour);
    drawLine(window, triangle.v1(), triangle.v2(), colour);
    drawLine(window, triangle.v2(), triangle.v0(), colour);
}

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

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle t, Colour c, float depthArr[HEIGHT][WIDTH])
{
    sortTrianglePointsVertical(t);
    CanvasPoint p = interpolatePoint(t.v0(), t.v1(), t.v2());
    CanvasTriangle t1(t.v0(), p, t.v1());
    CanvasTriangle t2(t.v2(), p, t.v1());
    fillTriangle(window, t1, c, depthArr, 1);
    fillTriangle(window, t2, c, depthArr, -1);
}

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

glm::vec3 calculateNormal(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2)
{
    glm::vec3 e0 = v1 - v0;
    glm::vec3 e1 = v2 - v0;
    return glm::normalize(glm::cross(e0, e1));
}

glm::vec3 projectCanvasPointOntoVertex(CanvasPoint p, float focalLength, glm::vec3 camera, glm::mat3 cameraOrientation)
{
    float z = focalLength;
    float x = ((p.x - WIDTH / 2) / focalLength) * -1 * z;
    float y = ((p.y - HEIGHT / 2) / focalLength) * z;
    x /= 320;
    y /= 320;
    glm::vec3 out = camera - (cameraOrientation * glm::vec3(x, y, z));
    // return cameraOrientation * (camera - glm::vec3(x, y, z));
    return out;
}

vertexNormalPair readObjFile(std::string filePath, std::unordered_map<std::string, Colour> pallette, float scale)
{
    std::vector<ModelTriangle> out;
    std::vector<objectVertex> vertices;
    std::vector<glm::vec3> vertexNormals;
    Colour c = WHITE;
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
            objectVertex vertex(
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
                // ignore points[1] == texture vertex
                std::vector<std::string> points = split(lineSplit[i], '/');
                int v = std::stoi(points[0]) - 1;
                t.vertices[i - 1] = vertices[v].vertex;
                t.vertexIdx[i - 1] = v;
                vertices[v].count += 1;
            }
            glm::vec3 normal = calculateNormal(t.vertices[0], t.vertices[1], t.vertices[2]);
            t.normal = normal;
            t.colour = c;
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
        else if (operation == "usemtl")
        {
            std::string colourName = lineSplit[1];
            if (pallette.find(colourName) == pallette.end())
            {
                c = WHITE;
            }
            else
            {
                c = pallette[colourName];
            }
        }
    }
    file.close();
    return vertexNormalPair{out, vertexNormals};
}

std::unordered_map<std::string, Colour> readMtlFile(std::string filePath)
{
    std::unordered_map<std::string, Colour> pallette;
    std::ifstream file(filePath);
    std::string line;
    std::string operation;
    std::string colourName;

    float r, g, b;

    while (getline(file, line))
    {
        std::vector<std::string> lineSplit = split(line, ' ');
        operation = lineSplit[0];

        if (operation == "newmtl")
        {
            colourName = lineSplit[1];

            if (getline(file, line))
            {
                std::vector<std::string> lineSplit = split(line, ' ');
                r = std::stof(lineSplit[1]);
                g = std::stof(lineSplit[2]);
                b = std::stof(lineSplit[3]);
                Colour c(colourName, round(r * 255), round(g * 255), round(b * 255));
                pallette[colourName] = c;
            }
        }
    }
    file.close();
    return pallette;
}

glm::mat3 x_rotation(float angle)
{
    glm::mat3 rotation(1.0, 0.0, 0.0,
                       0.0, std::cos(angle), std::sin(angle),
                       0.0, -1.0 * std::sin(angle), std::cos(angle));
    return rotation;
}

glm::mat3 y_rotation(float angle)
{
    glm::mat3 rotation(std::cos(angle), 0.0, -1.0 * std::sin(angle),
                       0.0, 1.0, 0.0,
                       std::sin(angle), 0.0, std::cos(angle));
    return rotation;
}

bool notOnTheScreen(CanvasPoint p)
{
    if ((p.x >= WIDTH || p.x < 0) || (p.y >= HEIGHT || p.y < 0) || (p.depth >= 0))
        return true;
    else
        return false;
}

glm::mat3 lookAt(glm::vec3 v, glm::vec3 c)
{
    glm::vec3 forward = glm::normalize((c - v));
    glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0, 1, 0), forward));
    glm::vec3 up = glm::normalize(glm::cross(forward, right));
    return glm::mat3(right, up, forward);
}

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

bool hasShadow(std::vector<ModelTriangle> modelPoints, glm::vec3 light, glm::vec3 surfacePoint, int index, int x, int y)
{
    glm::vec3 d = glm::normalize(light - surfacePoint);
    float distance = glm::distance(surfacePoint, light);
    float margin = 0.01;
    RayTriangleIntersection t = getClosestValidIntersection(surfacePoint, d, modelPoints, index, margin);
    if (t.triangleIndex < modelPoints.size() && t.distanceFromCamera < distance)
        return true;
    else
        return false;
}

float calculateProximity(glm::vec3 light, glm::vec3 point)
{
    float distance = glm::distance(light, point);
    float diffuse = 10 / (3 * 3.14 * distance * distance);
    if (diffuse >= 1.0)
        return 1.0;
    else
        return diffuse;
}

float calculateAngleOfIncidence(glm::vec3 point, glm::vec3 light, glm::vec3 normal)
{
    glm::vec3 direction = glm::normalize(light - point);
    float angle = glm::dot(direction, normal);
    if (angle <= 0.0)
        return 0.0f;
    else
        return angle;
}

float useDiffuseLighting(glm::vec3 point, glm::vec3 light, glm::vec3 normal)
{
    float brightness = calculateProximity(light, point);
    brightness *= calculateAngleOfIncidence(point, light, normal);
    return brightness;
}

float useSpecularLighting(glm::vec3 point, glm::vec3 light, glm::vec3 camera, glm::vec3 normal)
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

glm::vec3 calculateLightingColour(glm::vec3 point, glm::vec3 light, glm::vec3 camera, glm::vec3 normal, RayTriangleIntersection r, std::vector<ModelTriangle> modelPoints, int x, int y)
{
    float diffuseBrightness = std::max(useDiffuseLighting(point, light, normal), 0.0f);
    float ambientBrightness = 0.5f;
    float specularBrightness = std::max(useSpecularLighting(point, light, camera, normal), 0.0f);
    if (hasShadow(modelPoints, light, point, r.triangleIndex, x, y))
    {
        glm::vec3 c = returnLightingColour(r.intersectedTriangle.colour, diffuseBrightness, 0.0, ambientBrightness);
        return c;
    }
    else
    {
        glm::vec3 c = returnLightingColour(r.intersectedTriangle.colour, diffuseBrightness, specularBrightness, ambientBrightness);
        return c;
    }
}

Colour phongShading(glm::vec3 point, glm::vec3 light, glm::vec3 camera, glm::vec3 n0, glm::vec3 n1, glm::vec3 n2, RayTriangleIntersection r, std::vector<ModelTriangle> modelPoints, int x, int y)
{
    ModelTriangle t = r.intersectedTriangle;
    glm::vec3 e0 = t.vertices[1] - t.vertices[0];
    glm::vec3 e1 = t.vertices[2] - t.vertices[0];

    float area = 0.5f * glm::length(glm::cross(e0, e1));
    float u = glm::length(0.5f * glm::cross((point - t.vertices[0]), (point - t.vertices[1]))) / area;
    float v = glm::length(0.5f * glm::cross((point - t.vertices[0]), (point - t.vertices[2]))) / area;
    float w = 1.0f - (v + u);

    glm::vec3 normalAtPoint = glm::normalize(w * n0 + v * n1 + u * n2);

    glm::vec3 c = calculateLightingColour(point, light, camera, normalAtPoint, r, modelPoints, x, y);
    return Colour(c.r, c.g, c.b);
}

Colour gouraudShading(glm::vec3 point, glm::vec3 light, glm::vec3 camera, glm::vec3 n0, glm::vec3 n1, glm::vec3 n2, RayTriangleIntersection r, std::vector<ModelTriangle> modelPoints)
{
    ModelTriangle t = r.intersectedTriangle;
    float area = 0.5f * glm::length(glm::cross((t.vertices[1] - t.vertices[0]), (t.vertices[2] - t.vertices[0])));
    float u = glm::length(0.5f * glm::cross((point - t.vertices[0]), (point - t.vertices[1]))) / area; // region between v1 and v0
    float v = glm::length(0.5f * glm::cross((point - t.vertices[0]), (point - t.vertices[2]))) / area; // region between v2 and v0
    float w = (float)1 - v - u;                                                                        // region between v2 and v1
    glm::vec3 v0 = calculateLightingColour(r.intersectedTriangle.vertices[0], light, camera, n0, r, modelPoints, 0, 0);
    glm::vec3 v1 = calculateLightingColour(r.intersectedTriangle.vertices[1], light, camera, n1, r, modelPoints, 0, 0);
    glm::vec3 v2 = calculateLightingColour(r.intersectedTriangle.vertices[2], light, camera, n2, r, modelPoints, 0, 0);

    glm::vec3 c = w * v0 + v * v1 + u * v2;
    return Colour(c.r, c.g, c.b);
}

void drawUsingWireframes(DrawingWindow &window, std::unordered_map<std::string, Colour> pallette, std::vector<ModelTriangle> modelPoints, glm::vec3 &camera, float focalLength, glm::mat3 &cameraOrientation, bool orbit)
{
    for (size_t i = 0; i < modelPoints.size(); i++)
    {
        ModelTriangle t = modelPoints[i];
        CanvasTriangle cT(
            projectVertexOntoCanvasPoint(camera, focalLength, t.vertices[0], cameraOrientation),
            projectVertexOntoCanvasPoint(camera, focalLength, t.vertices[1], cameraOrientation),
            projectVertexOntoCanvasPoint(camera, focalLength, t.vertices[2], cameraOrientation));
        if (!(notOnTheScreen(cT.v0()) && notOnTheScreen(cT.v1()) && notOnTheScreen(cT.v2())))
            drawStrokedTriangle(window, cT, t.colour);
    }
    if (orbit)
    {
        camera = y_rotation(0.01) * camera;
        cameraOrientation = lookAt(glm::vec3(0, 0, 0), camera);
    }
}

void drawUsingRasterisation(DrawingWindow &window, std::unordered_map<std::string, Colour> pallette, std::vector<ModelTriangle> modelPoints, glm::vec3 &camera, float focalLength, glm::mat3 &cameraOrientation, bool orbit)
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
            drawFilledTriangle(window, cT, t.colour, depthArr);
    }

    if (orbit)
    {
        camera = y_rotation(0.01) * camera;
        cameraOrientation = lookAt(glm::vec3(0, 0, 0), camera);
    }
}

void drawUsingRayTracing(DrawingWindow &window, std::unordered_map<std::string, Colour> pallette, std::vector<ModelTriangle> modelPoints, std::vector<glm::vec3> vertexNormals, glm::vec3 &camera, float focalLength, glm::vec3 light, glm::mat3 cameraOrientation)
{
    float margin = 0.0;
    for (int y = 0; y < HEIGHT; y++)
    {
        for (int x = 0; x < WIDTH; x++)
        {
            glm::vec3 p = projectCanvasPointOntoVertex(CanvasPoint(x, y), focalLength, camera, cameraOrientation);
            glm::vec3 d = glm::normalize(p - camera);
            RayTriangleIntersection r = getClosestValidIntersection(camera, d, modelPoints, -1, margin);
            glm::vec3 point = camera + r.distanceFromCamera * d;
            // glm::vec3 c = calculateLightingColour(point, light, camera, r.intersectedTriangle.normal, r, modelPoints, x, y);

            // Colour c = gouraudShading(point, light, camera, vertexNormals[r.intersectedTriangle.vertexIdx[0]], vertexNormals[r.intersectedTriangle.vertexIdx[1]], vertexNormals[r.intersectedTriangle.vertexIdx[2]], r, modelPoints);
            Colour c = phongShading(point, light, camera, vertexNormals[r.intersectedTriangle.vertexIdx[0]], vertexNormals[r.intersectedTriangle.vertexIdx[1]], vertexNormals[r.intersectedTriangle.vertexIdx[2]], r, modelPoints, x, y);
            if (r.triangleIndex < modelPoints.size())
                window.setPixelColour(x, y, returnColour(c));
            else
                window.setPixelColour(x, y, returnColour(BLACK));
        }
    }
}

void draw(DrawingWindow &window, std::unordered_map<std::string, Colour> pallette, std::vector<ModelTriangle> modelPoints, std::vector<glm::vec3> vertexNormals, glm::vec3 &camera, float focalLength, glm::mat3 &cameraOrientation, bool orbit, glm::vec3 light, int renderType)
{
    window.clearPixels();
    if (renderType == 1)
        drawUsingWireframes(window, pallette, modelPoints, camera, focalLength, cameraOrientation, orbit);
    else if (renderType == 2)
        drawUsingRasterisation(window, pallette, modelPoints, camera, focalLength, cameraOrientation, orbit);
    else if (renderType == 3)
        drawUsingRayTracing(window, pallette, modelPoints, vertexNormals, camera, focalLength, light, cameraOrientation);
    CanvasPoint p = projectVertexOntoCanvasPoint(camera, focalLength, light, cameraOrientation);
    window.setPixelColour(p.x, p.y, returnColour(BLUE));
}

void handleEvent(SDL_Event event, DrawingWindow &window, glm::vec3 &camera, glm::mat3 &cameraOrientation, bool &orbit, int &renderType, glm::vec3 &light)
{
    glm::vec3 x_translate(0.2, 0.0, 0.0);
    glm::vec3 y_translate(0.0, 0.2, 0.0);
    glm::vec3 z_translate(0.0, 0.0, 0.2);

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
            cameraOrientation = lookAt(glm::vec3(0, 0, 0), camera);
        }
        else if (event.key.keysym.sym == SDLK_d)
        {
            std::cout << "+5 deg y-axis" << std::endl;
            camera = y_rotation(angle) * camera;
            cameraOrientation = lookAt(glm::vec3(0, 0, 0), camera);
        }
        else if (event.key.keysym.sym == SDLK_s)
        {
            std::cout << "-5 deg x-axis" << std::endl;
            if (camera.z >= 0)
                camera = x_rotation(angle) * camera;
            else
                camera = x_rotation(-1 * angle) * camera;
            cameraOrientation = lookAt(glm::vec3(0, 0, 0), camera);
        }
        else if (event.key.keysym.sym == SDLK_w)
        {
            std::cout << "+5 deg x-axis" << std::endl;
            if (camera.z >= 0)
                camera = x_rotation(-1 * angle) * camera;
            else
                camera = x_rotation(angle) * camera;
            cameraOrientation = lookAt(glm::vec3(0, 0, 0), camera);
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
            light += glm::vec3(0.0, 0.1, 0.0);
        }
        else if (event.key.keysym.sym == SDLK_g)
        {
            std::cout << "Move Light down" << std::endl;
            light += glm::vec3(0.0, -0.1, 0.0);
        }
        else if (event.key.keysym.sym == SDLK_f)
        {
            std::cout << "Move Light left" << std::endl;
            light += glm::vec3(-0.1, 0.0, 0.0);
        }
        else if (event.key.keysym.sym == SDLK_h)
        {
            std::cout << "Move Light right" << std::endl;
            light += glm::vec3(0.1, 0.0, 0.0);
        }
        else if (event.key.keysym.sym == SDLK_v)
        {
            std::cout << "Move Light out" << std::endl;
            light += glm::vec3(0.0, 0.0, 0.1);
        }
        else if (event.key.keysym.sym == SDLK_b)
        {
            std::cout << "Move Light in" << std::endl;
            light += glm::vec3(0.0, 0.0, -0.1);
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
    std::unordered_map<std::string, Colour> pallette = readMtlFile("../../models/cornell-box.mtl");
    vertexNormalPair v = readObjFile("../../models/Custom-box.obj", pallette, 2.8);
    std::vector<ModelTriangle> modelPoints = v.modelPoints;
    std::vector<glm::vec3> vertexNormals = v.vertexNormals;
    glm::vec3 camera(0.0, 0.0, 4.0);
    glm::vec3 light(0.0, 0.55, 1.0);
    glm::mat3 cameraOrientation = glm::mat3(1.0, 0.0, 0.0,
                                            0.0, 1.0, 0.0,
                                            0.0, 0.0, 1.0);
    float focalLength = 2.0;
    bool orbit = false;
    int renderType = 2;

    std::cout << modelPoints[0].colour << std::endl;

    while (true)
    {
        // We MUST poll for events - otherwise the window will freeze !
        if (window.pollForInputEvents(event))
            handleEvent(event, window, camera, cameraOrientation, orbit, renderType, light);
        draw(window, pallette, modelPoints, vertexNormals, camera, focalLength, cameraOrientation, orbit, light, renderType);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }
}