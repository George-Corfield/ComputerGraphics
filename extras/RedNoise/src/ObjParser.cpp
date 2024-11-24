#include <ModelTriangle.h>
#include <Utils.h>
#include <vector>
#include <glm/glm.hpp>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <Colour.h>

#define WHITE Colour("White", 255, 255, 255)

struct objectVertex
{
    glm::vec3 vertex;
    int count;

    objectVertex(const glm::vec3 &v, int c) : vertex(v), count(c) {}
};

struct vertexNormalPair
{
    std::vector<ModelTriangle> modelPoints;
    std::vector<glm::vec3> vertexNormals;
    std::vector<objectVertex> vertices;
};

glm::vec3 calculateNormal(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2)
{
    glm::vec3 e0 = v1 - v0;
    glm::vec3 e1 = v2 - v0;
    return glm::normalize(glm::cross(e0, e1));
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

vertexNormalPair readObjFile(std::string filePath, std::unordered_map<std::string, Colour> pallette)
{
    std::vector<ModelTriangle> out;
    std::vector<objectVertex> vertices;
    std::vector<glm::vec3> vertexNormals;
    Colour c;
    glm::mat3 xRotation = x_rotation(-1.5708f);
    glm::vec3 translation(0.65f, -2.85f, 1.4f);
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
            vertex.vertex /= 113.0f / 2.5f;
            vertex.vertex = xRotation * vertex.vertex;
            vertex.vertex += translation;
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
            }
            glm::vec3 normal = calculateNormal(t.vertices[0], t.vertices[1], t.vertices[2]);
            t.normal = normal;
            t.colour = c;
            t.shadingType = "phong";
            out.push_back(t);
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
    // std::cout << vertices[12].vertex.x << ";" << vertices[12].vertex.y << ";" << vertices[12].vertex.z << std::endl;
    // std::cout << vertices[13].vertex.x << ";" << vertices[13].vertex.y << ";" << vertices[13].vertex.z << std::endl;
    // std::cout << vertices[1].vertex.x << ";" << vertices[1].vertex.y << ";" << vertices[1].vertex.z << std::endl;
    // std::cout << vertices[108].vertex.x << ";" << vertices[108].vertex.y << ";" << vertices[108].vertex.z << std::endl;
    // std::cout << vertices[31].vertex.x << ";" << vertices[31].vertex.y << ";" << vertices[31].vertex.z << std::endl;
    // std::cout << vertices[2].vertex.x << ";" << vertices[2].vertex.y << ";" << vertices[2].vertex.z << std::endl;

    return vertexNormalPair{out, vertexNormals, vertices};
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

void writeToNewFile(std::string newFileName, vertexNormalPair v, int offset)
{
    std::vector<ModelTriangle> t = v.modelPoints;
    std::vector<glm::vec3> vertexNormals = v.vertexNormals;
    std::vector<objectVertex> vertices = v.vertices;

    std::fstream file(newFileName, std::fstream::app);
    std::string line;
    std::string operation;
    std::string colourName = "None";
    file << "\n";
    file << "o bunny_low\n";

    for (int i = 0; i < vertices.size(); i++)
    {
        line = "vn " + std::to_string(vertexNormals[i].x) + " " + std::to_string(vertexNormals[i].y) + " " + std::to_string(vertexNormals[i].z) + "\n";
        file << line;
        line = "v " + std::to_string(vertices[i].vertex.x) + " " + std::to_string(vertices[i].vertex.y) + " " + std::to_string(vertices[i].vertex.z) + "\n";
        file << line;
    }

    for (int i = 0; i < t.size(); i++)
    {
        if (colourName != t[i].colour.name)
        {
            file << "usemtl " + t[i].colour.name + "\n";
            file << "s phong\n";
            colourName = t[i].colour.name;
        }
        line = "f " + std::to_string(t[i].vertexIdx[0] + 1 + offset) + "//" + std::to_string(t[i].vertexIdx[0] + 1 + offset) + " " + std::to_string(t[i].vertexIdx[1] + 1 + offset) + "//" + std::to_string(t[i].vertexIdx[1] + 1 + offset) + " " + std::to_string(t[i].vertexIdx[2] + 1 + offset) + "//" + std::to_string(t[i].vertexIdx[2] + 1 + offset) + "\n";
        file << line;
    }
    file.close();
}

int main(int argc, char *argv[])
{
    std::unordered_map<std::string, Colour> pallette = readMtlFile("../../models/cornell-box.mtl");
    vertexNormalPair v = readObjFile("../../models/bunny-low.obj", pallette);
    std::vector<ModelTriangle> modelPoints = v.modelPoints;
    std::vector<glm::vec3> vertexNormals = v.vertexNormals;

    writeToNewFile("../../models/test.obj", v, 98);
}