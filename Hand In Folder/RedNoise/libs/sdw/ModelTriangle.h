#pragma once

#include <glm/glm.hpp>
#include <string>
#include <array>
#include "Colour.h"

struct ModelTriangle
{
	std::array<glm::vec3, 3> vertices{};
	std::array<glm::vec2, 3> texturePoints{};
	Colour colour{};
	glm::vec3 normal{};
	std::array<int, 3> vertexIdx{};
	std::string shadingType;

	ModelTriangle();
	ModelTriangle(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, Colour trigColour);
	friend std::ostream &operator<<(std::ostream &os, const ModelTriangle &triangle);
};
