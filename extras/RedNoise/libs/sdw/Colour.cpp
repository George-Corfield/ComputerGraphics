#include "Colour.h"
#include <utility>
#include "TextureMap.h"

Colour::Colour() = default;
Colour::Colour(int r, int g, int b) : red(r), green(g), blue(b), texture(-1) {}
Colour::Colour(std::string n, int r, int g, int b) : name(std::move(n)),
													 red(r), green(g), blue(b),
													 texture(-1) {}

std::ostream &operator<<(std::ostream &os, const Colour &colour)
{
	os << colour.name << " ["
	   << colour.red << ", "
	   << colour.green << ", "
	   << colour.blue << "]";
	return os;
}
