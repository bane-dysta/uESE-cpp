#ifndef UESE_SURFACE_HPP
#define UESE_SURFACE_HPP

#include <array>
#include <cstddef>
#include <vector>

namespace uese {

// Lengths are in angstrom, areas in angstrom^2, volume in angstrom^3.
struct SurfaceAtom {
    std::array<double, 3> center;
    double radius;
};

struct SurfaceElement {
    std::size_t atom;  // Zero-based index of the parent atom.
    std::array<double, 3> center;
    std::array<double, 3> normal;
    double area;
};

struct Surface {
    std::vector<SurfaceElement> elements;
    std::vector<double> atomic_areas;
    double area = 0.0;
    double volume = 0.0;
};

// Independent implementation of the published GEPOL93-style vdW-surface
// quadrature: 60 parent triangles, 4^(ndiv-1) samples per parent, removal
// of buried samples, and regrouping into at most 60 elements per sphere.
// No solvent probe, radius scaling, or auxiliary spheres are added here.
// ndiv must be in [1,5]. Throws std::invalid_argument for invalid inputs.
// The finite quadrature for volume uses the given coordinate origin,
// matching the original GEPOL volume convention. See the project notes.
Surface make_surface(const std::vector<SurfaceAtom>& atoms, int ndiv = 3);

}  // namespace uese

#endif
