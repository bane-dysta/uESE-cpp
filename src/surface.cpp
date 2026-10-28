#include "uese/surface.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace uese {
namespace {

using Vec = std::array<double, 3>;
using Face = std::array<std::size_t, 3>;
constexpr double pi = 3.141592653589793238462643383279502884;

Vec add(const Vec& a, const Vec& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Vec sub(const Vec& a, const Vec& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vec scale(const Vec& a, double s) {
    return {a[0] * s, a[1] * s, a[2] * s};
}

double dot(const Vec& a, const Vec& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double length(const Vec& a) {
    return std::hypot(a[0], a[1], a[2]);
}

Vec cross(const Vec& a, const Vec& b) {
    return {a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

bool finite(const Vec& a) {
    return std::isfinite(a[0]) && std::isfinite(a[1]) && std::isfinite(a[2]);
}

Vec unit(const Vec& v) {
    const double n = length(v);
    if (!(n > 0.0) || !std::isfinite(n)) {
        throw std::runtime_error("Surface geometry has a zero or non-finite vector.");
    }
    return scale(v, 1.0 / n);
}

// For the two small, fixed convex vertex sets used here, testing all triples
// is simpler than introducing a general convex-hull dependency. A face's
// plane leaves every remaining vertex on one side. No topology table or
// code from the original GEPOL implementation is included.
std::vector<Face> hull_faces(const std::vector<Vec>& vertices) {
    std::vector<Face> faces;
    for (std::size_t a = 0; a < vertices.size(); ++a) {
        for (std::size_t b = a + 1; b < vertices.size(); ++b) {
            for (std::size_t c = b + 1; c < vertices.size(); ++c) {
                const Vec n = cross(sub(vertices[b], vertices[a]),
                                    sub(vertices[c], vertices[a]));
                const double norm = length(n);
                if (norm < 1.0e-12) continue;
                bool positive = false;
                bool negative = false;
                const double tolerance = 1.0e-10 * norm;
                for (const auto& v : vertices) {
                    const double side = dot(n, sub(v, vertices[a]));
                    positive = positive || side > tolerance;
                    negative = negative || side < -tolerance;
                    if (positive && negative) break;
                }
                if (!(positive && negative)) faces.push_back({a, b, c});
            }
        }
    }
    return faces;
}

struct Polyhedron {
    std::vector<Vec> vertices;
    std::vector<Face> faces;
};

Polyhedron parent_polyhedron() {
    // An icosahedron is defined by its poles and two staggered pentagons.
    // Its 20 normalized face centers are the dual dodecahedron vertices.
    // Their union has 32 unit vertices and a 60-triangle convex hull: the
    // spherical pentakis-dodecahedron used by GEPOL. The upper pentagon
    // begins on the +x meridian; the lower one is rotated by pi/5.
    Polyhedron p;
    p.vertices.push_back({0.0, 0.0, 1.0});
    const double h = 1.0 / std::sqrt(5.0);
    const double r = std::sqrt(1.0 - h * h);
    for (int ring = 0; ring < 2; ++ring) {
        const double z = ring == 0 ? h : -h;
        const double phase = ring == 0 ? 0.0 : pi / 5.0;
        for (int j = 0; j < 5; ++j) {
            const double phi = phase + 2.0 * pi * static_cast<double>(j) / 5.0;
            p.vertices.push_back({r * std::cos(phi), r * std::sin(phi), z});
        }
    }
    p.vertices.push_back({0.0, 0.0, -1.0});
    const auto icosahedral_faces = hull_faces(p.vertices);
    if (icosahedral_faces.size() != 20) {
        throw std::runtime_error("Failed to construct the icosahedron surface grid.");
    }
    for (const auto& f : icosahedral_faces) {
        p.vertices.push_back(unit(add(add(p.vertices[f[0]], p.vertices[f[1]]),
                                      p.vertices[f[2]])));
    }
    p.faces = hull_faces(p.vertices);
    if (p.vertices.size() != 32 || p.faces.size() != 60) {
        throw std::runtime_error("Failed to construct the 60-triangle surface grid.");
    }
    return p;
}

void triangle_samples(const Vec& a, const Vec& b, const Vec& c,
                      int depth, std::vector<Vec>& samples) {
    if (depth == 0) {
        samples.push_back(unit(add(add(a, b), c)));
        return;
    }
    const Vec ab = unit(add(a, b));
    const Vec bc = unit(add(b, c));
    const Vec ca = unit(add(c, a));
    triangle_samples(a, ab, ca, depth - 1, samples);
    triangle_samples(ab, b, bc, depth - 1, samples);
    triangle_samples(ca, bc, c, depth - 1, samples);
    triangle_samples(ab, bc, ca, depth - 1, samples);
}

std::vector<std::vector<Vec>> grid(int ndiv) {
    const auto p = parent_polyhedron();
    std::vector<std::vector<Vec>> samples;
    samples.reserve(p.faces.size());
    std::size_t per_parent = 1;
    for (int d = 1; d < ndiv; ++d) per_parent *= 4;
    for (const auto& f : p.faces) {
        samples.emplace_back();
        samples.back().reserve(per_parent);
        triangle_samples(p.vertices[f[0]], p.vertices[f[1]], p.vertices[f[2]],
                         ndiv - 1, samples.back());
    }
    return samples;
}

void validate(const std::vector<SurfaceAtom>& atoms, int ndiv) {
    if (ndiv < 1 || ndiv > 5) {
        throw std::invalid_argument("Surface ndiv must be an integer from 1 to 5.");
    }
    if (atoms.empty()) throw std::invalid_argument("Surface requires at least one atom.");
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const auto& a = atoms[i];
        if (!finite(a.center) || !std::isfinite(a.radius) || !(a.radius > 0.0)) {
            throw std::invalid_argument("Surface atom " + std::to_string(i + 1) +
                                        " has invalid coordinates or radius.");
        }
        if (!std::isfinite(4.0 * pi * a.radius * a.radius)) {
            throw std::invalid_argument("Surface radius is too large for finite area.");
        }
        for (std::size_t j = 0; j < i; ++j) {
            const double d = length(sub(a.center, atoms[j].center));
            if (!std::isfinite(d)) {
                throw std::invalid_argument("Surface coordinate range is too large.");
            }
            if (d <= 1.0e-12) {
                throw std::invalid_argument("Surface atoms " + std::to_string(j + 1) +
                                            " and " + std::to_string(i + 1) +
                                            " have coincident centers.");
            }
        }
    }
}

}  // namespace

Surface make_surface(const std::vector<SurfaceAtom>& atoms, int ndiv) {
    validate(atoms, ndiv);
    const auto samples = grid(ndiv);
    const double sample_count = 60.0 * static_cast<double>(samples.front().size());
    Surface result;
    result.atomic_areas.resize(atoms.size(), 0.0);
    if (atoms.size() <= std::numeric_limits<std::size_t>::max() / 60) {
        result.elements.reserve(atoms.size() * 60);
    }
    for (std::size_t a = 0; a < atoms.size(); ++a) {
        const auto& atom = atoms[a];
        std::vector<std::size_t> neighbours;
        for (std::size_t b = 0; b < atoms.size(); ++b) {
            if (a != b && length(sub(atom.center, atoms[b].center)) <
                              atom.radius + atoms[b].radius) {
                neighbours.push_back(b);
            }
        }
        const double sample_area = 4.0 * pi * atom.radius * atom.radius / sample_count;
        for (const auto& parent_samples : samples) {
            Vec exposed_sum{0.0, 0.0, 0.0};
            std::size_t exposed_count = 0;
            for (const auto& normal : parent_samples) {
                const Vec point = add(atom.center, scale(normal, atom.radius));
                if (!finite(point)) {
                    throw std::invalid_argument("Surface point exceeds finite coordinate range.");
                }
                bool buried = false;
                for (const auto b : neighbours) {
                    if (length(sub(point, atoms[b].center)) < atoms[b].radius) {
                        buried = true;
                        break;
                    }
                }
                if (!buried) {
                    exposed_sum = add(exposed_sum, normal);
                    ++exposed_count;
                }
            }
            if (exposed_count == 0) continue;
            // Each parent is represented by one point, on the parent sphere.
            // The same equal-area child quadrature is used for covered and
            // uncovered triangles; exact spherical triangle areas would
            // define a different discretization from the 2021 reference.
            const Vec normal = unit(exposed_sum);
            const Vec center = add(atom.center, scale(normal, atom.radius));
            const double area = sample_area * static_cast<double>(exposed_count);
            result.elements.push_back({a, center, normal, area});
            result.atomic_areas[a] += area;
            result.area += area;
            result.volume += dot(center, normal) * area / 3.0;
        }
    }
    if (result.elements.empty()) {
        throw std::runtime_error("No exposed surface elements were generated.");
    }
    if (!std::isfinite(result.area) || !std::isfinite(result.volume)) {
        throw std::runtime_error("Non-finite surface area or volume.");
    }
    return result;
}

}  // namespace uese
