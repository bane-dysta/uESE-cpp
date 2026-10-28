#pragma once

#include "uese/io.hpp"
#include "uese/parameters.hpp"
#include "uese/surface.hpp"

#include <map>
#include <iosfwd>
#include <string>
#include <vector>

namespace uese {

// CODATA conversion for two unit charges a distance of 1 angstrom apart.
// This is a physical unit conversion, not a fitted model parameter.
constexpr double coulomb_kcal_angstrom = 332.063713299;

struct Options {
    int ndiv = 3;
    bool electrostatic_only = false;
    bool center_coordinates = false;
    std::size_t max_surface_points = 6000;
    std::map<int, double> radius_overrides;
};

struct EnergyTerms {
    double electrostatic = 0.0;
    double atomic_area = 0.0;
    double solvent_shift = 0.0;
    double linear_charge = 0.0;
    double nonlinear_charge = 0.0;
    double volume = 0.0;
    double short_range = 0.0;
    double correction() const;
    double total() const;
};

struct Result {
    Molecule molecule;
    Solvent solvent;
    Options options;
    Surface surface;
    std::vector<double> radii;
    std::vector<double> surface_charges; // conductor-limit q_i, not multiplied by f(epsilon)
    std::vector<double> atomic_surface_charges;
    std::vector<double> atomic_sigma;
    EnergyTerms energy;
    double dielectric_factor = 0.0;
    double conductor_energy = 0.0;
    double relative_residual = 0.0;
    double surface_charge_sum = 0.0;
    double sigma_min = 0.0, sigma_max = 0.0;
    bool cosmo_matrix_positive_definite = true;
    std::vector<std::string> warnings;
};

// Exact paper expressions (4)-(12), using a single fixed-charge linear solve.
Result calculate(Molecule molecule, Solvent solvent, const Options& options = {});
void write_text(const Result& result, std::ostream& stream);
void write_json(const Result& result, std::ostream& stream);
void write_surface_csv(const Result& result, const std::string& path);

} // namespace uese
