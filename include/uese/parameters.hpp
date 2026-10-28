#pragma once

#include <optional>
#include <string>
#include <vector>

namespace uese {

struct ElementParameters {
    int atomic_number;
    double kappa; // kcal mol^-1 angstrom^-2
    double p;     // kcal mol^-1, with charges in units of e
    double g;     // kcal mol^-1, with charges in units of e
};

struct ClassParameters {
    char code;
    double zeta;  // kcal mol^-1 angstrom^-3
    double src_a, src_b, src_c;
    std::vector<ElementParameters> elements;
};

struct Solvent {
    std::string name;
    std::string aliases; // '|' separated
    char solvent_class;
    double epsilon;
    std::optional<double> xi; // kcal mol^-1 angstrom^-2; only C and D use this
    std::string xi_source;
};

const ClassParameters& class_parameters(char code);
const ElementParameters& element_parameters(char code, int atomic_number);
double default_radius(int atomic_number);
std::string normalize_solvent_name(std::string name);
const std::vector<Solvent>& solvents();
Solvent find_solvent(const std::string& name);
std::string class_description(char code);

} // namespace uese
