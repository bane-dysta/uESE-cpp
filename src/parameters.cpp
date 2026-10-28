#include "uese/parameters.hpp"
#include "uese/io.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace uese {

// Vyboishchikov & Voityuk, J. Comput. Chem. 42 (2021), 1184-1194,
// DOI 10.1002/jcc.26531, Table 1, printed p.1187.
// A blank table entry is deliberately absent; it is NOT interpreted as zero.
const ClassParameters& class_parameters(char code) {
    static const ClassParameters a{'A', 0.0, 68500.0, 87700.0, -8.87, {
        {1, .024, -16.5, 0}, {6, -.003, 35.4, 0}, {7, 0, 65.8, 0},
        {8, .012, 14.2, 0}, {9, .054, 14.5, 0}, {15, .030, 0, 0},
        {16, -.022, 29.2, 0}, {17, -.006, 28.3, 0},
        {35, -.016, 29.8, 0}, {53, -.025, 23.0, 0}}};
    static const ClassParameters b{'B', -.070, 186000.0, 78500.0, -2.36, {
        {1, .045, -14.6, 0}, {6, .030, 5.5, 0}, {7, 0, 73.0, 0},
        {8, .072, 17.4, 0}, {9, .077, 3.0, 0}, {15, .110, 10.0, 0},
        {16, 0, 24.0, 0}, {17, .019, 24.0, 0},
        {35, 0, 25.6, 0}, {53, -.008, 26.0, 0}}};
    static const ClassParameters c{'C', -.236, 0, 0, 0, {
        {1, .151, -17.21, 7.19}, {6, .170, 12.89, 9.27},
        {7, .239, 77.5, 0}, {8, .052, 14.56, 21.82},
        {9, .200, 0, 12.0}, {16, .122, 7.6, 21.35},
        {17, .086, 0, 32.1}, {35, .083, 0, 32.7}}};
    static const ClassParameters d{'D', -.166, 0, 0, 0, {
        {1, .106, -6.15, 0}, {6, .115, 1.54, 0}, {7, -.092, 23.1, 0},
        {8, .056, 9.52, 0}, {9, .100, 6.17, 0}, {15, .600, -70.0, 0},
        {16, .050, 12.5, 0}, {17, .087, 1.94, 0},
        {35, .072, 4.95, 0}, {53, .077, -4.33, 0}}};
    switch (code) {
        case 'A': return a;
        case 'B': return b;
        case 'C': return c;
        case 'D': return d;
        default: throw std::runtime_error("Solvent class must be A, B, C, or D");
    }
}

const ElementParameters& element_parameters(char code, int atomic_number) {
    const auto& params = class_parameters(code);
    for (const auto& e : params.elements)
        if (e.atomic_number == atomic_number) return e;
    throw std::runtime_error("No published Table 1 correction parameters for " +
        element_symbol(atomic_number) + " in solvent class " + code +
        ". Missing parameters are not replaced by zero. "
        "Use --electrostatic-only only if the COSMO contribution alone is wanted.");
}

double default_radius(int atomic_number) {
    // Unscaled Bondi radii used for the published parameterized elements.
    switch (atomic_number) {
        case 1: return 1.20;
        case 6: return 1.70;
        case 7: return 1.55;
        case 8: return 1.52;
        case 9: return 1.47;
        case 15: return 1.80;
        case 16: return 1.80;
        case 17: return 1.75;
        case 35: return 1.85;
        case 53: return 1.98;
        // These two radii were verified in the author's v1.2 program output.
        // Their presence only enables electrostatic-only runs, not uESE totals.
        case 5: return 1.92;
        case 79: return 1.66;
        default: throw std::runtime_error("No verified default radius for " +
            element_symbol(atomic_number) + ". Supply --radius Element=value explicitly.");
    }
}

std::string normalize_solvent_name(std::string name) {
    std::string result;
    for (unsigned char c : name) {
        if (std::isalnum(c)) result += static_cast<char>(std::tolower(c));
    }
    return result;
}

const std::vector<Solvent>& solvents() {
    static const std::vector<Solvent> table = {
#include "solvents_data.inc"
    };
    return table;
}

Solvent find_solvent(const std::string& name) {
    const auto wanted = normalize_solvent_name(name);
    for (const auto& solvent : solvents()) {
        if (normalize_solvent_name(solvent.name) == wanted) return solvent;
        std::istringstream stream(solvent.aliases);
        for (std::string alias; std::getline(stream, alias, '|');)
            if (!alias.empty() && normalize_solvent_name(alias) == wanted) return solvent;
    }
    throw std::runtime_error("Unknown solvent '" + name + "'. Use --list-solvents, or "
        "--class A|B|C|D --epsilon VALUE --xi VALUE for an explicit custom solvent.");
}

std::string class_description(char code) {
    switch (code) {
        case 'A': return "water";
        case 'B': return "polar protic";
        case 'C': return "polar aprotic";
        case 'D': return "nonpolar";
        default: throw std::runtime_error("Invalid solvent class");
    }
}

} // namespace uese
