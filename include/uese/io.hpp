#pragma once

#include <array>
#include <string>
#include <vector>

namespace uese {

// All coordinates in the public API are in angstrom; charges are in units of e.
struct Atom {
    int atomic_number = 0;
    std::array<double, 3> position{};
    double charge = 0.0;
};

struct Molecule {
    std::string title;
    std::vector<Atom> atoms;
    double total_charge() const;
};

enum class InputUnits { automatic, angstrom, bohr };

constexpr double bohr_to_angstrom = 0.529177210903;

std::string element_symbol(int atomic_number);
int element_number(const std::string& token);
// Accepts the original "Atomic coordinates" / "Ground state charges" format,
// extended XYZ (symbol x y z charge), and a headerless five-column charge file.
// Automatic units default to angstrom. Explicit declarations include a whole
// "units=bohr" line/comment, a unit-bearing coordinate header, or a whole
// "atomic_units" line. An explicit InputUnits argument overrides declarations.
Molecule read_molecule(const std::string& path,
                       InputUnits units = InputUnits::automatic);

} // namespace uese
