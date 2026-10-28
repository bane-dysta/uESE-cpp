#include "uese/io.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace uese {
namespace {
const std::vector<std::string> symbols = {
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
}

std::vector<std::string> tokens(std::string line) {
    const auto comment = line.find_first_of("#!");
    if (comment != std::string::npos) line.resize(comment);
    std::istringstream stream(line);
    std::vector<std::string> result;
    for (std::string word; stream >> word;) result.push_back(word);
    return result;
}

double number(std::string token, const std::string& context) {
    std::replace(token.begin(), token.end(), 'D', 'E');
    std::replace(token.begin(), token.end(), 'd', 'e');
    try {
        std::size_t used = 0;
        const double x = std::stod(token, &used);
        if (used != token.size() || !std::isfinite(x)) throw std::invalid_argument("number");
        return x;
    } catch (const std::exception&) {
        throw std::runtime_error(context + ": invalid finite number '" + token + "'");
    }
}

std::string compact(std::string value) {
    value = lower(trim(value));
    value.erase(std::remove_if(value.begin(), value.end(),
        [](unsigned char c) { return std::isspace(c) != 0; }), value.end());
    return value;
}

// Unit declarations must occupy a complete directive, header suffix, or
// comment. Ordinary prose mentioning atomic units must never rescale input.
bool unit_value(const std::string& text, InputUnits& unit) {
    const auto t = compact(text);
    if (t == "bohr" || t == "bohrs" || t == "au" || t == "a.u." ||
        t == "atomic_units" || t == "atomic-units" || t == "atomicunits") {
        unit = InputUnits::bohr;
        return true;
    }
    if (t == "angstrom" || t == "angstroms") {
        unit = InputUnits::angstrom;
        return true;
    }
    return false;
}

bool unit_directive(const std::string& text, InputUnits& unit) {
    const auto t = compact(text);
    if (t == "atomic_units" || t == "atomic-units" || t == "atomicunits") {
        unit = InputUnits::bohr;
        return true;
    }
    if (t.find("units=") == 0) return unit_value(t.substr(6), unit);
    if (t.find("units:") == 0) return unit_value(t.substr(6), unit);
    if (t.find("coordinatesin") == 0) return unit_value(t.substr(13), unit);
    return false;
}

std::string code_part(const std::string& line) {
    return trim(line.substr(0, line.find_first_of("#!")));
}

bool coordinate_header(const std::string& text, InputUnits& header_unit) {
    const std::string prefix = "atomic coordinates";
    const auto t = lower(trim(text));
    if (t.compare(0, prefix.size(), prefix) != 0) return false;
    auto suffix = trim(t.substr(prefix.size()));
    if (!suffix.empty() && suffix.front() == ':') suffix = trim(suffix.substr(1));
    header_unit = InputUnits::automatic;
    if (suffix.empty()) return true;
    if (suffix.size() >= 2 &&
        ((suffix.front() == '(' && suffix.back() == ')') ||
         (suffix.front() == '[' && suffix.back() == ']'))) {
        suffix = trim(suffix.substr(1, suffix.size() - 2));
    }
    if (suffix.find("in ") == 0) suffix = trim(suffix.substr(3));
    return unit_value(suffix, header_unit) || unit_directive(suffix, header_unit);
}

bool charge_header(const std::string& text) {
    const auto t = lower(trim(text));
    return t == "ground state charges" || t == "ground state charges:";
}

bool units_only_line(const std::string& line) {
    InputUnits ignored = InputUnits::automatic;
    return unit_directive(code_part(line), ignored);
}

Atom read_atom(const std::vector<std::string>& row, bool with_charge,
               std::size_t line_number) {
    const std::string context = "line " + std::to_string(line_number);
    if (row.size() != (with_charge ? 5U : 4U)) {
        throw std::runtime_error(context + (with_charge
            ? ": expected element x y z charge (five columns)"
            : ": expected atomic-number/symbol x y z (four columns)"));
    }
    Atom atom;
    atom.atomic_number = element_number(row[0]);
    for (std::size_t k = 0; k < 3; ++k) atom.position[k] = number(row[k + 1], context);
    if (with_charge) atom.charge = number(row[4], context);
    return atom;
}
} // namespace

double Molecule::total_charge() const {
    double result = 0.0;
    for (const auto& atom : atoms) result += atom.charge;
    return result;
}

std::string element_symbol(int atomic_number) {
    if (atomic_number < 1 || atomic_number >= static_cast<int>(symbols.size()))
        throw std::runtime_error("Atomic number must be between 1 and 118");
    return symbols[static_cast<std::size_t>(atomic_number)];
}

int element_number(const std::string& token) {
    if (!token.empty() && std::all_of(token.begin(), token.end(),
                    [](unsigned char c) { return std::isdigit(c) != 0; })) {
        try {
            const int z = std::stoi(token);
            if (z >= 1 && z <= 118) return z;
        } catch (const std::exception&) {}
    }
    const auto value = lower(token);
    for (std::size_t i = 1; i < symbols.size(); ++i)
        if (lower(symbols[i]) == value) return static_cast<int>(i);
    throw std::runtime_error("Unknown element/atomic number '" + token + "'");
}

Molecule read_molecule(const std::string& path, InputUnits units) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot open input file '" + path + "'");
    std::vector<std::string> lines;
    for (std::string line; std::getline(file, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (lines.empty() && line.compare(0, 3, "\xef\xbb\xbf") == 0) line.erase(0, 3);
        lines.push_back(line);
    }
    InputUnits declared_units = InputUnits::automatic;
    auto note_units = [&](InputUnits found) {
        if (found == InputUnits::automatic) return;
        if (units == InputUnits::automatic && declared_units != InputUnits::automatic &&
            found != declared_units) {
            throw std::runtime_error("Conflicting coordinate unit declarations");
        }
        declared_units = found;
    };
    for (const auto& line : lines) {
        InputUnits found = InputUnits::automatic;
        const auto code = code_part(line);
        if (unit_directive(code, found)) note_units(found);
        if (coordinate_header(code, found)) note_units(found);
        const auto comment = line.find_first_of("#!");
        if (comment != std::string::npos &&
            unit_directive(line.substr(comment + 1), found)) note_units(found);
    }
    Molecule result;
    std::size_t first = 0;
    while (first < lines.size() &&
        (tokens(lines[first]).empty() || units_only_line(lines[first]))) ++first;
    if (first == lines.size()) throw std::runtime_error("Input file is empty");
    const auto first_row = tokens(lines[first]);
    const bool xyz_header = first_row.size() == 1 &&
        std::all_of(first_row[0].begin(), first_row[0].end(),
            [](unsigned char c) { return std::isdigit(c) != 0; });
    std::size_t coordinate_line = lines.size(), charge_line = lines.size();
    // XYZ's second line is an arbitrary title, even when it happens to equal
    // an original-format block heading. Select XYZ before scanning headings.
    if (!xyz_header) {
        for (std::size_t i = 0; i < lines.size(); ++i) {
            const auto t = code_part(lines[i]);
            InputUnits ignored = InputUnits::automatic;
            if (coordinate_header(t, ignored)) {
                if (coordinate_line != lines.size())
                    throw std::runtime_error("Multiple coordinate blocks are not supported");
                coordinate_line = i;
            }
            if (charge_header(t)) {
                if (charge_line != lines.size())
                    throw std::runtime_error("Multiple charge blocks are not supported");
                charge_line = i;
            }
        }
    }
    if (coordinate_line != lines.size()) {
        if (charge_line == lines.size() || charge_line <= coordinate_line)
            throw std::runtime_error("Expected 'Ground state charges' after coordinates");
        if (coordinate_line != 0) result.title = lines[0];
        for (std::size_t i = coordinate_line + 1; i < charge_line; ++i) {
            const auto row = tokens(lines[i]);
            if (row.empty() || units_only_line(lines[i])) continue;
            result.atoms.push_back(read_atom(row, false, i + 1));
        }
        std::vector<double> charges;
        for (std::size_t i = charge_line + 1; i < lines.size(); ++i) {
            const auto row = tokens(lines[i]);
            if (units_only_line(lines[i])) continue;
            for (const auto& value : row)
                charges.push_back(number(value, "charge block line " + std::to_string(i + 1)));
        }
        if (charges.size() != result.atoms.size())
            throw std::runtime_error("Charge count (" + std::to_string(charges.size()) +
                ") differs from atom count (" + std::to_string(result.atoms.size()) + ")");
        for (std::size_t i = 0; i < charges.size(); ++i) result.atoms[i].charge = charges[i];
    } else {
        if (xyz_header) {
            std::size_t count = 0;
            try { count = static_cast<std::size_t>(std::stoull(first_row[0])); }
            catch (const std::exception&) { throw std::runtime_error("Invalid XYZ atom count"); }
            if (count > 100000 || first + 1 >= lines.size())
                throw std::runtime_error("Invalid XYZ header or unreasonable atom count");
            result.title = lines[first + 1];
            for (std::size_t i = first + 2; i < lines.size(); ++i) {
                const auto row = tokens(lines[i]);
                if (!row.empty() && !units_only_line(lines[i]))
                    result.atoms.push_back(read_atom(row, true, i + 1));
            }
            if (result.atoms.size() != count)
                throw std::runtime_error("XYZ header atom count does not match coordinate records");
        } else {
            for (std::size_t i = first; i < lines.size(); ++i) {
                const auto row = tokens(lines[i]);
                if (!row.empty() && !units_only_line(lines[i]))
                    result.atoms.push_back(read_atom(row, true, i + 1));
            }
        }
    }
    if (result.atoms.empty()) throw std::runtime_error("Input contains no atoms");
    if (result.atoms.size() > 100000) throw std::runtime_error("Too many atoms in input");
    const double scale = (units == InputUnits::bohr ||
        (units == InputUnits::automatic && declared_units == InputUnits::bohr))
        ? bohr_to_angstrom : 1.0;
    for (auto& atom : result.atoms)
        for (double& coordinate : atom.position) coordinate *= scale;
    return result;
}

} // namespace uese
