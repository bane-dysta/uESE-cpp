#include "uese/model.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

namespace {
void help() {
    std::cout << R"(uese-cpp 1.0: independent, dependency-free C++17 uESE implementation

Usage:
  uese INPUT -solvent water [OPTIONS]
  uese INPUT --class C --epsilon 35.6881 --xi VALUE [OPTIONS]
  uese --list-solvents

Input:
  Original uESE/Multiwfn .inp (coordinates plus gas-phase CM5 charges), or
  extended XYZ: atom count, comment, then symbol x y z charge on each line.
  Coordinates default to angstrom; charge values are in units of e.

Options:
  -solvent, --solvent NAME  Named solvent (default: water)
  --class A|B|C|D           Explicit custom solvent class
  --epsilon VALUE          Explicit dielectric constant (>= 1)
  --xi VALUE               Explicit solvent shift, kcal/mol/angstrom^2
                           Required for custom C/D; zero is an explicit variant.
  --units angstrom|bohr     Override the coordinate units in the input
  --ndiv N                 GEPOL-style subsampling depth 1..5 (reference: 3)
  --format text|json       Output format (default: text)
  --surface FILE.csv       Export exposed tesserae and conductor charges
  --center                 Explicitly center coordinates before cavity integration
  --electrostatic-only     Calculate only COSMO; no solvation total is reported
  --radius ELEMENT=VALUE   Explicit radius override in angstrom; repeatable
  --max-points N           Dense-solver size guard (default: 6000)
  --list-solvents           List solvent names, classes, epsilon, and xi provenance
  --help, -h               Show this help
  --version                Show version

Model reference: DOI 10.1002/jcc.26531 (2021), equations (4)-(12), Table 1.
This program consumes CM5 charges; it does not perform a quantum-chemical
calculation or compute CM5 charges from an XYZ geometry alone.
Missing published elemental parameters are errors, not zero defaults.
)";
}

double real_number(const std::string& text, const std::string& flag) {
    try {
        std::size_t used = 0;
        double result = std::stod(text, &used);
        if (used != text.size() || !std::isfinite(result)) throw std::invalid_argument("number");
        return result;
    } catch (const std::exception&) {
        throw std::runtime_error(flag + " requires a finite number, received '" + text + "'");
    }
}

std::size_t positive_integer(const std::string& text, const std::string& flag) {
    if (text.empty() || !std::all_of(text.begin(), text.end(),
        [](unsigned char c) { return std::isdigit(c) != 0; }))
        throw std::runtime_error(flag + " requires a positive integer");
    try {
        const auto n = std::stoull(text);
        if (n == 0 || n > 100000000U) throw std::out_of_range("integer");
        return static_cast<std::size_t>(n);
    } catch (const std::exception&) { throw std::runtime_error(flag + " is out of range"); }
}

void list_solvents() {
    std::cout << "name\tclass\tepsilon\txi_kcal_mol_A2\txi_source\taliases\n"
              << std::setprecision(12);
    for (const auto& s : uese::solvents()) {
        std::cout << s.name << '\t' << s.solvent_class << '\t' << s.epsilon << '\t';
        if (s.xi) std::cout << *s.xi; else std::cout << "UNAVAILABLE";
        std::cout << '\t' << s.xi_source << '\t' << s.aliases << '\n';
    }
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 1) { help(); return 0; }
        std::string input, solvent_name = "water", format = "text", surface_path;
        std::optional<char> solvent_class;
        std::optional<double> epsilon, xi;
        bool named_solvent = false, do_list = false;
        uese::Options options;
        uese::InputUnits units = uese::InputUnits::automatic;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            auto value = [&]() -> std::string {
                if (i + 1 >= argc) throw std::runtime_error("Missing value after " + arg);
                return argv[++i];
            };
            if (arg == "--help" || arg == "-h") { help(); return 0; }
            if (arg == "--version") { std::cout << "uese-cpp 1.1 (2021 equations, SI Table S1 verified)\n"; return 0; }
            if (arg == "--list-solvents") do_list = true;
            else if (arg == "-solvent" || arg == "--solvent") { solvent_name = value(); named_solvent = true; }
            else if (arg == "--class") {
                auto v = value();
                if (v.size() != 1) throw std::runtime_error("--class must be A, B, C, or D");
                const char c = static_cast<char>(std::toupper(static_cast<unsigned char>(v[0])));
                uese::class_parameters(c);
                solvent_class = c;
            } else if (arg == "--epsilon") epsilon = real_number(value(), arg);
            else if (arg == "--xi") xi = real_number(value(), arg);
            else if (arg == "--units") {
                const auto v = value();
                if (v == "angstrom") units = uese::InputUnits::angstrom;
                else if (v == "bohr") units = uese::InputUnits::bohr;
                else throw std::runtime_error("--units must be angstrom or bohr");
            } else if (arg == "--ndiv") {
                const auto n = positive_integer(value(), arg);
                if (n > 5) throw std::runtime_error("--ndiv must be in [1,5]");
                options.ndiv = static_cast<int>(n);
            } else if (arg == "--format") {
                format = value();
                if (format != "text" && format != "json")
                    throw std::runtime_error("--format must be text or json");
            } else if (arg == "--surface") surface_path = value();
            else if (arg == "--center") options.center_coordinates = true;
            else if (arg == "--electrostatic-only") options.electrostatic_only = true;
            else if (arg == "--max-points") options.max_surface_points = positive_integer(value(), arg);
            else if (arg == "--radius") {
                const auto v = value();
                const auto split = v.find('=');
                if (split == std::string::npos) throw std::runtime_error("--radius requires Element=value");
                const auto z = uese::element_number(v.substr(0, split));
                const double r = real_number(v.substr(split + 1), arg);
                if (r <= 0.0) throw std::runtime_error("A radius must be positive");
                options.radius_overrides[z] = r;
            } else if (!arg.empty() && arg[0] == '-') {
                throw std::runtime_error("Unknown option '" + arg + "'");
            } else {
                if (!input.empty()) throw std::runtime_error("Only one input file can be given");
                input = arg;
            }
        }
        if (do_list) { list_solvents(); return 0; }
        if (input.empty()) throw std::runtime_error("An input file is required; see --help");
        uese::Solvent solvent;
        if (solvent_class && !named_solvent) {
            if (!epsilon) throw std::runtime_error("A custom solvent requires --epsilon");
            solvent = {"custom", "", *solvent_class, *epsilon, xi, "user supplied"};
            if ((*solvent_class == 'A' || *solvent_class == 'B') && !xi) solvent.xi = 0.0;
        } else {
            solvent = uese::find_solvent(solvent_name);
            if (solvent_class) solvent.solvent_class = *solvent_class;
            if (epsilon) solvent.epsilon = *epsilon;
            if (xi) { solvent.xi = *xi; solvent.xi_source = "user override"; }
            if (solvent_class && !xi && (solvent.solvent_class == 'C' || solvent.solvent_class == 'D'))
                throw std::runtime_error("Changing a named solvent's class requires an explicit --xi");
        }
        auto result = uese::calculate(uese::read_molecule(input, units), solvent, options);
        if (!surface_path.empty()) {
            if (std::filesystem::weakly_canonical(surface_path) == std::filesystem::weakly_canonical(input) ||
                (std::filesystem::exists(surface_path) && std::filesystem::equivalent(surface_path, input)))
                throw std::runtime_error("Surface output must not replace the input file");
            uese::write_surface_csv(result, surface_path);
        }
        if (format == "json") uese::write_json(result, std::cout);
        else uese::write_text(result, std::cout);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "uese: error: " << e.what() << '\n';
        return 1;
    }
}
