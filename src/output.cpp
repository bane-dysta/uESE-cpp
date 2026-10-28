#include "uese/model.hpp"

#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>

namespace uese {
namespace {
void quoted(std::ostream& out, const std::string& text) {
    out << '"';
    for (unsigned char c : text) {
        switch (c) {
            case '"': out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (c < 0x20) {
                    const char* hex = "0123456789abcdef";
                    out << "\\u00" << hex[c >> 4] << hex[c & 15];
                } else out << static_cast<char>(c);
        }
    }
    out << '"';
}
} // namespace

void write_text(const Result& r, std::ostream& out) {
    out << "uese-cpp 1.0 -- independent implementation of the 2021 uESE equations\n"
        << "Reference: Vyboishchikov & Voityuk, DOI 10.1002/jcc.26531\n\n"
        << "Calculation: " << (r.options.electrostatic_only ? "COSMO electrostatic contribution only"
                                                          : "uESE solvation free energy") << '\n'
        << "Solvent: " << r.solvent.name << "   Class: " << r.solvent.solvent_class
        << " (" << class_description(r.solvent.solvent_class) << ")\n"
        << std::fixed << std::setprecision(8)
        << "Dielectric constant: " << r.solvent.epsilon << '\n'
        << "Dielectric scaling factor: " << r.dielectric_factor << '\n';
    if (r.solvent.solvent_class == 'C' || r.solvent.solvent_class == 'D') {
        out << "Solvent shift xi [kcal/mol/angstrom^2]: ";
        if (r.solvent.xi) out << *r.solvent.xi; else out << "not used / unavailable";
        out << "\nSolvent shift source: " << r.solvent.xi_source << '\n';
    }
    out << "Atoms: " << r.molecule.atoms.size() << '\n'
        << "Input CM5 charge sum [e]: " << r.molecule.total_charge() << '\n'
        << "NDIV: " << r.options.ndiv << '\n'
        << "Surface points: " << r.surface.elements.size() << '\n'
        << "Cavity area [angstrom^2]: " << r.surface.area << '\n'
        << "Cavity volume [angstrom^3]: " << r.surface.volume << '\n'
        << "Conductor surface charge sum [e]: " << r.surface_charge_sum << '\n'
        << std::scientific << "Linear-equation relative residual: " << r.relative_residual
        << std::fixed << "\n\nEnergy components [kcal/mol]:\n"
        << "  COSMO electrostatic: " << std::setw(17) << r.energy.electrostatic << '\n';
    if (!r.options.electrostatic_only) {
        out << "  Atomic area term:   " << std::setw(17) << r.energy.atomic_area << '\n'
            << "  Solvent shift term: " << std::setw(17) << r.energy.solvent_shift << '\n'
            << "  Linear charge term:" << std::setw(17) << r.energy.linear_charge << '\n'
            << "  Charge q^2/|q| term:" << std::setw(16) << r.energy.nonlinear_charge << '\n'
            << "  Volume term:        " << std::setw(17) << r.energy.volume << '\n'
            << "  Short-range term:   " << std::setw(17) << r.energy.short_range << '\n'
            << "  Correction total:   " << std::setw(17) << r.energy.correction() << '\n'
            << "  TOTAL SOLVATION:    " << std::setw(17) << r.energy.total() << '\n';
    }
    out << "\nAtomic surface data (conductor charges, before dielectric scaling):\n"
        << " atom element       Q_input          radius            area            q_surface           sigma\n";
    for (std::size_t a = 0; a < r.molecule.atoms.size(); ++a) {
        out << std::setw(5) << a + 1 << std::setw(7) << element_symbol(r.molecule.atoms[a].atomic_number)
            << std::setw(16) << r.molecule.atoms[a].charge << std::setw(16) << r.radii[a]
            << std::setw(16) << r.surface.atomic_areas[a]
            << std::setw(21) << r.atomic_surface_charges[a];
        if (r.surface.atomic_areas[a] > 0.0) out << std::setw(16) << r.atomic_sigma[a];
        else out << std::setw(16) << "buried";
        out << '\n';
    }
    for (const auto& warning : r.warnings) out << "\nWARNING: " << warning << '\n';
}

void write_json(const Result& r, std::ostream& out) {
    out << std::setprecision(15) << "{\n  \"implementation\": \"uese-cpp 1.0\",\n"
        << "  \"model\": \"" << (r.options.electrostatic_only ? "COSMO-only" : "uESE-2021") << "\",\n"
        << "  \"doi\": \"10.1002/jcc.26531\",\n  \"title\": "; quoted(out, r.molecule.title);
    out << ",\n  \"units\": {\"length\":\"angstrom\",\"area\":\"angstrom^2\","
        << "\"volume\":\"angstrom^3\",\"charge\":\"e\",\"energy\":\"kcal/mol\"},\n"
        << "  \"solvent\": {\"name\":"; quoted(out, r.solvent.name);
    out << ",\"class\":\"" << r.solvent.solvent_class << "\",\"epsilon\":" << r.solvent.epsilon
        << ",\"xi\":";
    if (r.solvent.xi) out << *r.solvent.xi; else out << "null";
    out << ",\"xi_source\":"; quoted(out, r.solvent.xi_source); out << "},\n";
    out << "  \"settings\": {\"ndiv\":" << r.options.ndiv << ",\"centered\":"
        << (r.options.center_coordinates ? "true" : "false") << ",\"radius_overrides\":{ ";
    bool first = true;
    for (const auto& [z, radius] : r.options.radius_overrides) {
        if (!first) out << ',';
        quoted(out, element_symbol(z)); out << ':' << radius; first = false;
    }
    out << "}},\n  \"n_atoms\":" << r.molecule.atoms.size()
        << ",\n  \"n_surface_points\":" << r.surface.elements.size()
        << ",\n  \"solute_charge\":" << r.molecule.total_charge()
        << ",\n  \"surface_charge\":" << r.surface_charge_sum
        << ",\n  \"surface_area\":" << r.surface.area
        << ",\n  \"cavity_volume\":" << r.surface.volume
        << ",\n  \"dielectric_factor\":" << r.dielectric_factor
        << ",\n  \"energy\":{\"electrostatic\":" << r.energy.electrostatic;
    if (r.options.electrostatic_only) {
        out << ",\"correction\":null,\"total\":null";
    } else {
        out << ",\"atomic_area\":" << r.energy.atomic_area
            << ",\"solvent_shift\":" << r.energy.solvent_shift
            << ",\"linear_charge\":" << r.energy.linear_charge
            << ",\"nonlinear_charge\":" << r.energy.nonlinear_charge
            << ",\"volume\":" << r.energy.volume
            << ",\"short_range\":" << r.energy.short_range
            << ",\"correction\":" << r.energy.correction()
            << ",\"total\":" << r.energy.total();
    }
    out << "},\n  \"numerics\":{\"relative_residual\":" << r.relative_residual
        << ",\"positive_definite\":" << (r.cosmo_matrix_positive_definite ? "true" : "false")
        << ",\"conductor_energy\":" << r.conductor_energy
        << ",\"sigma_min\":" << r.sigma_min << ",\"sigma_max\":" << r.sigma_max
        << ",\"coulomb_kcal_angstrom\":" << coulomb_kcal_angstrom << "},\n  \"atoms\":[\n";
    for (std::size_t a = 0; a < r.molecule.atoms.size(); ++a) {
        const auto& atom = r.molecule.atoms[a];
        out << "    {\"index\":" << a + 1 << ",\"element\":";
        quoted(out, element_symbol(atom.atomic_number));
        out << ",\"xyz\":[" << atom.position[0] << ',' << atom.position[1] << ',' << atom.position[2]
            << "],\"input_charge\":" << atom.charge << ",\"radius\":" << r.radii[a]
            << ",\"surface_area\":" << r.surface.atomic_areas[a]
            << ",\"surface_charge\":" << r.atomic_surface_charges[a] << ",\"sigma\":";
        if (r.surface.atomic_areas[a] > 0.0) out << r.atomic_sigma[a]; else out << "null";
        out << '}' << (a + 1 == r.molecule.atoms.size() ? "\n" : ",\n");
    }
    out << "  ],\n  \"warnings\":[";
    for (std::size_t i = 0; i < r.warnings.size(); ++i) {
        if (i != 0) out << ',';
        quoted(out, r.warnings[i]);
    }
    out << "]\n}\n";
}

void write_surface_csv(const Result& r, const std::string& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write surface CSV '" + path + "'");
    out << "surface_index,atom_index,element,x_angstrom,y_angstrom,z_angstrom,"
        << "nx,ny,nz,area_angstrom2,q_conductor_e,sigma_conductor_e_per_angstrom2\n"
        << std::setprecision(15);
    for (std::size_t i = 0; i < r.surface.elements.size(); ++i) {
        const auto& e = r.surface.elements[i];
        out << i + 1 << ',' << e.atom + 1 << ','
            << element_symbol(r.molecule.atoms[e.atom].atomic_number);
        for (double x : e.center) out << ',' << x;
        for (double x : e.normal) out << ',' << x;
        out << ',' << e.area << ',' << r.surface_charges[i] << ','
            << r.surface_charges[i] / e.area << '\n';
    }
    if (!out) throw std::runtime_error("Failed while writing surface CSV '" + path + "'");
}

} // namespace uese
