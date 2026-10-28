#include "uese/model.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace uese {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;

double distance(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return std::hypot(std::hypot(a[0] - b[0], a[1] - b[1]), a[2] - b[2]);
}

double kernel(const SurfaceElement& a, const SurfaceElement& b, bool diagonal) {
    if (diagonal) return 1.0694 * std::sqrt(4.0 * pi / a.area);
    const double d = distance(a.center, b.center);
    if (d < 1e-10)
        throw std::runtime_error("Coincident surface centers: the COSMO matrix is singular");
    return 1.0 / d;
}

std::vector<double> make_matrix(const Surface& surface) {
    const std::size_t n = surface.elements.size();
    std::vector<double> a(n * n);
    for (std::size_t i = 0; i < n; ++i) {
        a[i * n + i] = kernel(surface.elements[i], surface.elements[i], true);
        for (std::size_t j = 0; j < i; ++j) {
            const double v = kernel(surface.elements[i], surface.elements[j], false);
            a[i * n + j] = v;
            a[j * n + i] = v;
        }
    }
    return a;
}

bool cholesky_solve(std::vector<double>& a, const std::vector<double>& rhs,
                    std::vector<double>& x) {
    const std::size_t n = rhs.size();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double value = a[i * n + j];
            for (std::size_t k = 0; k < j; ++k)
                value -= a[i * n + k] * a[j * n + k];
            if (i == j) {
                if (!(value > 1e-13 * a[i * n + i]) || !std::isfinite(value)) return false;
                a[i * n + i] = std::sqrt(value);
            } else {
                a[i * n + j] = value / a[j * n + j];
            }
        }
    }
    x = rhs;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < i; ++j) x[i] -= a[i * n + j] * x[j];
        x[i] /= a[i * n + i];
    }
    for (std::size_t ii = n; ii-- > 0;) {
        for (std::size_t j = ii + 1; j < n; ++j) x[ii] -= a[j * n + ii] * x[j];
        x[ii] /= a[ii * n + ii];
    }
    return true;
}

// A coarse collocation matrix can be indefinite. LU is a numerical fallback,
// not a change to the published electrostatic kernel or a hidden regularizer.
std::vector<double> lu_solve(std::vector<double>& a, std::vector<double> rhs) {
    const std::size_t n = rhs.size();
    std::vector<double> scale(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            scale[i] = std::max(scale[i], std::abs(a[i * n + j]));
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        double best = 0.0;
        for (std::size_t i = k; i < n; ++i) {
            const double score = scale[i] == 0.0 ? 0.0 : std::abs(a[i * n + k]) / scale[i];
            if (score > best) { best = score; pivot = i; }
        }
        if (best < 1e-13)
            throw std::runtime_error("Singular or numerically ill-conditioned COSMO matrix");
        if (pivot != k) {
            for (std::size_t j = k; j < n; ++j) std::swap(a[k * n + j], a[pivot * n + j]);
            std::swap(rhs[k], rhs[pivot]);
            std::swap(scale[k], scale[pivot]);
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            const double f = a[i * n + k] / a[k * n + k];
            a[i * n + k] = 0.0;
            if (f == 0.0) continue;
            for (std::size_t j = k + 1; j < n; ++j) a[i * n + j] -= f * a[k * n + j];
            rhs[i] -= f * rhs[k];
        }
    }
    for (std::size_t i = n; i-- > 0;) {
        for (std::size_t j = i + 1; j < n; ++j) rhs[i] -= a[i * n + j] * rhs[j];
        rhs[i] /= a[i * n + i];
        if (!std::isfinite(rhs[i])) throw std::runtime_error("Non-finite surface charge from solver");
    }
    return rhs;
}

std::string decimal(double value) {
    std::ostringstream stream;
    stream << std::setprecision(8) << value;
    return stream.str();
}
} // namespace

double EnergyTerms::correction() const {
    return atomic_area + solvent_shift + linear_charge + nonlinear_charge + volume + short_range;
}

double EnergyTerms::total() const { return electrostatic + correction(); }

Result calculate(Molecule molecule, Solvent solvent, const Options& options) {
    const auto& cp = class_parameters(solvent.solvent_class);
    if (!std::isfinite(solvent.epsilon) || solvent.epsilon < 1.0)
        throw std::runtime_error("The solvent dielectric constant must be finite and >= 1");
    const bool uses_shift = solvent.solvent_class == 'C' || solvent.solvent_class == 'D';
    if (!options.electrostatic_only && uses_shift && !solvent.xi)
        throw std::runtime_error("No verified xi_solv for '" + solvent.name +
            "'. Supply --xi explicitly from the paper's Table S1. "
            "--xi 0 requests the explicitly labelled zero-shift variant.");
    if (solvent.xi && !std::isfinite(*solvent.xi))
        throw std::runtime_error("xi_solv must be finite");
    if (!uses_shift && solvent.xi && std::abs(*solvent.xi) > 0)
        throw std::runtime_error("Classes A and B have no xi_solv term; nonzero --xi is invalid");
    if (molecule.atoms.empty()) throw std::runtime_error("The molecule contains no atoms");
    for (const auto& atom : molecule.atoms) {
        element_symbol(atom.atomic_number);
        if (!std::isfinite(atom.charge)) throw std::runtime_error("Atomic charge must be finite");
        for (double x : atom.position)
            if (!std::isfinite(x)) throw std::runtime_error("Atomic coordinates must be finite");
        if (!options.electrostatic_only) element_parameters(solvent.solvent_class, atom.atomic_number);
    }
    Result result;
    result.options = options;
    result.solvent = std::move(solvent);
    if (options.center_coordinates) {
        std::array<double, 3> center{};
        for (const auto& atom : molecule.atoms)
            for (std::size_t k = 0; k < 3; ++k) center[k] += atom.position[k];
        for (double& c : center) c /= static_cast<double>(molecule.atoms.size());
        for (auto& atom : molecule.atoms)
            for (std::size_t k = 0; k < 3; ++k) atom.position[k] -= center[k];
        result.warnings.push_back("Coordinates were explicitly centered. The finite-quadrature "
            "volume may differ from an uncentered reference run.");
    }
    result.molecule = std::move(molecule);
    if (options.ndiv != 3)
        result.warnings.push_back("NDIV differs from the reference setting 3; changing it does "
            "not constitute a re-fit of the published parameters or full COSMO-grid convergence.");
    if (!options.radius_overrides.empty())
        result.warnings.push_back("User-specified radii override the published unscaled-radius model.");
    if (std::abs(result.molecule.total_charge() - std::round(result.molecule.total_charge())) > 0.005)
        result.warnings.push_back("The input charge sum is " + decimal(result.molecule.total_charge()) +
            " e. Input charges have not been rounded or renormalized.");
    std::vector<SurfaceAtom> spheres;
    for (const auto& atom : result.molecule.atoms) {
        const auto override = options.radius_overrides.find(atom.atomic_number);
        const double radius = override == options.radius_overrides.end()
                            ? default_radius(atom.atomic_number) : override->second;
        result.radii.push_back(radius);
        spheres.push_back({atom.position, radius});
    }
    result.surface = make_surface(spheres, options.ndiv);
    const std::size_t n = result.surface.elements.size();
    if (n == 0) throw std::runtime_error("No exposed surface elements were generated");
    if (options.max_surface_points == 0 || n > options.max_surface_points)
        throw std::runtime_error("Generated " + std::to_string(n) +
            " surface points; the dense solver limit is " + std::to_string(options.max_surface_points) +
            ". Its storage is O(M^2) and factorization is O(M^3). "
            "Increase --max-points explicitly if sufficient resources are available.");
    if (n > std::numeric_limits<std::size_t>::max() / n / sizeof(double))
        throw std::runtime_error("COSMO matrix size overflows addressable memory");
    if (!options.electrostatic_only && cp.zeta != 0.0 && result.surface.volume <= 0.0)
        throw std::runtime_error("The finite-quadrature cavity volume is nonpositive. "
            "Check coordinates, or explicitly use --center for a far-translated molecule.");
    std::vector<double> potential(n, 0.0), rhs(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (const auto& atom : result.molecule.atoms) {
            const double r = distance(result.surface.elements[i].center, atom.position);
            if (r < 1e-10) throw std::runtime_error("A surface element coincides with an atomic center");
            potential[i] += atom.charge / r;
        }
        rhs[i] = -potential[i];
    }
    auto matrix = make_matrix(result.surface);
    result.cosmo_matrix_positive_definite = cholesky_solve(matrix, rhs, result.surface_charges);
    if (!result.cosmo_matrix_positive_definite) {
        matrix = make_matrix(result.surface);
        result.surface_charges = lu_solve(matrix, rhs);
        result.warnings.push_back("The discrete COSMO matrix is not positive definite. "
            "A pivoted LU solve was used without altering the published kernel; "
            "inspect surface discretization before relying on this result.");
    }
    // Check the original linear equations independently of the factorization.
    double residual = 0.0, matrix_norm = 0.0, rhs_norm = 0.0, charge_norm = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double row = potential[i], row_norm = 0.0;
        for (std::size_t j = 0; j < n; ++j) {
            const double aij = kernel(result.surface.elements[i], result.surface.elements[j], i == j);
            row += aij * result.surface_charges[j];
            row_norm += std::abs(aij);
        }
        residual = std::max(residual, std::abs(row));
        matrix_norm = std::max(matrix_norm, row_norm);
        rhs_norm = std::max(rhs_norm, std::abs(potential[i]));
        charge_norm = std::max(charge_norm, std::abs(result.surface_charges[i]));
    }
    const double denom = matrix_norm * charge_norm + rhs_norm;
    result.relative_residual = denom == 0.0 ? residual : residual / denom;
    if (!std::isfinite(result.relative_residual) || result.relative_residual > 1e-9)
        throw std::runtime_error("COSMO solve failed the linear-equation residual check");
    result.conductor_energy = 0.5 * coulomb_kcal_angstrom *
        std::inner_product(result.surface_charges.begin(), result.surface_charges.end(),
                           potential.begin(), 0.0);
    if (!std::isfinite(result.conductor_energy) || result.conductor_energy > 1e-7)
        throw std::runtime_error("Nonphysical positive/non-finite conductor energy; "
            "the surface discretization is numerically unstable");
    result.dielectric_factor = (result.solvent.epsilon - 1.0) /
        (result.solvent.epsilon + (result.solvent.solvent_class == 'D' ? 0.5 : 0.0));
    result.energy.electrostatic = result.dielectric_factor * result.conductor_energy;
    result.surface_charge_sum = std::accumulate(result.surface_charges.begin(),
                                                result.surface_charges.end(), 0.0);
    const double screening_error = std::abs(result.surface_charge_sum + result.molecule.total_charge());
    if (screening_error > 0.05)
        result.warnings.push_back("Conductor surface charge differs from minus the solute "
            "charge by " + decimal(screening_error) + " e. Charges were not rescaled.");
    const std::size_t atoms = result.molecule.atoms.size();
    result.atomic_surface_charges.assign(atoms, 0.0);
    result.atomic_sigma.assign(atoms, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        result.atomic_surface_charges[result.surface.elements[i].atom] += result.surface_charges[i];
    result.sigma_min = std::numeric_limits<double>::infinity();
    result.sigma_max = -std::numeric_limits<double>::infinity();
    for (std::size_t a = 0; a < atoms; ++a) {
        const double area = result.surface.atomic_areas[a];
        const double q = result.atomic_surface_charges[a];
        if (area > 0.0) {
            result.atomic_sigma[a] = q / area;
            result.sigma_min = std::min(result.sigma_min, result.atomic_sigma[a]);
            result.sigma_max = std::max(result.sigma_max, result.atomic_sigma[a]);
        }
        // Fully buried atoms have no q_A/S_A and do not enter SRC extrema.
        if (!options.electrostatic_only) {
            const auto& p = element_parameters(result.solvent.solvent_class,
                                               result.molecule.atoms[a].atomic_number);
            result.energy.atomic_area += p.kappa * area;
            result.energy.nonlinear_charge += p.p * (cp.code == 'D' ? std::abs(q) : q * q);
            if (cp.code == 'C') result.energy.linear_charge += p.g * q;
        }
    }
    if (!options.electrostatic_only) {
        result.energy.volume = cp.zeta * result.surface.volume;
        if (uses_shift) result.energy.solvent_shift = *result.solvent.xi * result.surface.area;
        if (cp.code == 'A' || cp.code == 'B') {
            const double arg = cp.src_a * std::pow(std::abs(result.sigma_min), 3) +
                               cp.src_b * std::pow(std::abs(result.sigma_max), 3);
            result.energy.short_range = cp.src_c * std::erf(arg);
        }
    }
    if (!std::isfinite(result.energy.total())) throw std::runtime_error("Non-finite total energy");
    return result;
}

} // namespace uese
