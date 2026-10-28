#!/usr/bin/env python3
"""Independent reference regressions and physical/numerical invariants.

Python is needed only to run these tests, never to run the C++ program.
The references are actual author-program outputs, not C++-generated gold files.
"""
import csv
import json
import io
import re
import math
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
EXE = Path(sys.argv[1] if len(sys.argv) > 1 else ROOT / "uese").resolve()
checks = 0


def require(condition, description):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(description)


def near(actual, expected, tolerance, description):
    require(abs(actual - expected) <= tolerance,
            f"{description}: got {actual:.15g}, expected {expected:.15g}, tolerance {tolerance}")


def run(path, *args):
    p = subprocess.run([str(EXE), str(path), *args, "--format", "json"],
                       text=True, capture_output=True)
    if p.returncode:
        raise AssertionError(f"Unexpected failure: {p.args}\n{p.stderr}")
    return json.loads(p.stdout)


def failure(path, args, message):
    p = subprocess.run([str(EXE), str(path), *args], text=True, capture_output=True)
    require(p.returncode != 0 and message.lower() in p.stderr.lower(),
            f"Expected rejection containing {message!r}, got {p.returncode}, {p.stderr!r}")


def xyz(path, atoms, title="CM5 numerical test", scale=1.0, charge_scale=1.0):
    text = f"{len(atoms)}\n{title}\n"
    for a in atoms:
        coords = " ".join(format(v * scale, ".17g") for v in a["xyz"])
        text += f"{a['element']} {coords} {a['input_charge'] * charge_scale:.17g}\n"
    path.write_text(text)
    return path


def main():
    pho = ROOT / "examples/PhO_minus.inp"
    nh3 = ROOT / "examples/NH3.inp"
    old = run(pho, "-solvent", "acetonitrile")
    require(old["n_surface_points"] == 387, "2021 reference has 387 surface points")
    near(old["surface_area"], 115.714, 0.0005, "2021 cavity area")
    near(old["cavity_volume"], 89.302, 0.0005, "2021 cavity volume")
    near(old["energy"]["electrostatic"], -66.946, 0.005, "2021 electrostatic energy")
    near(old["energy"]["correction"], 12.393, 0.005, "2021 correction")
    near(old["energy"]["total"], -54.554, 0.005, "2021 total")

    maximum = {"electrostatic": 0.0, "correction": 0.0, "total": 0.0}
    with (ROOT / "tests/reference_2025_Pho_92.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    require(len(rows) == 92, "complete 92-solvent reference catalog")
    for ref in rows:
        result = run(pho, "-solvent", ref["name"])
        require(result["solvent"]["class"] == ref["class"], f"class of {ref['name']}")
        near(result["solvent"]["epsilon"], float(ref["epsilon"]), 0.000051,
             f"epsilon of {ref['name']}")
        for term in maximum:
            delta = abs(result["energy"][term] - float(ref[term]))
            maximum[term] = max(maximum[term], delta)
            near(result["energy"][term], float(ref[term]), 0.01,
                 f"2025 PhO reference {ref['name']} / {term}")
        require(result["numerics"]["relative_residual"] < 1e-12,
                f"original-equation residual in {ref['name']}")

    # An independent neutral molecule catches ionic-only accidental calibration.
    nh3_refs = {
        "water": (-4.42, 0.30, -4.12),
        "methanol": (-4.34, 1.39, -2.95),
        "acetonitrile": (-4.35, 3.36, -0.99),
        "benzene": (-2.05, 0.91, -1.14),
    }
    for solvent, values in nh3_refs.items():
        r = run(nh3, "-solvent", solvent)
        for term, expected in zip(("electrostatic", "correction", "total"), values):
            near(r["energy"][term], expected, 0.01, f"2025 NH3 {solvent}/{term}")
        near(r["solute_charge"], 0.0001, 1e-12, "do not silently renormalize CM5 charges")

    listed = subprocess.run([str(EXE), "--list-solvents"], capture_output=True,
                            text=True, check=True).stdout.strip().splitlines()
    require(len(listed) == 93, "92 built-in solvents plus TSV header")

    # Independent reproduction of the original 2021 SI, Table S1, page 2.
    # The reference fixture was transcribed from the paper, not generated
    # from the C++ source or from the numerical reference output.
    with (ROOT / "data/table_s1_xi.tsv").open() as f:
        si_rows = list(csv.DictReader(f, delimiter="\t"))
    require(len(si_rows) == 77, "2021 SI Table S1 has 77 xi entries")
    require(sum(row["class"] == "C" for row in si_rows) == 20,
            "2021 SI Table S1 contains 20 C-class solvents")
    require(sum(row["class"] == "D" for row in si_rows) == 57,
            "2021 SI Table S1 contains 57 D-class solvents")
    listed_rows = list(csv.DictReader(io.StringIO("\n".join(listed)), delimiter="\t"))
    require(len(listed_rows) == 92, "C++ table has 92 solvents")
    lookup = {s["name"]: s for s in listed_rows}
    require(len(lookup) == 92, "all C++ solvent names are unique")
    for row in si_rows:
        item = lookup.get(row["cpp_name"])
        require(item is not None, f"2021 SI solvent missing: {row['si_name']}")
        require(item["class"] == row["class"], f"2021 SI class mismatch: {row['si_name']}")
        near(float(item["xi_kcal_mol_A2"]), float(row["xi_kcal_mol_A2"]),
             1e-12, f"2021 SI xi mismatch: {row['si_name']}")
        require("Table S1" in item["xi_source"],
                f"2021 SI source provenance missing: {row['si_name']}")
        resolved = run(pho, "--solvent", row["si_name"])
        require(resolved["solvent"]["name"] == row["cpp_name"],
                f"2021 SI printed solvent name not understood by C++: {row['si_name']}")
    for si_name, author_name in (("diphenyl ether", "phenylether"),
                                 ("tert-butylbenzene", "tbutylbenzene")):
        from_si = run(pho, "--solvent", si_name)
        from_author = run(pho, "--solvent", author_name)
        near(from_si["energy"]["total"], from_author["energy"]["total"],
             1e-12, f"2021 SI solvent-name alias: {si_name}")
    near(run(pho, "--solvent", "MeCN")["energy"]["total"], old["energy"]["total"],
         1e-12, "official solvent synonym")
    near(run(pho, "--class", "C", "--epsilon", "35.6881", "--xi", "0.0021")
         ["energy"]["total"], old["energy"]["total"], 1e-12,
         "explicit manual parameters reproduce named solvent")
    no_xi = run(pho, "--solvent", "acetonitrile", "--xi", "0")
    near(old["energy"]["total"] - no_xi["energy"]["total"],
         0.0021 * old["surface_area"], 1e-11, "xi only contributes xi times area")

    with tempfile.TemporaryDirectory(prefix="uese-tests-") as tmp:
        tmp = Path(tmp)
        # Cover every published element/class parameter against the independent
        # author executable, using synthetic overlapping tetrahedra. These are
        # numerical fixtures, not chemically meaningful solvation benchmarks.
        with (ROOT / "tests/reference_2025_elements_76.tsv").open() as f:
            element_refs = list(csv.DictReader(f, delimiter="\t"))
        require(len(element_refs) == 76, "all 76 element/class/charge-sign reference cases")
        for ref in element_refs:
            a = 1.5 * float(ref["author_radius"]) / (2 * math.sqrt(2))
            coords = [(a, a, a), (-a, -a, a), (-a, a, -a), (a, -a, -a)]
            p = tmp / "element_tetrahedron.xyzq"
            p.write_text("4\nSynthetic parameter coverage, not a chemical prediction\n" +
                         "".join(f"{ref['element']} {x:.8f} {y:.8f} {z:.8f} {ref['charge_per_atom']}\n"
                                 for x, y, z in coords))
            r = run(p, "--solvent", ref["solvent"])
            require(r["n_surface_points"] == int(ref["author_surface_points"]),
                    f"element surface count {ref['element']}/{ref['class']}")
            near(r["atoms"][0]["radius"], float(ref["author_radius"]), 1e-12,
                 f"element radius {ref['element']}")
            for term in ("electrostatic", "correction", "total"):
                near(r["energy"][term], float(ref["author_" + term]), 0.01,
                     f"element {ref['element']} class {ref['class']} q={ref['charge_per_atom']} {term}")
        baseline = run(nh3, "--electrostatic-only")
        atoms = baseline["atoms"]
        q2 = run(xyz(tmp / "doubled.xyzq", atoms, charge_scale=2.0), "--electrostatic-only")
        near(q2["energy"]["electrostatic"], 4 * baseline["energy"]["electrostatic"],
             1e-10, "electrostatic homogeneity E(2Q)=4E(Q)")
        q0 = run(xyz(tmp / "zero.xyzq", atoms, charge_scale=0.0), "--electrostatic-only")
        near(q0["energy"]["electrostatic"], 0, 1e-14, "zero charge has zero electrostatic energy")
        eps1 = run(nh3, "--epsilon", "1", "--electrostatic-only")
        near(eps1["energy"]["electrostatic"], 0, 1e-14, "epsilon=1 electrostatics vanish")
        require(baseline["energy"]["total"] is None, "COSMO-only must not report a uESE total")
        for solvent in ("water", "methanol", "acetonitrile", "benzene"):
            r = run(nh3, "--solvent", solvent)
            near(r["surface_charge"], baseline["surface_charge"], 1e-12,
                 "unscaled conductor surface charges must be solvent independent")

        # Unit parsing, arbitrary XYZ titles, D exponents, and count errors.
        bohr = xyz(tmp / "bohr.xyzq", atoms, "units=bohr", 1 / 0.529177210903)
        near(run(bohr, "--electrostatic-only")["energy"]["electrostatic"],
             baseline["energy"]["electrostatic"], 1e-10, "bohr input conversion")
        title = xyz(tmp / "title.xyzq", atoms, "Atomic coordinates are just this XYZ title")
        near(run(title, "--electrostatic-only")["energy"]["electrostatic"],
             baseline["energy"]["electrostatic"], 1e-10, "XYZ title is not an inp header")
        prose = xyz(tmp / "prose.xyzq", atoms,
                    "Coordinates in angstrom; atomic units are used only for charges.")
        near(run(prose, "--electrostatic-only")["energy"]["electrostatic"],
             baseline["energy"]["electrostatic"], 1e-10, "prose must not change coordinate units")
        unchanged = prose.read_bytes()
        failure(prose, ["--surface", str(prose.parent / "." / prose.name)],
                "must not replace the input")
        require(prose.read_bytes() == unchanged, "surface output cannot overwrite the input")
        fortran = tmp / "fortran.inp"
        fortran.write_text("Example\nAtomic coordinates\n6 0D0 0D0 0D0\n"
                           "Ground state charges\n1D0\n")
        sphere = run(fortran, "--electrostatic-only")
        radius = 1.7
        near(sphere["surface_area"], 4 * math.pi * radius**2, 1e-10, "single-sphere analytic area")
        near(sphere["cavity_volume"], 4 * math.pi * radius**3 / 3, 1e-10,
             "single-sphere analytic volume")
        born = -332.063713299 / (2 * radius)
        require(abs(sphere["numerics"]["conductor_energy"] / born - 1) < 2e-4,
                "single-sphere Born limit within the documented finite-mesh error")

        surface_csv = tmp / "surface.csv"
        exported = run(pho, "--surface", str(surface_csv))
        with surface_csv.open() as f:
            patches = list(csv.DictReader(f))
        require(len(patches) == 387, "CSV exports every tessera")
        near(sum(float(p["area_angstrom2"]) for p in patches), exported["surface_area"],
             1e-9, "CSV surface-area sum")
        near(sum(float(p["q_conductor_e"]) for p in patches), exported["surface_charge"],
             1e-12, "CSV charge sum")
        buried = tmp / "buried.xyzq"
        buried.write_text("2\nBuried hydrogen test\nC 0 0 0 -0.3\nH 0.1 0 0 0.3\n")
        b = run(buried)
        require(b["atoms"][1]["surface_area"] == 0 and b["atoms"][1]["sigma"] is None,
                "buried atoms must not introduce 0/0 charge-density extrema")
        for name, content, args, message in [
            ("bad_count.inp", "Atomic coordinates\n6 0 0 0\nGround state charges\n0 0\n",
             [], "Charge count"),
            ("nan.xyzq", "1\nBad\nC nan 0 0 0\n", [], "invalid finite number"),
            ("duplicate.xyzq", "2\nBad\nC 0 0 0 0\nH 0 0 0 0\n", [], "Coincident"),
            ("P.xyzq", "1\nPhosphorus\nP 0 0 0 0\n", ["--solvent", "acetonitrile"],
             "No published Table 1"),
            ("I.xyzq", "1\nIodine\nI 0 0 0 0\n", ["--solvent", "acetonitrile"],
             "No published Table 1"),
        ]:
            p = tmp / name
            p.write_text(content)
            failure(p, args, message)

    failure(pho, ["--ndiv", "0"], "out of range")
    failure(pho, ["--ndiv", "6"], "[1,5]")
    failure(pho, ["--max-points", "1"], "dense solver limit")
    failure(pho, ["--epsilon", "nan"], "finite number")
    failure(pho, ["--solvent", "unknown-solvent"], "Unknown solvent")
    failure(pho, ["--class", "C", "--epsilon", "20"], "No verified xi_solv")
    failure(pho, ["--xi", "0.1"], "no xi_solv term")
    failure(ROOT / "examples/B2H6.inp", [], "No published Table 1")
    failure(ROOT / "examples/Au20.inp", [], "No published Table 1")

    print(f"PASS: {checks} assertions")
    print("Author 2025 PhO / 92 solvents, maximum absolute differences (kcal/mol):")
    for key, value in maximum.items():
        print(f"  {key}: {value:.12f}")
    print("2021 paper Table 1 (separately reviewed), SI Table S1 (77 entries), 2021 printed")
    print("benchmark, neutral NH3, geometry, solver invariants, units, CSV, and rejections passed.")


if __name__ == "__main__":
    main()
