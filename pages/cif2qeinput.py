import streamlit as st
from ase.io import read
from ase.data import atomic_masses, chemical_symbols
from ase.geometry import cell_to_cellpar
import spglib

# Helper
def ase_symbol_to_z(symbol):
    return chemical_symbols.index(symbol)
    
def determine_ibrav(cell):
    a, b, c, alpha, beta, gamma = cell_to_cellpar(cell)
    a, b, c = round(a, 3), round(b, 3), round(c, 3)
    alpha, beta, gamma = round(alpha, 3), round(beta, 3), round(gamma, 3)

    if a == b == c and alpha == beta == gamma == 90.0:
        return 1  # cubic
    elif a == b != c and alpha == beta == gamma == 90.0:
        return 6  # tetragonal
    elif a != b != c and alpha == beta == gamma == 90.0:
        return 8  # orthorhombic
    elif a == b != c and alpha == beta == 90.0 and gamma == 120.0:
        return 4  # hexagonal
    elif a == b == c and alpha == beta == gamma != 90.0:
        return 5  # rhombohedral
    elif alpha == gamma == 90.0 and beta != 90.0:
        return 12  # monoclinic (base-centered)
    elif alpha != 90.0 or beta != 90.0 or gamma != 90.0:
        return 14  # triclinic
    else:
        return 0  # unknown

def get_spacegroup_info(atoms):
    lattice = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    dataset = spglib.get_symmetry_dataset((lattice, positions, numbers))
    return dataset['international'], dataset['number']

def generate_qe_input(atoms, pseudo_map, calc_type, ecutwfc, ecutrho, raw_ibrav, kpoints):
    species = sorted(set(atoms.get_chemical_symbols()))
    missing = [el for el in species if el not in pseudo_map]
    if missing:
        st.error(f"Missing pseudopotentials for: {', '.join(missing)}")
        return None, None

    lines = []
    lines.append("&control")
    lines.append(f"  calculation = '{calc_type}',")
    lines.append("  prefix = 'calc',")
    lines.append("  outdir = './tmp',")
    lines.append("/\n")

    lines.append("&system")
    lines.append(f"  ibrav = {raw_ibrav},")
    lines.append(f"  nat = {len(atoms)},")
    lines.append(f"  ntyp = {len(species)},")
    lines.append(f"  ecutwfc = {ecutwfc},")
    lines.append(f"  ecutrho = {ecutrho},")
    lines.append("/\n")

    lines.append("&electrons")
    lines.append("  conv_thr = 1.0d-6")
    lines.append("/\n")

    lines.append("CELL_PARAMETERS angstrom")
    for vec in atoms.get_cell():
        lines.append("  {:.10f} {:.10f} {:.10f}".format(*vec))

    lines.append("\nATOMIC_SPECIES")
    for el in species:
        z = ase_symbol_to_z(el)
        mass = atomic_masses[z]
        lines.append(f"{el}  {mass:.4f}  {pseudo_map[el].name}")

    lines.append("\nATOMIC_POSITIONS angstrom")
    for atom in atoms:
        pos = atom.position
        lines.append(f"{atom.symbol}  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}")

    lines.append("\nK_POINTS automatic")
    lines.append("  " + " ".join(kpoints) + " 0 0 0")

    qe_input_str = "\n".join(lines)
    run_script = f"""#!/bin/bash
pw.x < qe_input.in > output.log
"""

    return qe_input_str, run_script

# ------------------- STREAMLIT APP -------------------

st.set_page_config("CIF to QE", layout="centered")

# --- INIT SESSION STATE ---
if "page" not in st.session_state:
    st.session_state.page = 1
if "atoms" not in st.session_state:
    st.session_state.atoms = None
if "species" not in st.session_state:
    st.session_state.species = []

# --- PAGE 1: CIF Upload ---
if st.session_state.page == 1:
    st.title("📦 CIF to Quantum ESPRESSO Converter")
    st.header("Step 1: Upload a CIF File")

    cif_file = st.file_uploader("Upload CIF file", type=["cif"])

    if cif_file:
        try:
            atoms = read(cif_file)
            species = sorted(set(atoms.get_chemical_symbols()))
            cell = atoms.get_cell()
            a, b, c, alpha, beta, gamma = cell_to_cellpar(cell)
            
            spg_name, spg_num = get_spacegroup_info(atoms)
            st.success(f"🧬 Space Group: **{spg_name} ({spg_num})**")

            raw_ibrav = determine_ibrav(cell)
            st.markdown(f"**Lattice Parameters:** a = {a:.3f}, b = {b:.3f}, c = {c:.3f} Å")
            st.markdown(f"**Angles:** α = {alpha:.2f}°, β = {beta:.2f}°, γ = {gamma:.2f}°")

            force_ibrav0 = st.checkbox("🔁 Force `ibrav = 0` (include full CELL_PARAMETERS)", value=(raw_ibrav == 0))
            ibrav = 0 if force_ibrav0 else raw_ibrav
            st.info(f"Using `ibrav = {ibrav}`")
            
            st.success(f"✅ CIF parsed. Detected elements: {', '.join(species)}")

            if st.button("➡ Continue to Upload Pseudopotentials"):
                st.session_state.atoms = atoms
                st.session_state.species = species
                st.session_state.ibrav = ibrav
                st.session_state.page = 2
                st.experimental_rerun()
        except Exception as e:
            st.error(f"❌ Error reading CIF file:\n{e}")

# --- PAGE 2: Upload UPFs and Parameters ---
elif st.session_state.page == 2:
    st.title("🧪 QE Input Setup")
    st.header("Step 2: Upload Pseudopotentials per Element")

    atoms = st.session_state.atoms
    species = st.session_state.species
    ibrav = st.session_state.ibrav

    pseudo_map = {}
    for el in species:
        pseudo = st.file_uploader(f"Select .UPF for `{el}`", type=["upf", "UPF"], key=el)
        if pseudo:
            pseudo_map[el] = pseudo

    st.divider()
    st.subheader("⚙ QE Parameters")
    col1, col2, col3 = st.columns(3)
    with col1:
        calc_type = st.selectbox("Calculation Type", ["scf", "relax", "vc-relax", "nscf"])
    with col2:
        ecutwfc = st.number_input("ecutwfc (Ry)", min_value=10.0, value=40.0)
    with col3:
        ecutrho = st.number_input("ecutrho (Ry)", min_value=40.0, value=320.0)

    kpt_input = st.text_input("K-points Grid (e.g., 4 4 4)", "4 4 4")
    kpoints = kpt_input.strip().split()

    # Back button
    if st.button("⬅ Go Back"):
        st.session_state.page = 1
        st.experimental_rerun()

    if st.button("⚙ Generate QE Input File"):
        if len(pseudo_map) < len(species):
            st.error("Missing pseudopotentials for some elements.")
        elif len(kpoints) != 3 or not all(k.isdigit() for k in kpoints):
            st.error("Invalid k-point grid.")
        else:
            qe_input, run_script = generate_qe_input(
                atoms, pseudo_map, calc_type, ecutwfc, ecutrho, ibrav, kpoints
            )
            if qe_input:
                st.success("✅ QE input file generated!")

                st.subheader("📄 Input File Preview")
                st.code(qe_input, language="text")
                st.download_button("⬇ Download QE Input", qe_input, file_name="qe_input.in")

                st.subheader("📜 Run Script")
                st.code(run_script, language="bash")
                st.download_button("⬇ Download Run Script", run_script, file_name="pw.sh")
