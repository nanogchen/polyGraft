import streamlit as st
from stmol import showmol
import py3Dmol
import io
import zipfile
from datetime import datetime
import pandas as pd
import tempfile
import os
import sys

# 1. Get the absolute path to the directory containing 'src' and 'gui'
# This assumes your structure is: project_root/gui/app.py and project_root/src/polyGraft.py
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))

# 2. Add project root to sys.path so 'src' is findable
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# 3. Add 'src' itself to sys.path so polyGraft.py can find rtp_define.py
src_path = os.path.join(project_root, "src")
if src_path not in sys.path:
    sys.path.insert(0, src_path)

# NOW perform your imports
from src.polyGraft import polyGraft
from src.polymer import Polymer
from src.crystal import Crystal
from src.atomsk import Atomsk

# --- Helper: Global Reset ---
def reset_all():
    # Clears all session data to start fresh
    for key in list(st.session_state.keys()):
        del st.session_state[key]
    st.rerun()

# --- Helper: Molecule Viewer ---
def render_molecule(file_obj, file_format, label):
    """Renders molecular files in 3D using py3Dmol."""
    if file_obj is not None:
        st.write(f"🔍 **{label} Preview**")
        try:
            string_data = file_obj.getvalue().decode("utf-8")
            view = py3Dmol.view(width=450, height=350)
            # GROMACS uses 'gro'; LAMMPS data often maps best to 'xyz' for basic atom viewing
            fmt_type = 'gro' if file_format == "GROMACS" else 'xyz' 
            view.addModel(string_data, fmt_type)
            view.setStyle({'stick': {'radius': 0.2}, 'sphere': {'radius': 0.4}})
            view.zoomTo()
            showmol(view, height=350, width=450)
        except Exception as e:
            st.error(f"Could not render {label}: {e}")

def render_lammps_preview(file_obj, atom_style_str):
    if file_obj is None:
        return

    # 1. Parse atom_style for indices
    cols = atom_style_str.split()
    try:
        ix, iy, iz = cols.index('x'), cols.index('y'), cols.index('z')
        it = cols.index('type')
    except ValueError:
        st.error(f"Mapping Error: Style '{atom_style_str}' is missing x, y, z, or type.")
        return

    # 2. Extract Data
    lines = file_obj.getvalue().decode("utf-8").splitlines()
    xyz_data = []
    reading = False
    
    for line in lines:
        clean_line = line.split('#')[0].strip() # Remove comments
        if not clean_line: continue
        
        # Detect Start of Atoms section
        if "Atoms" in clean_line:
            reading = True
            continue
        
        # Detect End of Atoms section
        if reading and any(k in clean_line for k in ["Bonds", "Velocities", "Angles", "Dihedrals"]):
            reading = False
            break
            
        if reading:
            parts = clean_line.split()
            # Ensure the line actually looks like atom data (starts with an ID number)
            if len(parts) >= len(cols) and parts[0].isdigit():
                try:
                    # Assign a color based on type (Type 1=H, 2=He, 6=C, etc.)
                    elements = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"]
                    elem = elements[int(parts[it]) % len(elements)]
                    xyz_data.append(f"{elem} {parts[ix]} {parts[iy]} {parts[iz]}")
                except (ValueError, IndexError):
                    continue

    # 3. Visualization
    if xyz_data:
        st.success(f"Parsed {len(xyz_data)} atoms successfully.")
        xyz_string = f"{len(xyz_data)}\nConverted LAMMPS\n" + "\n".join(xyz_data)
        
        view = py3Dmol.view(width=500, height=400)
        view.addModel(xyz_string, 'xyz')
        # OVITO style: Spacefill spheres, no connecting 'net' lines
        view.setStyle({'sphere': {'colorscheme': 'Jmol', 'radius': 0.8}})
        view.zoomTo()
        showmol(view, height=400, width=500)
    else:
        st.error("No atom data found. Please verify that your 'Atoms' section starts with the word 'Atoms' and contains numeric data rows.")
        # Debugging: show the first 20 lines so you can see why it failed
        with st.expander("View File Header (Debug)"):
            st.code("\n".join(lines[:30]))

def get_file_content(file_input):
    """Safely extracts content from UploadedFile or local file buffer."""
    if file_input is None:
        return None
    if hasattr(file_input, "getvalue"):
        return file_input.getvalue().decode("utf-8")
    if hasattr(file_input, "read"):
        data = file_input.read()
        return data.decode("utf-8") if isinstance(data, bytes) else data
    return None

def render_substrate(file_input, filename, atom_style_str="id type x y z"):
    """
    Renders substrate based on file extension (PDB, GRO, or LAMMPS).
    """
    content = get_file_content(file_input)
    if not content:
        st.warning("No content found for substrate.")
        return

    ext = filename.split('.')[-1].lower()
    
    view = py3Dmol.view(width=500, height=400)

    if ext in ['pdb', 'gro']:
        # Use py3Dmol's native parser for PDB and GRO
        view.addModel(content, ext)
        # Style: CPK colors for standard elements
        view.setStyle({'sphere': {'colorscheme': 'Jmol', 'radius': 0.7}})
        st.caption(f"Rendering {ext.upper()} Substrate")
    
    elif ext in ['data', 'lammps', 'txt']:
        # Fallback to our custom LAMMPS parser for .data files
        xyz_string = parse_lammps_to_xyz(content, atom_style_str)
        if xyz_string:
            view.addModel(xyz_string, 'xyz')
            view.setStyle({'sphere': {'colorscheme': 'Jmol', 'radius': 0.7}})
            st.caption("Rendering LAMMPS Substrate")
        else:
            st.error("Could not parse LAMMPS data.")
            return
    else:
        st.error(f"Unsupported format: {ext}")
        return

    view.zoomTo()
    showmol(view, height=400, width=500)

def parse_lammps_to_xyz(content, atom_style_str):
    """Helper to convert LAMMPS 'Atoms' section to XYZ format."""
    cols = atom_style_str.split()
    try:
        ix, iy, iz, it = cols.index('x'), cols.index('y'), cols.index('z'), cols.index('type')
    except ValueError:
        return None

    lines = content.splitlines()
    xyz_atoms = []
    reading = False
    for line in lines:
        line = line.strip().split('#')[0]
        if "Atoms" in line:
            reading = True; continue
        if reading and any(k in line for k in ["Bonds", "Velocities", "Angles"]):
            reading = False; break
        if reading:
            p = line.split()
            if len(p) >= len(cols) and p[0].isdigit():
                elems = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"]
                e = elems[int(p[it]) % len(elems)]
                xyz_atoms.append(f"{e} {p[ix]} {p[iy]} {p[iz]}")
    
    if not xyz_atoms: return None
    return f"{len(xyz_atoms)}\nConverted\n" + "\n".join(xyz_atoms)
    
# --- App Styling ---
st.set_page_config(page_title="PolyGraft Dashboard", layout="wide")

# --- Initialize Session State ---
if 'history' not in st.session_state:
    st.session_state.history = []
if 'output_ready' not in st.session_state:
    st.session_state.output_ready = None

# --- Header ---
head_col1, head_col2 = st.columns([3, 1])
with head_col1:
    st.title("🧪 PolyGraft: Molecular Grafting Tool")
with head_col2:
    st.button("🔄 Reset All Parameters", width='stretch', on_click=reset_all)

st.divider()

# --- Main Layout ---
col_left, col_right = st.columns([1, 1], gap="large")

with col_left:
    st.subheader("⚙️ Configuration inputs")
    
    # 1. Basics
    with st.expander("Step 1 & 2: Format & Resolution", expanded=True):
        fmt = st.radio("Data Format:", ["GROMACS", "LAMMPS"], horizontal=True)
        res = st.radio("Resolution:", ["Atomistic", "Coarse Grained"], horizontal=True)

        if res == "Coarse Grained":
            # st.markdown("---")
            units = "Å" # Default for Atomistic            
            units = st.radio(
                "Select CG Units:",
                ["Å (Real)", "LJ (Reduced)"],horizontal=True,
                help="Choose between real units (Angstroms) or Lennard-Jones reduced units.",
                key="cg_unit_radio"
            )
            # Clean the string for backend processing
            units = "Å" if "Å" in units else "LJ"

    # 2. Graft Logic
    with st.expander("Step 3 & 4: Grafting & Files", expanded=True):
        gtype = st.radio("Graft Type:", ["Unigraft", "Bigraft"], horizontal=True)
        uploaded_files = {}
        
        atom_style = "N/A"
        if gtype == "Unigraft":
            if fmt == "GROMACS":
                uploaded_files['gro'] = st.file_uploader("Upload .gro", type=['gro'])
                uploaded_files['itp'] = st.file_uploader("Upload .itp", type=['itp'])
            else:
                atom_style = st.text_input(
                    "LAMMPS Atom Style Columns:", 
                    value="id resid type charge x y z",
                    key="lammps_style_uni", # <--- Unique Key
                    help="Define the order of columns in your Atoms section. Default: id resid type charge x y z"
                )
                uploaded_files['data'] = st.file_uploader("Upload LAMMPS Data", type=['data', 'txt'])
        else:
            tabA, tabB = st.tabs(["🧬 Graft A", "🧬 Graft B"])
            with tabA:
                if fmt == "GROMACS":
                    uploaded_files['gro_a'] = st.file_uploader("Gro A", type=['gro'])
                    uploaded_files['itp_a'] = st.file_uploader("Itp A", type=['itp'])
                else:
                    atom_style = st.text_input(
                    "LAMMPS Atom Style Columns:", 
                    value="id resid type charge x y z",
                    key="lammps_style_bi1", # <--- Unique Key                    
                    help="Define the order of columns in your Atoms section. Default: id resid type charge x y z"
                    )
                    uploaded_files['data_a'] = st.file_uploader("Data A", type=['data', 'txt'])
            with tabB:
                if fmt == "GROMACS":
                    uploaded_files['gro_b'] = st.file_uploader("Gro B", type=['gro'])
                    uploaded_files['itp_b'] = st.file_uploader("Itp B", type=['itp'])
                else:
                    atom_style = st.text_input(
                    "LAMMPS Atom Style Columns:", 
                    value="id resid type charge x y z",
                    key="lammps_style_bi2", # <--- Unique Key
                    help="Define the order of columns in your Atoms section. Default: id resid type charge x y z"
                    )
                    uploaded_files['data_b'] = st.file_uploader("Data B", type=['data', 'txt'])

            # bigrafting type
            bi_gft_fmt = st.radio("bi-Grafting mode:", ["homogeneous", "random", "Janus"], horizontal=True)

    # 3. Substrate Logic
    with st.expander("Step 5 & 6: Substrate Geometry", expanded=True):
        substrate = st.selectbox("Substrate Type:", ["Slab", "Rod", "Pore", "Sphere"])
        # NEW: Lattice Configuration Knobs
        # st.markdown("---")
        dims = {}
        c1, c2, c3 = st.columns(3)
        if substrate == "Slab":
            dims['Lx'] = c1.number_input("Lx [Å]", 0.1, 100.0, 50.0)
            dims['Ly'] = c2.number_input("Ly [Å]", 0.1, 100.0, 50.0)
            dims['Lz'] = c3.number_input("Lz [Å]", 0.1, 100.0, 10.0)
        elif substrate == "Rod":
            dims['R'] = c1.number_input("Radius [Å]", 0.1, 50.0, 20.0)
            dims['L'] = c2.number_input("Length [Å]", 0.1, 200.0, 20.0)
        elif substrate == "Pore":
            dims['R_in'] = c1.number_input("Inner R [Å]", 0.1, 50.0, 10.0)
            # dims['R_out'] = c2.number_input("Outer R [Å]", 0.1, 50.0, 50.0)
            dims['L'] = c2.number_input("Length [Å]", 0.1, 100.0, 10.0)
        else: # Sphere
            dims['R'] = c1.number_input("Radius [Å]", 0.1, 50.0, 20.0)

        st.caption("⚛️ Lattice Parameters")
        col_lat1, col_lat2 = st.columns(2)
        with col_lat1:
            lattice_type = st.text_input(
                "Lattice Type", 
                value="fcc", 
                key="input_lattice_type"
            )
        with col_lat2:
            lattice_const = st.number_input(
                "Lattice Constant [Å]", 
                value=4.08, 
                format="%.2f", 
                key="input_lattice_const"
            )

        # Create Tabs for Uploading vs Generating
        tab_generate,tab_upload = st.tabs(["🏗️ Generate Substrate","📤 Upload Substrate Files"])

        with tab_generate:
            st.info("Will generate the substrate, press generate.")                        

            if st.button("🏗️ Generate Substrate Now", use_container_width=True, key="btn_gen_sub"):
                with st.spinner(f"Generating {lattice_type} {substrate}..."):
                    try:
                        lattice = Atomsk(lattice_type=lattice_type, 
                                            lattice_const=lattice_const, 
                                            element='Au')

                        # --- CALL BACKEND SUBSTRATE GENERATOR ---
                        if substrate == "Slab":
                            ofile="Auslab.pdb"
                            lattice.gen_slab(dims['Lx'],dims['Ly'],dims['Lz'],outFile=ofile)
                            nanoslab = Crystal("nanoslab", 'Au', dims['Lx'],dims['Ly'],dims['Lz'])
                            nanoslab.readPDB(ofile, guessing_bond=True, lattice_const=lattice_const)

                        elif substrate == "Rod":
                            ofile="Aurod.pdb"
                            lattice.gen_rod(dims['R'], dims['L'], outFile=ofile)
                            nanorod = Crystal("nanorod", 'Au', dims['R'], dims['L'])
                            nanorod.readPDB(ofile, guessing_bond=True, lattice_const=lattice_const)

                        elif substrate == "Pore":
                            ofile="Aupore.pdb"
                            lattice.gen_pore(dims['R_in'], dims['L'], outFile=ofile)
                            nanopore = Crystal("nanopore", 'Au', dims['R_in'], dims['L'])
                            nanopore.readPDB(ofile, guessing_bond=True, lattice_const=lattice_const)

                        else:
                            ofile="AuNP.pdb"
                            lattice.gen_particle(dims['R'], outFile=ofile)
                            nanoparticle = Crystal("nanoparticle", 'Au', dims['R'])
                            nanoparticle.readPDB(ofile, guessing_bond=True, lattice_const=lattice_const)
                        
                        # For demo, we simulate a successful generation
                        st.session_state['generated_sub_ready'] = True
                        st.success(f"Successfully generated {substrate} ({lattice_type})!")
                        
                        # You could automatically load this into the 'uploaded_files' dict
                        # or display a specialized preview here.
                    except Exception as e:
                        st.error(f"Generation failed: {e}")

                    #save the path
                    st.session_state['last_gen_sub_path'] = ofile

        with tab_upload:
            uploaded_files['substrate'] = st.file_uploader("Upload substrate .gro/.pdb", type=['gro','pdb'])

            if st.button("🏗️ Upload Substrate Now", use_container_width=True, key="btn_upload_sub"):
                with st.spinner(f"Generating {lattice_type} {substrate}..."):
                    try:

                        if uploaded_files is not None:
                            with tempfile.NamedTemporaryFile(delete=False, suffix=f"_{uploaded_files['substrate'].name}") as tmp:
                                tmp.write(uploaded_files['substrate'].getvalue())
                                tmp_path = tmp.name  # This is the string filename the backend wants

                            # --- CALL BACKEND SUBSTRATE GENERATOR ---
                            if substrate == "Slab":                            
                                # lattice.gen_slab(dims['Lx'],dims['Ly'],dims['Lz'],outFile="Auslab.pdb")
                                nanoslab = Crystal("nanoslab", 'Au', dims['Lx'],dims['Ly'],dims['Lz'])
                                nanoslab.readPDB(tmp_path, guessing_bond=True, lattice_const=lattice_const)

                            elif substrate == "Rod":
                                # lattice.gen_rod(dims['R'], dims['L'], outFile="Aurod.pdb")
                                nanorod = Crystal("nanorod", 'Au', dims['R'], dims['L'])
                                nanorod.readPDB(tmp_path, guessing_bond=True, lattice_const=lattice_const)

                            elif substrate == "Pore":
                                # lattice.gen_pore(dims['R_in'], dims['L'], outFile="Aupore.pdb")
                                nanopore = Crystal("nanopore", 'Au', dims['R_in'], dims['L'])
                                nanopore.readPDB(tmp_path, guessing_bond=True, lattice_const=lattice_const)

                            else:
                                # lattice.gen_particle(dims['R'], outFile="AuNP.pdb")
                                nanoparticle = Crystal("nanoparticle", 'Au', dims['R'])
                                nanoparticle.readPDB(tmp_path, guessing_bond=True, lattice_const=lattice_const)
                            
                            # For demo, we simulate a successful generation
                            # st.session_state['generated_sub_ready'] = True
                            st.success(f"Successfully loaded {substrate} ({lattice_type})!")

                    except Exception as e:
                        st.error(f"Generation failed: {e}")
            
                    #save the path
                    st.session_state['last_gen_sub_path'] = tmp_path

    # 4. Density & Run
    with st.expander("Step 7 & 8: Density & Output", expanded=True):
        density = st.slider("Grafting Density (chains/nm²):", 0.01, 5.0, 0.5)
        default_fname = f"{substrate}_{density:.2f}"
        out_name = st.text_input("Output Filename:", 
                                key="oname",
                                value=default_fname)
        
        st.divider()
        # if st.button("🚀 GENERATE SYSTEM", width='stretch', type="primary"):
        #     # --- START BACKEND INTEGRATION ---
        #     # Replace the lines below with: result = your_main_function(fmt, res, gtype, uploaded_files, dims, density)
        
        # --- Step 8: Finalize & Run ---
        if st.button("🚀 GENERATE SYSTEM", use_container_width=True, type="primary"):
            with st.spinner("Executing Molecular Grafting Algorithms..."):
                
                # 1. Create a temporary directory to handle files
                with tempfile.TemporaryDirectory() as tmpdir:
                    
                    # 2. Save Uploaded Files to the Temp Directory
                    # This allows your backend to read them via standard file paths
                    paths = {}
                    for key, file_obj in uploaded_files.items():
                        if file_obj:
                            temp_path = os.path.join(tmpdir, file_obj.name)
                            with open(temp_path, "wb") as f:
                                f.write(file_obj.getbuffer())
                            paths[key] = temp_path

                    try:
                        # 3. CALL YOUR BACKEND FUNCTION HERE
                        # Example: 
                        from my_backend import run_grafting
                        output_path = run_grafting(
                            fmt=fmt, 
                            res=res, 
                            graft_type=gtype, 
                            structure_path=paths.get('gro') or paths.get('data'),
                            density=density,
                            **dims  # Passes Lx, Ly, etc. from the dims dictionary
                        )
                        
                        # --- FOR DEMO: Simulated result generation ---
                        # In reality, your function would return a path or a string
                        result_content = f"Generated {res} {fmt} system\nSubstrate: {substrate}\nDensity: {density}"
                        
                        # 4. Store the result in Session State for the Download Button
                        st.session_state.output_ready = {
                            "name": out_name, 
                            "fmt": fmt,
                            "content": result_content.encode() # Convert string to bytes
                            #         "content_gro": b"Simulated GROMACS .gro file content",
                            #         "content_itp": b"Simulated GROMACS .itp file content",
                            #         "content_lammps": b"Simulated LAMMPS .data file content"
                        }
                        
                        # 5. Log to History
                        st.session_state.history.append({
                            "Time": datetime.now().strftime("%H:%M:%S"),
                            "Project": out_name,
                            "Format": fmt,
                            "Res": res,
                            "Substrate": substrate
                        })
                        
                        st.success("System generated successfully!")
                        
                    except Exception as e:
                        st.error(f"Backend Error: {str(e)}")

                    st.session_state['final_assembly_path'] = final_path
            
with col_right:
    st.subheader("🖼️ Molecular Preview")
    with st.expander("🧬 Graft Visualization", expanded=True):
        # This section reacts live to uploads in the left column
        if gtype == "Unigraft":
            f_main = uploaded_files.get('gro') or uploaded_files.get('data')
            if f_main: 
                if fmt == "GROMACS": render_molecule(f_main, fmt, "Primary Graft")
                else: render_lammps_preview(f_main, atom_style)
        else:
            f_a = uploaded_files.get('gro_a') or uploaded_files.get('data_a')
            f_b = uploaded_files.get('gro_b') or uploaded_files.get('data_b')
            if f_a: 
                if fmt == "GROMACS": render_molecule(f_a, fmt, "Graft A")
                else: render_lammps_preview(f_a, atom_style)
            if f_b: 
                st.divider()
                if fmt == "GROMACS": render_molecule(f_b, fmt, "Graft B")
                else: render_lammps_preview(f_b, atom_style)

            if not f_a and not f_b: st.info("Upload A/B structure files to preview.")

    # 2. Substrate Preview Section
    st.divider()
    with st.expander("🧱 Substrate Visualization", expanded=True):
        # Priority 1: Check for a freshly generated backend file
        gen_path = st.session_state.get('last_gen_sub_path')
        
        # Priority 2: Check for an uploaded substrate file
        uploaded_sub = uploaded_files.get('substrate')

        if gen_path and os.path.exists(gen_path):
            fname = os.path.basename(gen_path)
            with open(gen_path, 'rb') as f:
                render_substrate(f, fname, atom_style)
                
        elif uploaded_sub:
            # Streamlit UploadedFile has a .name attribute
            render_substrate(uploaded_sub, uploaded_sub.name, atom_style)
        else:
            st.info("No substrate generated/loaded.")

    # 3. poly-substrate Preview Section
    st.divider()
    with st.expander("🚀 Assembled System Visualization", expanded=True):
        pg_path = st.session_state.get('final_assembly_path')

        if pg_path:
            pass

        else:
            st.info("Not yet generated!")

# --- Download Area ---
if st.session_state.output_ready:
    st.divider()
    d_col1, d_col2 = st.columns([1, 1])
    with d_col1:
        out = st.session_state.output_ready
        if out['fmt'] == "GROMACS":
            zip_buf = io.BytesIO()
            with zipfile.ZipFile(zip_buf, "w") as z:
                z.writestr(f"{out['name']}.gro", out['content_gro'])
                z.writestr(f"{out['name']}.itp", out['content_itp'])
            st.download_button(
                "💾 DOWNLOAD ZIP (GRO/ITP)", 
                zip_buf.getvalue(), 
                f"{out['name']}.zip", 
                width='stretch'
            )
        else:
            st.download_button(
                "💾 DOWNLOAD DATA FILE", 
                out['content_lammps'], 
                f"{out['name']}.data", 
                width='stretch'
            )

# --- History Log ---
st.divider()
st.subheader("📜 Session History")
if st.session_state.history:
    st.dataframe(pd.DataFrame(st.session_state.history), width='stretch')
