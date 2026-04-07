import streamlit as st
from stmol import showmol
import py3Dmol
import io
import zipfile
from datetime import datetime
import pandas as pd
import tempfile
import os

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
    st.subheader("⚙️ Configuration Knobs")
    
    # 1. Basics
    with st.expander("Step 1 & 2: Format & Resolution", expanded=True):
        fmt = st.radio("Data Format:", ["GROMACS", "LAMMPS"], horizontal=True)
        res = st.radio("Resolution:", ["Atomistic", "Coarse Grained"], horizontal=True)

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

    # 3. Substrate Logic
    with st.expander("Step 5 & 6: Substrate Geometry", expanded=True):
        substrate = st.selectbox("Substrate Type:", ["Slab", "Rod", "Pore", "Sphere"])
        dims = {}
        c1, c2, c3 = st.columns(3)
        if substrate == "Slab":
            dims['Lx'] = c1.number_input("Lx [nm]", 0.1, 100.0, 10.0)
            dims['Ly'] = c2.number_input("Ly [nm]", 0.1, 100.0, 10.0)
            dims['Lz'] = c3.number_input("Lz [nm]", 0.1, 100.0, 5.0)
        elif substrate == "Rod":
            dims['R'] = c1.number_input("Radius [nm]", 0.1, 50.0, 2.0)
            dims['L'] = c2.number_input("Length [nm]", 0.1, 200.0, 20.0)
        elif substrate == "Pore":
            dims['R_in'] = c1.number_input("Inner R [nm]", 0.1, 50.0, 3.0)
            dims['R_out'] = c2.number_input("Outer R [nm]", 0.1, 50.0, 5.0)
            dims['L'] = c3.number_input("Length [nm]", 0.1, 100.0, 10.0)
        else: # Sphere
            dims['R'] = c1.number_input("Radius [nm]", 0.1, 50.0, 5.0)

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
            
with col_right:
    st.subheader("🖼️ Molecular Preview")
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
