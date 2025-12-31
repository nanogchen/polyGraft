import streamlit as st
from stmol import showmol
import py3Dmol
import io
import zipfile

# --- App Styling ---
st.set_page_config(page_title="Molecular Grafting Tool", layout="wide")

# Initialize Session State
if 'step' not in st.session_state:
    st.session_state.step = 1
if 'data' not in st.session_state:
    st.session_state.data = {}
if 'output_files' not in st.session_state:
        st.session_state.output_files = None    

# --- Navigation Logic ---
def next_step(): st.session_state.step += 1
def prev_step(): st.session_state.step -= 1

# --- Molecule Viewer Function ---
def render_molecule(file_obj, file_format, label):
    if file_obj is not None:
        st.caption(f"**Preview: {label}**")
        string_data = file_obj.getvalue().decode("utf-8")
        view = py3Dmol.view(width=400, height=400)
        # .gro is standard; .data often needs conversion, 
        # but 3Dmol tries its best with 'xyz' or 'pdb' logic
        fmt = 'gro' if file_format == "GROMACS" else 'xyz' 
        view.addModel(string_data, fmt)
        view.setStyle({'stick': {}, 'sphere': {'radius': 0.5}})
        view.zoomTo()
        showmol(view, height=400, width=400)

# --- Sidebar Progress Tracker ---
with st.sidebar:
    st.title("Step Tracker")
    steps_list = [
        "1. Data Format", "2. Resolution", "3. Graft Type", 
        "4. File Import", "5. Substrate Type", "6. Dimensions", 
        "7. Grafting Density", "8. Generate"
    ]
    
    for i, step_name in enumerate(steps_list, 1):
        if st.session_state.step == i:
            st.markdown(f"### **➡️ {step_name}**") 
        elif st.session_state.step > i:
            st.markdown(f"✅ {step_name}")
        else:
            st.markdown(f"⚪ {step_name}")
            
    st.divider()
    st.subheader("Current Selection Summary")
    for key, value in st.session_state.data.items():
        # if key not in ['files', 'dimensions']:
        if key not in ['files']:
            st.caption(f"**{key.replace('_',' ').capitalize()}:** {value}")

    if st.button("🔄 Reset All"):
        st.session_state.step = 1
        st.session_state.data = {}
        st.rerun()

st.title("🧪 Molecular Grafting Configuration")

# --- Step 1: Data Format ---
if st.session_state.step == 1:
    st.header("Step 1: Select Data Format")
    st.session_state.data['format'] = st.radio("Choose Format:", ["GROMACS", "LAMMPS"])
    st.button("Next ➡️", on_click=next_step)

# --- Step 2: Resolution ---
elif st.session_state.step == 2:
    st.header("Step 2: Select Resolution")
    is_gromacs = st.session_state.data.get('format') == "GROMACS"
    
    if is_gromacs:
        st.info("💡 GROMACS selected: Only **Atomistic** resolution is supported.")
        res_options = ["Atomistic"]
    else:
        res_options = ["Atomistic", "Coarse Grained"]

    st.session_state.data['resolution'] = st.radio("Choose Resolution:", res_options)
    col1, col2 = st.columns(2)
    with col1: st.button("⬅️ Back", on_click=prev_step)
    with col2: st.button("Next ➡️", on_click=next_step)

# --- Step 3: Graft Type ---
elif st.session_state.step == 3:
    st.header("Step 3: Select Graft Type")
    st.session_state.data['graft_type'] = st.radio("Choose Type:", ["Unigraft", "Bigraft"])
    col1, col2 = st.columns(2)
    with col1: st.button("⬅️ Back", on_click=prev_step)
    with col2: st.button("Next ➡️", on_click=next_step)

# --- Step 4: Refined File Import (Conditional Logic) ---
elif st.session_state.step == 4:
    fmt = st.session_state.data.get('format')
    gtype = st.session_state.data.get('graft_type')
    st.header(f"Step 4: Import {fmt} Files ({gtype})")
    
    files = {}    
    ready = False
    if gtype == "Unigraft":
        c1, c2 = st.columns([1, 1])
        with c1:
            if fmt == "GROMACS":
                files['gro'] = st.file_uploader("Upload .gro", type=['gro'])
                files['itp'] = st.file_uploader("Upload .itp", type=['itp'])
                ready = all([files['gro'], files['itp']])
            else:
                files['data'] = st.file_uploader("Upload LAMMPS data", type=['data'])
                ready = files['data'] is not None
        with c2:
            f_to_show = files.get('gro') or files.get('data')
            if f_to_show: render_molecule(f_to_show, fmt, "Single Graft")

    else: # Bigraft 
        input_col, view = st.columns([1, 1])
        with input_col:
            if fmt == "GROMACS":
                st.subheader("Graft A")
                files['gro_a'] = st.file_uploader("gro A", type=['gro'])
                files['itp_a'] = st.file_uploader("itp A", type=['itp'])
                st.divider()
                st.subheader("Graft B")
                files['gro_b'] = st.file_uploader("gro B", type=['gro'])
                files['itp_b'] = st.file_uploader("itp B", type=['itp'])
                ready = all([files['gro_a'], files['itp_a'], files['gro_b'], files['itp_b']])
            else:
                st.subheader("Graft A")
                files['data_a'] = st.file_uploader("data A", type=['data'])
                st.divider()
                st.subheader("Graft B")
                files['data_b'] = st.file_uploader("data B", type=['data'])
                ready = all([files['data_a'], files['data_b']])
        
        with view:
            f_a = files.get('gro_a') or files.get('data_a')
            if f_a: render_molecule(f_a, fmt, "Graft A")
            f_b = files.get('gro_b') or files.get('data_b')
            if f_b: render_molecule(f_b, fmt, "Graft B")

    st.session_state.data['files'] = files
    col1, col2 = st.columns(2)
    with col1: st.button("⬅️ Back", on_click=prev_step)
    with col2: 
        if ready: st.button("Next ➡️", on_click=next_step)
        else: st.warning("Please upload all required files to proceed.")

# --- Step 5: Substrate Type ---
elif st.session_state.step == 5:
    st.header("Step 5: Select Substrate Type")
    st.session_state.data['substrate'] = st.selectbox("Choose Substrate:", ["Slab", "Rod", "Pore", "Sphere"])
    col1, col2 = st.columns(2)
    with col1: st.button("⬅️ Back", on_click=prev_step)
    with col2: st.button("Next ➡️", on_click=next_step)

# --- Step 6: Substrate Dimensions ---
elif st.session_state.step == 6:
    st.header(f"Step 6: Dimensions for {st.session_state.data['substrate']}")
    sub = st.session_state.data['substrate']
    dims = {}
    # if sub == "Slab":
    #     col_x, col_y, col_z = st.columns(3)
    #     dims['L_x'] = col_x.number_input("Length (x) [nm]", min_value=1.0, value=10.0)
    #     dims['L_y'] = col_y.number_input("Width (y) [nm]", min_value=1.0, value=10.0)
    #     dims['L_z'] = col_z.number_input("Height (z) [nm]", min_value=1.0, value=5.0)
        
    # elif sub == "Rod":
    #     col_r, col_l = st.columns(2)
    #     dims['radius'] = col_r.number_input("Radius [nm]", min_value=1.0, value=2.0)
    #     dims['length'] = col_l.number_input("Length [nm]", min_value=1.0, value=10.0)
        
    # elif sub == "Pore":
    #     col_ir, col_or, col_l = st.columns(3)
    #     dims['inner_radius'] = col_ir.number_input("Inner Radius [nm]", min_value=1.0, value=3.0)
    #     dims['outer_radius'] = col_or.number_input("Outer Radius [nm]", min_value=4.0, value=6.0)
    #     dims['length'] = col_l.number_input("Length [nm]", min_value=1.0, value=15.0)
        
    # elif sub == "Sphere":
    #     dims['radius'] = st.number_input("Radius [nm]", min_value=1.0, value=5.0)    
    if sub == "Slab":
        dims['L_x'] = st.number_input("Length [nm]", min_value=1.0, value=5.0)
        dims['L_y'] = st.number_input("Width [nm]", min_value=1.0, value=5.0)
        dims['L_z'] = st.number_input("Height [nm]", min_value=1.0, value=5.0)
    elif sub == "Rod":
        dims['radius'] = st.number_input("Radius [nm]", min_value=1.0)
        dims['length'] = st.number_input("Length [nm]", min_value=1.0)
    elif sub == "Pore":
        dims['inner_R'] = st.number_input("Inner Radius [nm]", min_value=1.0)
        dims['outer_R'] = st.number_input("Outer Radius [nm]", min_value=1.0)
        dims['length'] = st.number_input("Length [nm]", min_value=1.0)
    elif sub == "Sphere":
        dims['radius'] = st.number_input("Radius [nm]", min_value=1.0)

    st.session_state.data['dimensions'] = dims
    col1, col2 = st.columns(2)
    with col1: st.button("⬅️ Back", on_click=prev_step)
    with col2: st.button("Next ➡️", on_click=next_step)

# --- Step 7: Grafting Density ---
elif st.session_state.step == 7:
    st.header("Step 7: Set Grafting Density")
    # st.session_state.data['density'] = st.slider("Density (chains/nm²)", 0.01, 5.0, 0.5)
    sigma = {}
    sigma['sigma'] = st.number_input("Density (chains/nm²)", min_value=0.05)
    st.session_state.data['grafting_density'] = sigma
    col1, col2 = st.columns(2)
    with col1: st.button("⬅️ Back", on_click=prev_step)
    with col2: st.button("Next ➡️", on_click=next_step)

# --- Step 8: Finalize ---
elif st.session_state.step == 8:
    st.header("Step 8: Generate")
    out_name = st.text_input("Project Name (for filenames)", "my_molecular_system")
    
    # We use a dictionary to store multiple generated files    

    if st.button("🚀 Run Backend Generation"):
        with st.spinner("Grafting molecules and preparing files..."):
            
            # --- BACKEND MOCKUP START ---
            # Replace these with your actual backend function outputs
            if st.session_state.data['format'] == "GROMACS":
                # Simulated file content
                gro_content = f"GROMACS Generated File\nName: {out_name}\n"
                itp_content = f"GROMACS Topology File\nName: {out_name}\n"
                
                # Create a ZIP in memory
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED, False) as zip_file:
                    zip_file.writestr(f"{out_name}.gro", gro_content)
                    zip_file.writestr(f"{out_name}.itp", itp_content)
                
                st.session_state.output_files = {
                    "data": zip_buffer.getvalue(),
                    "extension": "zip",
                    "mime": "application/zip",
                    "label": "💾 Download GROMACS Files (.zip)"
                }
                
            else: # LAMMPS
                data_content = f"LAMMPS Data File\nName: {out_name}\n"
                st.session_state.output_files = {
                    "data": data_content,
                    "extension": "data",
                    "mime": "text/plain",
                    "label": "💾 Download LAMMPS File (.data)"
                }
            # --- BACKEND MOCKUP END ---
            
            st.success("Generation Complete!")

    # Display the download button if files are ready
    if st.session_state.output_files:
        st.divider()
        st.download_button(
            label=st.session_state.output_files['label'],
            data=st.session_state.output_files['data'],
            file_name=f"{out_name}.{st.session_state.output_files['extension']}",
            mime=st.session_state.output_files['mime']
        )
    
    st.button("⬅️ Back", on_click=prev_step)
