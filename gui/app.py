import streamlit as st
import os,sys
import tempfile
import importlib
import zipfile
import json

# Add project root to sys.path so 'src' is findable
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))

# others
if project_root not in sys.path:
	sys.path.insert(0, project_root)

src_path = os.path.join(project_root, "src")
if src_path not in sys.path:
	sys.path.insert(0, src_path)

src_path = os.path.join(project_root, "examples")
if src_path not in sys.path:
	sys.path.insert(0, src_path)    

# -----------------------------------------------------------------------------
# Configuration and UI Setup
# -----------------------------------------------------------------------------
st.set_page_config(layout="wide", page_title="Polymer-Grafted Nanostructure Builder")

st.title("PolyGraft: Molecular Grafting Tool")

# Define the two main columns
# col1: Input parameters (left)
# col2: Visualization/Rendering (right)
col1, col2 = st.columns([1, 1], gap="large")

# -----------------------------------------------------------------------------
# Left Column: User Parameters
# -----------------------------------------------------------------------------
with col1:
	st.header("⚙️ Configuration inputs")
	
	# 1. Base Settings
	col_fmt, col_res = st.columns(2)
	with col_fmt:
		data_format = st.selectbox("Data Format", ["GROMACS", "LAMMPS"])
	with col_res:
		resolution = st.selectbox("Resolution", ["Atomistic", "Coarse-Grained"])
	
	# 2. Conditional Units (If Coarse-Grained)
	units = None
	if resolution == "Coarse-Grained":
		units = st.radio("Units", ["Å", "LJ"], horizontal=True)
		
	# 3. Grafting Type
	grafting_type = st.selectbox("Grafting Type", ["Unigraft", "Bigraft"])

	st.divider()

	# 4. File Uploaders based on Format and Grafting Type
	st.subheader("Polymer Input Files")
	
	uploaded_files = {}
	
	if data_format == "GROMACS":
		if grafting_type == "Unigraft":
			uploaded_files["gro_1"] = st.file_uploader("Upload Polymer .gro", type=["gro"], key="gro_uni")
			uploaded_files["itp_1"] = st.file_uploader("Upload Polymer .itp", type=["itp"], key="itp_uni")
		elif grafting_type == "Bigraft":
			col_a, col_b = st.columns(2)
			with col_a:
				uploaded_files["gro_1"] = st.file_uploader("Upload Polymer 1 .gro", type=["gro"], key="gro_bi1")
				uploaded_files["itp_1"] = st.file_uploader("Upload Polymer 1 .itp", type=["itp"], key="itp_bi1")
			with col_b:
				uploaded_files["gro_2"] = st.file_uploader("Upload Polymer 2 .gro", type=["gro"], key="gro_bi2")
				uploaded_files["itp_2"] = st.file_uploader("Upload Polymer 2 .itp", type=["itp"], key="itp_bi2")
				
	elif data_format == "LAMMPS":
		if grafting_type == "Unigraft":
			uploaded_files["data_1"] = st.file_uploader("Upload Polymer .data", type=["data"], key="data_uni")
		elif grafting_type == "Bigraft":
			col_a, col_b = st.columns(2)
			with col_a:
				uploaded_files["data_1"] = st.file_uploader("Upload Polymer 1 .data", type=["data"], key="data_bi1")
			with col_b:
				uploaded_files["data_2"] = st.file_uploader("Upload Polymer 2 .data", type=["data"], key="data_bi2")
				
		# 5. LAMMPS specific input
		if resolution == "Coarse-Grained":
			if units == "Å":		
				atom_style = st.text_input("LAMMPS atom_style", value="id resid type charge x y z",
						help="Define the order of columns in your Atoms section. Default: id resid type charge x y z")
			else:
				atom_style = st.text_input("LAMMPS atom_style", value="id resid type x y z",
						help="Define the order of columns in your Atoms section. Default: id resid type x y z")

	st.divider()

	# 6. Substrate Geometry
	st.subheader("Substrate Geometry")
	geometry = st.selectbox("Shape", ["Slab", "Rod", "Pore", "Sphere"])
	
	geom_params = {}
	if geometry == "Slab":
		col_x, col_y, col_z = st.columns(3)
		with col_x: geom_params['Lx'] = st.number_input("Lx (Å or σ in LJ unit)", value=50.0)
		with col_y: geom_params['Ly'] = st.number_input("Ly (Å or σ in LJ unit)", value=50.0)
		with col_z: geom_params['Lz'] = st.number_input("Lz (Å or σ in LJ unit)", value=10.0)
	elif geometry == "Rod":
		col_r, col_l = st.columns(2)
		with col_r: geom_params['R'] = st.number_input("Radius (Å or σ in LJ unit)", value=10.0)
		with col_l: geom_params['L'] = st.number_input("Length (Å or σ in LJ unit)", value=50.0)
	elif geometry == "Pore":
		col_r, col_l = st.columns(2)
		with col_r: geom_params['R_in'] = st.number_input("Inner Radius (Å or σ in LJ unit)", value=30.0)
		with col_l: geom_params['L'] = st.number_input("Length (Å or σ in LJ unit)", value=80.0)
	elif geometry == "Sphere":
		geom_params['Radius'] = st.number_input("Radius (Å or σ in LJ unit)", value=20.0)

	st.caption("⚛️ Lattice Parameters")
	col_lat1, col_lat2 = st.columns(2)
	with col_lat1:
		lattice_type = st.selectbox("Lattice Type", ["fcc", "bcc", "sc"])
	with col_lat2:
		if units == "Å":		
			lattice_const = st.number_input("Lattice Constant (Å) or nearest neighbors in LJ (use 1)", value=4.08)
		else:
			lattice_const = st.number_input("Lattice Constant (Å) or nearest neighbors in LJ (use 1)", value=1.0)		

	st.divider()

	# 7. Lattice & Grafting Parameters
	st.subheader("Grafting Density")
	grafting_density = st.number_input("Grafting Density (chains/Å² or chains/σ²)", value=0.0120, format="%.4f")

	# 8. Generation Trigger
	st.divider()
	generate_btn = st.button("Generate Hybrid Structure", type="primary", use_container_width=True)

# -----------------------------------------------------------------------------
# Backend Routing & Execution Logic
# -----------------------------------------------------------------------------
def get_backend_module(fmt, res, g_type):
	"""Maps the UI selections to the corresponding examples/ backend directory."""
	is_cg = (res == "Coarse-Grained")
	is_gmx = (fmt == "GROMACS")
	is_bi = (g_type == "Bigraft")
	
	# Resolve Directory
	if not is_cg and not is_bi and is_gmx: return "examples_gmx"
	if not is_cg and not is_bi and not is_gmx: return "examples_lmp"
	if not is_cg and is_bi and is_gmx: return "examples_bi_gmx"
	if not is_cg and is_bi and not is_gmx: return "examples_bi_lmp"
	
	if is_cg and not is_bi and is_gmx: return "examples_cg_gmx"
	if is_cg and not is_bi and not is_gmx: return "examples_cg_lmp"
	if is_cg and is_bi and is_gmx: return "examples_cg_bi_gmx"
	if is_cg and is_bi and not is_gmx: return "examples_cg_bi_lmp"

# -----------------------------------------------------------------------------
# Right Column: Rendering Preview
# -----------------------------------------------------------------------------
with col2:
	st.header("Preview & Rendering")
	
	# Placeholder for the 3D viewer (e.g., py3Dmol or NGLview)
	preview_container = st.container(border=True, height=650)
	
	with preview_container:
		if not generate_btn:
			st.info("Set your parameters on the left and click **Generate Hybrid Structure** to view the preview.")
		else:

			# Create a temporary directory that automatically cleans up after the block finishes
			with tempfile.TemporaryDirectory() as tmpdirname:
				preview_container.success("Processing hybrid structure in temporary directory...")
				
				# Helper function to write Streamlit UploadedFile to the temp disk
				def save_uploadedfile(uploaded_file):
					if uploaded_file is not None:
						file_path = os.path.join(tmpdirname, uploaded_file.name)
						with open(file_path, "wb") as f:
							f.write(uploaded_file.getbuffer())
						return file_path
					return None

				# 1. Initialize file path variables
				file_paths = {}

				# 2. Save the uploaded files based on the format and grafting type
				if data_format == "GROMACS":
					if grafting_type == "Unigraft":
						file_paths['gro'] = save_uploadedfile(uploaded_files["gro_1"])
						file_paths['itp'] = save_uploadedfile(uploaded_files["itp_1"])
					elif grafting_type == "Bigraft":
						file_paths['gro1'] = save_uploadedfile(uploaded_files["gro_1"])
						file_paths['itp1'] = save_uploadedfile(uploaded_files["itp_1"])
						file_paths['gro2'] = save_uploadedfile(uploaded_files["gro_2"])
						file_paths['itp2'] = save_uploadedfile(uploaded_files["itp_2"])
				
				elif data_format == "LAMMPS":
					if grafting_type == "Unigraft":
						file_paths['data'] = save_uploadedfile(uploaded_files["data_1"])
					elif grafting_type == "Bigraft":
						file_paths['data1'] = save_uploadedfile(uploaded_files["data_1"])
						file_paths['data2'] = save_uploadedfile(uploaded_files["data_2"])

				# Check if required files were actually uploaded before proceeding
				missing_files = [k for k, v in file_paths.items() if v is None]

				# HARD STOP IF FILES ARE MISSING
				if missing_files or not file_paths:
					st.error(f"Missing or lost file uploads: {missing_files}. Please re-upload and try again.")
					st.stop() # This completely prevents the backend from running!

				else:				
					# Resolves to examples.examples_gmx, examples.examples_cg_bi_lmp, etc.
					target_backend_module = get_backend_module(data_format, resolution, grafting_type) 
					# Target the specific gen_*.py script inside examples*
					script_name = target_backend_module.replace("examples", "gen") + ".py" # e.g., "gen_lmp.py"
					script_path = os.path.join(project_root, "examples", target_backend_module, script_name)
					
					st.write("### Execution Log")
					st.write(f"Files saved to temp dir: `{tmpdirname}`")
					
					# Display payload to ensure parameters are passing correctly
					payload = {
						"Backend Target": f"{target_backend_module}/{script_name}",
						# "Saved Files": file_paths,
						"Format": data_format,
						"Resolution": resolution,
						"Unit": units if resolution == "Coarse-Grained" else "N/A",
						"Grafting Type": grafting_type,
						"Atom Style": atom_style if data_format == "LAMMPS" else "N/A",
						"Geometry": geometry,
						"Geometry Settings": geom_params,
						"Lattice Type": lattice_type,
						"Lattice Constant": lattice_const,
						"Grafting Density": grafting_density
					}
					
					with st.expander("View Backend Parameter Payload"):
						st.json(payload)
										
					# 4. Import the module dynamically and execute
					try:
						if os.path.isfile(script_path):							
							module_name = f"examples.{target_backend_module}.{target_backend_module.replace('examples', 'gen')}"
							gen_code = importlib.import_module(module_name)							
							st.info(f"Executing backend script: `{module_name}`...")
							
							with st.spinner(f"Generating {geometry} Brush..."):
								# CALLING THE BACKEND FUNCTION:
								if data_format == "GROMACS": # Atomistic or cg
									output_files = gen_code.gen(
										gro_file=file_paths.get('gro'),
										itp_file=file_paths.get('itp'),
										geometry=geometry,
										geom_params=geom_params,
										lattice_type=lattice_type,
										lattice_constant=lattice_const,
										grafting_density=grafting_density,
										output_dir=tmpdirname
									)

								elif data_format == "LAMMPS": # Atomistic or cg
									output_files = gen_code.gen(
										data_file=file_paths.get('data'),
										atom_style=atom_style,
										geometry=geometry,
										geom_params=geom_params,
										lattice_type=lattice_type,
										lattice_constant=lattice_const,
										grafting_density=grafting_density,
										output_dir=tmpdirname
									)
							
							# Handle the outputs (whether it's a single string or a tuple/list of strings)
							if output_files:
								# Convert a single string to a list so we can process it uniformly
								if isinstance(output_files, str):
									output_files = [output_files]
								
								# Filter out any paths that don't actually exist on the disk
								valid_files = [filepath for filepath in output_files if filepath and os.path.exists(filepath)]
								
								if valid_files:
									# === NEW: SAVE THE PAYLOAD TO A JSON FILE ===
									payload_filename = "generation_parameters.json"
									payload_filepath = os.path.join(tmpdirname, payload_filename)
									
									with open(payload_filepath, "w") as pf:
										json.dump(payload, pf, indent=4)
									
									# Add the payload file to our list of valid files to be zipped
									valid_files.append(payload_filepath)

									# Create a ZIP file inside the temporary directory
									zip_filename = "hybrid_structure.zip"
									zip_filepath = os.path.join(tmpdirname, zip_filename)
									
									with zipfile.ZipFile(zip_filepath, 'w', zipfile.ZIP_DEFLATED) as zipf:
										for filepath in valid_files:
											# Get just the file name (e.g., "output.gro")
											filename = os.path.basename(filepath)
											zipf.write(filepath, arcname=filename)
									
									# Offer the ZIP file for download
									with open(zip_filepath, "rb") as zfile:
										st.download_button(
											label=f"📦 Download All Output Files (.zip)",
											data=zfile,
											file_name=zip_filename,
											mime="application/zip",
											use_container_width=True
										)
									
									st.success(f"{geometry} brush generation complete! ({len(valid_files)} files bundled in zip)")
								else:
									st.error("Backend executed, but none of the returned file paths could be found on disk.")
							else:
								st.error(f"{module_name}.py executed, but nothing was returned.")

						else:
							st.error(f"Backend script missing! Expected to find `{script_name}` at:\n`{script_path}`")
							st.stop() # Halt execution so it doesn't crash trying to import it				

					except ModuleNotFoundError:
						st.error(f"Backend script `{module_name}.py` not found in `examples/examples_gmx/`. Ensure the directory structure is correct.")
					except AttributeError:
						st.error(f"The function `build_slab` was not found inside `{module_name}.py`. Please check the actual function name in your backend script.")
					except Exception as e:
						st.error(f"An error occurred during generation: {e}")

			# Temporary files deleted automatically after this block
			st.caption("Disk cleanup: Temporary files have been removed from the server.")
			