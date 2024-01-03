# 
# Copyright (C) Guang Chen
# 
# This file is part of polyGraft program
#
# polyGraft is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# polyGraft is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

import tkinter as tk
from tkinter import filedialog
import os,sys

pady=5

class polyGraftGUI:

	def __init__(self, master):
		self.master = master
		self.master.title("polyGraft")
		self.master.geometry("300x600")  # Setting window size

		self.current_step = tk.IntVar()
		self.current_step.set(1)

		# call the first step
		# self.create_widgets()
		self.step1()

	def get_resolution(self):

		self.model_resolution = self.resolution.get()
		print(f"Step 1 selected value: {self.model_resolution}")

	def get_format(self):

		self.model_format = self.format.get()
		print(f"Step 2 selected value: {self.model_format}")

	def get_graft_type(self):

		self.model_graft_type = self.graft_type.get()
		print(f"Step 3 selected value: {self.model_graft_type}")

	def get_substrate_type(self):

		self.model_substrate_type = self.substrate_type.get()
		print(f"Step 5 selected value: {self.model_substrate_type}")

	def get_slab_values(self):

		self.slab_length = self.ent_slab_length.get()
		self.slab_width = self.ent_slab_width.get()
		self.slab_height = self.ent_slab_height.get()

	def get_grafting_density(self):

		self.grafting_density = self.ent_grafting_density.get()

	def get_fname(self):

		self.outfname = self.ent_fname.get()

	def open_file_lmp_uni(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("data files", "*.data"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.lmp_uni_data_file = file_path
			self.ent_file_lmp_uni.delete(0, tk.END)
			self.ent_file_lmp_uni.insert(0, file_path)			

	def open_file_lmp_bi_1st(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("data files", "*.data"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.lmp_bi_data_file_1st = file_path
			self.ent_file_lmp_bi_1st.delete(0, tk.END)
			self.ent_file_lmp_bi_1st.insert(0, file_path)	

	def open_file_lmp_bi_2nd(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("data files", "*.data"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.lmp_bi_data_file_2nd = file_path
			self.ent_file_lmp_bi_2nd.delete(0, tk.END)
			self.ent_file_lmp_bi_2nd.insert(0, file_path)

	def open_file_gmx_uni_gro(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("gro files", "*.gro"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.gmx_uni_gro_file = file_path
			self.ent_file_gmx_uni_gro.delete(0, tk.END)
			self.ent_file_gmx_uni_gro.insert(0, file_path)		

	def open_file_gmx_uni_itp(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("itp files", "*.itp"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.gmx_uni_itp_file = file_path
			self.ent_file_gmx_uni_itp.delete(0, tk.END)
			self.ent_file_gmx_uni_itp.insert(0, file_path)	

	def open_file_gmx_bi_gro_1st(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("gro files", "*.gro"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.gmx_bi_gro_file_1st = file_path
			self.ent_file_gmx_bi_gro_1st.delete(0, tk.END)
			self.ent_file_gmx_bi_gro_1st.insert(0, file_path)

	def open_file_gmx_bi_itp_1st(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("itp files", "*.itp"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.gmx_bi_itp_file_1st = file_path
			self.ent_file_gmx_bi_itp_1st.delete(0, tk.END)
			self.ent_file_gmx_bi_itp_1st.insert(0, file_path)

	def open_file_gmx_bi_gro_2nd(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("gro files", "*.gro"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.gmx_bi_gro_file_2nd = file_path
			self.ent_file_gmx_bi_gro_2nd.delete(0, tk.END)
			self.ent_file_gmx_bi_gro_2nd.insert(0, file_path)

	def open_file_gmx_bi_itp_2nd(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("itp files", "*.itp"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.gmx_bi_itp_file_2nd = file_path
			self.ent_file_gmx_bi_itp_2nd.delete(0, tk.END)
			self.ent_file_gmx_bi_itp_2nd.insert(0, file_path)

	def step1(self):		
		self.frm_step1 = tk.Frame(self.master)

		self.lbl_step1 = tk.Label(self.frm_step1, text="Step 1: select data format")
		self.lbl_step1.pack(pady=pady)

		self.step1_options = ["GROMACS", "LAMMPS"]
		self.format = tk.StringVar()
		self.format.set(self.step1_options[0])

		self.select_box = tk.OptionMenu(self.frm_step1, self.format, *self.step1_options)
		self.select_box.pack(pady=pady)

		self.btn_format = tk.Button(self.frm_step1, text="OK", command=self.get_format)
		self.btn_format.pack(pady=pady)

		self.btn_step1 = tk.Button(self.frm_step1, text="Next", command=self.step2)
		self.btn_step1.pack(side=tk.RIGHT, pady=pady)

		self.frm_step1.pack()

	def step2(self):
		self.frm_step1.pack_forget()
		self.frm_step2 = tk.Frame(self.master)

		self.lbl_step2 = tk.Label(self.frm_step2, text="Step2: select resolution")
		self.lbl_step2.pack(pady=pady)

		if self.model_format == "GROMACS":
			self.step2_options = ["Atomistic"]
		else:
			self.step2_options = ["Atomistic", "Coarse-grained"]

		self.resolution = tk.StringVar()
		self.resolution.set(self.step2_options[0])		 
		
		self.select_box = tk.OptionMenu(self.frm_step2, self.resolution, *self.step2_options)
		self.select_box.pack(pady=pady)

		self.btn_resolution = tk.Button(self.frm_step2, text="OK", command=self.get_resolution)
		self.btn_resolution.pack(pady=pady)

		self.btn_step2 = tk.Button(self.frm_step2, text="Next", command=self.step3)
		self.btn_step2.pack(side=tk.RIGHT, pady=pady)

		# prev_button = tk.Button(self.frm_step2, text="Previous", command=self.step1)
		# prev_button.pack(side=tk.LEFT, pady=pady)

		self.frm_step2.pack()

	def step3(self):
		self.frm_step2.pack_forget()
		self.frm_step3 = tk.Frame(self.master)

		self.lbl_step3 = tk.Label(self.frm_step3, text="Step 3: select grafts type")
		self.lbl_step3.pack(pady=pady)

		self.step3_options = ["unigraft", "bigraft"]
		self.graft_type = tk.StringVar()
		self.graft_type.set(self.step3_options[0])

		self.select_box = tk.OptionMenu(self.frm_step3, self.graft_type, *self.step3_options)
		self.select_box.pack(pady=pady)

		self.btn_graft_type = tk.Button(self.frm_step3, text="OK", command=self.get_graft_type)
		self.btn_graft_type.pack(pady=pady)

		self.btn_step3 = tk.Button(self.frm_step3, text="Next", command=self.step4)
		self.btn_step3.pack(side=tk.RIGHT, pady=pady)

		# prev_button = tk.Button(self.frm_step3, text="Previous", command=self.step2)
		# prev_button.pack(side=tk.LEFT, pady=pady)

		self.frm_step3.pack()
	
	def step4(self):
		self.frm_step3.pack_forget()
		self.frm_step4 = tk.Frame(self.master)

		self.lbl_step4 = tk.Label(self.frm_step4, text="Step 4: import graft data file")
		self.lbl_step4.pack(pady=pady)

		if self.model_graft_type == "unigraft":
			self.step4_1()

		elif self.model_graft_type == "bigraft":
			self.step4_2()

		else:
			print(f"Unknown setting for graft_type ({self.graft_type})")
			sys.exit(0)		

		self.btn_graft_type = tk.Button(self.frm_step4, text="OK")
		self.btn_graft_type.pack(pady=pady)

		self.btn_step4 = tk.Button(self.frm_step4, text="Next", command=self.step5)
		self.btn_step4.pack(side=tk.RIGHT, pady=pady)

		# prev_button = tk.Button(self.frm_step4, text="Previous", command=self.step3)
		# prev_button.pack(side=tk.LEFT, pady=pady)

		self.frm_step4.pack()

	def step4_1(self):
		# for uni-grafts
		if self.model_format == "LAMMPS":
			# import .data file for graft			

			# file box
			btn_lmp_uni = tk.Button(self.frm_step4, text="Select LAMMPS DATA File", command=self.open_file_lmp_uni)
			btn_lmp_uni.pack(pady=pady)
			self.ent_file_lmp_uni = tk.Entry(self.frm_step4, width=20)
			self.ent_file_lmp_uni.pack(pady=pady)

		elif self.model_format == "GROMACS":
			# import .gro and .itp file for graft

			# file box
			btn_gmx_uni_gro = tk.Button(self.frm_step4, text="Select GROMACS gro File", command=self.open_file_gmx_uni_gro)
			btn_gmx_uni_gro.pack(pady=pady)
			self.ent_file_gmx_uni_gro = tk.Entry(self.frm_step4, width=20)
			self.ent_file_gmx_uni_gro.pack(pady=pady)			

			btn_gmx_uni_itp = tk.Button(self.frm_step4, text="Select GROMACS itp File", command=self.open_file_gmx_uni_itp)
			btn_gmx_uni_itp.pack(pady=pady)
			self.ent_file_gmx_uni_itp = tk.Entry(self.frm_step4, width=20)
			self.ent_file_gmx_uni_itp.pack(pady=pady)

	def step4_2(self):
		# for bi-grafts
		if self.model_format == "LAMMPS":
			# import .data file for graft

			# for the first graft
			btn_lmp_bi = tk.Button(self.frm_step4, text="Select LAMMPS DATA File for the 1st graft", command=self.open_file_lmp_bi_1st)
			btn_lmp_bi.pack(pady=pady)
			self.ent_file_lmp_bi_1st = tk.Entry(self.frm_step4, width=20)
			self.ent_file_lmp_bi_1st.pack(pady=pady)
			
			# for the second graft
			btn_lmp_bi = tk.Button(self.frm_step4, text="Select LAMMPS DATA File for the 2nd graft", command=self.open_file_lmp_bi_2nd)
			btn_lmp_bi.pack(pady=pady)
			self.ent_file_lmp_bi_2nd = tk.Entry(self.frm_step4, width=20)
			self.ent_file_lmp_bi_2nd.pack(pady=pady)			

		elif self.model_format == "GROMACS":
			# import .gro and .itp file for graft

			# for the first graft
			btn_gmx_bi_gro = tk.Button(self.frm_step4, text="Select GROMACS gro File for the 1st graft", command=self.open_file_gmx_bi_gro_1st)
			btn_gmx_bi_gro.pack(pady=pady)
			self.ent_file_gmx_bi_gro_1st = tk.Entry(self.frm_step4, width=20)
			self.ent_file_gmx_bi_gro_1st.pack(pady=pady)			

			btn_gmx_bi_itp = tk.Button(self.frm_step4, text="Select GROMACS itp File for the 1st graft", command=self.open_file_gmx_bi_itp_1st)
			btn_gmx_bi_itp.pack(pady=pady)
			self.ent_file_gmx_bi_itp_1st = tk.Entry(self.frm_step4, width=20)
			self.ent_file_gmx_bi_itp_1st.pack(pady=pady)			

			# for the second graft
			btn_gmx_bi_gro = tk.Button(self.frm_step4, text="Select GROMACS gro File for the 2nd graft", command=self.open_file_gmx_bi_gro_2nd)
			btn_gmx_bi_gro.pack(pady=pady)
			self.ent_file_gmx_bi_gro_2nd = tk.Entry(self.frm_step4, width=20)
			self.ent_file_gmx_bi_gro_2nd.pack(pady=pady)			

			btn_gmx_bi_itp = tk.Button(self.frm_step4, text="Select GROMACS itp File for the 2nd graft", command=self.open_file_gmx_bi_itp_2nd)
			btn_gmx_bi_itp.pack(pady=pady)
			self.ent_file_gmx_bi_itp_2nd = tk.Entry(self.frm_step4, width=20)
			self.ent_file_gmx_bi_itp_2nd.pack(pady=pady)			

	def step5(self):
		# set substrate shapes
		self.frm_step4.pack_forget()
		self.frm_step5 = tk.Frame(self.master)

		self.lbl_step5 = tk.Label(self.frm_step5, text="Step 5: select substrate type")
		self.lbl_step5.pack(pady=pady)

		self.step5_options = ["slab", "rod", "pore", "sphere"]
		self.substrate_type = tk.StringVar()
		self.substrate_type.set(self.step5_options[0])

		self.select_box = tk.OptionMenu(self.frm_step5, self.substrate_type, *self.step5_options)
		self.select_box.pack(pady=pady)

		self.btn_substrate_type = tk.Button(self.frm_step5, text="OK", command=self.get_substrate_type)
		self.btn_substrate_type.pack(pady=pady)

		self.btn_step5 = tk.Button(self.frm_step5, text="Next", command=self.step6)
		self.btn_step5.pack(side=tk.RIGHT, pady=pady)

		# prev_button = tk.Button(self.frm_step5, text="Previous", command=self.step4)
		# prev_button.pack(side=tk.LEFT, pady=pady)

		self.frm_step5.pack()

	def step6(self):
		self.frm_step5.pack_forget()
		self.frm_step6 = tk.Frame(self.master)

		if self.model_substrate_type == "slab":
			self.step6_slab()

		elif self.model_substrate_type == "rod":
			self.step6_rod()

		elif self.model_substrate_type == "pore":
			self.step6_pore()

		elif self.model_substrate_type == "sphere":
			self.step6_sphere()

		else:
			print(f"Unknown setting for substrate_type ({self.model_substrate_type})")
			sys.exit(0)

		self.btn_step6 = tk.Button(self.frm_step6, text="Next", command=self.step7)
		self.btn_step6.pack(side=tk.RIGHT, pady=pady)

		# prev_button = tk.Button(self.frm_step6, text="Previous", command=self.step5)
		# prev_button.pack(side=tk.LEFT, pady=pady)

		self.frm_step6.pack()

	def step6_slab(self):		

		# set substrate sizes
		self.lbl_step6 = tk.Label(self.frm_step6, text="Step 6: set substrate sizes")
		self.lbl_step6.pack(pady=pady)

		lbl_text = tk.Label(self.frm_step6, text="Please input the sizes of the substrate")
		lbl_text.pack(pady=pady)

		# three entry
		frm_entry = tk.Frame(self.frm_step6)
		frm_entry.rowconfigure([0, 1, 2], minsize=2, weight=1)
		frm_entry.columnconfigure([0, 1], minsize=2, weight=1)

		lbl_slab_length = tk.Label(frm_entry, text="Length:")		
		self.ent_slab_length = tk.Entry(frm_entry, width=10)
		lbl_slab_length.grid(row=0, column=0, sticky="w", pady=pady)		
		self.ent_slab_length.grid(row=0, column=1, sticky="e", pady=pady)

		lbl_slab_width = tk.Label(frm_entry, text="Width:")
		self.ent_slab_width = tk.Entry(frm_entry, width=10)
		lbl_slab_width.grid(row=1, column=0, sticky="w", pady=pady)
		self.ent_slab_width.grid(row=1, column=1, sticky="e", pady=pady)

		lbl_slab_height = tk.Label(frm_entry, text="Height:")
		self.ent_slab_height = tk.Entry(frm_entry, width=10)
		lbl_slab_height.grid(row=2, column=0, sticky="w", pady=pady)
		self.ent_slab_height.grid(row=2, column=1, sticky="e", pady=pady)

		frm_entry.pack()

		self.btn_graft_type = tk.Button(self.frm_step6, text="OK", command=self.get_slab_values)
		self.btn_graft_type.pack(pady=pady)

	def step7(self):
		# set grafting density
		self.frm_step6.pack_forget()
		self.frm_step7 = tk.Frame(self.master)

		self.lbl_step6 = tk.Label(self.frm_step7, text="Step 7: set grafting density")
		self.lbl_step6.pack(pady=pady)

		lbl_text = tk.Label(self.frm_step7, text="Please input the grafting density")
		lbl_text.pack(pady=pady)

		# one entry
		lbl_grafting_density = tk.Label(self.frm_step7, text="Grafting density:")
		self.ent_grafting_density = tk.Entry(self.frm_step7, width=10)
		self.ent_grafting_density.pack(pady=pady)

		self.btn_graft_type = tk.Button(self.frm_step7, text="OK", command=self.get_grafting_density)
		self.btn_graft_type.pack(pady=pady)

		self.btn_step7 = tk.Button(self.frm_step7, text="Next", command=self.step8)
		self.btn_step7.pack(side=tk.RIGHT, pady=pady)

		# prev_button = tk.Button(self.frm_step7, text="Previous", command=self.step6)
		# prev_button.pack(side=tk.LEFT, pady=pady)

		self.frm_step7.pack()

	def step8(self):
		# set output file name
		self.frm_step7.pack_forget()
		self.frm_step8 = tk.Frame(self.master)

		self.lbl_step8 = tk.Label(self.frm_step8, text="Step 8: set output file name")
		self.lbl_step8.pack(pady=pady)

		lbl_text = tk.Label(self.frm_step8, text="Please input the output file name")
		lbl_text.pack(pady=pady)

		# one entry
		self.ent_fname = tk.Entry(self.frm_step8, width=10)
		self.ent_fname.pack(pady=pady)

		self.btn_fname = tk.Button(self.frm_step8, text="OK", command=self.get_fname)
		self.btn_fname.pack(pady=pady)

		# prev_button = tk.Button(self.frm_step8, text="Previous", command=self.step7)
		# prev_button.pack(side=tk.LEFT, pady=pady)

		self.btn_step8 = tk.Button(self.frm_step8, text="Submit", command=self.call_polyGraft)
		self.btn_step8.pack(pady=pady)

		self.frm_step8.pack()

	def call_polyGraft(self):
		print("Generating the structure and topology ...")

def main():
	root = tk.Tk()
	app = polyGraftGUI(root)
	root.mainloop()

if __name__ == "__main__":
	main()

