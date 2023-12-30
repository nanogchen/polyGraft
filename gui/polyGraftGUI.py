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

step_dict = {1: "select resolution",
			 2: "select data format",
			 3: "select grafts type",
			 4: "import graft file",
			 5: "substrate type",
			 6: "substrate sizes",
			 7: "grafting density",
			 8: "output"}

class polyGraftGUI:

	def __init__(self, master):
		self.master = master
		self.master.title("polyGraft")
		self.master.geometry("400x600")  # Setting window size

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

	def open_file_lmp_bi_1st(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("data files", "*.data"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.lmp_bi_data_file_1st = file_path

	def open_file_lmp_bi_2nd(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("data files", "*.data"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.lmp_bi_data_file_2nd = file_path

	def open_file_gmx_uni_gro(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("gro files", "*.gro"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.gmx_uni_gro_file = file_path

	def open_file_gmx_uni_itp(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("itp files", "*.itp"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.gmx_uni_itp_file = file_path

	def open_file_gmx_bi_gro_1st(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("gro files", "*.gro"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.gmx_bi_gro_file_1st = file_path

	def open_file_gmx_bi_itp_1st(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("itp files", "*.itp"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.gmx_bi_itp_file_1st = file_path

	def open_file_gmx_bi_gro_2nd(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("gro files", "*.gro"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.gmx_bi_gro_file_2nd = file_path

	def open_file_gmx_bi_itp_2nd(self):
		file_path = filedialog.askopenfilename(
			title="Select File", filetypes=(("itp files", "*.itp"), ("All files", "*.*"))
		)
		if file_path:
			print("Selected file:", file_path)
			# Do something with the selected file path, e.g., display it in a label or process the file
			self.gmx_bi_itp_file_2nd = file_path

	def step1(self):
		self.lbl_step1 = tk.Label(self.master, text="Step1: select resolution")
		self.lbl_step1.pack(pady=pady)

		self.step1_options = ["Atomistic", "Coarse grained"]
		self.resolution = tk.StringVar()
		self.resolution.set(self.step1_options[0])		 
		
		self.select_box = tk.OptionMenu(self.master, self.resolution, *self.step1_options)
		self.select_box.pack(pady=pady)

		self.btn_resolution = tk.Button(self.master, text="OK", command=self.get_resolution)
		self.btn_resolution.pack(pady=pady)

		self.btn_step1 = tk.Button(self.master, text="Next", command=self.step2)
		self.btn_step1.pack(pady=pady)

		self.update_windows()

	def step2(self):

		self.lbl_step2 = tk.Label(self.master, text="Step 2: select data format")
		self.lbl_step2.pack(pady=pady)

		self.step2_options = ["GROMACS", "LAMMPS"]
		self.format = tk.StringVar()
		self.format.set(self.step2_options[0])

		self.select_box = tk.OptionMenu(self.master, self.format, *self.step2_options)
		self.select_box.pack(pady=pady)

		self.btn_format = tk.Button(self.master, text="OK", command=self.get_format)
		self.btn_format.pack(pady=pady)

		self.btn_step2 = tk.Button(self.master, text="Next", command=self.step3)
		self.btn_step2.pack(pady=pady)

		prev_button = tk.Button(self.master, text="Previous", command=self.step1)
		prev_button.pack(pady=pady)

		self.update_windows()

	def step3(self):
		self.lbl_step3 = tk.Label(self.master, text="Step 3: select grafts type")
		self.lbl_step3.pack(pady=pady)

		self.step3_options = ["unigraft", "bigraft"]
		self.graft_type = tk.StringVar()
		self.graft_type.set(self.step3_options[0])

		self.select_box = tk.OptionMenu(self.master, self.graft_type, *self.step3_options)
		self.select_box.pack(pady=pady)

		self.btn_graft_type = tk.Button(self.master, text="OK", command=self.get_graft_type)
		self.btn_graft_type.pack(pady=pady)

		self.btn_step3 = tk.Button(self.master, text="Next", command=self.step4)
		self.btn_step3.pack(pady=pady)

		prev_button = tk.Button(self.master, text="Previous", command=self.step2)
		prev_button.pack(pady=pady)

		self.update_windows()
	
	def step4(self):

		self.lbl_step4 = tk.Label(self.master, text="Step 4: import graft data file")
		self.lbl_step4.pack(pady=pady)

		if self.model_graft_type == "unigraft":
			self.step4_1()

		elif self.model_graft_type == "bigraft":
			self.step4_2()

		else:
			print(f"Unknown setting for graft_type ({self.graft_type})")
			sys.exit(0)		

		self.btn_graft_type = tk.Button(self.master, text="OK")
		self.btn_graft_type.pack(pady=pady)

		self.btn_step4 = tk.Button(self.master, text="Next", command=self.step5)
		self.btn_step4.pack(pady=pady)

		prev_button = tk.Button(self.master, text="Previous", command=self.step3)
		prev_button.pack(pady=pady)

		self.update_windows()

	def step4_1(self):
		# for uni-grafts
		if self.model_format == "LAMMPS":
			# import .data file for graft			

			# file box
			btn_lmp_uni = tk.Button(self.master, text="Select LAMMPS DATA File", command=self.open_file_lmp_uni)
			btn_lmp_uni.pack(pady=pady)

		elif self.model_format == "GROMACS":
			# import .gro and .itp file for graft

			# file box
			btn_gmx_uni_gro = tk.Button(self.master, text="Select GROMACS gro File", command=self.open_file_gmx_uni_gro)
			btn_gmx_uni_gro.pack(pady=pady)

			btn_gmx_uni_itp = tk.Button(self.master, text="Select GROMACS itp File", command=self.open_file_gmx_uni_itp)
			btn_gmx_uni_itp.pack(pady=pady)

	def step4_2(self):
		# for bi-grafts
		if self.model_format == "LAMMPS":
			# import .data file for graft

			# for the first graft
			btn_lmp_bi = tk.Button(self.master, text="Select LAMMPS DATA File for the 1st graft", command=self.open_file_lmp_bi_1st)
			btn_lmp_bi.pack(pady=pady)

			# for the second graft
			btn_lmp_bi = tk.Button(self.master, text="Select LAMMPS DATA File for the 2nd graft", command=self.open_file_lmp_bi_2nd)
			btn_lmp_bi.pack(pady=pady)

		elif self.model_format == "GROMACS":
			# import .gro and .itp file for graft

			# for the first graft
			btn_gmx_bi_gro = tk.Button(self.master, text="Select GROMACS gro File for the 1st graft", command=self.open_file_gmx_bi_gro_1st)
			btn_gmx_bi_gro.pack(pady=pady)

			btn_gmx_bi_itp = tk.Button(self.master, text="Select GROMACS itp File for the 1st graft", command=self.open_file_gmx_bi_itp_1st)
			btn_gmx_bi_itp.pack(pady=pady)

			# for the second graft
			btn_gmx_bi_gro = tk.Button(self.master, text="Select GROMACS gro File for the 2nd graft", command=self.open_file_gmx_bi_gro_2nd)
			btn_gmx_bi_gro.pack(pady=pady)

			btn_gmx_bi_itp = tk.Button(self.master, text="Select GROMACS itp File for the 2nd graft", command=self.open_file_gmx_bi_itp_2nd)
			btn_gmx_bi_itp.pack(pady=pady)

	def step5(self):
		# set substrate shapes
		self.lbl_step5 = tk.Label(self.master, text="Step 5: select substrate type")
		self.lbl_step5.pack(pady=pady)

		self.step5_options = ["slab", "rod", "pore", "sphere"]
		self.substrate_type = tk.StringVar()
		self.substrate_type.set(self.step5_options[0])

		self.select_box = tk.OptionMenu(self.master, self.substrate_type, *self.step5_options)
		self.select_box.pack(pady=pady)

		self.btn_substrate_type = tk.Button(self.master, text="OK", command=self.get_substrate_type)
		self.btn_substrate_type.pack(pady=pady)

		self.btn_step5 = tk.Button(self.master, text="Next", command=self.step6)
		self.btn_step5.pack(pady=pady)

		prev_button = tk.Button(self.master, text="Previous", command=self.step4)
		prev_button.pack(pady=pady)

		self.update_windows()

	def step6(self):

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

		self.btn_step6 = tk.Button(self.master, text="Next", command=self.step7)
		self.btn_step6.pack(pady=pady)

		prev_button = tk.Button(self.master, text="Previous", command=self.step5)
		prev_button.pack(pady=pady)

		self.update_windows()

	def step6_slab(self):

		# set substrate sizes
		self.lbl_step6 = tk.Label(self.master, text="Step 6: set substrate sizes")
		self.lbl_step6.pack(pady=pady)

		lbl_text = tk.Label(self.master, text="Please input the sizes of the substrate (unit in Angstrom)")
		lbl_text.pack(pady=pady)

		# three entry
		lbl_slab_length = tk.Label(self.master, text="Length:")
		self.ent_slab_length = tk.Entry(self.master, width=10)
		self.ent_slab_length.pack(pady=pady)

		lbl_slab_width = tk.Label(self.master, text="Width:")
		self.ent_slab_width = tk.Entry(self.master, width=10)
		self.ent_slab_width.pack(pady=pady)

		lbl_slab_height = tk.Label(self.master, text="Height:")
		self.ent_slab_height = tk.Entry(self.master, width=10)
		self.ent_slab_height.pack(pady=pady)

		self.btn_graft_type = tk.Button(self.master, text="OK", command=self.get_slab_values)
		self.btn_graft_type.pack(pady=pady)

		self.btn_step6 = tk.Button(self.master, text="Next", command=self.step7)
		self.btn_step6.pack(pady=pady)

	def step7(self):
		# set grafting density
		self.lbl_step6 = tk.Label(self.master, text="Step 7: set grafting density")
		self.lbl_step6.pack(pady=pady)

		lbl_text = tk.Label(self.master, text="Please input the grafting density (unit in 1/Angstrom2)")
		lbl_text.pack(pady=pady)

		# one entry
		lbl_grafting_density = tk.Label(self.master, text="Grafting density:")
		self.ent_grafting_density = tk.Entry(self.master, width=10)
		self.ent_grafting_density.pack(pady=pady)

		self.btn_graft_type = tk.Button(self.master, text="OK", command=self.get_grafting_density)
		self.btn_graft_type.pack(pady=pady)

		self.btn_step7 = tk.Button(self.master, text="Next", command=self.step8)
		self.btn_step7.pack(pady=pady)

		prev_button = tk.Button(self.master, text="Previous", command=self.step6)
		prev_button.pack(pady=pady)

		self.update_windows()

	def step8(self):
		# set output file name
		self.lbl_step8 = tk.Label(self.master, text="Step 8: set output file name")
		self.lbl_step8.pack(pady=pady)

		lbl_text = tk.Label(self.master, text="Please input the output file name")
		lbl_text.pack(pady=pady)

		# one entry
		self.ent_fname = tk.Entry(self.master, width=10)
		self.ent_fname.pack(pady=pady)

		self.btn_fname = tk.Button(self.master, text="OK", command=self.get_fname)
		self.btn_fname.pack(pady=pady)

		prev_button = tk.Button(self.master, text="Previous", command=self.step7)
		prev_button.pack(pady=pady)

		self.btn_step8 = tk.Button(self.master, text="Submit", command=self.call_polyGraft)
		self.btn_step8.pack(pady=pady)

		self.update_windows()

	def prev_step(self):
		self.current_step.set(self.current_step.get() - 1)
		step = self.current_step.get()

		if step == 1:
			self.prev_button.config(state=tk.DISABLED)
		else:
			self.prev_button.config(state=tk.NORMAL)

		self.next_button.config(state=tk.NORMAL)

		self.update_windows()

	def call_polyGraft(self):
		print("Generating the structure and topology ...")

	def update_windows(self):
		step = self.current_step.get()

		if step == 2:
			self.lbl_step2.config(text="Step 2: select data format")

		elif step == 3:
			self.lbl_step3.config(text="Step 3: select grafts type")

		elif step == 4:
			self.lbl_step4.config(text="Step 4: import graft data file")

		elif step == 5:
			self.lbl_step5.config(text="Step 5: select substrate type")

		elif step == 6:
			self.lbl_step6.config(text="Step 6: set substrate sizes")

		elif step == 7:
			self.lbl_step7.config(text="Step 7: set grafting density")

		elif step == 8:
			self.lbl_step8.config(text="Step 8: set output file name")

def main():
	root = tk.Tk()
	app = polyGraftGUI(root)
	root.mainloop()

if __name__ == "__main__":
	main()

