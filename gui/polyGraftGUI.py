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
			 2: "select format",
			 3: "grafts type",
			 4: "substrate type",
			 5: "grafting density",
			 6: "output"}

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
	
	def step4(self):
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

	def step4_1(self):
		# for uni-grafts
		if self.model_format == "LAMMPS":
			# import .data file for graft
			self.lbl_step4_1 = tk.Label(self.master, text="Step 4: import graft data file")
			self.lbl_step4_1.pack(pady=pady)

			# file box
			btn_lmp_uni = tk.Button(self.master, text="Select LAMMPS DATA File", command=self.open_file_lmp_uni)
			btn_lmp_uni.pack(pady=pady)

		elif self.model_format == "GROMACS":
			# import .gro and .itp file for graft
			self.lbl_step4_2 = tk.Label(self.master, text="Step 4: import graft gro/itp file")
			self.lbl_step4_2.pack(pady=pady)

			# file box
			btn_gmx_uni_gro = tk.Button(self.master, text="Select GROMACS gro File", command=self.open_file_gmx_uni_gro)
			btn_gmx_uni_gro.pack(pady=pady)

			btn_gmx_uni_itp = tk.Button(self.master, text="Select GROMACS itp File", command=self.open_file_gmx_uni_itp)
			btn_gmx_uni_itp.pack(pady=pady)

	def step4_2(self):
		# for bi-grafts
		if self.model_format == "LAMMPS":
			# import .data file for graft
			self.lbl_step4_1 = tk.Label(self.master, text="Step 4: import graft data file")
			self.lbl_step4_1.pack(pady=pady)

			# for the first graft
			btn_lmp_bi = tk.Button(self.master, text="Select LAMMPS DATA File for the 1st graft", command=self.open_file_lmp_bi_1st)
			btn_lmp_bi.pack(pady=pady)

			# for the second graft
			btn_lmp_bi = tk.Button(self.master, text="Select LAMMPS DATA File for the 2nd graft", command=self.open_file_lmp_bi_2nd)
			btn_lmp_bi.pack(pady=pady)

		elif self.model_format == "GROMACS":
			# import .gro and .itp file for graft
			self.lbl_step4_2 = tk.Label(self.master, text="Step 4: import graft gro/itp file")
			self.lbl_step4_2.pack(pady=pady)

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
		self.select_box.pack(pady=5)

		self.btn_substrate_type = tk.Button(self.master, text="OK", command=self.get_substrate_type)
		self.btn_substrate_type.pack(pady=pady)

		self.btn_step5 = tk.Button(self.master, text="Next", command=self.step6)
		self.btn_step5.pack(pady=pady)

	def step6(self):
		pass
		
	def create_widgets(self):
		# step1: resolution settings
		self.label = tk.Label(self.master, text="Step 1: select resolution")
		self.label.pack(pady=10)

		self.options = ["Atomistic", "Coarse grained"]
		self.selection = tk.StringVar()
		self.selection.set(self.options[0])
		self.label_choice = tk.Label(self.master, text="Please select resolution:")
		self.label_choice.pack()

		self.select_box = tk.OptionMenu(self.master, self.selection, *self.options)
		self.select_box.pack(pady=5)

		self.next_button = tk.Button(self.master, text="Next", command=self.next_step)
		self.next_button.pack(pady=5)

		self.prev_button = tk.Button(self.master, text="Previous", command=self.prev_step, state=tk.DISABLED)
		self.prev_button.pack(pady=5)

		self.run_button = tk.Button(self.master, text="Run", command=self.run_process, state=tk.DISABLED)
		self.run_button.pack(pady=5)

	def next_step(self):
		self.current_step.set(self.current_step.get() + 1)
		step = self.current_step.get()

		if step == 2:

			self.label = tk.Label(self.master, text="Step 2: select format")
			self.label.pack(pady=10)

			self.options = ["GROMACS", "LAMMPS"]
			self.selection = tk.StringVar()
			self.selection.set(self.options[0])
			self.label_choice = tk.Label(self.master, text="Please select format:")
			self.label_choice.pack()

			self.select_box = tk.OptionMenu(self.master, self.selection, *self.options)
			self.select_box.pack(pady=5)

		elif step == 3:

			self.label = tk.Label(self.master, text="Step 3: select grafts type")
			self.label.pack(pady=10)

			self.options = ["unigraft", "bigraft"]
			self.selection = tk.StringVar()
			self.selection.set(self.options[0])
			self.label_choice = tk.Label(self.master, text="Please select grafts type:")
			self.label_choice.pack()

			self.select_box = tk.OptionMenu(self.master, self.selection, *self.options)
			self.select_box.pack(pady=5)

		elif step == 4:

			self.label = tk.Label(self.master, text="Step 4: select substrate type")
			self.label.pack(pady=10)

			self.options = ["slab", "rod", "pore", "sphere"]
			self.selection = tk.StringVar()
			self.selection.set(self.options[0])
			self.label_choice = tk.Label(self.master, text="Please select substrate type:")
			self.label_choice.pack()

			self.select_box = tk.OptionMenu(self.master, self.selection, *self.options)
			self.select_box.pack(pady=5)

		elif step == 5:
			self.label = tk.Label(self.master, text="Step 5: set substrate size")
			self.label.pack(pady=10)

			self.file_path = tk.StringVar()
			self.entry = tk.Entry(self.master, textvariable=self.file_path)
			self.entry.pack(pady=5)

			user_input = self.entry.get()
			selected_option = self.selection.get()
			self.label_choice = tk.Label(self.master, text="Please set substrate size:")
			self.label_choice.pack()

		elif step == 6:
			self.label = tk.Label(self.master, text="Step 6: set output file name")
			self.label.pack(pady=10)

			user_input = self.entry.get()
			selected_option = self.selection.get()
			self.label_choice = tk.Label(self.master, text="Please set output file name:")
			self.label_choice.pack()

			self.select_box = tk.OptionMenu(self.master, self.selection, *self.options)
			self.select_box.pack(pady=5)

		if step == 6:
			self.next_button.config(state=tk.DISABLED)
		else:
			self.next_button.config(state=tk.NORMAL)

		self.prev_button.config(state=tk.NORMAL)

		self.update_widgets()

	def prev_step(self):
		self.current_step.set(self.current_step.get() - 1)
		step = self.current_step.get()

		if step == 1:
			self.prev_button.config(state=tk.DISABLED)
		else:
			self.prev_button.config(state=tk.NORMAL)

		self.next_button.config(state=tk.NORMAL)

		self.update_widgets()

	def run_process(self):
		print("Executing the process...")

	def update_widgets(self):
		step = self.current_step.get()
		self.label.config(text=f"Step {step}: {step_dict[step]}")

		if step == 6:
			self.run_button.config(state=tk.NORMAL)
		else:
			self.run_button.config(state=tk.DISABLED)


def main():
	root = tk.Tk()
	app = polyGraftGUI(root)
	root.mainloop()

if __name__ == "__main__":
	main()

