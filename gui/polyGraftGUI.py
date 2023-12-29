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
		self.master.geometry("800x600")  # Setting window size

		self.current_step = tk.IntVar()
		self.current_step.set(1)

		self.create_widgets()

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

