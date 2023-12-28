import tkinter as tk
from tkinter import filedialog

class SixStepGUI:
    def __init__(self, master):
        self.master = master
        self.master.title("Six Step GUI")
        self.master.geometry("800x600")  # Setting window size

        self.current_step = tk.IntVar()
        self.current_step.set(1)

        self.create_widgets()

    def create_widgets(self):
        self.label = tk.Label(self.master, text="Step 1")
        self.label.pack(pady=10)

        self.file_path = tk.StringVar()
        self.entry = tk.Entry(self.master, textvariable=self.file_path)
        self.entry.pack(pady=5)

        self.select_file_button = tk.Button(self.master, text="Select File", command=self.select_file)
        self.select_file_button.pack()

        self.options = ["Choice 1", "Choice 2"]
        self.selection = tk.StringVar()
        self.selection.set(self.options[0])

        self.label_choice = tk.Label(self.master, text="Please select choice:")
        self.label_choice.pack()

        self.select_box = tk.OptionMenu(self.master, self.selection, *self.options)
        self.select_box.pack(pady=5)

        self.next_button = tk.Button(self.master, text="Next", command=self.next_step)
        self.next_button.pack(pady=5)

        self.prev_button = tk.Button(self.master, text="Previous", command=self.prev_step, state=tk.DISABLED)
        self.prev_button.pack(pady=5)

        self.run_button = tk.Button(self.master, text="Run", command=self.run_process, state=tk.DISABLED)
        self.run_button.pack(pady=5)

    def update_widgets(self):
        step = self.current_step.get()
        self.label.config(text=f"Step {step}")

        if step == 6:
            self.run_button.config(state=tk.NORMAL)
        else:
            self.run_button.config(state=tk.DISABLED)

    def select_file(self):
        file_path = filedialog.askopenfilename()
        self.file_path.set(file_path)

    def next_step(self):
        user_input = self.entry.get()
        selected_option = self.selection.get()
        print(f"Step {self.current_step.get()} - User input: {user_input}, Selected: {selected_option}")

        self.current_step.set(self.current_step.get() + 1)
        step = self.current_step.get()

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

def main():
    root = tk.Tk()
    app = SixStepGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
