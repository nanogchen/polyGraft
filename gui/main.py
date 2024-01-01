import tkinter as tk
from tkinter import filedialog

def next_step(step, next_step_widget):
    # Hide current step
    step.pack_forget()
    # Show next step
    next_step_widget.pack()

def back_step(step, prev_step_widget):
    # Hide current step
    step.pack_forget()
    # Show previous step
    prev_step_widget.pack()

def browse_file():
    filename = filedialog.askopenfilename(initialdir=".", title="Select a File")
    file_entry.delete(0, tk.END)
    file_entry.insert(0, filename)

root = tk.Tk()
root.title("polyGraft GUI")
root.geometry("400x600")  # Setting window size

# Step 1
step1 = tk.Frame(root)
label1 = tk.Label(step1, text="Step 1: Choose an option", width=40)
label1.pack(padx=20, pady=20)
choices = ["Option 1", "Option 2", "Option 3"]
choice_var = tk.StringVar()
choice_var.set(choices[0])
option_menu = tk.OptionMenu(step1, choice_var, *choices)
option_menu.pack()
button1_next = tk.Button(step1, text="Next", command=lambda: next_step(step1, step2))
button1_next.pack(anchor='e')

# Step 2
step2 = tk.Frame(root)
label2 = tk.Label(step2, text="Step 2: Select a file")
label2.pack(padx=20, pady=20)
file_entry = tk.Entry(step2, width=40)
file_entry.pack()
browse_button = tk.Button(step2, text="Browse", command=browse_file)
browse_button.pack()
button2_back = tk.Button(step2, text="Back", command=lambda: back_step(step2, step1))
button2_back.pack(side=tk.LEFT)
button2_next = tk.Button(step2, text="Next", command=lambda: next_step(step2, step3))
button2_next.pack(side=tk.RIGHT)

# Step 3
step3 = tk.Frame(root)
label3 = tk.Label(step3, text="Step 3: Enter some input")
label3.pack(padx=20, pady=20)
input_label = tk.Label(step3, text="Enter text:")
input_label.pack(anchor='w')
input_entry = tk.Entry(step3, width=40)
input_entry.pack()
button3_back = tk.Button(step3, text="Back", command=lambda: back_step(step3, step2))
button3_back.pack(anchor='w')

button3_finish = tk.Button(step3, text="Submit", command=root.destroy)
button3_finish.pack(padx=20, pady=20)

# Show initial step
step1.pack()

root.mainloop()
