# following tutorial at https://www.youtube.com/watch?v=itRLRfuL_PQ

import tkinter as tk
from PIL import Image, ImageTk
import cellanneal
import pandas as pd
from tkinter.filedialog import askopenfile

root = tk.Tk()

canvas = tk.Canvas(root, width=600, height=300)
canvas.grid(columnspan=3, rowspan=3)

# Logo
logo = Image.open('logo.png')
# convert pillow image to tkinter image
logo = ImageTk.PhotoImage(logo)
# place image inside label widget (both lines below are needed)
logo_label = tk.Label(image=logo)
logo_label.image = logo
# place the logo label inside the box using grid method
logo_label.grid(column=1, row=0)

# Instructions
instructions = tk.Label(
                root,
                text="Select bulk data to be deconvolved from your computer.",
                font=("Linux Libertine Mono O", 12))
instructions.grid(columnspan=3, column=0, row=1)

def open_file():
        browse_text.set("loading...")
        file = askopenfile(
                parent=root,
                mode='rb',
                title="Choose bulk data.",
                filetypes=[("csv file", "*.csv")])
        if file:
                print("hello.")
                bulk_df = pd.read_csv(file, index_col=0)
                print(bulk_df.head())
# Browse data
browse_text = tk.StringVar()
browse_btn = tk.Button(
                root,
                textvariable=browse_text,
                command=lambda: open_file(),
                font=("Linux Libertine Mono O", 12),
                highlightbackground="sea green",
                fg="black",
                height=2,
                width=15)
browse_text.set("Browse")
browse_btn.grid(column=1, row=2)

# add whitespace at the bottom of the app
canvas = tk.Canvas(root, width=600, height=250)
canvas.grid(columnspan=3)


root.mainloop()
