import tkinter as tk
from PIL import Image, ImageTk
import cellanneal
import pandas as pd
from tkinter.filedialog import askopenfile


class cellgui:

    def __init__(self, root):
        self.root = root
        root.title("cellanneal")

        # place holders for required data
        self.bulk_df = pd.DataFrame()
        self.celltype_df = pd.DataFrame()

        # parameters for the deconvolution procedure, set to default values
        self.bulk_min = 1e-5  # minimum expression in bulk
        self.bulk_max = 0.01  # maximum expression in bulk
        self.disp_min = 0.5  # minimum dispersion for highly variable genes
        self.maxiter = 1000  # maximum annealing iterations
        self.repeat = 10  # deconvolution repeats in case of repeatanneal

        # basic layout considerations
        self.canvas = tk.Canvas(root, width=800, height=600)
        self.canvas.grid(columnspan=7, rowspan=8)

        """ logo """
        logo = Image.open('logo.png')
        # convert pillow image to tkinter image
        logo = ImageTk.PhotoImage(logo)
        # place image inside label widget (both lines below are needed)
        logo_label = tk.Label(image=logo)
        logo_label.image = logo
        # place the logo label inside the box using grid method
        logo_label.grid(column=1, row=0)

        """ main section labels, structure """
        self.sec1_label = tk.Label(root, text="1) Import data.", font=('Helvetica', 14, 'bold'))
        self.sec2_label = tk.Label(root, text="2) Set parameters. \n[optional]", font=('Helvetica', 14, 'bold'))
        self.sec3_label = tk.Label(root, text="3) Run deconvolution.", font=('Helvetica', 14, 'bold'))
        self.sec1_label.grid(row=1, column=0, sticky='w')
        self.sec2_label.grid(row=4, column=0, sticky='w')
        self.sec3_label.grid(row=7, column=0, sticky='w')

        """ data import section """
        # import bulk data
        self.bulk_folder_path = tk.StringVar()
        # label to indicate that we want bulk data here
        self.bulk_data_label = tk.Label(
                                    root,
                                    text="Select bulk data (*.csv).")
        self.bulk_data_label.grid(row=1, column=1, columnspan=2, sticky=tk.W+tk.E)
        # path entry field
        self.bulk_data_entry = tk.Entry(
                                    root,
                                    textvariable=self.bulk_folder_path)
        self.bulk_data_entry.grid(row=1, column=3, columnspan=2, sticky=tk.W+tk.E)
        # file system browse button
        self.bulk_browse_button = tk.Button(
                                        root,
                                        text="Browse file system",
                                        command=lambda: self.import_bulk_data())
        self.bulk_browse_button.grid(row=1, column=5, columnspan=2, sticky=tk.W+tk.E)

        # import celltype data
        self.celltype_folder_path = tk.StringVar()
        # label to indicate that we want celltype data here
        self.celltype_data_label = tk.Label(root, text="Select celltype data (*.csv).")
        self.celltype_data_label.grid(row=2, column=1, columnspan=2, sticky=tk.W+tk.E)
        # path entry field
        self.celltype_data_entry = tk.Entry(root, textvariable=self.celltype_folder_path)
        self.celltype_data_entry.grid(row=2, column=3, columnspan=2, sticky=tk.W+tk.E)
        # file system browse button
        self.celltype_browse_button = tk.Button(
                                        root,
                                        text="Browse file system",
                                        command=lambda: self.import_celltype_data())
        self.celltype_browse_button.grid(row=2, column=5, columnspan=2, sticky=tk.W+tk.E)

        """ parameter section """


        """ run deconvolution section """

    # methods
    def import_bulk_data(self):
        file = askopenfile(
                parent=root,
                mode='rb',
                title="Choose a csv file.",
                filetypes=[("csv file", "*.csv")])
        if file:
            self.bulk_df = pd.read_csv(file, index_col=0)
            self.bulk_folder_path.set(file.name)
            print(self.bulk_df.head())

    def import_celltype_data(self):
        file = askopenfile(
                parent=root,
                mode='rb',
                title="Choose a csv file.",
                filetypes=[("csv file", "*.csv")])
        if file:
            self.celltype_df = pd.read_csv(file, index_col=0)
            self.celltype_folder_path.set(file.name)
            print(self.celltype_df.head())

    def validate_float(self, key_input):

        try:
            self.entered_number = int(new_text)
            return True
        except ValueError:
            return False


root = tk.Tk()
my_gui = cellgui(root)
root.mainloop()
