import tkinter as tk
from PIL import Image, ImageTk
import cellanneal
import pandas as pd
from tkinter.filedialog import askopenfile, askdirectory
from tkinter import messagebox

class cellgui:

    def __init__(self, root):
        self.root = root
        root.title("cellanneal")

        # place holders for required data
        self.bulk_df = pd.DataFrame()
        self.bulk_df_is_set = 0  # to check later whether data has been selected
        self.celltype_df = pd.DataFrame()
        self.celltype_df_is_set = 0

        # output path
        self.output_folder_path = tk.StringVar()
        self.output_path_is_set = 0

        # parameters for the deconvolution procedure, set to default values
        self.bulk_min = 1e-5  # minimum expression in bulk
        self.bulk_min_default = 1e-5  # default values are requried for reset functionality

        self.bulk_max = 0.01  # maximum expression in bulk
        self.bulk_max_default = 0.01

        self.disp_min = 0.5  # minimum dispersion for highly variable genes
        self.disp_min_default = 0.5

        self.maxiter = 1000  # maximum annealing iterations
        self.maxiter_default = 1000

        self.N_repeat = 10  # deconvolution repeats in case of repeatanneal
        self.N_repeat_default = 10

        # basic layout considerations
        self.canvas = tk.Canvas(root, width=800, height=600)
        self.canvas.grid(columnspan=7, rowspan=12)
        # for easier grid layout changes
        i_i = 1  # start of import section
        p_i = 5  # start of parameter section
        d_i = 11  # start of deconv section

        """ logo """
        logo = Image.open('logo.png')
        # convert pillow image to tkinter image
        logo = ImageTk.PhotoImage(logo)
        # place image inside label widget (both lines below are needed)
        logo_label = tk.Label(image=logo)
        logo_label.image = logo
        # place the logo label inside the box using grid method
        logo_label.grid(column=3, row=0)

        """ main section labels, structure """
        self.sec1_label = tk.Label(root, text="1) Select source data \nand output folder.", font=('Helvetica', 14, 'bold'))
        self.sec2_label = tk.Label(root, text="2) Set parameters. \n[optional]", font=('Helvetica', 14, 'bold'))
        self.sec3_label = tk.Label(root, text="3) Run deconvolution.", font=('Helvetica', 14, 'bold'))
        self.sec1_label.grid(row=i_i, column=0, sticky='w')
        self.sec2_label.grid(row=p_i, column=0, sticky='w')
        self.sec3_label.grid(row=d_i, column=0, sticky='w')

        """ data import section """
        # import bulk data
        self.bulk_folder_path = tk.StringVar()
        # label to indicate that we want bulk data here
        self.bulk_data_label = tk.Label(
                                    root,
                                    text="Select bulk data (*.csv).")
        self.bulk_data_label.grid(row=i_i, column=1, columnspan=2, sticky=tk.W)
        # path entry field
        self.bulk_data_entry = tk.Entry(
                                    root,
                                    textvariable=self.bulk_folder_path)
        self.bulk_data_entry.grid(row=i_i, column=3, columnspan=2, sticky=tk.W+tk.E)
        # file system browse button
        self.bulk_browse_button = tk.Button(
                                        root,
                                        text="Browse file system",
                                        command=lambda: self.import_bulk_data())
        self.bulk_browse_button.grid(row=i_i, column=5, columnspan=2, sticky=tk.W+tk.E)

        # import celltype data
        self.celltype_folder_path = tk.StringVar()
        # label to indicate that we want celltype data here
        self.celltype_data_label = tk.Label(root, text="Select celltype data (*.csv).")
        self.celltype_data_label.grid(row=i_i+1, column=1, columnspan=2, sticky=tk.W)
        # path entry field
        self.celltype_data_entry = tk.Entry(root, textvariable=self.celltype_folder_path)
        self.celltype_data_entry.grid(row=i_i+1, column=3, columnspan=2, sticky=tk.W+tk.E)
        # file system browse button
        self.celltype_browse_button = tk.Button(
                                        root,
                                        text="Browse file system",
                                        command=lambda: self.import_celltype_data())
        self.celltype_browse_button.grid(row=i_i+1, column=5, columnspan=2, sticky=tk.W+tk.E)

        # select output folder
        # label to indicate that we want celltype data here
        self.output_folder_label = tk.Label(root, text="Select folder to store results.")
        self.output_folder_label.grid(row=i_i+2, column=1, columnspan=2, sticky=tk.W)
        # path entry field
        self.output_folder_entry = tk.Entry(root, textvariable=self.output_folder_path)
        self.output_folder_entry.grid(row=i_i+2, column=3, columnspan=2, sticky=tk.W+tk.E)
        # file system browse button
        self.celltype_browse_button = tk.Button(
                                        root,
                                        text="Browse file system",
                                        command=lambda: self.select_output_folder())
        self.celltype_browse_button.grid(row=i_i+2, column=5, columnspan=2, sticky=tk.W+tk.E)

        """ parameter section """
        # for parameter bulk_min
        # title label and current value
        self.bulk_min_label = tk.Label(root, text="bulk_min", font="-weight bold")
        self.bulk_min_label.grid(row=p_i, column=1, columnspan=2, sticky=tk.W+tk.E)
        self.bulk_min_current_label = tk.Label(root, text="current value: {}".format(self.bulk_min), font="-slant italic")
        self.bulk_min_current_label.grid(row=p_i+1, column=1, columnspan=2, sticky=tk.W+tk.E)

        # entry field
        self.bulk_min_entry = tk.Entry(
                                    root,
                                    width=8)
        self.bulk_min_entry.grid(row=p_i+2, column=1, sticky=tk.E)

        # set button
        self.bulk_min_set_button = tk.Button(
                                        root,
                                        text='set',
                                        command=lambda: self.set_bulk_min(),
                                        width=8)
        self.bulk_min_set_button.grid(row=p_i+2, column=2, sticky=tk.W)

        # for parameter bulk_max
        # title label and current value
        self.bulk_max_label = tk.Label(root, text="bulk_max", font="-weight bold")
        self.bulk_max_label.grid(row=p_i, column=3, columnspan=2, sticky=tk.W+tk.E)
        self.bulk_max_current_label = tk.Label(root, text="current value: {}".format(self.bulk_max), font="-slant italic")
        self.bulk_max_current_label.grid(row=p_i+1, column=3, columnspan=2, sticky=tk.W+tk.E)

        # entry field
        self.bulk_max_entry = tk.Entry(
                                    root,
                                    width=8)
        self.bulk_max_entry.grid(row=p_i+2, column=3, sticky=tk.E)

        # set button
        self.bulk_max_set_button = tk.Button(
                                        root,
                                        text='set',
                                        command=lambda: self.set_bulk_max(),
                                        width=8)
        self.bulk_max_set_button.grid(row=p_i+2, column=4, sticky=tk.W)


        # for parameter disp_min
        # title label and current value
        self.disp_min_label = tk.Label(root, text="disp_min", font="-weight bold")
        self.disp_min_label.grid(row=p_i, column=5, columnspan=2, sticky=tk.W+tk.E)
        self.disp_min_current_label = tk.Label(root, text="current value: {}".format(self.disp_min), font="-slant italic")
        self.disp_min_current_label.grid(row=p_i+1, column=5, columnspan=2, sticky=tk.W+tk.E)

        # entry field
        self.disp_min_entry = tk.Entry(
                                    root,
                                    width=8)
        self.disp_min_entry.grid(row=p_i+2, column=5, sticky=tk.E)

        # set button
        self.disp_min_set_button = tk.Button(
                                        root,
                                        text='set',
                                        command=lambda: self.set_disp_min(),
                                        width=8)
        self.disp_min_set_button.grid(row=p_i+2, column=6, sticky=tk.W)

        # for parameter maxiter
        # title label and current value
        self.maxiter_label = tk.Label(root, text="maxiter", font="-weight bold")
        self.maxiter_label.grid(row=p_i+3, column=1, columnspan=2, sticky=tk.W+tk.E)
        self.maxiter_current_label = tk.Label(root, text="current value: {}".format(self.maxiter), font="-slant italic")
        self.maxiter_current_label.grid(row=p_i+4, column=1, columnspan=2, sticky=tk.W+tk.E)

        # entry field
        self.maxiter_entry = tk.Entry(
                                    root,
                                    width=8)
        self.maxiter_entry.grid(row=p_i+5, column=1, sticky=tk.E)

        # set button
        self.maxiter_set_button = tk.Button(
                                        root,
                                        text='set',
                                        command=lambda: self.set_maxiter(),
                                        width=8)
        self.maxiter_set_button.grid(row=p_i+5, column=2, sticky=tk.W)

        # for parameter N_repeat
        # title label and current value
        self.N_repeat_label = tk.Label(root, text="N_repeat", font="-weight bold")
        self.N_repeat_label.grid(row=p_i+3, column=3, columnspan=2, sticky=tk.W+tk.E)
        self.N_repeat_current_label = tk.Label(root, text="current value: {}".format(self.N_repeat), font="-slant italic")
        self.N_repeat_current_label.grid(row=p_i+4, column=3, columnspan=2, sticky=tk.W+tk.E)

        # entry field
        self.N_repeat_entry = tk.Entry(
                                    root,
                                    width=8)
        self.N_repeat_entry.grid(row=p_i+5, column=3, sticky=tk.E)

        # set button
        self.N_repeat_set_button = tk.Button(
                                        root,
                                        text='set',
                                        command=lambda: self.set_N_repeat(),
                                        width=8)
        self.N_repeat_set_button.grid(row=p_i+5, column=4, sticky=tk.W)


        """ run deconvolution section """

    # methods
    def import_bulk_data(self):
        file = askopenfile(
                parent=root,
                mode='rb',
                title="Choose a csv file.",
                filetypes=[("csv file", "*.csv")])
        if file:
            try:
                self.bulk_df = pd.read_csv(file, index_col=0)
                self.bulk_folder_path.set(file.name)
                print(file.name)
                self.bulk_df_is_set = 1
                print(self.bulk_df.head())
            except:
                messagebox.showerror("Import error", """Your bulk data file could not be imported. Please check the documentation for format requirements and look at the example bulk data file.""")

    def import_celltype_data(self):
        file = askopenfile(
                parent=root,
                mode='rb',
                title="Choose a csv file.",
                filetypes=[("csv file", "*.csv")])
        if file:
            try:
                self.celltype_df = pd.read_csv(file, index_col=0)
                self.celltype_folder_path.set(file.name)
                self.celltype_df_is_set = 1
                print(self.celltype_df.head())
            except:
                messagebox.showerror("Import error", """Your celltype data file could not be imported. Please check the documentation for format requirements and look at the example celltype data file.""")

    def select_output_folder(self):
        folder = askdirectory(
                parent=root,
                title="Choose a folder.")
        if folder:
            print(folder)
            self.output_folder_path.set(folder)
            self.output_path_is_set = 1
            print('1', self.output_folder_path.get())



    def set_bulk_min(self):
        # get input and check validity
        input = self.bulk_min_entry.get()
        # can it be interpreted as float?
        try:
            new_val = float(input)
        except ValueError:
            messagebox.showerror("Input error", """Please provdide a numerical value between 0 (inclusive) and 1 (exclusive). Examples: "0", "0.01", "2e-5".""")
            return 0
        # check if value is between 0 and 1 and smaller than bulk_max
        if (0 <= new_val < 1) and (new_val < self.bulk_max):
            # update bulk_min
            self.bulk_min = new_val
            # update displa of current bulk_min
            self.bulk_min_current_label['text'] = "current value: {}".format(self.bulk_min)
        elif not (0 <= new_val < 1):
            messagebox.showerror("Input error", """Please provdide a numerical value between 0 (inclusive) and 1 (exclusive). Examples: "0", "0.01", "2e-5".""")
        elif not (new_val < self.bulk_max):
            messagebox.showerror("Input error", """Your value must be smaller than bulk_max.""")

    def set_bulk_max(self):
        # get input and check validity
        input = self.bulk_max_entry.get()
        # can it be interpreted as float?
        try:
            new_val = float(input)
        except ValueError:
            messagebox.showerror("Input error", """Please provdide a numerical value between 0 (exclusive) and 1 (inclusive). Examples: "0.01", "2e-5", "1".""")
            return 0
        # check if value is between 0 and 1 and if it is larger than bulk_min
        if (0 < new_val <= 1) and (new_val > self.bulk_min):
            # update bulk_min
            self.bulk_max = new_val
            # update display of current bulk_min
            self.bulk_max_current_label['text'] = "current value: {}".format(self.bulk_max)
        elif not (0 < new_val <= 1):
            messagebox.showerror("Input error", """Please provdide a numerical value between 0 (exclusive) and 1 (inclusive). Examples: "0.01", "2e-5", "1".""")
        elif not (new_val > self.bulk_min):
            messagebox.showerror("Input error", """Your value must be larger than bulk_min.""")

    def set_disp_min(self):
        # get input and check validity
        input = self.disp_min_entry.get()
        # can it be interpreted as float?
        try:
            new_val = float(input)
        except ValueError:
            messagebox.showerror("Input error", """Please provdide a positive numerical value. Examples: "0.5", "5e-1", "1".""")
            return 0
        # check if value is >0
        if new_val >= 0:
            # update disp_min
            self.disp_min = new_val
            # update display of current bulk_min
            self.disp_min_current_label['text'] = "current value: {}".format(self.disp_min)
        else:
            messagebox.showerror("Input error", """Please provdide a positive numerical value. Examples: "0.5", "5e-1", "1".""")

    def set_maxiter(self):
        # get input and check validity
        input = self.maxiter_entry.get()
        print(input)
        # can it be interpreted as float?
        try:
            new_val = int(input)
        except ValueError:
            messagebox.showerror("Input error", """Please provdide a positive integer. Examples: "50", "587", "1000".""")
            return 0
        # check if value is >0
        if new_val >= 1:
            # update disp_min
            self.maxiter = new_val
            # update display of current bulk_min
            self.maxiter_current_label['text'] = "current value: {}".format(self.maxiter)
        else:
            messagebox.showerror("Input error", """Please provdide a positive integer. Examples: "50", "587", "1000".""")

    def set_N_repeat(self):
        # get input and check validity
        input = self.N_repeat_entry.get()
        print(input)
        # can it be interpreted as float?
        try:
            new_val = int(input)
        except ValueError:
            messagebox.showerror("Input error", """Please provdide a positive integer. Examples: "5", "10", "21".""")
            return 0
        # check if value is >0
        if new_val >= 1:
            # update disp_min
            self.N_repeat = new_val
            # update display of current bulk_min
            self.N_repeat_current_label['text'] = "current value: {}".format(self.N_repeat)
        else:
            messagebox.showerror("Input error", """Please provdide a positive integer. Examples: "5", "10", "21".""")


root = tk.Tk()
my_gui = cellgui(root)
root.mainloop()
