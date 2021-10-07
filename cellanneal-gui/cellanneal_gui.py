import tkinter as tk
from PIL import Image, ImageTk
from pandas import DataFrame, read_excel, read_csv
from tkinter.filedialog import askopenfile, askdirectory
from tkinter import messagebox
from pathlib import Path
from cellanneal import cellanneal_pipe
import openpyxl  # for xlsx import
import xlrd  # for xls import
import sys
from os import path
import time


def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = path.abspath(".")

    return path.join(base_path, relative_path)


class cellgui:

    def __init__(self, root):
        self.root = root
        root.title("cellanneal")

        # place holders for required data
        self.bulk_df = DataFrame()
        self.bulk_data_path = tk.StringVar()
        self.bulk_df_is_set = 0  # check later whether data has been selected
        self.celltype_df = DataFrame()
        self.celltype_data_path = tk.StringVar()
        self.celltype_df_is_set = 0

        # output path
        self.output_path = tk.StringVar()
        self.output_path_is_set = 0

        # parameters for the deconvolution procedure, set to default values
        self.bulk_min = 1e-5  # minimum expression in bulk
        self._bulk_min_default = 1e-5

        self.bulk_max = 0.01  # maximum expression in bulk
        self._bulk_max_default = 0.01

        self.disp_min = 0.5  # minimum dispersion for highly variable genes
        self._disp_min_default = 0.5

        self.maxiter = 1000  # maximum annealing iterations
        self._maxiter_default = 1000

        # basic layout considerations
        self.canvas = tk.Canvas(root, width=700, height=500)
        self.canvas.grid(columnspan=6, rowspan=15)
        # for easier grid layout changes
        i_i = 4  # start of import section
        p_i = 9  # start of parameter section
        d_i = 14  # start of deconv section

        """ logo and welcome section """
        self.logo = Image.open(resource_path('logo_orange.png'))
        # convert pillow image to tkinter image
        self.logo = ImageTk.PhotoImage(self.logo)
        # place image inside label widget (both lines below are needed)
        self.logo_label = tk.Label(image=self.logo)
        self.logo_label.image = self.logo
        # place the self.logo label inside the box using grid method
        self.logo_label.grid(column=1, columnspan=1, sticky=tk.W, row=0,
                             padx=5, pady=5, rowspan=2)

        self.welcome_label = tk.Label(
            text='Welcome to cellanneal, the user-friendly bulk deconvolution software.',
            wraplength=300, font="-weight bold")
        self.welcome_label.grid(column=2, columnspan=4, sticky=tk.W+tk.E, row=0,
                                pady=10)

        # a button that allows to show a small documentation
        self.instruct_button = tk.Button(root,
                                         text="Instructions",
                                         command=lambda: self.show_instructions(),
                                         highlightbackground='#f47a60',
                                         highlightthickness=0.1,)
        self.instruct_button.grid(row=1, column=3,
                                  columnspan=2, sticky=tk.W+tk.E,
                                  padx=5, pady=1,
                                  ipadx=10, ipady=10)
        # for buttons belows
        self.ca_button = Image.open(resource_path('cellanneal_button.png'))
        self.ca_button = ImageTk.PhotoImage(self.ca_button)
        self.ca_label = tk.Label(image=self.ca_button)
        self.ca_label.image = self.ca_button

        # spacers for tidier look
        spacer1 = tk.Label(root, text="")
        spacer1.grid(row=i_i-2, column=0, padx=10, pady=2)

        spacer2 = tk.Label(root, text="")
        spacer2.grid(row=p_i-2, column=6, padx=10, pady=2)


        """ data import section """
        # import bulk data
        # title
        self.param_label = tk.Label(root, text="deconvolution data", font="-weight bold")
        self.param_label.grid(row=i_i-1, column=1, columnspan=1, sticky=tk.W)
        # label to indicate that we want bulk data here
        self.bulk_data_label = tk.Label(
                                    root,
                                    text="Select mixture data (.csv, .txt, .xlsx).")
        self.bulk_data_label.grid(row=i_i, column=1, columnspan=1, sticky=tk.W)
        # path entry field
        self.bulk_data_entry = tk.Entry(
                                    root,
                                    textvariable=self.bulk_data_path,
                                    state='readonly')
        self.bulk_data_entry.grid(row=i_i, column=2, columnspan=2, sticky=tk.W+tk.E)
        # file system browse button
        self.bulk_browse_button = tk.Button(
                                        root,
                                        text="Browse file system",
                                        command=lambda: self.import_bulk_data())
        self.bulk_browse_button.grid(row=i_i, column=4, columnspan=2, sticky=tk.W+tk.E)

        # import celltype data
        # label to indicate that we want celltype data here
        self.celltype_data_label = tk.Label(root, text="Select signature data (.csv, .txt, .xlsx).")
        self.celltype_data_label.grid(row=i_i+1, column=1, columnspan=2, sticky=tk.W)
        # path entry field
        self.celltype_data_entry = tk.Entry(
                                        root,
                                        textvariable=self.celltype_data_path,
                                        state='readonly')
        self.celltype_data_entry.grid(row=i_i+1, column=2, columnspan=2, sticky=tk.W+tk.E)
        # file system browse button
        self.celltype_browse_button = tk.Button(
                                        root,
                                        text="Browse file system",
                                        command=lambda: self.import_celltype_data())
        self.celltype_browse_button.grid(row=i_i+1, column=4, columnspan=2, sticky=tk.W+tk.E)

        # select output folder
        # label to indicate that we want celltype data here
        self.output_folder_label = tk.Label(root, text="Select folder to store results.")
        self.output_folder_label.grid(row=i_i+2, column=1, columnspan=2, sticky=tk.W)
        # path entry field
        self.output_folder_entry = tk.Entry(root,
                                            textvariable=self.output_path,
                                            state='readonly')
        self.output_folder_entry.grid(row=i_i+2, column=2, columnspan=2, sticky=tk.W+tk.E)
        # file system browse button
        self.celltype_browse_button = tk.Button(
                                        root,
                                        text="Browse file system",
                                        command=lambda: self.select_output_folder())
        self.celltype_browse_button.grid(row=i_i+2, column=4, columnspan=2, sticky=tk.W+tk.E)

        """ parameter section """
        # title
        self.param_label = tk.Label(root, text="deconvolution parameters", font="-weight bold")
        self.param_label.grid(row=p_i-1, column=1, columnspan=3, sticky=tk.W)
        # list parameter settings with set of labels
        self.bulk_min_param_label = tk.Label(
            root,
            text="minimum expression in mixture")
        self.bulk_min_param_label.grid(row=p_i, column=1, sticky=tk.W)

        self.bulk_max_param_label = tk.Label(
            root,
            text="maximum expression in mixture")
        self.bulk_max_param_label.grid(row=p_i+1, column=1, sticky=tk.W)

        self.disp_min_param_label = tk.Label(
            root,
            text="minimum scaled dispersion")
        self.disp_min_param_label.grid(row=p_i+2, column=1, sticky=tk.W)

        self.maxiter_param_label = tk.Label(
            root,
            text="maximum number of iterations")
        self.maxiter_param_label.grid(row=p_i+3, column=1, sticky=tk.W)

        # corresponding set of values
        self.bulk_min_value_label = tk.Label(
            root,
            text="{}".format(self.bulk_min))
        self.bulk_min_value_label.grid(row=p_i, column=2, sticky=tk.W)

        self.bulk_max_value_label = tk.Label(
            root,
            text="{}".format(self.bulk_max))
        self.bulk_max_value_label.grid(row=p_i+1, column=2, sticky=tk.W)

        self.disp_min_value_label = tk.Label(
            root,
            text="{}".format(self.disp_min))
        self.disp_min_value_label.grid(row=p_i+2, column=2, sticky=tk.W)

        self.maxiter_value_label = tk.Label(
            root,
            text="{}".format(self.maxiter))
        self.maxiter_value_label.grid(row=p_i+3, column=2, sticky=tk.W)

        # button for editing parameters which spawns new window
        self.parameter_change_button = tk.Button(root,
                                                 text='Change\nparameters',
                                                 command=lambda: self.open_param_window(),
                                                 highlightbackground='#f47a60',
                                                 highlightthickness=0.1,)
        self.parameter_change_button.grid(row=p_i,
                                          rowspan=2,
                                          column=3,
                                          columnspan=2,
                                          padx=0, pady=5,
                                          ipadx=17, ipady=5)

        # button for resetting all parameters to default
        self.default_button = tk.Button(root,
                                        text='Reset to\ndefault values',
                                        command=lambda: self.reset_default_params(),
                                        highlightbackground='#f47a60',
                                        highlightthickness=0.1,)
        self.default_button.grid(row=p_i+2,
                                 rowspan=2,
                                 column=3,
                                 columnspan=2,
                                 padx=0, pady=5,
                                 ipadx=8, ipady=5)




        """ run deconvolution section """
        # make button for cellanneal
        self.cellanneal_button = tk.Button(
                                        root,
                                        text='cellanneal',
                                        font="-weight bold ",
                                        image=self.ca_button,
                                        highlightbackground='#f47a60',
                                        highlightthickness=0.1,
                                        command=lambda: self.cellanneal(),
                                        width=20)
        self.cellanneal_button.grid(
                                row=d_i,
                                column=1,
                                columnspan=5,
                                sticky=tk.W+tk.E,
                                padx=10, pady=50)

    # methods
    def show_instructions(self):
        readme_window = tk.Toplevel(root)
        readme_window.title("Instructions")
        readme_window.canvas = tk.Canvas(root, width=700, height=600)
        # write text into this window
        wl = 700
        self.readme_text0 = tk.Label(readme_window, wraplength=wl,  justify=tk.LEFT, text=""" """)
        self.readme_text0.grid(row=13, column=0, padx=10)

        self.readme_text1 = tk.Label(readme_window, wraplength=wl, text="Welcome to cellanneal.", font="-weight bold")
        self.readme_text1.grid(row=1, column=0, padx=10)

        self.readme_text2 = tk.Label(readme_window, wraplength=wl, text="For the full documentation and example files please see https://github.com/LiBuchauer/cellanneal.")
        self.readme_text2.grid(row=2, column=0, padx=10)

        self.readme_text3 = tk.Label(readme_window, wraplength=wl, text="\n\n1) Deconvolution data.\n", font="-weight bold")
        self.readme_text3.grid(row=3, column=0, padx=10)

        self.readme_text4 = tk.Label(readme_window, wraplength=wl, justify=tk.LEFT, text="cellanneal accepts text files (*.csv and *.txt) as well as excel files (*.xlsx and *.xls) as inputs for both mixture and signature data provided that they are formatted as specified in the examples. \nSpecifically, gene names need to appear in the first column for both mixture and signature data files, and sample names (for mixture data file) or cell type names (for signature data file) need to appear in the first row.")
        self.readme_text4.grid(row=4, column=0, padx=10)

        self.readme_text5 = tk.Label(readme_window, wraplength=wl, text="\n\n2) Deconvolution parameters.\n", font="-weight bold")
        self.readme_text5.grid(row=5, column=0, padx=10)

        self.readme_text6 = tk.Label(readme_window, wraplength=wl, justify=tk.LEFT, text="The first three parameters govern the set of genes underlying the deconvolution process for each sample; the fourth parameter (iteration number) specifies for how long to run the optimisation process. \nMinimum expression in mixture - minimum required expression level in the mixture sample (where total expression is normalised to sum up to 1) for a gene to be considered, default=1e-5. Allowed values are in the range [0, 1) but must be smaller than the maximum allowed expression. This parameter allows to exclude lowly expressed and potentially noisy genes. \nMaximum expression in mixture - maximum allowed expression level in the mixture sample (where total expression is normalised to sum up to 1) for a gene to be considered, default=0.01. Allowed values are in the range (0, 1] but must be larger than the minimum allowed expression. This parameter allows to exclude very highly expressed, potential contaminant, genes. \nMinimum scaled dispersion - minimum scaled dispersion (variance/mean) over cell types for a gene to be considered, default=0.5. The value indicates the number of standard deviations which the dispersion of a specific gene lies above or below the mean when compared to genes of similar expression. All numerical values are allowed, but reasonable values for most cases lie between 0 and 3 as this parameter is used to select genes which vary strongly across cell types in the signature file. \nMaximum number of iterations - the maximum number of iterations through the logical chain of the underlying optimisation algorithm, scipy’s dual annealing (scipy.optimize.dual_annealing). Default=1000, after which typical problems have converged. Problems with a very high number of celltypes may require a higher number of iterations.")
        self.readme_text6.grid(row=6, column=0, padx=10)

        self.readme_text11 = tk.Label(readme_window, wraplength=wl, text="\n\n3) Running deconvolution with cellanneal.\n", font="-weight bold")
        self.readme_text11.grid(row=11, column=0, padx=10)

        self.readme_text12 = tk.Label(readme_window, wraplength=wl,  justify=tk.LEFT, text="""Click on "run cellanneal" and watch the progress in the accompanying console. The main window freezes while deconvolution is in progress. Result tables and figures are stored to the user-specified output folder.""")
        self.readme_text12.grid(row=12, column=0, padx=10)

        self.readme_text13 = tk.Label(readme_window, wraplength=wl,  justify=tk.LEFT, text=""" """)
        self.readme_text13.grid(row=13, column=0, padx=10)

    def import_bulk_data(self):
        file = askopenfile(
                parent=root,
                mode='rb',
                title="Choose a mixture data file.",
                filetypes=[("tabular data files", ".csv .txt .xlsx .xls")])
        if file:
            self.bulk_data_path.set(file.name)
            self.bulk_df_is_set = 1

    def import_celltype_data(self):
        file = askopenfile(
                parent=root,
                mode='rb',
                title="Choose a signature data file.",
                filetypes=[("tabular data files", ".csv .txt .xlsx")])
        if file:
            self.celltype_data_path.set(file.name)
            self.celltype_df_is_set = 1

    def select_output_folder(self):
        folder = askdirectory(
                parent=root,
                title="Choose a folder.")
        if folder:
            self.output_path.set(folder)
            self.output_path_is_set = 1

    def open_param_window(self):
        par_window = tk.Toplevel(root)
        par_window.title("Set parameters")
        par_window.canvas = tk.Canvas(root, width=150, height=400)
        # grab attention for this window until it is closed
        par_window.grab_set()

        # empty first column
        spacer1 = tk.Label(par_window, text="")
        spacer1.grid(row=0, column=0, padx=5, pady=0)

        spacer2 = tk.Label(par_window, text="")
        spacer2.grid(row=17, column=3, padx=5, pady=0)

        # for parameter bulk_min
        # title label and current value
        self.bulk_min_label = tk.Label(par_window, text="minimum expression in mixture", font="-weight bold -size 13")
        self.bulk_min_label.grid(row=1, column=1, columnspan=2, sticky=tk.W+tk.E)

        # entry field
        self.bulk_min_entry = tk.Entry(
                                    par_window,
                                    width=8)
        self.bulk_min_entry.grid(row=3, column=1, sticky=tk.E)
        self.bulk_min_entry.insert(tk.END, self.bulk_min)
        # set button
        self.bulk_min_set_button = tk.Button(
                                        par_window,
                                        text='set',
                                        command=lambda: self.set_bulk_min(),
                                        width=8)
        self.bulk_min_set_button.grid(row=3, column=2, sticky=tk.W)

        # for parameter bulk_max
        # title label and current value
        self.bulk_max_label = tk.Label(par_window, text="maximum expression in mixture", font="-weight bold")
        self.bulk_max_label.grid(row=4, column=1, columnspan=2, sticky=tk.W+tk.E)

        # entry field
        self.bulk_max_entry = tk.Entry(
                                    par_window,
                                    width=8)
        self.bulk_max_entry.grid(row=6, column=1, sticky=tk.E)
        self.bulk_max_entry.insert(tk.END, self.bulk_max)

        # set button
        self.bulk_max_set_button = tk.Button(
                                        par_window,
                                        text='set',
                                        command=lambda: self.set_bulk_max(),
                                        width=8)
        self.bulk_max_set_button.grid(row=6, column=2, sticky=tk.W)


        # for parameter disp_min
        # title label and current value
        self.disp_min_label = tk.Label(par_window, text="minimum scaled dispersion", font="-weight bold")
        self.disp_min_label.grid(row=7, column=1, columnspan=2, sticky=tk.W+tk.E)

        # entry field
        self.disp_min_entry = tk.Entry(
                                    par_window,
                                    width=8)
        self.disp_min_entry.grid(row=9, column=1, sticky=tk.E)
        self.disp_min_entry.insert(tk.END, self.disp_min)

        # set button
        self.disp_min_set_button = tk.Button(
                                        par_window,
                                        text='set',
                                        command=lambda: self.set_disp_min(),
                                        width=8)
        self.disp_min_set_button.grid(row=9, column=2, sticky=tk.W)

        # for parameter maxiter
        # title label and current value
        self.maxiter_label = tk.Label(par_window, text="maximum number of iterations", font="-weight bold")
        self.maxiter_label.grid(row=10, column=1, columnspan=2, sticky=tk.W+tk.E)

        # entry field
        self.maxiter_entry = tk.Entry(
                                    par_window,
                                    width=8)
        self.maxiter_entry.grid(row=12, column=1, sticky=tk.E)
        self.maxiter_entry.insert(tk.END, self.maxiter)

        # set button
        self.maxiter_set_button = tk.Button(
                                        par_window,
                                        text='set',
                                        command=lambda: self.set_maxiter(),
                                        width=8)
        self.maxiter_set_button.grid(row=12, column=2, sticky=tk.W)

        # buton for wuitting and return to main window
        self.return_button = tk.Button(par_window,
                                       text='Return to main window',
                                       command=par_window.destroy,
                                       highlightbackground='#f47a60',
                                       highlightthickness=0.1,)
        self.return_button.grid(row=16, column=1,
                                columnspan=2, sticky=tk.W+tk.E, padx=5, pady=10,
                                ipadx=5, ipady=5)

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
            # update display of current bulk_min
            self.bulk_min_value_label['text'] = "{}".format(self.bulk_min)
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
            self.bulk_max_value_label['text'] = "{}".format(self.bulk_max)
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
            # update disp_min
            self.disp_min = new_val
            # update display of current bulk_min
            self.disp_min_value_label['text'] = "{}".format(self.disp_min)
        except ValueError:
            messagebox.showerror("Input error", """Please provdide a positive value. Examples: "0", "0.5", "2.7".""")
            return 0

    def set_maxiter(self):
        # get input and check validity
        input = self.maxiter_entry.get()
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
            self.maxiter_value_label['text'] = "{}".format(self.maxiter)
        else:
            messagebox.showerror("Input error", """Please provdide a positive integer. Examples: "50", "587", "1000".""")

    def reset_default_params(self):
        self.bulk_min = self._bulk_min_default
        self.bulk_min_value_label['text'] = "{}".format(self.bulk_min)

        self.bulk_max = self._bulk_max_default
        self.bulk_max_value_label['text'] = "{}".format(self.bulk_max)

        self.disp_min = self._disp_min_default
        self.disp_min_value_label['text'] = "{}".format(self.disp_min)

        self.maxiter = self._maxiter_default
        self.maxiter_value_label['text'] = "{}".format(self.maxiter)

    def cellanneal(self):
        # check if input and output is set
        if self.bulk_df_is_set == 0:
            messagebox.showerror("Data error", """Please select a mixture data file in section 1).""")
            return 0
        if self.celltype_df_is_set == 0:
            messagebox.showerror("Data error", """Please select a signature data file in section 1).""")
            return 0
        if self.output_path_is_set == 0:
            messagebox.showerror("Data error", """Please select a folder for storing results in section 1).""")
            return 0

        print("\n+++ Welcome to cellanneal! +++")
        print("{}\n".format(time.ctime()))

        # import bulk and celltype data
        # bulk data
        print("\n+++ Importing mixture data ... +++ \n")
        try:
            bfile = self.bulk_data_path.get()
            # depending on extension, use different import function
            if bfile.split(".")[-1] in ["csv", 'txt']:
                self.bulk_df = read_csv(bfile, index_col=0,
                                        sep=None)
            elif bfile.split(".")[-1] in ["xlsx"]:
                self.bulk_df = read_excel(bfile, index_col=0, engine='openpyxl')
            elif bfile.split(".")[-1] in ["xls"]:
                self.bulk_df = read_excel(bfile, index_col=0, engine='xlrd')
            else:
                raise ImportError
            # here, in order to make further course case insensitive,
            # change all gene names to uppercase only
            self.bulk_df.index = self.bulk_df.index.str.upper()
            # also, if there are duplicate genes, the are summed here
            self.bulk_df = self.bulk_df.groupby(self.bulk_df.index).sum()
            # finally, if there are nan's after import, set them to 0 to
            # avoid further issues
            self.bulk_df = self.bulk_df.fillna(0)

        except:
            messagebox.showerror("Import error", """Your mixture data file could not be imported. Please check the documentation for format requirements and look at the example mixture data files.""")
            print("+++ Aborted. +++")
            return 0

        # celltype data
        print("\n+++ Importing signature data ... +++ \n")
        try:
            # depending on extension, use different import function
            cfile = self.celltype_data_path.get()
            if cfile.split(".")[-1] in ["csv", 'txt']:
                self.celltype_df = read_csv(cfile, index_col=0,
                                            sep=None)
            elif cfile.split(".")[-1] in ["xlsx"]:
                self.celltype_df = read_excel(cfile, index_col=0, engine='openpyxl')
            elif cfile.split(".")[-1] in ["xls"]:
                self.celltype_df = read_excel(cfile, index_col=0, engine='xlrd')
            else:
                raise ImportError
            # here, in order to make further course case insensitive,
            # change all gene names to uppercase only
            self.celltype_df.index = self.celltype_df.index.str.upper()
            # also, if there are duplicate genes, the are summed here
            self.celltype_df = self.celltype_df.groupby(self.celltype_df.index).sum()
            # finally, if there are nan's after import, set them to 0 to
            # avoid further issues
            self.celltype_df = self.celltype_df.fillna(0)

        except:
            messagebox.showerror("Import error", """Your signature data file could not be imported. Please check the documentation for format requirements and look at the example signature data files.""")
            print("+++ Aborted. +++")
            return 0

        cellanneal_pipe(
            Path(self.celltype_data_path.get()),  # path object!
            self.celltype_df,
            Path(self.bulk_data_path.get()),  # path object!
            self.bulk_df,
            self.disp_min,
            self.bulk_min,
            self.bulk_max,
            self.maxiter,
            Path(self.output_path.get()))  # path object!

root = tk.Tk()
ca_gui = cellgui(root)
root.mainloop()
