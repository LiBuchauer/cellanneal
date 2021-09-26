import tkinter as tk
from PIL import Image, ImageTk
import pandas as pd
from tkinter.filedialog import askopenfile, askdirectory
import tkinter.scrolledtext as scrolledtext
from tkinter import messagebox
from pathlib import Path
from cellanneal import cellanneal_pipe, repeatanneal_pipe

import sys
import os
from subprocess import Popen, PIPE
from threading import Thread

class cellgui:

    def __init__(self, root):
        self.root = root
        root.title("cellanneal")

        # place holders for required data
        self.bulk_df = pd.DataFrame()
        self.bulk_folder_path = tk.StringVar()
        self.bulk_df_is_set = 0  # check later whether data has been selected
        self.celltype_df = pd.DataFrame()
        self.celltype_folder_path = tk.StringVar()
        self.celltype_df_is_set = 0

        # output path
        self.output_path = tk.StringVar()
        self.output_path_is_set = 0

        # parameters for the deconvolution procedure, set to default values
        self.bulk_min = 1e-5  # minimum expression in bulk
        self.bulk_min_default = 1e-5

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
        self.canvas.grid(columnspan=7, rowspan=14)
        # for easier grid layout changes
        i_i = 1  # start of import section
        p_i = 5  # start of parameter section
        d_i = 11  # start of deconv section
        t_i = 13  # start of text display section

        """ self.logo """
        self.logo = Image.open('logo_orange.png')
        # convert pillow image to tkinter image
        self.logo = ImageTk.PhotoImage(self.logo)
        # place image inside label widget (both lines below are needed)
        self.logo_label = tk.Label(image=self.logo)
        self.logo_label.image = self.logo
        # place the self.logo label inside the box using grid method
        self.logo_label.grid(column=3, row=0)

        # for buttons belows
        self.ca_button = Image.open('cellanneal_button.png')
        self.ca_button = ImageTk.PhotoImage(self.ca_button)
        self.ca_label = tk.Label(image=self.ca_button)
        self.ca_label.image = self.ca_button

        self.ra_button = Image.open('repeatanneal_button.png')
        self.ra_button = ImageTk.PhotoImage(self.ra_button)
        self.ra_label = tk.Label(image=self.ra_button)
        self.ra_label.image = self.ra_button

        self.running_img = Image.open('running.png')
        self.running_img = ImageTk.PhotoImage(self.running_img)
        self.running_label = tk.Label(image=self.running_img)
        self.running_label.image = self.running_img

        """ main section labels, structure """
        self.sec1_label = tk.Label(root, text="1) Select source data \nand output folder.", font=('Helvetica', 14, 'bold'))
        self.sec2_label = tk.Label(root, text="2) Set parameters. \n[optional]", font=('Helvetica', 14, 'bold'))
        self.sec3_label = tk.Label(root, text="3) Run deconvolution.", font=('Helvetica', 14, 'bold'))
        self.sec1_label.grid(row=i_i, column=0, sticky='w')
        self.sec2_label.grid(row=p_i, column=0, sticky='w')
        self.sec3_label.grid(row=d_i, column=0, sticky='w')

        """ data import section """
        # import bulk data
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
        self.output_folder_entry = tk.Entry(root, textvariable=self.output_path)
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
        # make button for cellanneal
        self.cellanneal_button = tk.Button(
                                        root,
                                        text='cellanneal',
                                        font="-weight bold ",
                                        image=self.ca_button,
                                        highlightbackground='#f47a60',
                                        command=lambda: self.cellanneal())
        self.cellanneal_button.grid(
                                row=d_i,
                                column=1,
                                columnspan=3,
                                sticky=tk.W+tk.E,
                                padx=10, pady=10)

        # make button for cellanneal
        self.repeatanneal_button = tk.Button(
                                        root,
                                        text='repeatanneal',
                                        font="-weight bold ",
                                        image=self.ra_button,
                                        command=lambda: self.repeatanneal(),
                                        highlightbackground='#f47a60')
        self.repeatanneal_button.grid(
                                    row=d_i,
                                    column=4,
                                    columnspan=3,
                                    sticky=tk.W+tk.E,
                                    padx=10, pady=10)

        """ progress display section """
        self.progress_label = tk.Label(root, text="Progress updates", font="-weight bold")
        self.progress_label.grid(
                            row=t_i,
                            column=1,
                            columnspan=6,
                            sticky=tk.W+tk.E)
        self.progress_text = scrolledtext.ScrolledText(root, height=10, width=50)
        self.progress_text.grid(
                                row=t_i+1,
                                column=1,
                                columnspan=6,
                                sticky=tk.W+tk.E)
        self.progress_text.insert(tk.END, 'would be great to see progress updates here')
        # Create a buffer for the stdout
        self.stdout_data = ""
        # A tkinter loop that will show `self.stdout_data` on the screen
        self.show_stdout()

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
                self.bulk_df_is_set = 1
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
            except:
                messagebox.showerror("Import error", """Your celltype data file could not be imported. Please check the documentation for format requirements and look at the example celltype data file.""")

    def select_output_folder(self):
        folder = askdirectory(
                parent=root,
                title="Choose a folder.")
        if folder:
            self.output_path.set(folder)
            self.output_path_is_set = 1

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

    def cellanneal(self):
        # check if input and output is set
        if self.bulk_df_is_set == 0:
            messagebox.showerror("Data error", """Please select a bulk data file in section 1).""")
            return 0
        if self.celltype_df_is_set == 0:
            messagebox.showerror("Data error", """Please select a celltype data file in section 1).""")
            return 0
        if self.output_path_is_set == 0:
            messagebox.showerror("Data error", """Please select a folder for storing results in section 1).""")
            return 0

        # the progress text box should be emptied when a new round is started
        self.progress_text.config(state=tk.NORMAL)
        self.progress_text.delete('1.0', tk.END)

        # start subprocess
        self.subprocess = Popen([sys.executable, "-u",
                                 'cellanneal_pipeline_script.py',
                                 self.celltype_folder_path.get(),
                                 self.bulk_folder_path.get(),
                                 str(self.disp_min),
                                 str(self.bulk_min),
                                 str(self.bulk_max),
                                 str(self.maxiter),
                                 str(self.output_path.get())], stdout=PIPE)

        # Create a new thread that will read stdout and write the data to
        # `self.stdout_buffer`
        thread = Thread(
                    target=self.read_output,
                    args=(self.subprocess.stdout, ))
        thread.start()

    def read_output(self, pipe):
        """Read subprocess' output and store it in `self.stdout_data`."""
        while True:
            data = os.read(pipe.fileno(), 1 << 20)
            # Windows uses: "\r\n" instead of "\n" for new lines.
            data = data.replace(b"\r\n", b"\n")
            if data:
                self.stdout_data += data.decode()
            else:  # clean up
                self.root.after(5000, self.stop)  # stop in 5 seconds
                return None

    def show_stdout(self):
        """Read `self.stdout_data` and put the data in the GUI."""
        data = self.stdout_data
        # after having grabbed it, empty the string
        self.stdout_data = ""
        self.progress_text.insert(tk.END, data)
        self.progress_text.see(tk.END)
        self.root.after(100, self.show_stdout)

    def stop(self, stopping=[]):
        """Stop subprocess"""
        if stopping:
            return # avoid killing subprocess more than once
        stopping.append(True)

        self.subprocess.terminate() # tell the subprocess to exit

        # kill subprocess if it hasn't exited after a countdown
        def kill_after(countdown):
            if self.subprocess.poll() is None: # subprocess hasn't exited yet
                countdown -= 1
                if countdown < 0: # do kill
                    self.subprocess.kill() # more likely to kill on *nix
                else:
                    self.root.after(1000, kill_after, countdown)
                    return # continue countdown in a second

            self.subprocess.stdout.close()  # close fd
            self.subprocess.wait()          # wait for the subprocess' exit

        kill_after(countdown=5)



    def repeatanneal(self):
        # check if input and output is set
        if self.bulk_df_is_set == 0:
            messagebox.showerror("Data error", """Please select a bulk data file in section 1).""")
            return 0
        if self.celltype_df_is_set == 0:
            messagebox.showerror("Data error", """Please select a celltype data file in section 1).""")
            return 0
        if self.output_path_is_set == 0:
            messagebox.showerror("Data error", """Please select a folder for storing results in section 1).""")
            return 0

        # the progress text box should be emptied when a new round is started
        self.progress_text.config(state=tk.NORMAL)
        self.progress_text.delete('1.0', tk.END)

        # start subprocess
        self.subprocess = Popen([sys.executable, "-u",
                                 'repeatanneal_pipeline_script.py',
                                 self.celltype_folder_path.get(),
                                 self.bulk_folder_path.get(),
                                 str(self.disp_min),
                                 str(self.bulk_min),
                                 str(self.bulk_max),
                                 str(self.maxiter),
                                 str(self.N_repeat),
                                 str(self.output_path.get())], stdout=PIPE)

        # Create a new thread that will read stdout and write the data to
        # `self.stdout_buffer`
        thread = Thread(
                    target=self.read_output,
                    args=(self.subprocess.stdout, ))
        thread.start()

    def quit(self):
        self.subprocess.kill() # exit subprocess if GUI is closed (zombie!)
        self.root.destroy()

root = tk.Tk()
ca_gui = cellgui(root)
# make sure all processes are closed when user clicks X
root.protocol("WM_DELETE_WINDOW", ca_gui.quit)
root.mainloop()
