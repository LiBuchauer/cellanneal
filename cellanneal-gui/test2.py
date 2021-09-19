from tkinter import *
from tkinter import ttk
from tkinter import filedialog
gui = Tk()
gui.geometry("400x400")
gui.title("FC")

class FolderSelect(Frame):
    def __init__(self,parent=None,folderDescription="",**kw):
        Frame.__init__(self,master=parent,**kw)
        self.folderPath = StringVar()
        self.lblName = Label(self, text=folderDescription)
        self.lblName.grid(row=0,column=0)
        self.entPath = Entry(self, textvariable=self.folderPath)
        self.entPath.grid(row=0,column=1)
        self.btnFind = ttk.Button(self, text="Browse Folder",command=self.setFolderPath)
        self.btnFind.grid(row=0,column=2)
    def setFolderPath(self):
        folder_selected = filedialog.askdirectory()
        self.folderPath.set(folder_selected)
    @property
    def folder_path(self):
        return self.folderPath.get()

def doStuff():
    folder1 = directory1Select.folder_path
    folder2 = directory2Select.folder_path
    folder3 = directory3Select.folder_path
    print("Doing stuff with folder", folder1, folder2, folder3)

folderPath = StringVar()

directory1Select = FolderSelect(gui,"Select Folder 1")
directory1Select.grid(row=0)

directory2Select = FolderSelect(gui,"Select Folder 2")
directory2Select.grid(row=1)

directory3Select = FolderSelect(gui,"Select Folder 3")
directory3Select.grid(row=2)


c = ttk.Button(gui, text="find", command=doStuff)
c.grid(row=4,column=0)
gui.mainloop()
