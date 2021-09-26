import logging
import os
import sys
from subprocess import Popen, PIPE, STDOUT
from threading import Thread

try:
    import tkinter as tk
except ImportError: # Python 3
    import tkinter as tk

info = logging.getLogger(__name__).info

# define dummy subprocess to generate some output
cmd = [sys.executable or "python", "-u", "-c", """
import itertools, time

for i in itertools.count():
    print(i)
    time.sleep(0.5)
"""]

cmd = [sys.executable, "-u",
 'cellanneal_pipeline_script.py',
'/Users/lisa/Desktop/ca_gui/sc_ref_human_liver_patient5.csv',
'/Users/lisa/Desktop/ca_gui/bulk_human_liver_tumor_patient5.csv',
'0.5',
'1e-05',
'0.01',
'20',
'/Users/lisa/Desktop/ca_gui']


class ShowProcessOutputDemo:
    def __init__(self, root):
        """Start subprocess, make GUI widgets."""
        self.root = root

        # start subprocess using a button
        tk.Button(root, text="Start subprocess", command=self.start).pack()
        # stop subprocess using a button
        tk.Button(root, text="Stop subprocess", command=self.stop).pack()

        self.label = tk.Label(root)  # put subprocess output here
        self.label.pack()

        # Create a buffer for the stdout
        self.stdout_data = ""

        # A tkinter loop that will show `self.stdout_data` on the screen
        self.show_stdout()

    def start(self):
        self.proc = Popen(cmd, stdout=PIPE, stderr=STDOUT)
        # Create a new thread that will read stdout and write the data to
        # `self.stdout_buffer`
        thread = Thread(target=self.read_output, args=(self.proc.stdout, ))
        thread.start()

    def read_output(self, pipe):
        """Read subprocess' output and store it in `self.stdout_data`."""
        while True:
            data = os.read(pipe.fileno(), 1 << 20)
            # Windows uses: "\r\n" instead of "\n" for new lines.
            data = data.replace(b"\r\n", b"\n")
            if data:
                info("got: %r", data)
                self.stdout_data += data.decode()
            else:  # clean up
                info("eof")
                self.root.after(5000, self.stop) # stop in 5 seconds
                return None

    def show_stdout(self):
        """Read `self.stdout_data` and put the data in the GUI."""
        print(self.stdout_data.strip("\n"))
        self.label.config(text=self.stdout_data.strip("\n"))
        self.root.after(100, self.show_stdout)

    def stop(self, stopping=[]):
        """Stop subprocess and quit GUI."""
        if stopping:
            return # avoid killing subprocess more than once
        stopping.append(True)

        info("stopping")
        self.proc.terminate() # tell the subprocess to exit

        # kill subprocess if it hasn't exited after a countdown
        def kill_after(countdown):
            if self.proc.poll() is None: # subprocess hasn't exited yet
                countdown -= 1
                if countdown < 0: # do kill
                    info("killing")
                    self.proc.kill() # more likely to kill on *nix
                else:
                    self.root.after(1000, kill_after, countdown)
                    return # continue countdown in a second

            self.proc.stdout.close()  # close fd
            self.proc.wait()          # wait for the subprocess' exit
            self.root.destroy()       # exit GUI
        kill_after(countdown=5)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
root = tk.Tk()
app = ShowProcessOutputDemo(root)
root.protocol("WM_DELETE_WINDOW", app.stop) # exit subprocess if GUI is closed
root.mainloop()
info("exited")
