from tkinter import *
import socket

client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

HOST = "0.0.0.0"  # The server's hostname or IP address
PORT = 49468  # The port used by the server

client_socket.connect((HOST, PORT))


class Window(Frame):

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.master = master

        # widget can take all window
        self.pack(fill=BOTH, expand=1)

        text = Label(self, text="ChiTech interactive window", fg="black")
        text.place(x=0, y=0)

        self.power_label = Label(text="Power :", fg="black")
        self.power_label.place(x=250, y=30)

        self.period_label = Label(text="Period:", fg="black")
        self.period_label.place(x=250, y=60)

        self.rho_label = Label(text="Rho   :", fg="black")
        self.rho_label.place(x=250, y=90)

        self.tmp_label = Label(text="Temp  :", fg="black")
        self.tmp_label.place(x=250, y=120)

        self.physics_counter = 0
        self.physics_counter_label = Label(text="0", fg="black")
        self.physics_counter_label.place(x=250, y=0)
        self.physics_ping()

        query_button = Button(self, text="Query", command=click_query_button)
        query_button.place(x=0, y=30)

        exit_button = Button(self, text="Exit", command=click_exit_button)
        exit_button.place(x=0, y=60)

        add_rho_button = Button(self, text="Add Rho", command=click_add_rho_button)
        add_rho_button.place(x=0, y=90)

        rem_rho_button = Button(self, text="Remove Rho", command=click_rem_rho_button)
        rem_rho_button.place(x=0, y=120)

        slow_add_rho_button = Button(self, text="Add Rho", command=click_slow_add_rho_button)
        slow_add_rho_button.place(x=120, y=90)

        slow_rem_rho_button = Button(self, text="Remove Rho", command=click_slow_rem_rho_button)
        slow_rem_rho_button.place(x=120, y=120)

        export_vtk_button = Button(self, text="Export VTK", command=click_export_vtk_button)
        export_vtk_button.place(x=0, y=150)

        stop_button = Button(self, text="STOP", command=click_stop_button)
        stop_button.place(x=120, y=150)

    def physics_ping(self):
        self.physics_counter += 1
        self.physics_counter_label.configure(text=str(self.physics_counter))

        message = b"GET / HTTP/1.1\r\n"
        message += b"\r\n"
        message += b'return reactor_power;'
        message += b'return chi.PostProcessorGetValue("period(s)");'
        message += b'return rho;'
        message += b'return avg_temp'
        message += b'\r\n'
        client_socket.sendall(message)
        data = client_socket.recv(1024)
        results = interpret_chi_message(data.decode('utf-8'))

        try:
            power = float(results[0])
            period = float(results[1])
            rho = float(results[2])
            temp = float(results[3])

            self.power_label.configure(text=f"Power : {power:.3e}")
            self.period_label.configure(text=f"Period: {period:+.4f}")
            self.rho_label.configure(text=f"Rho   : {rho:.3f}")
            self.tmp_label.configure(text=f"Temp  : {temp:.2f}")

        except:
            print("oopsie: " + results.__str__())

        self.after(int(1000.0 / 30), self.physics_ping)


def click_exit_button():
    exit()


def click_query_button():
    message = b"GET / HTTP/1.1\r\n"
    message += b"\r\n"
    message += b'chiLog(LOG_0, "Hello from python");'
    message += b'return chi.PostProcessorGetValue("population");'
    message += b'return reactor_power;'
    message += b'return chi.PostProcessorGetValue("period(s)");'
    message += b'rho_added=rho_added + 0.0001'
    message += b'\r\n'
    client_socket.sendall(message)
    data = client_socket.recv(1024)
    interpret_chi_message(data.decode('utf-8'), verbose=True)


def click_add_rho_button():
    message = b"POST / HTTP/1.1\r\n"
    message += b"\r\n"
    message += b'rho_added=rho_added + 2.1'
    message += b'\r\n'
    client_socket.sendall(message)
    data = client_socket.recv(1024)
    interpret_chi_message(data.decode('utf-8'), verbose=True)


def click_rem_rho_button():
    message = b"POST / HTTP/1.1\r\n"
    message += b"\r\n"
    message += b'rho_added=rho_added - 2.1'
    message += b'\r\n'
    client_socket.sendall(message)
    data = client_socket.recv(1024)
    interpret_chi_message(data.decode('utf-8'), verbose=True)

def click_slow_add_rho_button():
    message = b"POST / HTTP/1.1\r\n"
    message += b"\r\n"
    message += b'rho_added=rho_added + 0.01'
    message += b'\r\n'
    client_socket.sendall(message)
    data = client_socket.recv(1024)
    interpret_chi_message(data.decode('utf-8'), verbose=True)


def click_slow_rem_rho_button():
    message = b"POST / HTTP/1.1\r\n"
    message += b"\r\n"
    message += b'rho_added=rho_added - 0.01'
    message += b'\r\n'
    client_socket.sendall(message)
    data = client_socket.recv(1024)
    interpret_chi_message(data.decode('utf-8'), verbose=True)

def click_export_vtk_button():
    message = b"POST / HTTP/1.1\r\n"
    message += b"\r\n"
    message += b'chiExportMultiFieldFunctionToVTK(ffs_timestamp1, "ZTimeStamp1")'
    message += b'chiExportMultiFieldFunctionToVTK(ffs_timestamp2, "ZTimeStamp2")'
    message += b'\r\n'
    client_socket.sendall(message)
    data = client_socket.recv(1024)
    interpret_chi_message(data.decode('utf-8'), verbose=True)

def click_stop_button():
    message = b"POST / HTTP/1.1\r\n"
    message += b"\r\n"
    message += b'alive = false'
    message += b'\r\n'
    client_socket.sendall(message)
    data = client_socket.recv(1024)
    interpret_chi_message(data.decode('utf-8'), verbose=True)

def interpret_chi_message(message: str, verbose: bool = False):
    lines = message.splitlines()
    data_lines = []
    data_active = False
    for line in lines:
        if len(line) == 0:
            data_active = True
            continue
        if data_active:
            data_lines.append(line)

    results = []
    for line in data_lines:
        line_results = line.split(";")
        for result in line_results:
            if len(result) != 0:
                results.append(result)

    if verbose:
        print(results)
    return results


scalars_from_chi_tech = {""}

root = Tk()
app = Window(root)
root.wm_title("Tkinter button")
root.geometry("600x200")
root.mainloop()
