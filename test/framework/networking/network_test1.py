from tkinter import *
import socket

client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

HOST = "0.0.0.0"  # The server's hostname or IP address
PORT = 49468  # The port used by the server

client_socket.connect((HOST, PORT))

from tkinter import *
import socket

client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

HOST = "0.0.0.0"  # The server's hostname or IP address
PORT = 49468  # The port used by the server

client_socket.connect((HOST, PORT))


def make_chi_post(commands: list):
    message = b"POST / HTTP/1.1\r\n"
    message += b"\r\n"
    for command in commands:
        command_str: str = command
        message += command_str.encode('utf-8')
    message += b'\r\n'

    return message


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


class Window(Frame):

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.master = master

        # widget can take all window
        self.pack(fill=BOTH, expand=1)

        text = Label(self, text="ChiTech interactive window", fg="black")
        text.place(x=0, y=0)

        self.physics_counter = 0
        self.physics_counter_label = Label(text="0", fg="black")
        self.physics_counter_label.place(x=250, y=0)
        self.physics_ping()

        self.time_label = Label(self, text="ChiTech time=", fg="black")
        self.time_label.place(x=250, y=30)

        stop_button = Button(self, text="STOP", command=click_stop_button)
        stop_button.place(x=120, y=150)

        self.input1 = Text(self, width=20, height=5, fg="black", bg="white", bd=1)
        self.input1.place(x=0, y=30)

        post_button = Button(self, text="Post", command=click_post_button)
        post_button.place(x=0, y=100)

    def physics_ping(self):
        self.physics_counter += 1
        self.physics_counter_label.configure(text=str(self.physics_counter))

        message = make_chi_post([f'return chiProgramTime()'])
        client_socket.sendall(message)
        data = client_socket.recv(1024)
        results = interpret_chi_message(data.decode('utf-8'))

        try:
            time = float(results[0])
            self.time_label.configure(text="ChiTech time=" + str(time))
        except:
            print("oopsie")

        self.after(int(1000.0 / 30), self.physics_ping)


def click_stop_button():
    message = make_chi_post(['alive = false'])
    client_socket.sendall(message)
    data = client_socket.recv(1024)
    interpret_chi_message(data.decode('utf-8'), verbose=True)


def click_post_button():
    txt = app.input1.get(1.0, "end-1c")
    message = make_chi_post([txt])
    client_socket.sendall(message)
    data = client_socket.recv(1024)
    results = interpret_chi_message(data.decode('utf-8'))

    for result in results:
        if result == "{}":
            continue
        print(result)


root = Tk()
app = Window(root)
root.wm_title("Tkinter button")
root.geometry("600x200")
root.mainloop()
