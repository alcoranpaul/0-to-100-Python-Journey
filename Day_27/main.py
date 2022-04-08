from tkinter import *
from tkinter import ttk
from Unit7 import *

window = Tk()
window.title("My First HUI")
window.minsize(width=500, height=500)

myNotebook = ttk.Notebook(window)
myNotebook.pack()


# Tab 1
## Functions

def calculateBinomial():
    if tab1_x.get() != "" or tab1_n.get() != "" or tab1_p.get() != "":
        bi_x = int(tab1_x.get())
        bi_n = int(tab1_n.get())
        bi_p = float(tab1_p.get())
        result = binomial(bi_x, bi_n, bi_p)
        tab1_results["text"] = f"Result: {result:.4f}"


def clear():
    tab1_p.delete(0, END)
    tab1_x.delete(0, END)
    tab1_n.delete(0, END)
    tab1_results["text"] = "Result: "


myTab1 = Frame(myNotebook, width=500, height=500)
myTab1.pack(fill="both", expand=1)

# textFrame = Frame(myTab1, width=500, height=100, highlightbackground='red', highlightthickness=3)
# textFrame.pack()

myLabel = Label(myTab1, text="x: number of successes seeking", font=20, anchor="w")
myLabel.pack()
tab1_x = Entry(myTab1, width=10, font=15)
tab1_x.pack()

myLabel = Label(myTab1, text="n: total population", font=20, anchor="w")
myLabel.pack()
tab1_n = Entry(myTab1, width=10, font=15)
tab1_n.pack()

myLabel = Label(myTab1, text="p: probability of success", font=20, anchor="w")
myLabel.pack()
tab1_p = Entry(myTab1, width=10, font=15)
tab1_p.pack()

tab1_results = Label(myTab1, text="Result: ")
tab1_results.pack()

myButton = Button(myTab1, text="Calculate", command=calculateBinomial, font=8)
myButton.place(x=400, y=210)

resetButton = Button(myTab1, text="Clear all", command=clear, font=8)
resetButton.place(x=300, y=210)

myNotebook.add(myTab1, text="Binomial")

# Tab 2
myTab2 = Frame(myNotebook, width=500, height=500, bg='red')
myTab2.pack(fill="both", expand=1)
myNotebook.add(myTab2, text="Geometric")

window.mainloop()
