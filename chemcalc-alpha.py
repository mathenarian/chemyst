#!/bin/env python3
#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>

# Author: Bruno Henc hencb@openmailbox.org
import peq, tkinter
window = tkinter.Tk()
def callback():
    s = ent.get()
    cheq = peq.chem_eq(s)
    text.insert(tkinter.CURRENT, str(cheq)+"\n")
lbl = tkinter.Label(window, text= "Enter equation:")
ent = tkinter.Entry(window, width= "300", xscrollcommand=True)
btn = tkinter.Button(window, text="Calculate!", command=callback)
text = tkinter.Text(window)
window.geometry("500x500")
lbl.pack()
ent.pack()
btn.pack()
text.pack()
window.mainloop()
