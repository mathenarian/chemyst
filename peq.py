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

def findlines(start, end, l):
    si = se = 0
    for i in range(len(l)):
        if l[i].find(start) != -1:
            si = i
            break
    for j in range(si, len(l)):
        if l[j].find(end) != -1:
            se = j
            break
    r = []
    for i in range(si, se+1):
        r.append(l[i])
    postr = []
    for j in range(se, len(l)):
        postr.append(l[j])
    return r, postr

def grabtable(name="table.txt"):
    url = "http://en.wikipedia.org/wiki/Extended_periodic_table_(large_version)"
    rawtable = []
    import urllib.request
    local_filename, headers = urllib.request.urlretrieve(url)
    with open(local_filename, encoding="utf-8") as cfile:
        clist = list(cfile)
        ps = '<h2><span class="mw-headline" id="Preparation">Preparation</span>'
        pe = "<h3>"
        r = findlines(ps, pe, clist)[0]
        left = r
        with open(file=name, mode="w", encoding="utf-8") as dfile:
            while True:
                r, left = findlines('<td title=\"', '/td>', left)
                symbol = ""
                name = ""
                mmass = ""
                st = "<td title=\""
                for i in r:
                    if i.find(st) != -1:
                        name = i[i.find(st)+len(st):i.find(":"):]
                    if i.find('<a href=') != -1 and i.find("</big>") != -1:
                        symbol = i[i.find("<b>")+len("<b>"):i.find("</b>"):]
                    if i.find('class="noprint"') != -1:
                        ti = i[i.find("<small>")+len("<small>"):i.find("<br />"):]
                        dfile.write(symbol+" "+name+" "+mmas+"\n")
                if " ".join(r).find("Unseptbium") != -1: 
                    break
    return 0

def gettable(name="table.txt"):
    with open(name) as cfile:
        names = []
        molmass = []
        elements = []
        for i in cfile:
            l = i.split()
            names += [l[0]]
            elements += [l[1]]
            molmass += [l[2]]
    return names, elements, molmass

class matrix:
    def __init__(self, rl, pl):
        self.solutions = [1] + [0 for v in range(len(rl)+len(pl)-1)]
        self.matrix = [[1] + [0 for v in range(len(rl)+len(pl)-1)]]
        genf = []
        for i in rl:
            for j in range(len(i.formula)):
                if i.formula[j] != 0 and j not in genf:
                    genf += [j]
                    temp = [rl[k].formula[j] for k in range(len(rl))] \
                        +[-1*pl[l].formula[j] for l in range(len(pl))]
                    if self.matrix.count(temp) == 0 and len(self.matrix) < len(self.solutions):
                        self.matrix.append(temp)
        self.mtx = self.matrix
        matrix.solve(self)
    def eliminate_rows(self):
        order = {}
        while len(list(order.keys())) < len(self.mtx):
            for i in range(len(self.mtx)):
                for p in range(len(self.mtx)):
                    if self.mtx[p][i] != 0:
                        inorder = False
                        for r in order.keys():
                            if p == order[r]:
                                inorder = True
                                break
                        if inorder:
                            continue
                        div = self.mtx[p][i]
                        order[i] = p
                        for j in range(len(self.mtx)):
                            if p == j:
                                continue
                            koef = abs(self.mtx[j][i]) / abs(div)

                            if div > 0:
                                c = 1
                            else:
                                c = -1
                            if self.mtx[j][i] < 0:
                                self.solutions[j] += c*koef*self.solutions[p]
                                for k in range(len(self.mtx[j])):
                                    self.mtx[j][k] = self.mtx[j][k] +  \
                                                    c * koef * self.mtx[p][k]
    
                            elif self.mtx[j][i] > 0:
                                self.solutions[j] -= c*koef*self.solutions[p]
                                for l in range(len(self.mtx[j])):
                                    self.mtx[j][l] = self.mtx[j][l] - \
                                                     c * koef * self.mtx[p][l]
        rmtx = []
        rsol = []
        keys = list(order.keys())
        keys.sort()
        #print(order)
        for o in keys:
            rmtx.append(self.mtx[order[o]])
            rsol.append(self.solutions[order[o]])
        self.mtx = rmtx
        self.solutions = rsol
        #return self.mtx, self.solutions
    def fix_minus_diagonal(self):
        for i in range(len(self.mtx)):
            if self.mtx[i][i] < 0:
                c = abs(self.mtx[i][i])
                self.solutions[i] *= -1 / c
                for j in range(len(self.mtx[i])):
                    self.mtx[i][j] *= -1 / c
            else:
                c = abs(self.mtx[i][i])
                self.solutions[i] *= 1 / c
                for k in range(len(self.mtx[i])):
                    self.mtx[i][k] *= 1 / c
        #return self.mtx, self.solutions
    def solve(self): #Gauss-Jordan elimination
        matrix.print(self)
        matrix.eliminate_rows(self)
        matrix.fix_minus_diagonal(self)
        matrix.fix_sol(self)
        matrix.print(self)
        #return self.solutions    
    def fix_sol(self):
        mini = abs(self.solutions[0])
        FoundFloat = False
        for i in self.solutions:
            if round(i) != i:
                FoundFloat = True
                if abs(i) < mini:
                    mini = abs(i)
        if FoundFloat:
            for j in range(len(self.solutions)):
                self.solutions[j] /= mini
        #return self.solutions

    def fix_minus_diagonal(self):
        for i in range(len(self.mtx)):
            if self.mtx[i][i] < 0:
                c = abs(self.mtx[i][i])
                self.solutions[i] *= -1 / c
                for j in range(len(self.mtx[i])):
                  self.mtx[i][j] *= -1 / c
            else:   
                c = abs(self.mtx[i][i])
                self.solutions[i] *= 1 / c
                for k in range(len(self.mtx[i])):
                    self.mtx[i][k] *= 1 / c
        #return self.mtx, self.solutions  
    def fdifzero(l):
        k = []
        for i in range(len(l)):
            if l[i] != 0:
               k += [i]
        return k

    def maxnz(l):
        for i in l:
            if i != 0:
                maxi = i
        return maxi

    def ninmtx(self, index):
        count = 0
        for i in self.mtx:
            if index in i:
                count += 1
        return count

    def print(self):
        print(self.mtx)
        print(self.solutions)
        print(len(self.mtx), len(self.solutions))
        for o in range(len(self.mtx)):
            for p in range(len(self.mtx[o])):
                print("{}".format(self.mtx[o][p]), end = "   ")
            print("|", self.solutions[o])
        print()
        return ""
    def __repr__(self):
        s = ""
        for o in range(len(self.mtx)):
            for p in range(len(self.mtx[o])):
                s += "{}".format(self.mtx[o][p]) + "\t" 
        
            s += "| {}".format(self.solutions[o]) + "\n"
        s += "\n"
        return s
    def __str__(self):
        s = ""
        for o in range(len(self.mtx)):
            for p in range(len(self.mtx[o])):
                s += "{}".format(self.mtx[o][p]) + "\t" 
        
            s += "| {}".format(self.solutions[o]) + "\n"
        s += "\n"
        return s
    def nindict(inx, dl):
        count = 0
        for i in dl.keys():
            if inx in dl[i]:
                count += 1
        return count

class chem_eq:
    def __init__(self, s, sepsides="=", sepitems="+" ):
        self.s = s
        self.sepsides = sepsides
        self.sepitems = sepitems
        #if sepsides not in s:
            #return ["Error"], ["Error"]
        l = s.split(sepsides)
        self.rl = []
        self.pl = []
        self.rlk = []
        self.plk = []
        for i in l[0].split(sepitems):
            if i.isdigit():
                continue
            else:
                self.rl += [ch_formula(i)]
                self.rlk += [1]
        for j in l[1].split(sepitems):
            if j.isdigit():
                continue
            else:
                self.pl += [ch_formula(j)]
                self.plk += [1]
        #print(self.rl)
        elr = chem_eq.fnumofel(self.rl)
        elp = chem_eq.fnumofel(self.pl)
        if elr != elp:
            for il in list(set(elr).symmetric_difference(set(elp))):
                print(elements[il]+ " not found on one side of reaction")
            #return -5
        self.matrix = matrix(self.rl, self.pl)
        self.solutions = self.matrix.solutions
        #print(self.solutions)
    def fnumofel(l):
        afound = []
        for i in l:
            found = []
            for j in range(len(i.formula)):
                if i.formula[j] != 0:
                    found += [j]
            for k in found:
                fflag= False
                for ii in afound:
                    if ii ==  k:
                        fflag = True
                if not fflag:
                    afound += [k]
        afound.sort()
        return afound

    def genstreq(self):
        rl = self.rl
        pl = self.pl
        sol = self.solutions
        preset = 2
        sub = ""
        for j in range(len(rl)):
            if round(sol[j]) == sol[j]:
                if sol[j] == 1:
                    sub += "{} + ".format(rl[j])
                else:
                    sub += "{} {} + ".format(round(sol[j]), rl[j])
            else:
                sub += "{} {} + ".format(round(sol[j]), rl[j])
        sub = sub.rstrip(" + ")
        sub += " > "
        for k in range(len(pl)):
            if round(sol[len(rl)+k]) == sol[len(rl)+k]:
                if sol[len(rl)+k] == 1:
                    sub += "{} ".format(pl[k])
                else:
                    sub += "{} {} + ".format(round(sol[len(rl)+k]), pl[k])
            else:
                sub += "{} {} + ".format(round(sol[len(rl)+k]), pl[k])
        return sub.rstrip(" + ")

    def gen_masses(self):
        rl = self.rl
        pl = self.pl
        sol = self.solutions
        s = []
        cs = []
        masses = [rl[i].molmass() for i in range(len(rl))]
        masses += [pl[j].molmass() for j in range(len(pl))]
        for i in range(len(rl)):
            s += ["{} :\t{}".format(rl[i], masses[i])]
            cs += ["{} {} :\t{}".format(sol[i], rl[i], sol[i]*masses[i])]
        for j in range(len(pl)):
            s += ["{} :\t{}".format(pl[j], masses[len(rl)+j])]
            cs += ["{} {}:\t{}".format(sol[len(rl)+j], pl[j], \
                     sol[len(rl)+j]*masses[len(rl)+j])]
        return s, cs

    def __repr__(self): #TODO:Balanced eq
        s = "\n\n{}\n\n\n".format(chem_eq.genstreq(self))
        temp, tempp =  gen_masses(rl, pl, sol)
        for i in range(len(temp)):
            s += "{}\t\t{}\n".format(temp[i], tempp[i])
        return s

    def __str__(self): #TODO:Balanced eq
        s = "\n\n{}\n\n\n".format(chem_eq.genstreq(self))
        temp, tempp =  chem_eq.gen_masses(self)
        for i in range(len(temp)):
            s += "{}\t\t{}\n".format(temp[i], tempp[i])
        return s

class ch_formula:
    def __init__(self, s):
        self.original = s
        s = s.replace('[', '(')
        s = s.replace(']', ')')
        s = s.replace('{', '(')
        s = s.replace('}', ')')
        global names
        l = []
        self.formula  = [0 for u in range(len(names))]
        for i in s.split("("):
            l += i.split(")")
        coefl = []
        for j in range(len(l)-1, -1, -1):
            if l[j].isdigit():
                coefl  += [float(l[j])]    
                continue
            else:
                if len(coefl) == 0:
                    coef = 1
                else:
                    coef = 1
                    for k in coefl:
                        coef *= k
                    coefl.pop(len(coefl)-1)
                templ = ch_formula.parse_normal(l[j])
                for i in range(len(templ)):
                    self.formula[i] += templ[i]*coef
        self.molformula = ch_formula.mol_formula(self)
       

    def parse_normal(s):
        global names
        l = [0 for u in range(len(names))]
        for j in range(len(names)-1, -1, -1):
            inx = s.find(names[j])
            sc = s
            while inx != -1:
                num = ""
                for i in range(inx+len(names[j]), len(sc)):
                    if sc[i].isdigit():
                        num += sc[i]
                    elif sc[i] == ".":
                        continue
                    else:
                        break
                if num == "":  
                    num = "1"
                l[j] += float(num)
                sc = sc[inx+len(num)::]
                inx = sc.find(names[j])
            s=s.replace(names[j], "|")
        return l
    def mol_formula(self, l=[]):
        if len(l) == 0:
            l = self.formula
        decimals=5
        global names
        s = ""
        for i in range(len(l)):
            if l[i] != 0:
                if round(l[i]) == l[i]:
                    if l[i] == 1:
                        s += "{}".format(names[i])

                    else:
                        s += "{}{}".format(names[i], int(l[i]))
                else:
                    s += "{}{}".format(names[i], l[i])
        return s
    def molmass(self, l = []):
        global molmass
        r = []
        mass = 0
        if len(l) == 0:
            l = self.formula
        for i in range(len(l)):
            s = molmass[i].lstrip("[")
            s = s.rstrip("]").strip()
            if s.find("("):
                s = s[:len(s)-3:]
            if s.find("?") != -1:
                continue
            if l[i] != 0:
                mass += float(s)*l[i]
        return mass           
    def __add__(self):
        pass
    def __str__(self):
        return self.molformula
    def __repr__(self):
        return self.original
def callback(ent, text):
    s = ent.get()
    if s == "":
        return None
    cheq = chem_eq(s)
    text.insert(tkinter.CURRENT, str(cheq)+"\n")
def main():
    window = tkinter.Tk()
    lbl = tkinter.Label(window, text= "Enter equation:")
    ent = tkinter.Entry(window, width= "300", xscrollcommand=True)
    text = tkinter.Text(window)
    btn = tkinter.Button(window, text="Calculate!", command=callback(ent, text))
    window.geometry("500x500")
    lbl.pack()
    ent.pack()
    btn.pack()
    text.pack()
    window.mainloop()    
if __name__ == "__main__":
    tablename = "table.txt"
    import tkinter, urllib
    try:
        fname = open(tablename)
        fname.close()
    except IOError:
        try:
            grabtable(tablename)
        except urllib.error.URLError:
            print("Cannot establish internet connection")
            print("Please move the periodic table file manually")
            print("and move it to the directory from where you execute")
            print("the program and rename it to table.txt")
        else:
            names, elements, molmass = gettable(tablename)
    else:
        names, elements, molmass = gettable(tablename)
        main()
    #"C12H22O11 + KClO3 > KCl + CO2 + H2O"
    #rl, pl, sol = p_parseeq("C12H22O11 + KClO3 > KCl + CO2 + H2O")
    #rl, pl, sol = p_parseeq("H2O2 + KI + H2SO4 > I2 + K2SO4 + H2O")

    #rl, pl, sol = p_parseeq("10 KNO3 + 3 S + 8 C > 2 K2CO3 + 3 K2SO4 + 6 CO2 + 5 N2 ")
    #"N2 + H2 > NH3")
    #rl, pl, sol = p_parseeq("KNO3 + S +  C > K2CO3 + K2SO4 + CO2 + N2") 
    #cli_print_eq(rl, pl, sol)
else:
    tablename = "table.txt"
    import tkinter, urllib
    try:
        fname = open(tablename)
        fname.close()
    except IOError:
        try:
            grabtable(tablename)
        except urllib.error.URLError:
            print("Cannot establish internet connection")
            print("Please move the periodic table file manually")
            print("and move it to the directory from where you execute")
            print("the program and rename it to table.txt")
        else:
            names, elements, molmass = gettable(tablename)
    else:
        names, elements, molmass = gettable(tablename)      
