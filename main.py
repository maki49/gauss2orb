import argparse
import numpy as np
import matplotlib.pyplot as plt

def gauss(coef, beta, r):
    return coef*np.exp(-abs(beta)*r**2)

class Chi:
    l=0 # angular momentum
    coef=[] # coefficients
    exp=[] # exponents
    
    def nprim (self):
        return len(self.coef)
    
    def __init__(self, l, coef, exp):
        self.l=l
        self.coef=coef
        self.exp=exp
        
    def __str__(self):
        return " l="+str(self.l)+",\n coef="+str(self.coef)+",\n exp="+str(self.exp)

    def __call__(self, r):
        return np.sum([gauss(self.coef[i], self.exp[i], r) for i in range(self.nprim())])
    
    

def read_fhiaims(file):
    orb_dict={}
    elem=""
    with open(file, 'r') as f:
        lines=f.readlines()
    for iline in range(len(lines)):
        if lines[iline].strip(): # skip blank lines
            l=lines[iline].split()
            if l[0] == "gaussian":
                if lines[iline-1].split()[0]=="#":  # new element
                    elem=lines[iline-1].split()[1]
                    orb_dict[elem]=[]
                angmom=int(l[1])
                nprim=int(l[2])
                coef=[]
                exp=[]
                if(nprim>1):
                    for i in range(nprim):
                        iline+=1
                        l=lines[iline].split()
                        coef.append(float(l[0]))
                        exp.append(float(l[1]))
                    orb_dict[elem].append(Chi(angmom, coef, exp) )
                else:
                    orb_dict[elem].append(Chi(angmom, [1.0], [float(l[3])]))
    return orb_dict

def plot_orbital(orb_dict,  rcut, dr, family):
    r=np.arange(0.0, rcut, dr)
    for elem, chis in orb_dict.items():
        plt.figure()
        plt.title(elem)
        plt.xlabel("r (a.u.)")
        for chi in chis:
            if(chi.nprim()<6):# plot only valence orbitals. Note: 6 is only suitable for 6-31G basis set
                y=np.array([chi(r[i]) for i in range(len(r))])  # use __call__ function to evaluate value on the grid
                plt.plot(r, y, label="l="+str(chi.l)+", nprim="+str(chi.nprim()))
        plt.legend()
        plt.savefig(elem+".png", bbox_inches='tight')
    
    
def count_l(chis, l):
    count=0
    for chi in chis:
        if chi.l==l:
            count+=1
    return count
    
def convert_abacus(orb_dict, rcut, dr, family ):
    r=np.arange(0.0, rcut, dr)
    for elem, chis in orb_dict.items():
        with open(elem+"-"+family+"-"+str(rcut)+"au.orb", 'w') as f:
            f.write("-------------------------------------------------------------------------------\n")
            f.write("Element                     "+elem+"\n")
            f.write("Radius Cutoff(a.u.)         "+str(rcut)+"\n")
            f.write("Lmax                        "+str(chis[-1].l)+"\n")
            f.write("Number of Sorbital-->       "+str(count_l(chis, 0))+"\n")
            if(chis[-1].l>=1):
                f.write("Number of Porbital-->       "+str(count_l(chis, 1))+"\n")
            if(chis[-1].l>=2):
                f.write("Number of Dorbital-->       "+str(count_l(chis, 2))+"\n")
            f.write("-------------------------------------------------------------------------------\n")
            f.write("SUMMARY  END\n\n")
            f.write("Mesh                        "+str(len(r))+"\n")
            f.write("dr                          "+str(dr)+"\n")
            nchi_in_l=0
            for ichi in range(len(chis)):
                chi=chis[ichi]
                if(chi.nprim()<int(family[0])):# output only valence orbitals
                    f.write("                Type                   L                   N\n")
                    if(ichi>0 and chis[ichi-1].l!=chi.l):
                        nchi_in_l=0
                    f.write("                   0                   "+str(chi.l)+"                 "+str(nchi_in_l)+"\n")
                    nchi_in_l+=1
                    y=np.array([chi(r[i]) for i in range(len(r))])  # use __call__ function to evaluate value on the grid
                    for i in range(len(r)):
                        f.write("{:10.12e} ".format(y[i]))
                        if((i+1)%4==0 or i==len(r)-1):
                            f.write("\n")
                    
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str, help='Input filename.')
    parser.add_argument('--family', type=str, help='Basis family.', default="6-31G")
    parser.add_argument('--rcut', type=float, help='Cutoff radius.',default=10.0)
    parser.add_argument('--dr', type=float, help='Radial step.', default=0.01)
    args=parser.parse_args()
    orb_dict=read_fhiaims(args.filename)
    # plot_orbital(orb_dict, args.rcut, args.dr, args.family)
    convert_abacus(orb_dict, args.rcut, args.dr, args.family)
    
    