"""
Author: Reuben Allen
Date: 6/27/2022

Given an amino acid sequence, this simple program will return
the molecular formula for the protein and the masses of various
isotopically-enriched variants for that protein. This program will
only handle the 20 canonical amino acids in single-letter abbreviation
format."""

# import python modules
import numpy as np

# dictionary giving the molecular formula of the amino acids
# format: "single-letter code": # of C, # of H, # of O, # of N, # of S
# this assumes the amino acid has formed two peptide bonds
amine_dict = {
    "G": [2,3,1,1,0],
    "A": [3,5,1,1,0],
    "V": [5,9,1,1,0],
    "L": [6,11,1,1,0],
    "I": [6,11,1,1,0],
    "F": [9,9,1,1,0],
    "Y": [9,9,2,1,0],
    "W": [11,10,1,2,0],
    "S": [3,5,2,1,0],
    "T": [4,7,2,1,0],
    "P": [5,7,1,1,0],
    "Q": [5,8,2,2,0],
    "C": [3,5,1,1,1],
    "M": [5,9,1,1,1],
    "N": [4,6,2,2,0],
    "D": [4,5,3,1,0],
    "E": [5,7,3,1,0],
    "K": [6,12,1,2,0],
    "R": [6,12,1,4,0],
    "H": [6,7,1,3,0]
    }

# dictionary with the most abundant isotope for each element [index,isotope]
isotope_dict = {
    "C": [0,13.00335],
    "H": [1,2.014102],
    "O": [2,17.99916],
    "N": [3,15.00011],
    "S": [4,33.96787],
}


# numpy array with mean and modal mass of isotopic distribution
# order from left to right is C,H,O,N,S
# top is mean, bottom is most abundant isotope
mass_arr = np.array([[12.0108,1.00795,15.9994,14.0067,32.066],[12,1.007825,15.99491,14.00307,31.97207]])

def check_str(amine_dict):
    """This function fetches user inputs and checks them before proceeding."""
    isotope_list = []

    while True:
        # get sequence
        print('Enter the amino acid sequence as a continuous string of single letter abbreviations:')
        sequence = input()
        sequence_check = 0

        for key in amine_dict:
            sequence_num = sequence.count(key) # will return the number of occurences of a particular amino acid in the sequence
            sequence_check += sequence_num

        if sequence_check != len(sequence): # 
            print("Invalid Characters Entered. Please Try Again!")
            continue

        else:
            break
        # get isotopes
    while True:
        print('Enter the atomic symbol of the desired isotope (multiple labels should be separate by ";"):')
        print('Key: "N" - Nitrogen 15')
        print('     "O" - Oxygen 18')
        print('     "H" - Deuterium')
        print('     "C" - Carbon 13')
        print('     "S" - Sulfur 34')
        print('     Custom: (H,C,O,N,S),MW')
        isotope_str= input()
        if isotope_str == '':
            break
        else:
            try:
                temp = isotope_str.split(";")
                for i in temp:
                    entry = i.split(",")
                    symbol = entry[0] # atomic symbol

                    if isotope_dict.get(symbol) == None:
                        print('Invalid Atomic Symbol. Please Try Again!')
                        raise ValueError

                    if len(entry) > 1:   
                        isotope = float(entry[1]) # isotope
                        if isotope <= 0:
                            raise ValueError
                                
                    else:
                        isotope = None
                    isotope_list.append([symbol,isotope])
                break

            except ValueError:
                print('Invalid Isotope Entered. Check that the mass is a positive real number')
                continue
    return sequence, isotope_list

# class to store polypeptide info
class Peptide:
    def __init__(self,sequence,isotopes):
        self.sequence = sequence
        self.isotopes = isotopes
        self.formula = self.get_formula()
        self.mass = self.calculate_mass()
        self.isotope_mass = self.enrich()

    def get_formula(self):
        formula = np.zeros(5) # array for counting the number of atoms
        for i in self.sequence:
            formula = np.add(formula,np.array(amine_dict[i]))
        return formula.astype(np.int64)

    def calculate_mass(self):
        mass = np.array([18.015,18.01056]) # starting mass to account for loss of water
        aa_mass = np.multiply(self.formula,mass_arr) # get mass of all atoms from formula
        mass = np.add(np.sum(aa_mass,axis=1).reshape((1,2)),mass)
        return mass

    def enrich(self):
        mass = self.mass
        if len(self.isotopes) > 0:
            for i in self.isotopes: # iterate through isotope labels
                symbol = i[0] # get atomic symbol
                key = isotope_dict[symbol][0]
                num_atoms = self.formula[key] # number of atoms of that element in protein
                mass = np.subtract(self.mass,np.multiply(num_atoms,mass_arr[:,key])) # remove this component from mass

                if i[1] == None: # preset isotope
                    isotope = isotope_dict[symbol][1]
                    mass = np.add(mass,num_atoms*isotope)
                else:
                    isotope = i[1]
                    mass = np.add(mass,num_atoms*isotope)
        return mass
    
    def summarize(self):
        print("")
        print("Molecular Formula: C-{0} H-{1} O-{2} N-{3} S-{4}".format(self.formula[0],self.formula[1],
        self.formula[2],self.formula[3],self.formula[4]))
        print("Average Mass (amu)")
        print("  Normal: {0}".format(round(self.mass[0,0],5)))
        print("  Enriched: {0}".format(round(self.isotope_mass[0,0],5)))
        print("Monoisotopic Mass (amu)")
        print("  Normal: {0}".format(round(self.mass[0,1],5)))
        print("  Enriched: {0}".format(round(self.isotope_mass[0,1],5)))
        print("")

# main code
if __name__ == "__main__":
    sequence,isotopes = check_str(amine_dict) # get user input
    polypeptide = Peptide(sequence,isotopes)
    polypeptide.summarize()
