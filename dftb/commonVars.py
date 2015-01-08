#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np

###################### Chemistry ######################
valence = {'H' : 1, 'He' : 1, 'Li' : 4, 'Be' : 4, 'B' : 4, 'C' : 4, 'N' : 4, 'O' : 4, 'F' : 4, 'Ne' : 4, 'Na' : 9, 'Mg' : 9, 'Al' : 9, 'Si' : 9, 'P' : 9, 'S' : 9, 'Cl' : 9, 'Ar' : 9, 'K' : 9, 'Ca' : 9, 'Sc' : 9, 'Ti' : 9, 'V' : 9, 'Cr' : 9, 'Mn' : 9, 'Fe' : 9, 'Co' : 9, 'Ni' : 9, 'Cu' : 9, 'Zn' : 9, 'Ga' : 9, 'Ge' : 9, 'As' : 9, 'Se' : 9, 'Br' : 9, 'Kr' : 9, 'Rb' : 9, 'Sr' : 9, 'Y' : 9, 'Zr' : 9, 'Nb' : 9, 'Mo' : 9, 'Tc' : 9, 'Ru' : 9, 'Rh' : 9, 'Pd' : 9, 'Ag' : 9, 'Cd' : 9, 'In' : 9, 'Sn' : 9, 'Sb' : 9, 'Te' : 9, 'I' : 9, 'Xe' : 9, 'Cs' : 9, 'Ba' : 9}
valence_names = ['s1','p1','p2','p3','d1','d2','d3','d4','d5']

###################### Constants ######################
deg2rad = np.pi / 180.
np.set_printoptions(precision=8)
qp = np.dtype('f16') # quadruple precision
dp = np.dtype('f8') # double precision


if __name__ == '__main__':
  print("""
****************
Common variables
S.Nenon 2012
****************

Edit file to see/edit variables and constants

""")

