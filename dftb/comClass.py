#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys
from dftb.parse import parse

#========================================
# Gestion des fichiers com
#========================================
class comFile():
  '''
  this class permits to manage and convert .com files
  '''
  def __init__(self,filename):
    '''
    cree l'objet comFile et initialise les attributs:
	atom_names: 	['atom1', 'atom2', ... , 'atomN']
	n_atoms: 	nombre d'atomes
	geometry: 	[ [1, 'atom1', x1, y1, z1] , ... , [N, 'atomN', xN, yN, zN] ]
    '''
    
    self.atom_names = [] # ['atom1', 'atom2', ... , 'atomN']
    self.atomNumbers = ['xxx','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Uuu','Uub','Uut','Uuq','Uup','Uuh','Uus','Uuo']
    
    # open file
    try:
      comhandle = open(filename,'r')
    except:
      sys.exit(filename+' not found')
    comfile = comhandle.readlines()
    comhandle.close()
    
    self.isPeriodic = False
    self.lattice = []
    self.geometry = []
    # read geometry
    self.n_atoms = 1
    whitelines = 0
    for i in comfile:
      line = parse(i)
      if len(line) == 0:
        whitelines+=1
      if '%' in i or '#' in i or len(line) < 3 or whitelines < 2:
        continue
      elif whitelines == 3:
        break
      else:
        try:      
          self.geometry.append([self.n_atoms,self.atomNumbers[int(line[0])],float(line[1]),float(line[2]),float(line[3])])
          self.atom_names.append(self.atomNumbers[int(line[0])])
        except:
          if line[0] == "Tv":
            self.isPeriodic = True
            self.lattice.append([float(line[1]),float(line[2]),float(line[3])])
            self.n_atoms -= 1
            continue
          else:
            self.geometry.append([self.n_atoms,line[0],float(line[1]),float(line[2]),float(line[3])])
            self.atom_names.append(line[0])
        self.n_atoms += 1
    self.atom_names = list(set(self.atom_names))
    self.n_atoms = len(self.geometry)
     ##################################### End of attributes #############################################################
    
  def writeGen(self,outname):
    self.atomline = ''
    for i in self.atom_names:
      self.atomline +=i.rjust(3)
    with open(outname,'w') as fout:
      if self.isPeriodic:
        fout.write('%5d  S\n'%self.n_atoms)
      else:
        fout.write('%5d  C\n'%self.n_atoms)
      fout.write(self.atomline+'\n')
      for i in range(self.n_atoms):
        try:
          fout.write('%5d%2d %19.10E %19.10E %19.10E\n'%(self.geometry[i][0],self.atom_names.index(self.geometry[i][1])+1,self.geometry[i][2],self.geometry[i][3],self.geometry[i][4]))
        except:
          print(self.geometry[i])
      if self.isPeriodic:
        fout.write('%20.10E %19.10E %19.10E\n'%(0.0,0.0,0.0))
        for i in range(3):  
          fout.write('%20.10E %19.10E %19.10E\n'%(self.lattice[i][0],self.lattice[i][1],self.lattice[i][2]))
    print(outname+' written')


if __name__ == '__main__':
  sys.exit('class comFile():'+comFile.__doc__)
