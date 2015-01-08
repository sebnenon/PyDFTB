#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys

#========================================
# Parse les fichiers eigenvec.out
#========================================
class eigvecFile():
  '''
  this class permits to manage eigenvec.out files
  '''
  def __init__(self,bandfile='band.out'):
    '''
    
    '''
    try:
      eigFile = open('eigenvec.out','r')
    except:
      sys.exit('eigenvec.out must be present in the folder')
    
    gen = genFile()
    band = bandFile(bandfile)
    
    self.n_atoms = gen.n_atoms
    self.n_eigvecs = len(band.getEigvals())
    self.n_kpt = band.n_kpt
    self.atomlist = gen.getAtomList()
    self.atomtypes = gen.getAtomTypes()
    self.eigvals = band.getEigvals()
    
    atom_names = gen.getAtomTypes()
    self.max_oa = 0
    for i in atom_names:
      if valence[i] > self.max_oa:
	self.max_oa = valence[i]
	
    del(band)
    del(gen)
    
    self.eigenvecs = np.zeros((self.n_atoms,self.n_eigvecs,self.n_kpt,self.max_oa))
    
    print('parsing eigenvec.out')
    at_cpt,eig_cpt,kp_cpt,oa_cpt = 0,0,0,0
    ################################################
    eigFile.readline()
    eigFile.readline()
    eigFile.readline()
    print(eigFile.readline())
    for kp_cpt in range(self.n_kpt):
      print("kpoint %d/%d"%(kp_cpt+1,self.n_kpt))
      for eig_cpt in range(self.n_eigvecs):
	for at_cpt in range(self.n_atoms):
	  for oa_cpt in range(valence[self.atomlist[at_cpt][1]]):
	    line = eigFile.readline()
	    try:
	      self.eigenvecs[at_cpt,eig_cpt,kp_cpt,oa_cpt] += float(parse(line)[-1])
	    except:
	      print('Error: kp %d eig %d atome %d oa %d: "%s" total valence:%d'%(kp_cpt+1,eig_cpt+1,at_cpt+1,oa_cpt+1,line,valence[self.atomlist[at_cpt][1]]))
	  eigFile.readline()
      	eigFile.readline()
	eigFile.readline()
    print('done')
    eigFile.close()
    
  
  def contrib(self,name):
    '''
    projette les contributions sur les atomes et oa
    '''
    
    #recupere une liste [ [num,valence,nom], ... ] des atomes traites
    self.indexes = []
    treated_atomtypes = []
    for i in range(self.n_atoms):
      self.indexes.append([self.atomlist[i][0] - 1,valence[self.atomlist[i][1]],self.atomlist[i][1]])
      if self.atomlist[i-1][1] not in treated_atomtypes:
	treated_atomtypes.append(self.atomlist[i-1][1])
    
    
    #ecriture des entetes
    self.outfileoa = '#Eigenvalue'.ljust(15)
    self.outfileatom = '#Eigenvalue'.ljust(15)
    for atom_cpt in self.indexes:
      self.outfileatom += ("%s_%d"%(self.atomlist[atom_cpt[0]][1],self.atomlist[atom_cpt[0]][0])).rjust(15)
      if self.atomlist[atom_cpt[0]][1] not in self.outfileoa:
	for oa_cpt in range(valence[self.atomlist[atom_cpt[0]][1]]):
	  self.outfileoa += ("%s(%s)"%(self.atomlist[atom_cpt[0]][1],valence_names[oa_cpt])).rjust(15)
    self.outfileoa += '\n'
    self.outfileatom += '\n'
    
    #somme sur les atomes
    for eig_cpt in range(self.n_eigvecs):
      lineatom = ('%.6f'%self.eigvals[eig_cpt]).ljust(15)
      for atom in self.indexes:
	atom_cpt = atom[0]
	atomsum = 0.0
	for oa_cpt in range(atom[1]):
	  for k_cpt in range(self.n_kpt):
	    atomsum += self.eigenvecs[atom_cpt,eig_cpt,k_cpt,oa_cpt]
	lineatom += ('%.6f'%atomsum).rjust(15)
      self.outfileatom += "%s\n"%lineatom

    #somme sur les oa
    #phase 1 somme sur les atomes de la liste et les points k
    bufferoa = np.zeros((len(treated_atomtypes),self.n_eigvecs,self.max_oa))
    for eig_cpt in range(self.n_eigvecs):
      for atom in self.indexes:
	atom_cpt = atom[0]
	atom_name = atom[2]
	for oa_cpt in range(valence[atom_name]):
	  for k_cpt in range(self.n_kpt):
	    bufferoa[treated_atomtypes.index(atom_name),eig_cpt,oa_cpt] += self.eigenvecs[atom_cpt,eig_cpt,k_cpt,oa_cpt]
    #phase 2 creation du fichier
    for eig_cpt in range(self.n_eigvecs):
      lineoa = ('%.6f'%self.eigvals[eig_cpt]).ljust(15)
      for atom in treated_atomtypes:
	atom_cpt = treated_atomtypes.index(atom)
	for oa_cpt in range(valence[atom]):
	  lineoa += ('%.6f'%bufferoa[atom_cpt,eig_cpt,oa_cpt]).rjust(15)
	  
      self.outfileoa += "%s\n"%lineoa
      
    oaf = open('%s_oapdos.txt'%name,'w')
    oaf.write(self.outfileoa)
    oaf.close()
    
    atf = open('%s_atpdos.txt'%name,'w')
    atf.write(self.outfileatom)
    atf.close()


if __name__ == '__main__':
  sys.exit('class eigvecFile():'+eigvecFile.__doc__)
