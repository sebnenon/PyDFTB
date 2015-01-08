#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys

#========================================
# Parse les fichiers band.out
# pour recuperer le nombre de points k, de
# niveaux occupes et les niveaux occupes
#========================================
class bandFile():
  '''
  this class permits to manage *band.out files
  '''
  def __init__(self,fichier='band.out'):
    '''
    parse le fichier band out
    '''
    try:
      bandhandle = open(fichier,'r')
    except:
      sys.exit('%s must be present in the folder'%fichier)
    bandfile = bandhandle.readlines()
    bandhandle.close()
    
    self.n_kpt = 0
    self.bufflines = []
    self.bandlines = []
    n_eig = 0
    for i in [parse(j) for j in bandfile]:
      if i == []:
	continue
      elif 'KPT' in i[0]:
	self.n_kpt += 1
	n_eig = 0
	if self.n_kpt >= 2:
	  self.bandlines.append(self.bufflines)
	  del(self.bufflines[:])
	continue
      else:
	self.bufflines.append(i)
	if float(i[-1]) == 0.0:
	  continue
	else:
	  n_eig += 1
	  self.n_occup_eigvecs = n_eig
    self.bandlines.append(self.bufflines)
      
  def getBands(self):
    '''
    recupere tout le fichier dans un tableau kpt*eigval,contrib,occup
    '''
    return self.bandlines
    
  def getEigvals(self):
    '''
    recupere les niveaux
    '''
    return [float(i[0]) for i in self.bandlines[:][0]]
    
  def getOccEigvals(self):
    '''
    recupere les niveaux occupes
    '''
    return [float(i[0]) for i in self.bandlines[:][0] if float(i[-1]) > 1.0]


if __name__ == '__main__':
  sys.exit('class bandFile():'+bandFile.__doc__)
