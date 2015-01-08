#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys
from parse import parse
from geometry import translation,barycenter,colinearize,rotation

#========================================
# Gestion des fichiers gen
#========================================
class genFile():
  '''
  this class permits to manage and convert .gen files
  '''
  def __init__(self,filein):
    '''
    cree l'objet genFile et initialise les attributs:
	atom_names: 	['atom1', 'atom2', ... , 'atomN']
	n_atoms: 	nombre d'atomes
	geometry: 	[ [1, 'atom1', x1, y1, z1] , ... , [N, 'atomN', xN, yN, zN] ]
	lattice:	[[Ox, Oy, Oz], [v1x, v1y, v1z], [v2x, v2y, v2z], [v3x, v3y, v3z]]
    '''
    
    self.atom_names = [] # ['atom1', 'atom2', ... , 'atomN']
    self.filename = filein
    
    # open file
    genhandle = open(filein,'r')
    genfile = genhandle.readlines()
    while parse(genfile[-1].replace('\n','')) == []:
      del(genfile[-1])
    genhandle.close()
    
    # read header
    self.n_atoms = int(parse(genfile[0])[0])
    if 'C' in parse(genfile[0])[1]:
	    self.isPeriodic = False
    else:
	    self.isPeriodic = True
    self.atom_names = parse(genfile[1])
    
    self.lattice = []
    self.geometry = []
    if self.isPeriodic:
      # read geometry
      for i in genfile[2:-4]:
        line = parse(i)
        self.geometry.append([int(line[0]),self.atom_names[int(line[1])-1],float(line[2]),float(line[3]),float(line[4])])
      
      # read lattice vectors
      for i in genfile[-4:]:
        line = parse(i)
        self.lattice.append([float(line[0]),float(line[1]),float(line[2])])
    else:
      # read geometry
      for i in genfile[2:]:
        line = parse(i)
        self.geometry.append([int(line[0]),self.atom_names[int(line[1])-1],float(line[2]),float(line[3]),float(line[4])])  

     ##################################### End of attributes #############################################################
    
  def writeCom(self,fileout,kind=None):
    with open(fileout,'w') as fout:
      fout.write('%Mem=3900mb\n')
      fout.write('#p rhf/sto3g\n')
      fout.write('\n')
      fout.write('comment\n')
      fout.write('\n')
      fout.write(' 0 1 \n')
      for i in range(self.n_atoms):
        fout.write('%5s %20.10f %20.10f %20.10f\n'%(self.geometry[i][1],self.geometry[i][2],self.geometry[i][3],self.geometry[i][4]))
      if self.isPeriodic:
        for i in range(1,4):
          fout.write('%5s %20.10f %20.10f %20.10f\n'%('Tv',self.lattice[i][0],self.lattice[i][1],self.lattice[i][2]))
      fout.write('\n')
    print(fileout+' written')
    return 1
      
  def writePos(self,const=None):
    with open('POSCAR','w') as fout:
      fout.write('title\n') # title
      fout.write('   {:.14f}\n'.format(1.0)) # scaling factor
      if not self.isPeriodic:
          xxx=max([x[2] for x in self.geometry])*2.0
          yyy=max([x[3] for x in self.geometry])*2.0
          zzz=max([x[4] for x in self.geometry])*2.0
          self.lattice=[[0.0,0.0,0.0],[xxx,0.0,0.0],[0.0,yyy,0.0],[0.0,0.0,zzz]]
      for i in self.lattice[1:]:
            fout.write('{:>23.16f}{:>23.16f}{:>23.16f}\n'.format(i[0],i[1],i[2])) # lattice
      buffn=''
      buffd=''
      for i in self.atom_names:
      	buffn+='   {}'.format(i)
      	cpt=0
      	for j in self.geometry:
      		if j[1]==i:
      			cpt+=1
      	buffd+='{:>6d}'.format(cpt)
      buffn+=' \n'
      buffd+=' \n'
      fout.write(buffn) # atom names
      fout.write(buffd) # atom numbers
      fout.write('Selective dynamics\nDirect\n') # keywords
      for i in self.atom_names:
      	for j in self.geometry:
      		if j[1]==i:
      			fout.write('{:>20.16f}{:>20.16f}{:>20.16f}{:>4s}{:>4s}{:>4s}\n'.format(j[2]/self.lattice[1][0],j[3]/self.lattice[2][1],j[4]/self.lattice[3][2],'T','T','T')) # geometry
      fout.write('\n\n')
    print('POSCAR written')
    return 1
  
  # modif geometry
  def geo2bary(self):
    geo = []
    for i in self.geometry:
      geo.append([float(j) for j in i[2:]])
    geo2 = translation(geo,-1*barycenter(geo))
    newgeom = []
    cpt = 0
    for i in geo2:
      newgeom.append([self.geometry[cpt][0],self.geometry[cpt][1],i[0],i[1],i[2]])
      cpt+=1
    self.geometry = newgeom
    return 1
  
  def trans(self,vec1):
    geo = []
    for i in self.geometry:
      geo.append([float(j) for j in i[2:]])
    geo2 = translation(geo,vec1)
    newgeom = []
    cpt = 0
    for i in geo2:
      newgeom.append([self.geometry[cpt][0],self.geometry[cpt][1],i[0],i[1],i[2]])
      cpt+=1
    self.geometry = newgeom
    return 1
    
  def rot(self,vec1,vec2):
    geo = []
    for i in self.geometry:
      geo.append([float(j) for j in i[2:]])
    geo2 = colinearize(geo,vec1,vec2)
    newgeom = []
    cpt = 0
    for i in geo2:
      newgeom.append([self.geometry[cpt][0],self.geometry[cpt][1],i[0],i[1],i[2]])
      cpt+=1
    self.geometry = newgeom
    if self.isPeriodic:
        lat = []
        for i in self.lattice:
          lat.append([float(j) for j in i])
        lat2 = colinearize(lat,vec1,vec2)
        newlat = []
        cpt = 0
        for i in lat2:
          newlat.append([i[0],i[1],i[2]])
          cpt+=1
        self.lattice = newlat
    return 1
 
  def rotangle(self,anglex,angley,anglez):
    geo = []
    for i in self.geometry:
      geo.append([float(j) for j in i[2:]])
    geo2 = rotation(geo,anglex,angley,anglez,'d')
    newgeom = []
    cpt = 0
    for i in geo2:
      newgeom.append([self.geometry[cpt][0],self.geometry[cpt][1],i[0],i[1],i[2]])
      cpt+=1
    self.geometry = newgeom
    if self.isPeriodic:
        lat = []
        for i in self.lattice:
          lat.append([float(j) for j in i])
        lat2 = colinearize(lat,vec1,vec2)
        newlat = []
        cpt = 0
        for i in lat2:
          newlat.append([i[0],i[1],i[2]])
          cpt+=1
        self.lattice = newlat
    return 1
	
  
  # save geom
  def savegen(self,fich):
    fout = open(fich,'w')
    # set vars
    atomstr = ''
    for i in self.atom_names:
      atomstr+=i.rjust(3)
    if self.isPeriodic:
      calc = 'S'
      footer = ''
      for i in self.lattice:
        footer += ("%.10E"%i[0]).rjust(20)+("%.10E"%i[1]).rjust(20)+("%.10E\n"%i[2]).rjust(20)
    else:
      calc = 'C'
      footer = ''
    # writing header
    fout.write(str(self.n_atoms).rjust(5)+("%c"%(calc)).rjust(3)+"\n"+atomstr+"\n")
    # writing geom
    for i in self.geometry:
      fout.write(("%d"%i[0]).rjust(5)+("%d"%(self.atom_names.index(i[1])+1)).rjust(2)+("%.10E"%i[2]).rjust(20)+("%.10E"%i[3]).rjust(20)+("%.10E"%i[4]).rjust(20)+"\n")
    fout.write(footer)
    fout.close()
    print(fich+' written')
    return 1


if __name__ == '__main__':
  sys.exit('class genFile():'+genFile.__doc__)
