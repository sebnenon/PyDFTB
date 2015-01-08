#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys,os
import numpy as np
from commonVars import *

###################### Functions ######################
#---------------------Vector norm---------------------
def normalize(vec):
  """normalize(vec)
Returns the normalized vector vec
vec can be a 1xN list or array (whatever N)
"""
  norm = 0.
  for i in vec:
    norm += i*i
  return vec/np.sqrt(norm)

#---------------------Barycenter---------------------
def barycenter(geom):
  """barycenter(geom)
Returns the coordinates of the barycenter of a points distribution
"""
  bar = np.zeros(3)  
  nAt = 0
  for i in geom:
    nAt += 1
    bar = bar + i
  return bar/nAt

#---------------------Translation---------------------
def translation(Geomin,vector):
  """translation(Geomin,vector)
translates a set of cartesian points (Geomin) along a vector
"""
  # fill matrix with 0.0
  T = np.zeros((4,4))
  # set diagonal to 1
  for i in range(4):
    T[i][i] = 1.
  # set translation values
  for i in range(3):
    T[i][3] = vector[i]
  # translate
  Geomout = []
  for i in Geomin:
    j = i
    j.append(1.)
    Geomout.append(np.dot(T,j)[:3])
  return Geomout

#---------------------Quaternion to matrix---------------------
def quat2mat(quat):
  """quat2mat(quat):
Returns the rotation matrix associated to the quaternion quat.
quat can be a 1x4 list or array.
"""
  M = np.zeros((3,3))
  x = quat[0]
  y = quat[1]
  z = quat[2]
  w = quat[3]
  xx = x * x
  xy = x * y
  xz = x * z
  xw = x * w
  yy = y * y
  yz = y * z
  yw = y * w
  zz = z * z
  zw = z * w
  M[0][0]  = 1. - 2. * ( yy + zz )
  M[0][1]  = 2. * ( xy - zw )
  M[0][2]  = 2. * ( xz + yw )
  M[1][0]  = 2. * ( xy + zw )
  M[1][1]  = 1. - 2. * ( xx + zz )
  M[1][2]  = 2. * ( yz - xw )
  M[2][0]  = 2. * ( xz - yw )
  M[2][1]  = 2. * ( yz + xw )
  M[2][2] = 1. - 2. * ( xx + yy )
  return M
  
#---------------------Rotation---------------------
def rotation(point_list, x_rot,y_rot,z_rot,unit='d'):
  '''rotate(point_list, x_angle,y_angle,z_angle,unit):
Rotates a set of points around x, y, and z from specified angles (in specified unit (d or r), default: d) by the quaternion method.
Set of points must be a list of x,y,z cartesian coordinates.
Returns a list
  '''
  # manage input angles
  angles = np.array([x_rot,y_rot,z_rot])
  if unit == 'r':
    for i in angles:
      if i > 2*np.pi:
        sys.exit('Error: Are you sure angles are in radians?')
    pass
  elif unit == 'd':
    angles = deg2rad * angles
  else:
    sys.exit('Bad unit specification: should be d (degree) or r (radian)')
  
  # create input vectors list
  Geomin = []
  for i in point_list:
    Geomin.append(np.array(i))
  
  # generate rotation matrix
  c = np.cos(angles)
  s = np.sin(angles)
  Rx = np.array([[1,0,0],[0,c[0],s[0]],[0,-s[0],c[0]]])
  Ry = np.array([[c[1],0,s[1]],[0,1,0],[-s[1],0,c[1]]])
  Rz = np.array([[c[2],-s[2],0],[s[2],c[2],0],[0,0,1]])
  R = np.dot(Rx, np.dot(Ry, Rz))
  
   # apply the rotation
  Geomout = []
  for i in Geomin:
    Geomout.append(np.dot(R,i))
  
  geomout_list = []
  for i in Geomout: 
    geomout_list.append(list(i))
  return geomout_list

#---------------------Vectors alignment---------------------
def colinearize(point_list, vec1, vec2, debug=False):
  """colinearize(point_list, vec1, vec2):
Rotates a set of points so that vec1 becomes colinear to vec2.
Useful for example to aligne the dipole moment (vec1) of a system (point_list) on a direction (vec2).
Set of points must be a list of x,y,z cartesian coordinates.
Returns a list of new coordinates
  """
  # create input vectors list and format vectors
  Geomin = []
  for i in point_list:
    Geomin.append(np.array(i))
  V1 = normalize(np.array(vec1))
  V2 = normalize(np.array(vec2))
  if debug:
    print('Create and format vectors:\nV1:')
    print(V1)
    print('V2:')
    print(V2)
  
  # calculation of rotation axis and angle
  Axis = normalize(np.cross(V1,V2)) # cross product to determine rotation axis
  Angle = np.arccos(np.dot(V1,V2)) # angle between vectors
  if debug:
    print('\nRotation axis and angle:\nAxis:')
    print(Axis)
    print('Angle:')
    print(Angle,Angle*180./np.pi)
  
  # quaternions generation and normalisation
  Q_p = np.zeros(4)
  Q_m = np.zeros(4)
  sina = np.sin(Angle/2.)
  cosa = np.cos(Angle/2.)
  for i in range(3):
    Q_p[i] = Axis[i]*sina
  Q_p[3] = cosa
  Q_m = -1.*Q_p
  Q_m[3] = Q_p[3]
  Q_m = normalize(Q_m)
  Q_p = normalize(Q_p)
  if debug:
    print('\nQuaternions:\nQp:')
    print(Q_p)
    print('Qm:')
    print(Q_m)
  
  # conversion into rotation matrix
  Mp = np.zeros((3,3))
  Mm = np.zeros((3,3))
  Mp = quat2mat(Q_p)
  Mm = quat2mat(Q_m)
  if debug:
    print('\nMatrices:\nMp:')
    print(Mp)
    print('Mm:')
    print(Mm)
  
  # choice of good rotation matrix (for good angle)
  V2p = np.dot(Mp,V1)
  V2m = np.dot(Mm,V1)
  R = np.zeros((3,3))
  V2pc = sum(abs(V2-V2p))
  V2mc = sum(abs(V2-V2m))
  if V2pc > V2mc:
    R = Mm
    err = V2mc
  else:
    R = Mp
    err = V2pc
  if debug:
    print('\nV1:')
    print(V1)
    print('V2:')
    print(V2)
    print('\nV2p, diff:')
    print(V2p,V2pc)
    print('V2m, diff:')
    print(V2m,V2mc)
    print('\n chosen matrix:')
    print(R)
  
  # apply the rotation
  Geomout = []
  for i in Geomin:
    Geomout.append(np.dot(R,i))
  
  geomout_list = []
  for i in Geomout: 
    geomout_list.append(list(i))
  return geomout_list

#---------------------Vectors alignment---------------------
def splitangle(vec1, vec2, n, debug=False):
  """splitangle(vec1, vec2, n):
Generates n vectors between vec1 and vec2.
Returns a list of new coordinates
  """
  # create input vectors list and format vectors
  V1 = normalize(np.array(vec1))
  V2 = normalize(np.array(vec2))
  if debug:
    print('Create and format vectors:\nV1:')
    print(V1)
    print('V2:')
    print(V2)
  
  # calculation of rotation axis and angle
  Axis = normalize(np.cross(V1,V2)) # cross product to determine rotation axis
  Angle = np.arccos(np.dot(V1,V2)) # angle between vectors
  frac_Angle = Angle/float(n)
  if debug:
    print('\nRotation axis and angle:\nAxis:')
    print(Axis)
    print('Angle:')
    print(Angle,Angle*180./np.pi)
    print('Angle elementaire:')
    print(frac_Angle,frac_Angle*180./np.pi)
  
  # quaternions generation and normalisation
  Q_p = np.zeros(4)
  Q_m = np.zeros(4)
  sina = np.sin(frac_Angle/2.)
  cosa = np.cos(frac_Angle/2.)
  for i in range(3):
    Q_p[i] = Axis[i]*sina
  Q_p[3] = cosa
  Q_m = -1.*Q_p
  Q_m[3] = Q_p[3]
  Q_m = normalize(Q_m)
  Q_p = normalize(Q_p)
  if debug:
    print('\nQuaternions:\nQp:')
    print(Q_p)
    print('Qm:')
    print(Q_m)
  
  # conversion into rotation matrix
  Mp = np.zeros((3,3))
  Mm = np.zeros((3,3))
  Mp = quat2mat(Q_p)
  Mm = quat2mat(Q_m)
  if debug:
    print('\nMatrices:\nMp:')
    print(Mp)
    print('Mm:')
    print(Mm)
  
  # choice of good rotation matrix (for good angle)
  V2p = np.dot(Mp,V1)
  V2m = np.dot(Mm,V1)
  R = np.zeros((3,3))
  V2pc = sum(abs(V2-V2p))
  V2mc = sum(abs(V2-V2m))
  if V2pc > V2mc:
    R = Mm
    err = V2mc
  else:
    R = Mp
    err = V2pc
  if debug:
    print('\nV1:')
    print(V1)
    print('V2:')
    print(V2)
    print('\nV2p, diff:')
    print(V2p,V2pc)
    print('V2m, diff:')
    print(V2m,V2mc)
    print('\n chosen matrix:')
    print(R)
  
  # apply the rotation
  Geomout = []
  Geomout.append(vec1)
  for i in range(n):
    Geomout.append(np.dot(R,Geomout[i]))
  
  geomout_list = []
  for i in Geomout: 
    geomout_list.append(list(i))
  return geomout_list

#---------------------Tools---------------------

def isPtInCell(pl,v1l,v2l):
	"""isPtInCell(pl,v1l,v2l)
takes lists pl, v1l and v2l in the form [x,y] and returns true if p is in the cell
defined by v1 and v2
	"""
	# list to arrays
	v1=np.array(v1l)
	v2=np.array(v2l)
	p=np.array(pl)
	# column vectors
	vc1=v1[:,np.newaxis]
	vc2=v2[:,np.newaxis]
	vp=p[:,np.newaxis]
	# order vectors to get a positive surface
	mx1=np.concatenate((vc1,vc2),axis=1)
	mx2=np.concatenate((vc2,vc1),axis=1)
	if np.linalg.det(mx1)<0:
		vc3=vc2
		vc2=vc1
		vc1=vc3
	vtot=np.linalg.det(np.concatenate((vc1,vc2),axis=1))
	# prepares the translated point
	vp2=vp-(vc1+vc2)
	# computes the surfaces
	mx1=np.concatenate((vc1,vp),axis=1)
	vmx1=np.linalg.det(mx1)
	mx2=np.concatenate((vp,vc2),axis=1)
	vmx2=np.linalg.det(mx2)
	mx3=np.concatenate((vp2,vc2),axis=1)
	vmx3=np.linalg.det(mx3)
	mx4=np.concatenate((vc1,vp2),axis=1)
	vmx4=np.linalg.det(mx4)
	# if point not in the box, or translated point in the box
	if vmx1<0 or vmx2<0 or vmx3>=0 or vmx4>=0:
		return False
	else:
		return True

###################### Main ######################
if __name__ == '__main__':
  print("""
**************
Geometry tools
S.Nenon 2014
**************

Function list:

""")
# print(vecprod.__doc__)
  print(translation.__doc__)
  print(quat2mat.__doc__)
  print(barycenter.__doc__)
  print(rotation.__doc__)
  print(colinearize.__doc__)
  print(normalize.__doc__)
  print(isPtInCell.__doc__)







