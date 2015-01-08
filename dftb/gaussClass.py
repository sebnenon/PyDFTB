#!/usr/bin/env python
#-*- coding: utf-8 -*-

import sys
import numpy as np


#========================================
# Classe importee de dosplot 
# ((c)Balint Aradi 2006)
#========================================
class SmearingFunction(object):
  """Abstract class for the dos"""
  pass


class Gauss(SmearingFunction):
  """Gaussian smearing"""

  epsilon = 1e-8

  def __init__(self, mult,coefficient, exponent, center):
    self._coeff = coefficient
    self._exp = exponent
    self._center = center
    self._mult = mult
    tmp = np.log(self._coeff / self.epsilon)
    if tmp > 30:
      error("Coefficient for Gaussian too big (%f)" % self._exp)
    self._tol = tmp / self._exp
    self._exp *= -1

  def __call__(self, xx):
    tmp = (xx - self._center)**2
    if tmp < self._tol:
      return self._coeff *self._mult* np.exp(self._exp * tmp)
    else:
      return 0.0

if __name__ == '__main__':
  sys.exit('Smearing functions classes')

  
      

