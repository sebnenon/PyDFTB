#!/usr/bin/env python
#-*- coding: utf-8 -*-

import re

#========================================
# Permet de parser une chaine
#========================================
def parse(chaine):
  '''
  Parses a chain by eliminating all \t \n \r symbols
  Returns a list of words.
  Example:
    parse('This is\t a test sentence\n')
      returns
    ['This','is','a','test','sentence']
  '''
  regex = re.compile('\s')
  chaine = re.sub(regex,' ',chaine)
  while chaine.find('  ') !=  -1:
    chaine = chaine.replace('  ', ' ')
  chaine = chaine.lstrip(' ').rstrip(' ')
  return chaine.split()

if __name__ == '__main__':
  import sys
  if len(sys.argv) < 2:
    sys.exit('Nothing to parse'+parse.__doc__)
  else:
    print(parse(sys.argv[1:]))
