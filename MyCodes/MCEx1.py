# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 14:18:15 2014

@author: meza
"""

#!/usr/bin/env python

'''
 Monte Carlo Program to calculate <r>, <r2>, <E> in harmonic potential
 Kalju Kahn, kalju@chem.ucsb.edu, UC Santa Barbara
 Help and Inspiration by Prof. Horia Metiu is hereby acknowledged

 Usage: MC_Harmonic.py param_file temperature num_points r_initial stepsize
 param_file is a plain text file with three entries: ID k r_eq.  For example:
 CO      2743.0  1.1283 

 Copyright 2005 Kalju Kahn
    
 Copying and distribution of this file, with or without modification,
 are permitted in any medium without royalty provided the copyright
 notice and this notice are preserved.
'''

import sys
from math import exp			# Math exp is much faster than SciPy exp
from random import random		# Mersienne Twister as the core generator

# Check if all the input is provided 
if len(sys.argv ) !=6 :
   sys.exit ( "Usage: MC_Harmonic.py param_file temperature num_points r_initial stepsize" )

#
# Potential parameters come from the external file
#
try:
   file = open( sys.argv[ 1 ] )
except IOError:
   print >> sys.stderr, "Parameter file could not be opened"
   sys.exit( 1 )

parameters = file.readlines()
for parameter in parameters:
    fields = parameter.split()
    mol  = (fields [0] )
    k    = float (fields [1] )  	# Harmonic force constant
    r_eq = float (fields [2] )  	# Potential minimum distance
#
# Monte Carlo run parameters come from the command line
#
temp    = float(sys.argv[ 2 ])			# Temperature
points  = int(sys.argv[ 3 ])    		# Number of MC attempts
r_ini   = float(sys.argv[ 4 ])  		# Initial distance in MC
step    = float(sys.argv[ 5 ])  		# Step size in MC

R      = 1.9872156e-3    # Gas constant  kcal/mol/degree
beta   = 1 / ( R * temp)
print '################################################'
print '#            HARMONIC POTENTIAL                #'
print '################################################'
print ' '
print 'Harmonic force constant k: ' ,k,    'kcal/mol/A**2'
print 'Potential minimum:         ' ,r_eq, 'Ang'
print ' '
#
# Hamiltonian: Harmonic potential
#
def hamiltonian(r):
    pot = 0.5 * k * (r - r_eq)**2 
    return pot
#
# Monte Carlo code	
#
acc = 0
mc_dst = 0.0
mc_dst2 = 0.0
mc_ener = 0.0

r = r_ini

for point in range(1,points):
   r_new = r + step * ( random() - 0.5 ) * 2.0
   v     = hamiltonian(r)
   v_new = hamiltonian(r_new)
   if v_new < v:				# Downhill move always accepted
      r = r_new
      v = v_new
      acc = acc + 1
   else:					# Uphill move accepted with luck
      A_mov = exp( -beta *(v_new - v) )
      if A_mov > random():
         r = r_new
         v = v_new
         acc = acc + 1
						# Update regardless of acceptance!
   mc_dst  = mc_dst + r
   mc_dst2 = mc_dst2 + r*r
   mc_ener = mc_ener + v

#
# Print out averages
#
mc_dst_av  = mc_dst / float(points)
mc_dst2_av = mc_dst2 / float(points)
mc_ener_av = mc_ener / float(points)

print 'Acceptance is:', 100 * acc / float(points), '%'
print ' '
print 'Monte Carlo Averages:'
print 'Mean distance, Ang:                 ', mc_dst_av
print 'Mean squared distance, Ang**2       ', mc_dst2_av
print 'Mean MC square displacemet, Ang**2  ', mc_dst2_av - mc_dst_av*mc_dst_av
print 'Mean MC potential energy, kcal/mol: ', mc_ener_av
