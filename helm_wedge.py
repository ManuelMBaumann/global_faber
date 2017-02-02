from __future__ import print_function, division

import matplotlib
matplotlib.use('agg')

from nutils import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import mmwrite


@log.title
def makeplots( domain, geom, Lx, Lz, value, name, title, ndigits=0, index=None, clim=None, lineOn=False):
    
  points, colors = domain.elem_eval( [ geom, value ], ischeme='bezier3', separate=True )
  
  with plot.PyPlot( name, ndigits=ndigits, figsize=(5,6), index=index ) as plt:
    plt.mesh( points, colors, triangulate='bezier', edgecolors='none' )
    plt.title(title)
    plt.xlabel('x [m]')
    plt.ylabel('z [m]')
  
    plt.xticks([0, Lx/2.0, Lx], ['0', '300', '600'])
    plt.yticks([-0, -0.4*Lz, -0.8*Lz, -Lz], ['0', '400', '800', '1000'])
  
    if clim is not None:
      plt.clim(*clim)
    plt.colorbar()
    
    if lineOn:
      # Only for wedge problem
      plt.plot([0, 600],[-400, -500],'k')
      plt.plot([0, 600],[-800, -600],'k')
            
    
def makevtk(domain, geom, c, sol, p_basis, name):
  
  Nom = sol.shape[0]
  vtk_geom, vtk_c = domain.simplex.elem_eval( [ geom, c ], ischeme='vtk', separate=True )
  with plot.VTKFile( name ) as vtk:
      vtk.unstructuredgrid( vtk_geom )
      vtk.pointdataarray( 'c', vtk_c )
      for i in range(0,Nom):
          pres = p_basis.dot( sol[i,:] ).real
          vtk_pres = domain.simplex.elem_eval( pres, ischeme='vtk', separate=True )
          vtk.pointdataarray( 'pres_'+str(i), vtk_pres )


def makespyplot( matrix, name, imgtype=None ):

  if not scipy.sparse.isspmatrix( matrix ):
      matrix = matrix.toscipy()

  with plot.PyPlot( name, ndigits=0, imgtype=imgtype ) as plt:
    plt.spy( matrix, markersize=0.8, color='black')
    plt.title( name+', nnz = '+str(matrix.nnz) )

      
def point_eval(func, domain, geom, point):
  domain = domain[tuple(slice(0, p) if p > 0 else slice(None) for p in point)]
  for p in point:
    domain = domain.boundary['right' if p > 0 else 'left']
  return numpy.asarray(domain.integrate( func, geometry=geom, ischeme='gauss2' ))

  
def helm_mat(c, ndims, nx, ny, nz, p_basis, domain, geom):

  # define Helmholtz eqn.
  laplace    = function.outer( p_basis.grad(geom), p_basis.grad(geom)  ).sum(-1)
  mass       = function.outer( p_basis, p_basis ) / c**2
  sommerfeld = function.outer( p_basis, p_basis ) / c
    
  # Build matrices
  if ndims == 2:
      sommerfeld_boundary = 'left,right,bottom,top'
      source_position = nx//2, nz
  else:
      sommerfeld_boundary = 'left,right,bottom,top,front,back'
      source_position = nx//2, ny//2, nz
      
  K, M = domain.integrate( [laplace, mass], geometry=geom, ischeme='gauss2' )
  C = domain.boundary[sommerfeld_boundary].integrate( sommerfeld, geometry=geom, ischeme='gauss2' )
    
  # Define RHS
  rhs = point_eval(p_basis, domain, geom, source_position)
 
  return K, C, M, rhs

def main( ndims=2, degree=1, dx = 10.0, dy = 10.0, dz = 10.0, 
          freq=[16.0,32.0], storing=False, plots=True):
      
  # domain size
  Lx = 600.0
  Ly = 600.0
  Lz = 1000.0
  
  # parameters
  freq = np.array(freq)
  om   = 2.0*np.pi*freq
  Nom  = len(om)
  
  # define physical params:
  c0 = 2000.0
  c1 = 1500.0
  c2 = 3000.0
  
  if ndims == 2:
      nx = int(np.round(Lx/dx))+1
      ny = 1
      nz = int(np.round(Lz/dz))+1
      verts_x = np.linspace( 0, Lx, nx )
      verts_z = np.linspace( -Lz, 0, nz )
      verts = [verts_x, verts_z]
  else:
      nx = int(np.round(Lx/dx))+1
      ny = int(np.round(Ly/dy))+1
      nz = int(np.round(Lz/dz))+1
      verts_x = np.linspace( 0, Lx, nx )
      verts_y = np.linspace( 0, Ly, ny )
      verts_z = np.linspace( -Lz, 0, nz )
      verts = [verts_x, verts_y, verts_z]

  domain, geom = mesh.rectilinear(verts)
  p_basis      = domain.splinefunc( degree=degree )
  nelems       = int( len(p_basis)/ndims )
  l_p_basis    = len( p_basis )  
  
  print( 'problem size: '+str(nx-1+degree)+' x '+str(ny-1+degree)+' x '+str(nz-1+degree) )
  print( 'dofs        : '+str(l_p_basis) )
  
  # wedge problem:
  c = function.select(
      [function.greater(geom[-1]+400+geom[0]/6, 0), function.greater(geom[-1]+800-geom[0]/3, 0)],
      [c0, c1], c2)
  
  # Create discretization matrices with nutils
  K, C, M, rhs = helm_mat(c, ndims, nx, ny, nz, p_basis, domain, geom)

  if storing:
      mmwrite('matlab_io/K.mtx', K.toscipy())
      mmwrite('matlab_io/C.mtx', C.toscipy())
      mmwrite('matlab_io/M.mtx', M.toscipy())
      np.savetxt('matlab_io/b.txt', rhs)
     
  #Solve in python
  print('Use python solver')
  sol = np.zeros((Nom, l_p_basis), dtype=complex)
  for k in range(0,Nom):
      matrix = K + 1j*om[k]*C - om[k]**2*M
      sol[k,:] = scipy.sparse.linalg.spsolve( matrix.toscipy().tocsc(), rhs)

      
  if plots:
      if ndims==2:
          makeplots( domain, geom, Lx, Lz, c, 'c', 'Velocity profile [m/s]' )
          for k in range(0,Nom):
              pres = p_basis.dot( sol[k,:] ).real
              makeplots( domain, geom, Lx, Lz, pres, 'pres'+str(k), 'Pressure at {} Hz'.format(freq[k]), lineOn=True )
      else:
          makevtk(domain, geom, c, sol, p_basis, 'helm_wedge3d')
        
util.run( main )
