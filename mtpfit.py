#!/usr/bin/env python3

import sys
import numpy as np

class GaussianCube():
    '''
    GaussianCube Class:
    Includes a bunch of methods to manipulate cube data
    '''
    def __init__(self,fname=None):
        if fname != None:
            try:
                self.read_cube(fname)
            except IOError as e:
                print( "File used as input: %s" % fname )
                print( "File error ({0}): {1}".format(e.errno, e.strerror))
                self.terminate_code()
        else:
            self.default_values()
        return None

    def terminate_code(self):
        print( "Code terminating now")
        exit()
        return None

    def default_values(self):
        self.natoms=0
        self.comment1=0
        self.comment2=0
        self.origin=np.array([0,0,0])
        self.NX=0
        self.NY=0
        self.NZ=0
        self.X=0
        self.Y=0
        self.Z=0
        self.atoms=['0']
        self.atomsXYZ=[0,0,0]
        self.data=[0]
        return None

    def read_cube(self,fname):
        """
        Method to read cube file. Just needs the filename
        """

        with open(fname, 'r') as fin:
            self.filename = fname
            self.comment1 = fin.readline() #Save 1st comment
            self.comment2 = fin.readline() #Save 2nd comment
            nOrigin = fin.readline().split() # Number of Atoms and Origin
            self.natoms = int(nOrigin[0]) #Number of Atoms
            self.origin = np.array([float(nOrigin[1]),float(nOrigin[2]),float(nOrigin[3])]) #Position of Origin
            nVoxel = fin.readline().split() #Number of Voxels
            self.NX = int(nVoxel[0])
            self.X = np.array([float(nVoxel[1]),float(nVoxel[2]),float(nVoxel[3])])
            nVoxel = fin.readline().split() #
            self.NY = int(nVoxel[0])
            self.Y = np.array([float(nVoxel[1]),float(nVoxel[2]),float(nVoxel[3])])
            nVoxel = fin.readline().split() #
            self.NZ = int(nVoxel[0])
            self.Z = np.array([float(nVoxel[1]),float(nVoxel[2]),float(nVoxel[3])])
            self.atoms = []
            self.atomsXYZ = []
            for atom in range(self.natoms):
                line= fin.readline().split()
                self.atoms.append(line[0])
                self.atomsXYZ.append(list(map(float,[line[2], line[3], line[4]])))
            self.data = np.zeros((self.NX,self.NY,self.NZ))
            i= int(0)
            for s in fin:
                for v in s.split():
                    self.data[int(i/(self.NY*self.NZ)), int((i/self.NZ)%self.NY), int(i%self.NZ)] = float(v)
                    i+=1
            # if i != self.NX*self.NY*self.NZ: raise NameError, "FSCK!"
        return None

    def write_cube(self,fname,comment='written cube'):
        '''
        Write out a Gaussian Cube file
        '''
        try:
            with open(fname,'w') as fout:
                if len(comment.split('\n')) != 2:
                    print( 'Comment line NEEDS to be two lines!')
                    self.terminate_code()
                fout.write('%s\n' % comment)
                fout.write("%4d %.6f %.6f %.6f\n" % (self.natoms, self.origin[0], self.origin[1], self.origin[2]))
                fout.write("%4d %.6f %.6f %.6f\n" % (self.NX, self.X[0], self.X[1], self.X[2]))
                fout.write("%4d %.6f %.6f %.6f\n" % (self.NY, self.Y[0], self.Y[1], self.Y[2]))
                fout.write("%4d %.6f %.6f %.6f\n" % (self.NZ, self.Z[0], self.Z[1], self.Z[2]))
                for atom,xyz in zip(self.atoms,self.atomsXYZ):
                    fout.write("%s %d %6.3f %6.3f %6.3f\n" % (atom, 0, xyz[0], xyz[1], xyz[2]))
                for ix in range(self.NX):
                   for iy in range(self.NY):
                       for iz in range(self.NZ):
                           fout.write("%.5e " % self.data[ix,iy,iz]),
                           if (iz % 6 == 5): fout.write('\n')
                       fout.write("\n")
        except IOError as e:
            print( "File used as output does not work: %s" % fname)
            print( "File error ({0}): {1}".format(e.errno, e.strerror))
            self.terminate_code()
        return None

    def to_grid_vals(self):
        grid_xyz = np.zeros((self.data.size,3))
        grid_val = np.zeros((self.data.size))
        i = 0
        for ix in range(self.NX):
            for iy in range(self.NY):
                for iz in range(self.NZ):
                    grid_xyz[i,:] = self.origin + ix*self.X + iy*self.Y + iz*self.Z
                    grid_val[i]   = self.data[ix,iy,iz]
                    i += 1
        return grid_xyz, grid_val


lam  = 1e-6 #regularization term (keeps coefficients small), 0 for no regularization
hartreetokcal = 627.509
pot_cube_file = ''
dens_cube_file = ''
mtpl_file = "fitted-mtpl.dat"
lmax = 5
qtot = 0.0
  
f = open(mtpl_file,'w')  

#############
# Read command line input

def usage():
    print ("Usage: python3 fit.py -pot [potcube] -dens [denscube] [-lmax [lmax]]",\
          "[-qtot [qtot]]")

for i in range(len(sys.argv)):
  if sys.argv[i] == '-pot':
    pot_cube_file = sys.argv[i+1]
  elif sys.argv[i] == '-dens':
    dens_cube_file = sys.argv[i+1]
  elif sys.argv[i] == '-lmax':
    lmax = int(sys.argv[i+1])
  elif sys.argv[i] == '-qtot':
    qtot = float(sys.argv[i+1])
  elif sys.argv[i] == '-h':
    print("Usage: python fit.py -pun [file] [-par [parfile]] [-h]")


if pot_cube_file == '' or dens_cube_file == '':
  usage()
  raise Exception('Incorrect arguments to mtpfit.py') 

#read cube files
dens_cube = GaussianCube(dens_cube_file)
pot_cube = GaussianCube(pot_cube_file)

#extract grid values
dens_xyz, dens_val = dens_cube.to_grid_vals()
esp_xyz, esp_val = pot_cube.to_grid_vals()

#atom coordinates
atoms_xyz = np.asarray(pot_cube.atomsXYZ)

#number of atoms
Natom = atoms_xyz.shape[0]

#number of spherical harmonics coefficients
L = 0
for l in range(lmax+1):
    L += 2*l+1

def extract_grid_vals(xyz, dens, esp, lower_bound=1.0e-3, upper_bound=3.162277660e-4):
    idx = np.arange(dens.size)
    idx = idx[dens < lower_bound]
    xyz  = xyz[idx]
    dens = dens[idx]
    esp  = esp[idx]
    idx = np.arange(dens.size)
    idx = idx[dens > upper_bound]
    xyz  = xyz[idx]
    dens = dens[idx]
    esp  = esp[idx]
    return xyz, esp

grid_xyz, grid_val = extract_grid_vals(dens_xyz, dens_val, esp_val)

#gridsize
Ngrid = grid_val.size

'''print grid as xyz
print(str(grid_val.size)+"\n")
for xyz, val in zip(grid_xyz, grid_val):
    if val < 0:
        print("O " + str(xyz[0]) + " " + str(xyz[1]) + " " + str(xyz[2]) + " " + str(val))
    else:
        print("N " + str(xyz[0]) + " " + str(xyz[1]) + " " + str(xyz[2]) + " " + str(val))
quit()
#'''

'''
l,m indices of spherical harmonics
cartesian vector between multipole center and evaluation point
The functions are taken from the appendix of "The Theory of Intermolecular Forces" by A. J. Stone
'''
def multipole_esp(l, m, rvec):
    x = 0
    y = 1
    z = 2 
    rmag = np.linalg.norm(rvec)
    r = rvec/rmag
    #appropriate power law
    Rpow = 1/rmag**(l+1)
    if l == 0: #monopole
        return Rpow
    elif l == 1: #dipole
        if m == 0: #10
            return Rpow * r[z]
        elif m == 1: #11c
            return Rpow * r[x]
        elif m == -1: #11s
            return Rpow * r[y]
        else: 
            print('m =', m, ' is not supported for l =', l)
            quit()
    elif l == 2: #quadrupole
        if m == 0: #20
            return Rpow * 0.5 * (3*r[z]**2 - 1)
        elif m == 1: #21c
            return Rpow * np.sqrt(3) * r[x] * r[z]
        elif m == -1: #21s
            return Rpow * np.sqrt(3) * r[y] * r[z]
        elif m == 2: #22c
            return Rpow * 0.5 * np.sqrt(3) * (r[x]**2 - r[y]**2)
        elif m == -2: #22s
            return Rpow * np.sqrt(3) * r[x] * r[y]
        else: 
            print('m =', m, ' is not supported for l =', l)
            quit()
    elif l == 3: #octopole
        if m == 0: #30
            return Rpow * 0.5 * (5*r[z]**3 - 3*r[z])
        elif m == 1: #31c
            return Rpow * 0.25 * np.sqrt(6) * r[x] * (5*r[z]**2 - 1)
        elif m == -1: #31s
            return Rpow * 0.25 * np.sqrt(6) * r[y] * (5*r[z]**2 - 1)
        elif m == 2: #32c
            return Rpow * 0.5 * np.sqrt(15) * r[z] * (r[x]**2 - r[y]**2)
        elif m == -2: #32s
            return Rpow * np.sqrt(15) * r[x] * r[y] * r[z]
        elif m == 3: #33c
            return Rpow * 0.25 * np.sqrt(10) * r[x] * (r[x]**2 - 3*r[y]**2)
        elif m == -3: #33s
            return Rpow * 0.25 * np.sqrt(10) * r[y] * (3*r[x]**2 - r[y]**2)
        else: 
            print('m =', m, ' is not supported for l =', l)
            quit()
    elif l == 4: #hexadecapole
        if m == 0: #40
            return Rpow * 0.125 * (35*r[z]**4 - 30*r[z]**2 + 3)
        elif m == 1: #41c
            return Rpow * 0.25 * np.sqrt(10) * (7*r[x]*r[z]**3 - 3*r[x]*r[z])
        elif m == -1: #41s
            return Rpow * 0.25 * np.sqrt(10) * (7*r[y]*r[z]**3 - 3*r[y]*r[z]) 
        elif m == 2: #42c
            return Rpow * 0.25 * np.sqrt(5) * (7*r[z]**2 - 1) * (r[x]**2 - r[y]**2)
        elif m == -2: #42s
            return Rpow * 0.5 * np.sqrt(5) * (7*r[z]**2 - 1) * r[x] * r[y]
        elif m == 3: #43c
            return Rpow * 0.25 * np.sqrt(70) * r[x] * r[z] * (r[x]**2 - 3*r[y]**2)
        elif m == -3: #43s
            return Rpow * 0.25 * np.sqrt(70) * r[y] * r[z] * (3*r[x]**2 - r[y]**2)
        elif m == 4: #44c
            return Rpow * 0.125 * np.sqrt(35) * (r[x]**4 - 6*r[x]**2*r[y]**2 + r[y]**4)
        elif m == -4: #44s
            return Rpow * 0.5 * np.sqrt(35) * r[x] * r[y] * (r[x]**2 - r[y]**2)
        else: 
            print('m =', m, ' is not supported for l =', l)
            quit()
    elif l == 5: #ditriantapole
        if m == 0: #50
            return Rpow * 0.125 * (63*r[z]**5 - 70*r[z]**3 + 15*r[z])
        elif m == 1: #51c
            return Rpow * 0.125 * np.sqrt(15) * (21*r[x]*r[z]**4 - 14*r[x]*r[z]**2 + r[x])
        elif m == -1: #51s
            return Rpow * 0.125 * np.sqrt(15) * (21*r[y]*r[z]**4 - 14*r[y]*r[z]**2 + r[y])
        elif m == 2: #52c
            return Rpow * 0.25 * np.sqrt(105) * (3*r[x]**2*r[z]**3 - 3*r[y]**2*r[z]**3 - r[x]**2*r[z] + r[y]**2*r[z])
        elif m == -2: #52s
            return Rpow * 0.5 * np.sqrt(105) * (3*r[x]*r[y]*r[z]**3 - r[x]*r[y]*r[z])
        elif m == 3: #53c
            return Rpow * 0.0625 * np.sqrt(70) * (9*r[x]**3*r[z]**2 - 27*r[x]*r[y]**2*r[z]**2 - r[x]**3 + 3*r[x]*r[y]**2)
        elif m == -3: #53s
            return Rpow * 0.0625 * np.sqrt(70) * (27*r[x]**2*r[y]*r[z]**2 - 9*r[y]**3*r[z]**2 - 3*r[x]**2*r[y] + r[y]**3)
        elif m == 4: #54c
            return Rpow * 0.375 * np.sqrt(35) * (r[x]**4*r[z] - 6*r[x]**2*r[y]**2*r[z] + r[y]**4*r[z])
        elif m == -4: #54s
            return Rpow * 1.5 * np.sqrt(35) * (r[x]**3*r[y]*r[z] - r[x]*r[y]**3*r[z])
        elif m == 5: #55c
            return Rpow * 0.1875 * np.sqrt(14) * (r[x]**5 - 10*r[x]**3*r[y]**2 + 5*r[x]*r[y]**4)
        elif m == -5: #55s
            return Rpow * 0.1875 * np.sqrt(14) * (5*r[x]**4*r[y] - 10*r[x]**2*r[y]**3 + r[y]**5)
        else: 
            print('m =', m, ' is not supported for l =', l)
            quit()
    else:
        print('l =', l, ' is not supported')
        quit()

#build design matrix
A = np.zeros((Ngrid, Natom*L))
for n in range(Ngrid):
    i = 0
    for l in range(lmax+1):
        for a in range(Natom):
            r = grid_xyz[n] - atoms_xyz[a]
            for m in range(-l,l+1):
                A[n,i] = multipole_esp(l, m, r)
                i += 1

#build constraint matrix (so total charge is conserved)
B = np.zeros((1,Natom*L))
B[:,:Natom] = 1
d = np.zeros((1))
d[:] = qtot

def lse(A, b, B, d, lam=0.0):
    """
    Equality-contrained least squares.
    The following algorithm minimizes ||Ax - b|| subject to the
    constrain Bx = d.
    Parameters
    ----------
    A : array-like, shape=[m, n]
    B : array-like, shape=[p, n]
    b : array-like, shape=[m]
    d : array-like, shape=[p]
    lam : regularization
    Reference
    ---------
    Matrix Computations, Golub & van Loan, algorithm 12.1.2
    Examples
    --------
    >>> A = np.array([[0, 1], [2, 3], [3, 4.5]])
    >>> b = np.array([1, 1])
    >>> # equality constrain: ||x|| = 1.
    >>> B = np.ones((1, 3))
    >>> d = np.ones(1)
    >>> lse(A.T, b, B, d)
    array([-0.5,  3.5, -2. ])
    """
    from scipy import linalg
    A, b, B, d = map(np.asanyarray, (A, b, B, d))
    p = B.shape[0]
    Q, R = linalg.qr(B.T)
    y = linalg.solve_triangular(R[:p, :p].T, d)
    A = np.dot(A, Q)
    if lam == 0: #unregularized
        z = linalg.lstsq(A[:, p:], b - np.dot(A[:, :p], y))[0].ravel()
    else: #regularized
        z = linalg.lstsq(A[:, p:].T.dot(A[:, p:]) + lam*np.eye(A[:,p:].shape[1]), A[:, p:].T.dot(b - np.dot(A[:, :p], y)))[0].ravel()
    return np.dot(Q[:, :p], y) + np.dot(Q[:, p:], z)


x = lse(A, grid_val, B, d, lam)

#write results to mtpl_file (passed as argument above)
i = 0
rmse = np.sqrt(np.mean((np.matmul(A,x)-grid_val)**2))
maxe = np.amax(np.absolute(np.matmul(A,x)-grid_val))
maxmep = np.amax(np.absolute(grid_val)) * hartreetokcal
meanabsmep = np.mean(np.absolute(grid_val)) * hartreetokcal
print("Max MEP = "+str(maxmep)+"\n")
print("Mean abs MEP = "+str(meanabsmep)+"\n")
f.write("Qtot" + str(np.sum(x[:Natom])) + '\n')
f.write("# RMSE: " + str(rmse*hartreetokcal) + ", Max. Err: "+str(maxe*hartreetokcal)+" kcal/mol Qa(l,m) [a = atom index]\n")
f.write(str(lmax) + '\n')
for l in range(lmax+1):
    for a in range(Natom):
        for m in range(-l,l+1):
            f.write('Q'+str(a+1)+'('+str(l)+','+str(m)+')' + '     ' + str(x[i]) + '\n')
            i += 1


