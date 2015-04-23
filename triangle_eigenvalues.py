# this python module uses numpy and a simple finite elements pde solver to find discretized eigenvalues and eigenfunctions of a prescribed euclidean triangle

from scipy import linalg, sparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

from operator import itemgetter
import os
import time

class Triangle:
	
	def __init__(self,defining_matrix):
		#print defining_matrix
		self.T = defining_matrix/np.sign(np.linalg.det(defining_matrix))
		#print self.T
		self.meshes = dict()
		self.eigs = dict()
		x1 = defining_matrix[0][0]
		x2 = defining_matrix[0][1]
		y1 = defining_matrix[1][0]
		y2 = defining_matrix[1][1]
		self.name = str(self.sides())

	# the namespace is organized as follows:
	# ~/triangles/[name]/[mesh fineness]/[eigenvalues.npy or eigenvectors.npy]
	# (at some future time, will add /[neumann or dirichlet] between fineness and eigenfiles)
	def ensureDirectory(self):
		# check if folders exist
		base_dir = 'triangles'
		if not os.path.exists(base_dir):
			os.makedirs(base_dir)
		this_dir = os.path.join(base_dir,self.name)
		if not os.path.exists(this_dir):
			os.makedirs(this_dir)
		for num_meshes in self.meshes:
			cur_mesh_dir = os.path.join(this_dir,str(num_meshes))
			if not os.path.exists(cur_mesh_dir):
				os.makedirs(cur_mesh_dir)
		return


	# save all currently loaded eigenvalues, eigenvectors to file
	def saveToFile(self):
		self.ensureDirectory()
		for num_mesh in self.meshes:
			mesh_dir = os.path.join('triangles',self.name,str(num_mesh))
			eigvalfile = open(mesh_dir + '/eigenvalues.npy', 'w+')
			eigvecfile = open(mesh_dir + '/eigenvectors.npy', 'w+')
			np.save(eigvalfile,self.eigs[num_mesh][0])
			np.save(eigvecfile,self.eigs[num_mesh][1])
			eigvalfile.close()
			eigvecfile.close()

	def sides(self):
		v0 = self.T[0,:]
		v1 = self.T[1,:]
		v2 = v1 - v0
		print v0,v1,v2
		l0 = np.sqrt(v0.dot(v0))
		l1 = np.sqrt(v1.dot(v1))
		l2 = np.sqrt(v2.dot(v2))
		print l0,l1,l2
		sidelist = [l0,l1,l2]
		sidelist.sort()
		print sidelist
		return sidelist
		
		

	# load a stuff from file
	def loadEigs(self,num):
		print num
		# the number of mesh points is the last pathname component of the current directory
		num_mesh = int(num)
		cur_dir = os.path.join('triangles',self.name,num)
		# if the directory exists and contains the desired files, load
		if os.path.exists(cur_dir) and 'eigenvalues.npy' in os.listdir(cur_dir) and 'eigenvectors.npy' in os.listdir(cur_dir):
			print 'files found!'
			print 'loading files ...'
			eigvalfile = open(cur_dir + '/eigenvalues.npy', 'r')
			eigvecfile = open(cur_dir + '/eigenvectors.npy', 'r')
			
			eigvals = np.load(eigvalfile)
			eigvecs = np.load(eigvecfile)
			
			eigvalfile.close()
			eigvecfile.close()
			
			self.addMesh(num_mesh)
			self.eigs[num_mesh] = (eigvals,eigvecs)

	def addMesh(self,num_part,dirichlet=False):
		self.meshes[num_part] = [createLapl(self.T,num_part+1,dirichlet), createL2InnerProduct(self.T,num_part+1,dirichlet)]


	def computeEig(self,num_part,dirichlet=False):
		if num_part in self.meshes:
			A,B = self.meshes[num_part]
		else:
			self.addMesh(num_part,dirichlet)
			A,B = self.meshes[num_part]
		#print 'A = ' + str(A)
		#print 'B = ' + str(B)
		#print 'computing eigenvalues of ' + str(A.shape) + ' generalized eigenvalue problem ... '
		start_time = time.time()
		l1,l2 = linalg.eigh(A,B)

		elapsed_time = time.time() - start_time
		print elapsed_time

		idx = l1.argsort()
		l1 = l1[idx]
		l2 = l2[:,idx]

		self.eigs[num_part] = (l1,l2)

	def returnEigs(self,num_part,n):
		if num_part not in self.eigs:
			self.computeEig(num_part)
		elist = list(self.eigs[num_part][0])
		trunc_eig = [x for x in elist if x < n]
		return trunc_eig

	def returnEigFunc(self,num_part,n):
		if num_part not in self.eigs:
			self.computeEig(num_part)
		X,Y = xyAxes(num_part,self.T)
		Z = self.eigs[num_part][1][:,n]
		return X,Y,Z
	
	def returnEigFuncFigure(self,num_part,n):
		if num_part not in self.eigs:
			self.computeEig(num_part)
		X,Y = xyAxes(num_part,self.T)
		Z = self.eigs[num_part][1][:,n]
		triang = tri.Triangulation(X,Y)
		return X,Y,Z,triang

	def plotMesh(self,num_part):
		if num_part not in self.meshes:
			self.addMesh(num_part)
		X = []
		Y = []
		for j in range(num_part+1):
			for k in range(j+1):
				v = self.T.dot(np.array([j-k,k]))/num_part
				X.append(v[0])
				Y.append(v[1])
		triang = tri.Triangulation(X,Y)
		plt.triplot(triang)
		plt.show()

	def plotEigvect(self,num_part,n,graph=False,contour=True):
		if num_part not in self.eigs:
			self.computeEig(num_part)
		ntheigpair = [self.eigs[num_part][0][n], self.eigs[num_part][1][:,n]]

		Z = ntheigpair[1]

		X,Y = xyAxes(num_part,self.T)

		triang = tri.Triangulation(X,Y)
		if graph:
			fig = plt.figure()
			ax = fig.add_subplot(111,projection='3d')
			ax.plot_trisurf(X,Y,Z, cmap=plt.get_cmap('copper'),linewidth=0.1)
			
		elif contour:
			plt.tricontourf(triang,Z,100,cmap=plt.get_cmap('copper'))
			#plt.colorbar()
		else:
			plt.tricontour(triang,Z,100,color='black')

		plt.show()

	def plotNodal(self,num_part,nlist,dimx,dimy):
		if num_part not in self.eigs:
			self.computeEig(num_part)

		m = len(nlist)

		for x in range(dimx):
			print x
			for y in range(dimy):
				s = x*dimx + y
				if s < m:
					print "  " + str(s)
					n = nlist[s]
					print "  " + str((x*dimx+1)*100 + (y+1)*10 + 1)
					plt.subplot((x*dimx+1)*100 + (y+1)*10 + 1)
					ntheigpair = [self.eigs[num_part][0][n], self.eigs[num_part][1][:,n]]

					Z = ntheigpair[1]

					X,Y = xyAxes(num_part,self.T)

					triang = tri.Triangulation(X,Y)

					plt.tricontour(triang,Z,0,color='black')

					oX = []
					oY = []
					for j in range(2):
						for k in range(j+1):
							v = self.T.dot(np.array([j-k,k]))
							oX.append(v[0])
							oY.append(v[1])
					otriang = tri.Triangulation(oX,oY)
					plt.triplot(otriang)

					plt.axis('equal')

		plt.show()

	def checkAgainstOnes(self,num_part):
		if num_part not in self.eigs:
			self.computeEig(num_part)
		A = self.meshes[num_part][0]
		o = np.ones( P(num_part+1) )
		print A.dot(o)

	# L2 norm of the nth eigenvector
	def l2norm(self,num_part,n):
		v = self.eigs[num_part][1][:,n]
		B = self.meshes[num_part][1]
		return v.transpose().dot(B.dot(v))
	
	def checkEigvect(self,num_part,n):
		v = self.eigs[num_part][1][:,n]
		e = self.eigs[num_part][0][n]
		A,B = self.meshes[num_part]
		print 'eigval: ' + str(e)
		return linalg.norm(A.dot(v) - e*B.dot(v))
		

# create or load a triangle from file
# input: 2x2 numpy array
def loadTriangle( tri ):
	# check if the triangle exists first!
	# don't want to do unnecessary work!
	x1 = tri[0][0]
	x2 = tri[0][1]
	y1 = tri[1][0]
	y2 = tri[1][1]
	# the name of the directory containing info, if it exists
	name = str(x1) + '-' + str(y1) + '_' + str(x2) + '-' + str(y2)
	# path to directory:
	tripath = 'triangles/' + name 
	print 'path to triangle:' + tripath

	# create triangle object to return
	T = Triangle( tri )

	# check if it exists
	if os.path.exists( tripath ):
		for p in os.listdir(tripath):
			print 'has: ' + p
			T.loadEigs(p)
	
	
	return T


def xyAxes(num_part,tri):
	X = []
	Y = []
	for j in range(num_part+1):
		for k in range(j+1):
			v = tri.dot(np.array([j-k,k]))/num_part
			X.append(v[0])
			Y.append(v[1])
	return X,Y

def createAdjacencyMatrix( tri, num_vert ):
	adjmat = np.array([[1.]])
	for i in range(num_vert):
		if i == num_vert - 1:
			adjmat = addValues(adjmat,tri,kind='adj',last=True)
		else:
			adjmat = addValues(adjmat,tri,kind='adj')
	return adjmat

def createL2InnerProduct( tri, num_vert, dirichlet=False ):
	
	A = np.linalg.det(tri)/2.
	
	n_tot = P(num_vert)
	
	C = A/(12.*(num_vert-1)**2)
	
	ipmat = np.zeros( (n_tot,n_tot) )
	
	# we go diagonal-by-diagonal
	# i indexes the diagonals
	for n in range( num_vert ):
		# always work with i+1 because the 0th diagonal is trivial and the (i+1)th = num_vert-th diagonal is the last one

		# the 0th diagonal is trivial
		if n > 0:
			# now we work down the diagonal
			# (this is nice because the ith diagonal has n+1 things in it)
			# note that we will not have any edges along the diagonal with first vertex at n, so we index edges by 0,...,n-1
			for i in range( n ):
				
				# handle the two triangles touching P(n)+i and P(n)+(i+1), i.e., this vertex and the next one along the diagonal
				# we add the relevant coeff matrix to laplmat

				# first, the element down and left
				# diagonal
				ipmat[P(n)+i][P(n)+i] += 2*C
				ipmat[P(n)+(i+1)][P(n)+(i+1)] += 2*C
				ipmat[P(n-1)+i][P(n-1)+i] += 2*C
				# off-diagonal
				ipmat[P(n)+i][P(n)+(i+1)] += C
				ipmat[P(n)+(i+1)][P(n)+i] += C
				ipmat[P(n)+i][P(n-1)+i] += C
				ipmat[P(n-1)+i][P(n)+i] += C
				ipmat[P(n)+(i+1)][P(n-1)+i] += C
				ipmat[P(n-1)+i][P(n)+(i+1)] += C
				
				
				# this will handle the 'exterior' elements, i.e., up and to the right; we want it to switch OFF when on the last diagonal
				if n+1 < num_vert:
					# diagonal
					ipmat[P(n)+i][P(n)+i] += 2*C
					ipmat[P(n)+(i+1)][P(n)+(i+1)] += 2*C
					ipmat[P(n+1)+(i+1)][P(n+1)+(i+1)] += 2*C
					# off-diagonal
					ipmat[P(n)+i][P(n)+(i+1)] += C
					ipmat[P(n)+(i+1)][P(n)+i] += C
					ipmat[P(n)+i][P(n+1)+(i+1)] += C
					ipmat[P(n+1)+(i+1)][P(n)+i] += C
					ipmat[P(n)+(i+1)][P(n+1)+(i+1)] += C
					ipmat[P(n+1)+(i+1)][P(n)+(i+1)] += C


	# if this is for a dirichlet problem, we need to set the boundary rows/columns to [ 0 0 ... 0 1 0 ... 0 ]
	if dirichlet:
		# the boundary indices are: 0, P(i), P(i-1), P(n-1) through P(n)
		for n in range( num_vert ):

			tmp = np.zeros( n_tot )
			
			# handle exterior diagonal
			if n+1 == num_vert:
				for k in range(P(n),P(n+1)):
					tmp[k] = 1.
					ipmat[k,:] = tmp
					ipmat[:,k] = tmp.transpose()

			else:
				# handle left side
				tmp[P(n)] = 1.
				ipmat[P(n),:] = tmp
				ipmat[:,P(n)] = tmp.transpose()

				# reset
				tmp[P(n)] = 0.

				# handle bottom
				tmp[P(n)+1] = 1.
				ipmat[P(n)+1,:] = tmp
				ipmat[:,P(n)+1] = tmp.transpose()			

			

	return ipmat

def createLapl( tri,num_vert, dirichlet=False ):
		
	# compute quantities for the current iteration
	D = (np.linalg.det(tri))**2
	A = np.linalg.det(tri)/2.

	u = tri[:,0]
	v = tri[:,1]
	w = v - u

	gugu = u.dot(u)/(4*A)
	gugv = -u.dot(v)/(4*A)
	gugw = -gugv - gugu
	gvgv = v.dot(v)/(4*A)
	gvgw = -gugv - gvgv
	gwgw = w.dot(w)/(4*A)

	coeff_matrix = np.array( [ [gwgw, gvgw, gugw], [gvgw, gvgv, gugv], [gugw, gugv, gugu] ] )
	
	n_tot = P(num_vert)	
	
	laplmat = np.zeros( (n_tot,n_tot) )
	
	# we go diagonal-by-diagonal
	# i indexes the diagonals
	for n in range( num_vert ):
		# always work with i+1 because the 0th diagonal is trivial and the (i+1)th = num_vert-th diagonal is the last one

		# the 0th diagonal is trivial
		if n > 0:
			# now we work down the diagonal
			# (this is nice because the ith diagonal has n+1 things in it)
			# note that we will not have any edges along the diagonal with first vertex at n, so we index edges by 0,...,n-1
			for i in range( n ):
				
				# handle the two triangles touching P(n)+i and P(n)+(i+1), i.e., this vertex and the next one along the diagonal
				# we add the relevant coeff matrix to laplmat

				# first, the element down and left
				laplmat[P(n)+i][P(n)+i] += gvgv
				laplmat[P(n)+(i+1)][P(n)+(i+1)] += gugu
				laplmat[P(n-1)+i][P(n-1)+i] += gwgw
				laplmat[P(n)+i][P(n)+(i+1)] += gugv
				laplmat[P(n)+(i+1)][P(n)+i] += gugv
				laplmat[P(n)+i][P(n-1)+i] += gvgw
				laplmat[P(n-1)+i][P(n)+i] += gvgw
				laplmat[P(n)+(i+1)][P(n-1)+i] += gugw
				laplmat[P(n-1)+i][P(n)+(i+1)] += gugw
				
				
				# this will handle the 'exterior' elements, i.e., up and to the right; we want it to switch OFF when on the last diagonal
				if n+1 < num_vert:
					laplmat[P(n)+i][P(n)+i] += gugu
					laplmat[P(n)+(i+1)][P(n)+(i+1)] += gvgv
					laplmat[P(n+1)+(i+1)][P(n+1)+(i+1)] += gwgw
					laplmat[P(n)+i][P(n)+(i+1)] += gugv
					laplmat[P(n)+(i+1)][P(n)+i] += gugv
					laplmat[P(n)+i][P(n+1)+(i+1)] += gugw
					laplmat[P(n+1)+(i+1)][P(n)+i] += gugw
					laplmat[P(n)+(i+1)][P(n+1)+(i+1)] += gvgw
					laplmat[P(n+1)+(i+1)][P(n)+(i+1)] += gvgw

	# if this is for a dirichlet problem, we need to set the boundary rows/columns to [ 0 0 ... 0 1 0 ... 0 ]
	if dirichlet:
		# the boundary indices are: 0, P(i), P(i-1), P(n-1) through P(n)
		for n in range( num_vert ):

			tmp = np.zeros( n_tot )
			
			# handle exterior diagonal
			if n+1 == num_vert:
				for k in range(P(n),P(n+1)):
					tmp[k] = 1.
					laplmat[k,:] = tmp
					laplmat[:,k] = tmp.transpose()

			else:
				# handle left side
				tmp[P(n)] = 1.
				laplmat[P(n),:] = tmp
				laplmat[:,P(n)] = tmp.transpose()

				# reset
				tmp[P(n)] = 0.

				# handle bottom
				tmp[P(n)+1] = 1.
				laplmat[P(n)+1,:] = tmp
				laplmat[:,P(n)+1] = tmp.transpose()			

			


	
	return laplmat

def P(n):
	return n*(n+1)/2


# this has been deprecated
def addValues( mat, tri, kind="adj", last=False ):

	x,y = mat.shape
	n = (np.sqrt(1+8*x)-1)/2
	newx = int(x+n+1)
	oldx = (n-1)*n/2

	new_adj_matrix = np.zeros( (newx,newx) )

	if kind == "adj":
		for i in range(x):
			for j in range(x
):
				new_adj_matrix[i][j] = mat[i][j]

		
		for k in range(x,newx):
			new_adj_matrix[k][k] = 1.
			
			if (k - x) > 0:
				new_adj_matrix[k][k-1] = 1.
			if (newx - k) > 1:
				new_adj_matrix[k][k+1] = 1.
		
			nind = oldx + k - x
			if k - x > 0:
				new_adj_matrix[k][nind-1] = 1.
				new_adj_matrix[nind-1][k] = 1.
			if nind < x:
				new_adj_matrix[k][nind] = 1.
				new_adj_matrix[nind][k] = 1.
		
		return new_adj_matrix
