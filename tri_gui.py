import numpy as np
from scipy.special import *
from triangle_eigenvalues import *
from Tkinter import *

import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

import matplotlib.pyplot as plt
import matplotlib.tri as tri

class TriangleGui:
	
	def __init__(self,root,canvas_width,canvas_height,base_refine,eiglevel):
		self.root=root
		self.canvas_width = canvas_width
		self.canvas_height = canvas_height
		self.base_refine = base_refine
		self.current_refine = base_refine
		self.eiglevel = eiglevel
		self.nodes = False

		# set up the first triangle
		self.current_tri = Triangle( [[0.5,0.],[0.5,np.sqrt(3.)/2.]] )
		self.current_tri.computeEig(self.current_refine)


		# now we build and pack the state data
		self.info_frame = Frame(master=root, width=50, height=canvas_height)
		self.info_frame.pack(side=BOTTOM)

		self.info_tri = Text(self.info_frame,height=3,width=50)
		self.info_tri.pack(side=TOP)
		
		self.info_refine = Text(self.info_frame,height=2,width=50)
		self.info_refine.pack(side=RIGHT)
		self.info_eiglevel = Text(self.info_frame,height=2,width=50)
		self.info_eiglevel.pack(side=LEFT)

		self.update_info()

		# build and pack other controls
		self.control_frame = Frame(master=root,width=canvas_width/4,height=canvas_height)
		
		self.control_frame.pack(side=TOP)

		self.b = Button(self.control_frame, text="QUIT", command = self.quit)
		self.b.pack(side=LEFT)

		self.c = Button(self.control_frame, text="Toggle Nodal", command=self.toggle)
		self.c.pack(side=LEFT)
		
		self.s = Button(self.control_frame, text="Save Image", command=self.save)
		self.s.pack(side=LEFT)
		
		self.eq = Button(self.control_frame, text="Standard Equilateral", command=self.standard_eq)
		self.eq.pack(side=LEFT)
		
		self.rt = Button(self.control_frame, text="Standard Right", command=self.standard_rt)
		self.rt.pack(side=LEFT)

		# set up the matplotlib nonsense
		self.fig, self.ax = mpl.pyplot.subplots()

		self.ax.clear()
		self.ax.autoscale(enable=False)
		self.ax.axis('off')
		
		X,Y,Z,triang = self.current_tri.returnEigFuncFigure(self.current_refine,self.eiglevel)
		triang = tri.Triangulation(X,Y)
		self.ax.tricontourf(triang,Z,100,cmap=plt.get_cmap('copper'))

		self.canvas = FigureCanvasTkAgg(self.fig,master=root)
		self.canvas.show()
		self.canvas.get_tk_widget().pack()

		self.canvas.mpl_connect("button_press_event", self.setVertex)
		self.canvas.mpl_connect("key_press_event", self.keyEvent)

		self.ax.set_xlim([0,1])
		self.ax.set_ylim([0,1])
		self.ax.autoscale(enable=False)

		self.fig.tight_layout()

		self.redraw()




	# command/control methods
	def save(self):
		self.ensureDirectory()
		if self.nodes:
			self.fig.savefig("pictures/("+self.current_tri.name+")-("+str(self.eiglevel)+","+str(self.current_refine)+")-"+"nodal.png")
		else:
			self.fig.savefig("pictures/("+self.current_tri.name+")-("+str(self.eiglevel)+","+str(self.current_refine)+").png")

	def quit(self):
		sys.exit()

	def toggle(self):
		self.nodes = not self.nodes
		self.redraw()

	def update_info(self):
		self.info_refine.delete(1.0,END)
		self.info_refine.insert(END,"Current side partition number: " + str(self.current_refine) + "\nCurrent mesh size: " + str(self.current_refine*(self.current_refine+1)/2))
		self.info_eiglevel.delete(1.0,END)
		self.info_eiglevel.insert(END,"Current index: " + str(self.eiglevel) + "\n  Eigenvalue: " + str(self.current_tri.eigs[self.current_refine][0][self.eiglevel]))
		self.info_tri.delete(1.0,END)
		self.info_tri.insert(END,'Current triangle:\n' + str(self.current_tri.T))

	def standard_eq(self):
		self.current_tri = Triangle( np.array([[1.0,0.5],[0.0,np.sqrt(3.)/2.]]) )
		self.redraw()

	def standard_rt(self):
		self.current_tri = Triangle( np.array([[0.6,0.],[0.,0.6]]) )
		self.redraw()

	# triangle manipulation methods
	def refine(self,event):
		if self.current_refine < :
			self.current_refine += 5
		self.redraw()

	def coarsen(self,event):
		if self.current_refine > 5:
			self.current_refine -= 5
		self.redraw()

	def energy_up(self,event):
		self.eiglevel += 1
		self.redraw()
	
	def energy_down(self,event):
		if self.eiglevel > 1:
			self.eiglevel -= 1
		self.redraw()

	def replot(self,event):
		self.redraw()
		
	def redraw(self):
		self.ax.clear()
		self.ax.autoscale(enable=False)
		
		X,Y,Z,triang = self.current_tri.returnEigFuncFigure(self.current_refine,self.eiglevel)
		triang = tri.Triangulation(X,Y)
		
		if Z[0] < 0:
			Z = list((-1)*np.array(Z))
		
		if self.nodes:
			mat = self.current_tri.T
			X = list(mat[0,:])
			Y = list(mat[1,:])
			X.append(0)
			Y.append(0)
			outline = tri.Triangulation(X,Y)
			self.ax.triplot(outline)
			self.ax.tricontour(triang,Z,0)
		else:
			self.ax.tricontourf(triang,Z,100,cmap=plt.get_cmap('copper'))

		self.canvas.draw()
		
		self.update_info()

	def keyEvent(self,event):
		if event.key == 'r':
			self.refine(event)
		if event.key == 'c':
			self.coarsen(event)
		if event.key == 'up':
			self.energy_up(event)
		if event.key == 'down':
			self.energy_down(event)
		if event.key == ' ':
			self.redraw()

	def setVertex(self,event):
		if event.button == 1:
			self.setNewX(event)
		if event.button == 3:
			self.setNewY(event)

	def setNewX(self,event):
		#print 'click!'
		point = (event.xdata,event.ydata)
		#print point
		v0,v1 = point
		mat = np.copy(self.current_tri.T)
		mat[:,0] = v0,v1
		if np.linalg.det(mat) > 0:
			self.current_tri = Triangle(mat)
			self.redraw()

	def setNewY(self,event):
		point = (event.xdata,event.ydata)
		#print point
		v0,v1 = point
		mat = np.copy(self.current_tri.T)
		mat[:,1] = v0,v1
		if np.linalg.det(mat) > 0:
			self.current_tri = Triangle(mat)
			self.redraw()

	# logic support methods
	def toCanvasCoords(self,point):
		#print 'converting!'
		return (self.canvas_width*point[0],self.canvas_height - self.canvas_height*point[1])

	def fromCanvasCoords(self,point):
		#print 'converting!'
		return (float(point[0])/self.canvas_width, 1-float(point[1])/self.canvas_height)

	def colormap(self,r):
		s = ''
		n = int(1+98*1./(1+np.exp(-r)))
		#print n
		s += 'gray' + str(n)
		#print s
		return s
		
	def ensureDirectory(self):
		# check if folders exist
		base_dir = 'pictures'
		if not os.path.exists(base_dir):
			os.makedirs(base_dir)
		return

if __name__ == '__main__':
	root = Tk()

	TriangleGui(root,600,600,10,1)

	root.mainloop()