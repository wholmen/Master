from pylab import *

class electronbasis():
	def __init__(self, N, rs, Nparticles):
		##################################################
		##
		##  Initialize basis:
		##  N = number of shells
		##  rs = parameter for volume
		##  Nparticles = Number of holes 
		##
		##################################################

		self.rs = rs
		self.states = []
		self.nstates = 0
		self.nparticles = Nparticles
		self.nshells = N - 1
		self.Nm = N + 1

		self.k_step = 2*(self.Nm + 1)
		Nm = N
		n = 0 # current shells
		ene_integer = 0
		while n <= self.nshells:
			is_shell = False
			for x in range(-Nm, Nm+1):
				for y in range(-Nm, Nm+1):
					for z in range(-Nm, Nm+1):
						e = x*x + y*y + z*z
						if e == ene_integer:
							is_shell = True
							self.nstates += 2
							self.states.append([e, x, y, z, 1])
							self.states.append([e, x, y, z,-1])

			if is_shell:
				n += 1
			ene_integer += 1

		self.L3 = (4*pi*self.nparticles*self.rs**3) / 3.0
		self.L2 = self.L3**(2/3.0)
		self.L = pow(self.L3, 1/3.0)

		for i in range(self.nstates):
			self.states[i][0] *= 2*(pi**2) / self.L**2

		self.states = array(self.states)

	def hfenergy(self, nParticles):
		# Calculate HF-energy (reference energy) for nParticles particles
		e0 = 0.0
		if nParticles <= self.nstates:
			for i in range(nParticles):
				e0 += self.h(i,i)
				for j in range(nParticles):
					if j != i:
						e0 += 0.5*self.v(i,j,i,j)
		else:
			# Safety for cases where nParticles exceeds size of basis
			print "Not enough basis states"

		return e0

		def h(self,p,q):
			# Return single particle energy
			return self.states[p,0]*(p==q)

		def v(self,p,q,r,s):
			# Two body interaction for electron gas
			val = 0
			terms = 0.0
			term1 = 0.0
			term2 = 0.0
			kdpl = self.kdplus(p,q,r,s)
			if kdpl != 0:
				val = 1.0 / self.L3
				