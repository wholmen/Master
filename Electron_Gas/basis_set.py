from pylab import *


class basis_set:
	def __init__(self, Nparticles, rs, Nshells):
		self.Nparticles = Nparticles
		self.rs = rs
		self.Nshells = Nshells

		self.states = []

		n = 0
		nx2ny2 = 0

		for shell in range(self.Nshells):
			# shell describes the value of nx^2 + ny^2

			for nx in range(-shell,shell+1):
				for ny in range(-shell,shell+1):
					for nz in range(-shell,shell+1):
						# Finding all possible combinations of nx and ny at this shell
						
						e = nx**2 + ny**2 + nz**2
						if e == shell:
							# Remove the combinations which either exceeds energy level
							# or duplicate lower levels 

							for ms in [-1,1]:
								# Including both spins

								s = [e,nx,ny,nz,ms]
								self.states.append(s)

		# States are initiated. Calculating energies and constants
		self.nstates = size(self.states,0)
		
		self.L3 = (4*pi*self.rs**3*self.Nparticles) / 3.0
		self.L2 = self.L3**(2/3.0)
		self.L = self.L3**(1/3.0)

		self.kstep = 2*pi / self.L

		for i in range(self.nstates):
			self.states[i][0] *= 2*pi**2/self.L**2




	def one_body(self,p,q):

		return 0

	def two_body(self,p,q,r,s):


		return 0


b = basis_set(10,1,2)
print size(b.states,0)