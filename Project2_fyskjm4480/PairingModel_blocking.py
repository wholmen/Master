from pylab import*

class basis_set:
	def __init__(self,Nparticles,Nshells,g,d):
		self.Nparticles = Nparticles;
		self.Nshells = Nshells
		self.g = g
		self.d = d

		self.states = []
		self.nstates = 0

		for p in range(Nshells):

			for ms in [-1,1]:
				self.states.append([p,ms])
				self.nstates += 1

		self.states = array(self.states)

	def ReferenceEnergy(self):
		Energy = 0.0

		if self.Nparticles < self.nstates:
			for p in range(self.Nparticles):
				Energy += self.OneBodyEnergy(p,p)

				for q in range(self.Nparticles):
					Energy += 0.5*self.TwoBodyEnergy(p,q,p,q)

		else:
			print "Calculations can not be done for Nparticles > Nstates."
		
		return Energy

	def OneBodyEnergy(self,a,b):
		p = 0.0
		if self.KD_pstate(a,b):
			p = self.states[a,0]*self.d

		return p

	def TwoBodyEnergy(self,a,b,c,d): 
		asym = 0.0

		if self.KD_pstate(a,b) * self.KD_pstate(c,d) == 1:		
			if self.KD_spin(a,b) == 0:
				if self.KD_spin(c,d) == 0:
					
					if self.KD_spin(a,c) == 1:
						asym = -self.g/2.0
					elif self.KD_spin(a,c) == 0:
						asym = self.g/2.0

		return asym

	def KD_integer(self,a,b):
		# Kroenecker delta for integers a and b
		return 1.0*(a == b)

	def KD_pstate(self,a,b):
		# Kroenecker delta checking if the p-level for a and b are equal
		return 1.0 * self.KD_integer(self.states[a,0], self.states[b,0])

	def KD_spin(self,a,b):
		# Kroenecker delta cheking if the spin is equal for a and b
		return 1.0 * self.KD_integer(self.states[a,1], self.states[b,1])









"""
Test of reference energy is computed correctly. It is.

Np = 4; Ns = 4; g = 0.6; d = 1.0; N=101
g = linspace(-1,1,N)

Ref0 = zeros(N); Ref1 = zeros(N);
i = 0
for G in g:
	b = basis_set(Np, Ns, G, d)

	Ref0[i] = 2 - G
	Ref1[i] = b.ReferenceEnergy()
	i += 1

	#print b.states

	#print b.ReferenceEnergy()

plot(g,Ref0,'r')
plot(g,Ref1,'b')
show()

"""