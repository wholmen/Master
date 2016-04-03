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

		self.states = array(self.states)

	def hf_energy(self):
		# Calculate the reference energy
		

	def one_body(self,p,q):
		# As the list states[p,0] holds 2*pi^2 / L * (nx^2 + ny^2 + nz^2), we do not need new calculations
		# p==q for delta_pq
		return self.states[p,0]*(p==q)

	def two_body(self,p,q,r,s):

		asym = 0

		if self.KDelta_sum(p,q,r,s) == 1:
			asym = 4*pi/self.L3

			if self.KDelta_spin(p,r)*self.KDelta_spin(q,s) == 1:
				asym1 = (1 - self.KDelta_k(p,r)) / self.Absolute_Distance(p,r)**2

			elif self.KDelta_spin(p,s)*self.KDelta_spin(q,r) == 1:
				asym2 = (1 - self.KDelta_k(p,s)) / self.Absolute_Distance(p,s)**2

		return asym*(asym1 - asym2)

	def KDelta_nr(self,a,b):
		# Kroenecker delta for integer comparison. Returns 0 if a and b are unequal
		return 1 * (a==b)
		
	def KDelta_array(self,a,b):
		# Kroenecker delta for array. Returns 0 if any term in a and b are unequal
		d = 1.0
		for i in range(len(a)):
			d *= (a[i]==b[i])
		return d

	def KDelta_k(self,p,q):
		# Kroenecker delta for two wave numbers kp and kq
		return self.KDelta_array( self.states[p,1:4], self.states[q,1:4])

	def KDelta_sum(self,p,q,r,s):
		# Kroenecker delta for sum of p+q and r+s. 
		return self.KDelta_array(self.states[p,1:4]+self.states[q,1:4], self.states[r,1:4]+self.states[s,1:4])

	def KDelta_spin(self,p,q):
		# Kroenecker delta for spin
		return 1 * (self.states[p,4]==self.states[q,4])

	def Absolute_Distance(self,p,q):
		return norm( self.states[p,1:4]-self.states[q,1:4] )





b = basis_set(10,1,2)
print size(b.states,0)