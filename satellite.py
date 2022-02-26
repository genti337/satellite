import numpy
import math
import matplotlib.pyplot as plt

# Satellite Class
class satellite:
	def __init__(self, perigee_alt, e, long_asc, inclination, arg_peri, log_file):
		# Earth's Gravitational Constant
		self.u_earth = 3.98574405096e14

		# Earth Mass
		self.mass_earth = 5.972e24

		# Earth Radius
		self.earth_radius = 6.357e6

		# Universal Gravitational Constant
		self.gm = 6.673e-11

		# Eccentricity of the Oribt
		self.e = e

		# Perigree of the Orbit
		self.rp = self.earth_radius + perigee_alt

		# Calculate the Orbital Angular Momentum
		self.h = math.sqrt(self.rp * self.u_earth * (1.0 + self.e))

		# Apogee of the Orbit
		self.ra = math.pow(self.h, 2) / (self.u_earth * (1.0 - self.e))

		# Orbit Longitude of Ascending Node
		self.long_asc = long_asc

		# Inclination of the Orbit
		self.inclination = inclination

		# Argument of Periaposis
		self.arg_peri = arg_peri

		# Rotation from Perifocal to Inertial Reference Frame
		x1 = math.cos(self.long_asc) * math.cos(self.arg_peri) - math.sin(self.long_asc) * math.cos(self.inclination) * math.sin(self.arg_peri)
		x2 = math.sin(self.long_asc) * math.cos(self.arg_peri) + math.cos(self.long_asc) * math.cos(self.inclination) * math.sin(self.arg_peri)
		x3 = math.sin(self.inclination) * math.sin(self.arg_peri)
		y1 = -math.cos(self.long_asc) * math.sin(self.arg_peri) - math.sin(self.long_asc) * math.cos(self.inclination) * math.cos(self.arg_peri)
		y2 = -math.sin(self.long_asc) * math.sin(self.arg_peri) + math.cos(self.long_asc) * math.cos(self.inclination) * math.cos(self.arg_peri)
		y3 = math.sin(self.inclination) * math.cos(self.arg_peri)
		z1 = math.sin(self.inclination) * math.sin(self.arg_peri)
		z2 = -math.sin(self.inclination) * math.cos(self.arg_peri)
		z3 = math.cos(self.inclination)
		self.T_eci2perifocal = numpy.array([[x1, x2, x3], [y1, y2, y3], [z1, z2, z3]])
		self.T_perifocal2eci = numpy.transpose(self.T_eci2perifocal)

		# Position/Velocity Vectors
		self.perifocalPos = [0.0, 0.0, 0.0]
		self.perifocalVel = [0.0, 0.0, 0.0]
		self.eciPos = [0.0, 0.0, 0.0]
		self.eciVel = [0.0, 0.0, 0.0]
		self.eciAcc = [0.0, 0.0, 0.0]

		# Log File
		self.log_file = log_file
		self.writer = open('Log_data/' + self.log_file, 'w')

	# Initialize the Satellite Orbital State
	def initialize(self, true_anom):
		# Current Radius Vector
		radius = (math.pow(self.h, 2) / self.u_earth) / (1.0 + self.e * math.cos(true_anom))

		# Position/Velocity wrt Perifocal Frame
		self.perifocalPos = [radius * math.cos(true_anom), radius * math.sin(true_anom), 0.0]
		self.perifocalVel = [-(self.u_earth / self.h) * math.sin(true_anom), (self.u_earth / self.h) * (self.e + math.cos(true_anom)), 0.0]	

		# Inertial Position/Velocity
		self.eciPos = numpy.dot(self.T_perifocal2eci, self.perifocalPos)
		self.eciVel = numpy.dot(self.T_perifocal2eci, self.perifocalVel)

		# Write the Log Header Information
		self.writer.write('time,eciPosX,eciPosY,eciPosZ,periPosP,periPosQ,periPosW\n')

	# Propagate the Satellite Motion
	def propagate(self, num_orbits, dt):
		# Orbital Period
		T = 2.0 * math.pi * math.sqrt(math.pow((self.rp + self.ra) / 2.0, 3) / self.u_earth)

		# Time to Propagate the Orbit
		prop_time = T * num_orbits

		time = 0.0
		while time < prop_time:
			# Gravitational Acceleration Magnitude
			grav_accel = (self.gm * self.mass_earth) / math.pow(numpy.linalg.norm(self.eciPos), 2)
	
			# Normalize the Position Vector
			eciPosNorm = [self.eciPos[0] / numpy.linalg.norm(self.eciPos),
			              self.eciPos[1] / numpy.linalg.norm(self.eciPos),
			              self.eciPos[2] / numpy.linalg.norm(self.eciPos)]
	
			i = 0
			while i < 3:
				# Inertial Acceleration
				self.eciAcc[i] = eciPosNorm[i] * -grav_accel
	
				# Integrate the Position/Velocity
				self.eciVel[i] = self.eciVel[i] + self.eciAcc[i] * dt
				self.eciPos[i] = self.eciPos[i] + self.eciVel[i] * dt
	
				i = i + 1

			# Inertial Position/Velocity
			self.perifocalPos = numpy.dot(self.T_eci2perifocal, self.eciPos)
			self.perifocalVel = numpy.dot(self.T_eci2perifocal, self.eciVel)

			# Write the Data to the Log File
			self.writer.write("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n" % (time, self.eciPos[0], self.eciPos[1], self.eciPos[2], self.perifocalPos[0], self.perifocalPos[1], self.perifocalPos[2]))

			# Increment the Time
			time = time + dt

# Instantiate instances of the Satellite Class

# Satellite 1
sat1 = satellite(413000.0, 0.0, 0.0, 0.0, 0.0, 'log_sat1.csv')
sat1.initialize(math.radians(0.0))
sat1.propagate(3.0, 0.5)

# Satellite 2
sat2 = satellite(413000.0, 0.1, 0.0, math.radians(150.0), 0.0, 'log_sat2.csv')
sat2.initialize(math.radians(0.0))
sat2.propagate(3.0, 0.5)

# Satellite 3
sat3 = satellite(435000.0, 0.2, 0.0, math.radians(30.0), 0.0, 'log_sat3.csv')
sat3.initialize(math.radians(0.0))
sat3.propagate(3.0, 0.5)

# Satellite 4
sat4 = satellite(413000.0, 0.2, 0.0, math.radians(10.0), math.radians(30.0), 'log_sat4.csv')
sat4.initialize(math.radians(0.0))
sat4.propagate(3.0, 0.5)

# Satellite 5
sat5 = satellite(450000.0, 0.2, 0.0, math.radians(10.0), math.radians(30.0), 'log_sat5.csv')
sat5.initialize(math.radians(0.0))
sat5.propagate(3.0, 0.5)
