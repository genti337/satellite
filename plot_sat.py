import argparse
import os
import matplotlib.pyplot as plt
import pandas

# Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument('--log_file', type=str, default='NA', help='Directory containing log files to plot')
parser.add_argument('--plot_peri', action='store_true', help='Plots the spacecraft motion in the Perifocal Frame')

args = parser.parse_args()

# Get a list of the log files
if args.log_file != 'NA':
	log_files = [args.log_file]
else:
	log_files = os.listdir(os.getcwd() + '/Log_data')

# Create Figure for 3D Plot
if args.plot_peri:
	ax = plt.axes()
else:
	ax = plt.axes(projection='3d')

# Read in the data from all of the log files
for log_file in log_files:
	# Create a Pandas Data Frame
	df = pandas.read_csv(os.getcwd() + '/Log_data/%s' % log_file)

	# Plot the Satellite Motion in the Perifocal Frame
	if args.plot_peri:
		ax.plot(df['periPosP'], df['periPosQ'])
		ax.set_xlabel('P (m)')
		ax.set_ylabel('Q (m)')
		ax.set_title('Position w.r.t. Perifocal Frame')
	# Plot the Satellite Motion in the Inertial Frame
	else:
		ax.plot3D(df['eciPosX'], df['eciPosY'], df['eciPosZ'])
		ax.set_xlabel('X (m)')
		ax.set_ylabel('Y (m)')
		ax.set_zlabel('Z (m)')
		ax.set_title('Position w.r.t. Inertial Frame')
                                                       
plt.show()
