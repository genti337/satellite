# Introduction

This tool can be used to help visualize a spacecraft's orbit around the Earth.  The desired orbit can be set from the orbital 
elements used to uniquely define an orbit.  The spacecraft's motion is then propagated by integrating the gravitational acceleration. 
It is also possible to propagate multiple orbit by creating additional instances of the satellite class.  This tool also includes a 
script to plot the spacecraft's inertial position in a 3D plot or the motion in the perifocal frame of reference.

My approach to the problem was to make it easy for the user to select the orbit to analyze through the orbital elements and let the 
script calculate the inertial position and velocity needed to initialize the spacecraft in that orbit.  I also made the satellite a 
a class object in case there was more than orbit that the user wanted to see simultaneously.  The data for each object is is logged 
in a comma separated file which is then processed by the plotting tool.  I chose to plot the spacecraft motion in 3D space for the 
inertial position because depending on the orbit inclination it may be difficult to visualization in the X-Y or Y-Z planes.  The option 
to view the motion in the perifocal frame was also added because it is easy to visualize the orbit in the orbital plane.

A couple things to note is that this model assumes the only force acting on the spacecraft is the gravitational force and the J2 effect due to 
the Earth's oblateness was not taken into account so this is an ideal orbit.  The integration method is also very simple and should be updated in the 
future to support more accurate results.  This is also a 3 degree of freedom model since the rotational motion is not simulated at this point.

# Installation
This repository uses PDM to manage the required python packages to run the scripts.  The pyproject.toml file has been added to the 
repository to track the needed packages and versions.  To install the required packages the following command can be run.

```bash
pdm install
```

# How to Run
There are two scripts provided with this toolset, the first is a class to propagate the satellite motion and the second is plotting tool 
to process the logged data.  The first step is using this tool is configure the satellite class instances in satellite.py with the following 
steps.

1. Create a new instance of the satellite class and set the desired orbit.  The arguments for the satellite constructor are as follows.

- Perigee Altitude (m)
- Orbit Eccentricity
- Longitude of Ascending Node (r)
- Inclination (r)
- Augment of Periapsis (r)
- File Name to Log Data

```python
sat1 = satellite(413000.0, 0.2, 0.0, math.radians(10.0), math.radians(30.0), 'log_sat4.csv')
```
2. Call the initialization routine with the desired true anomaly to start the satellite's motion.
```python
sat1.initialize(math.radians(0.0))
```
3. Call the propagation function with the desired number of orbits and the integration time step.
```python
sat1.propagate(3.0, 0.5)
```
4. Run the script in a terminal.  A log file will created for each instance of the of the satellite class in the Log_data directory.
```bash
pdm run python3.9 satellite.py
```

After the script to propagate the orbits for the satellites have been run, the data can be view using the plot_sat.py script.  The script includes 
a couple of options for plotting the data described below.

1. The script will create a 3D plot of the spacecraft's inertial position in 3D for all satellite instances.
```bash
pdm run python3.9 plot_sat.py
```
2. There is an option to plot the desired log file if desired using the --log_file argument.
```bash
pdm run python3.9 plot_sate.py --log_file=log_sat3.csv
```
3. Finally, the data can also be plotted in the Perifocal Frame using the --plot_peri argument.
```bash
pdm run python3.9 plot_sate.py --plot_peri
```
