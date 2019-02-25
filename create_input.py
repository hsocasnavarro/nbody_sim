#!/usr/bin/env python3
import struct, sys, numpy as np
import record_to_model as recmod

# Normalization values (position, time and mass are normalized
#   to these values for better numerical precision)
Norm_x=1.
Norm_t=1.
Norm_m=1.
# Box size
sizex=Norm_x
sizey=sizex
sizez=sizey
# smoothing parameter
eps=0.01*Norm_x
# Time step
dt=Norm_t
# Number of steps
nt=2000
# Whether to use periodic (1) or open (0) boundary conditions
periodic=1 # Periodic

# Set to an integer greater than one if you do not wish to record every
#   frame to disk. Example, frameskip=10 means record 1 out of every
#   10 frames
frameskip=1 

# Number of particles
NN=50000 
m=np.ones(NN) # set all masses to 1
x=np.random.uniform(-sizex,sizex,(NN,3)) # set positions to uniform distribution
v=np.zeros((NN,3)) # set all velocities to zero


record_length=NN+NN*3+NN*3 # Length of each frame (we record m, x and v at each step)

header_parameters=[NN, Norm_x, Norm_t, Norm_m, eps, dt, nt, frameskip, periodic, sizey/Norm_x, sizez/Norm_x]

if len(header_parameters) > record_length:
    print('Error. Size of record not enough to hold parameters ',NN)
    sys.exit(1)

# First record has parameter information
first_record=header_parameters+[0.]*(recl-len(headerparams))
f=open('input.bin','bw')
f.write(struct.pack(str(len(first_record))+'d',*first_record))

# Second record has the initial values for m,x,v
record=recmod.model_to_record(m,x,v) # Convert variables to record array
f.write(struct.pack(str(len(record))+'d',*record))

f.close()

