#!/usr/bin/env python3

import struct, sys, numpy as np
import record_to_model as recmod

filename='output.bin'
#read parameters
f=open(filename,'br')
rec0=struct.unpack('8d',f.read(8*8))
NN=int(rec0[0]+.5)
Norm_x=rec0[1]
Norm_t=rec0[2]
Norm_m=rec0[3]
eps=rec0[4]
dt=rec0[5]
nt=int(rec0[6]+.5)
frameskip=int(rec0[7]+.5)
recl=NN+NN*3+NN*3

# Point to first record with a frame
irec=1
f.seek(recl*8*irec)
rec0=struct.unpack(str(recl)+'d',f.read(8*recl))
(m,x,v)=recmod.record_to_model(NN,rec0)


f.close()

