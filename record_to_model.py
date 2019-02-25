import numpy as np
def record_to_model(NN,rec):
    m=np.zeros(NN)
    v=np.zeros([NN,3])
    x=np.zeros([NN,3])

    i0=0
    m=np.array(rec[i0:i0+NN])
    i0=i0+NN
    for idir in range(3):
        x[:,idir]=rec[i0:i0+NN]
        i0=i0+NN
    for idir in range(3):
        v[:,idir]=rec[i0:i0+NN]
        i0=i0+NN

    return (m,x,v)

def model_to_record(m,x,v):
    NN=len(m)
    recl=NN+NN*3+NN*3
    rec=np.zeros(recl)
    i0=0
    rec[i0:i0+NN]=m
    i0=i0+NN
    for idir in range(3):
        rec[i0:i0+NN]=x[:,idir]
        i0=i0+NN
    for idir in range(3):
        rec[i0:i0+NN]=v[:,idir]
        i0=i0+NN
    
    return (rec)
