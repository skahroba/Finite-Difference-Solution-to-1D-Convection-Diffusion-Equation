import numpy as np
from tridag import *

def ConvDiffsolver(numpar, physpar, bc):
    numpar.alfa3 = numpar.dt/(physpar.peclet*(numpar.dx**2))
    
    varspatial  = np.rec.array([(np.zeros([numpar.nx,1])), (np.zeros([numpar.nx,1]))], 
          dtype=[('c_old','float64'), ('cmod','float64'), ])
    
    pv = np.arange(0.0, numpar.tstop+numpar.dt, numpar.dt)
    vartemporal = np.rec.array([(pv), (np.zeros([len(pv)])), (np.zeros([len(pv)]))], 
          dtype=[('pv','float64'), ('cout','float64'),('time','float64') ])
    
    a = np.zeros([numpar.nx,3])
    b = np.zeros([numpar.nx,1])

    a[:,0]=-(numpar.alfa1+numpar.alfa3);
    a[:,1]=1+numpar.alfa1+2*numpar.alfa3;
    a[0,1]=1+numpar.alfa1+3*numpar.alfa3;
    a[numpar.nx-1,1]=1+numpar.alfa1+numpar.alfa3;
    a[:,2]=-numpar.alfa3;

    for j in range(0,numpar.nt):
        b=varspatial.cmod
        b[0,0]=b[0,0]+(numpar.alfa1+2*numpar.alfa3)*bc.bcL
        varspatial.cmod=tridag(a,b,numpar.nx)
        vartemporal.cout[j+1]=varspatial.cmod[numpar.nx-1,0]
        vartemporal.time[j+1]=vartemporal.time[j]+numpar.dt

    return varspatial, vartemporal

