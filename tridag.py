import numpy as np

def tridag(A,b,N):
    L = np.zeros([N,1])
    U = np.zeros([N,2])
    z = np.zeros([N,1])
    y = np.zeros([N,1])
    
    U[0,0]=A[0,1];
    U[0,1]=A[0,2];
    z[0,0]=b[0,0];
    
    for i in range(1,N):
        U[i,1]=A[i,2]
        L[i,0]=A[i,0]/U[i-1,0]
        U[i,0]=A[i,1]-L[i,0]*U[i-1,1]
        z[i,0]=b[i,0]-L[i,0]*z[i-1,0]
    
    
    y[N-1,0]=z[N-1,0]/U[N-1,0]
    
    for i in range(N-2,-1,-1):
        y[i,0]=(z[i,0]-U[i,1]*y[i+1,0])/U[i,0]
    
    return y
