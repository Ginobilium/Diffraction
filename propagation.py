def propogation(u0, d, N = 512, dL = 5.12e-3, lmb = 632.8e-9)
    #Parameter
    df = 1.0/dL
    k = np.pi*2.0/632.8e-9
    D= dL*dL/(N*lmb)
  
    #phase
    def phase(i,j):
        i -= N//2
        j -= N//2
        return (i*df)*(i*df)+(j*df)*(j*df)
    ph  = np.fromfunction(phase,shape=(N,N))

    #H
    H = np.exp(1.0j*k*d)*np.exp(-1.0j*lmb*np.pi*d*ph)

    #Result
    return (np.fft.ifft2(fftshift(H)*np.fft.fft2(u0)*dL*dL/(N*N))*N*N*df*df, d/D)
