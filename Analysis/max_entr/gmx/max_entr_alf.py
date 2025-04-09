def max_entr_alf(sysname,emapfile,matfile,dmapfile,ialffile,oalffile,pmin,kB,T):
    """Update the alphas in data-driven term according to maximum entropy principle."""
    """
        Inputs:
        emapfile is the inputting experimental Hi-C file.
        matfile is the inputting Binv file.
        dmapfile is the inputting map of simulated contact probability minus experimental one.
        ialffile is the inputting alpha.pyw file to be updated.
        oalffile is the outputting updated alpha.pyw file.
        kB is the Boltzmann constent.
        T is temperature.
    """

    """Import Modules"""
    import datetime
    import numpy as np
    import h5py
    from pyw import pyw
    from scipy.io import loadmat
    import matplotlib.pyplot as plt
    import seaborn as sns

    Pa=pyw(emapfile,'Map')[0]
    ### Even though there is only one column of data, Pa is still a hierarchical list. So "[0]" is
    ### essential.
    N=int((len(Pa))**0.5)
    P=np.zeros((N,N),dtype=float)
    for i in range(N):
        for j in range(N):
            P[i,j]=Pa[i*N+j]
    hipair=[]
    for i in range(N):
        for j in range(i):
            if P[i,j]>=pmin and abs(i-j)>2:
                hipair.append([i,j])

    # mat=loadmat(f'{matfile}')
    # B_=mat['B_']
    with h5py.File(matfile,'r') as f:
        B_=f['B_']
        B_=B_[:]
    delta_f_ave=np.array(pyw(dmapfile,'minus'))
    delta_f_ave=delta_f_ave.reshape(-1,1)
    ### For array format, neither "#.T" nor "np.transpose(#)" can turn it into a column vector. For
    ### this format, they all can do.
    delta_alf=kB*T*np.matmul(B_,delta_f_ave)
    delta_alf=delta_alf.reshape(1,-1)
    alf0_all=pyw(ialffile,'Alpha')
    alf0=[]
    for i in range(len(alf0_all[0])):
        j,k=divmod(i,N)
        if [j,k] in hipair:
            alf0+=[alf0_all[0][i]]
    alf0=np.array(alf0)
    alf=alf0+delta_alf
    Alf=np.zeros((N,N),dtype=float)
    # Alf+=2
    for i in range(N):
        for j in range(N):
            if [i,j] in hipair:
                k=hipair.index([i,j])
                # Alf[i,j]=alf[0][k]
                # Alf[j,i]=alf[0][k]
                if alf[0][k]<=10000 and alf[0][k]>=-10000:
                    Alf[i,j]=alf[0][k]
                    Alf[j,i]=alf[0][k]
                if alf[0][k]<-10000:
                    Alf[i,j]=-10000
                    Alf[j,i]=-10000
                if alf[0][k]>10000:
                    Alf[i,j]=10000
                    Alf[j,i]=10000

    """Output Data"""
    pyw=open(oalffile,'w')
    pyw.write('Alpha of data-driven term:\n')
    for i in range(N):
        for j in range(N):
            pyw.write(f'{Alf[i,j]}\n')

    """Plot"""
    opathlist=oalffile.split('/')[:len(oalffile.split('/'))-1]
    opath=''
    for i in range(len(opathlist)):
        opath+=opathlist[i]+'/'

    plt.figure()
    ax=sns.heatmap(Alf,center=0)
    ax.invert_yaxis()
    plt.xlabel('Index of particles')
    plt.ylabel('Index of particles')
    plt.savefig(opath+f'alf_{sysname}.png',format='png')
    plt.savefig(opath+f'alf_{sysname}.pdf',format='pdf')