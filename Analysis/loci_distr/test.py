for i in range(nsw):
    ax.fill_between(t[1:-9],d_mean[i,1:-9]-d_std[i,1:-9]/10,d_mean[i,1:-9]+d_std[i,1:-9]/10,color=color[i],alpha=0.1)
    for k in range(ncl):
        if (i,k)!=(0,1):
            ax.plot([min(t[1:-9]),max(t[1:-9])],[d0_mean[i,k],d0_mean[i,k]],'--',color=color0[k][i],linewidth=wth,label=label0[k][i])
    ax2.plot(t[1:-9],d_mean[i,1:-9],'-',color=color[i],linewidth=wth,label=sw[i])