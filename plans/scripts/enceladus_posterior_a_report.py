"""Run (a) posterior summary: reproducibility, quantiles, derived isostatic quantities,
correlations, zb-by-log10w bins, and the zb-w joint figure (zb_w_joint.png).
Usage: python plans/scripts/enceladus_posterior_a_report.py [result_dir]"""
import pickle,json,sys,numpy as np
D=(sys.argv[1].rstrip('/')+'/') if len(sys.argv)>1 else 'validation_reports/enceladus_isostasy/posterior_a/'
R=[pickle.load(open(D+f'rep{i}_result.pkl','rb')) for i in (1,2)]
a,b=R
print('names',a.param_names)
print('bit-identical samples',np.array_equal(a.samples,b.samples),'loglike',np.array_equal(a.log_likelihoods,b.log_likelihoods),'weights',np.array_equal(a.weights,b.weights))
w=a.weights/ a.weights.sum(); print('weight unique',np.unique(a.weights).size,'Kish ESS %.0f'%(1/np.sum(w**2)))
print('logZ',a.metadata['log_Z'],a.metadata['log_Z_err'],'skipped',a.metadata['tidal_solve_skipped'])
def q(x,w=w):
    i=np.argsort(x); c=np.cumsum(w[i]); return np.interp([.16,.5,.84],c,x[i])
X=a.samples; out={}
for j,n in enumerate(a.param_names):
    m=np.average(X[:,j],weights=w); s=np.sqrt(np.average((X[:,j]-m)**2,weights=w))
    out[n]=dict(q=q(X[:,j]).tolist(),sd=s); print('%-20s %s sd %.4g'%(n,np.round(q(X[:,j]),4),s))
iso=a.metadata['isostatic_derived']
for k,v in iso.items():
    v=np.asarray(v,float); f=np.isfinite(v).mean()
    if f>0: print('%-20s finite %.4f %s'%(k,f,np.round(q(np.where(np.isfinite(v),v,0)),5)))
print('max|mass_res|',np.nanmax(np.abs(iso['mass_residual_after'])))
C=np.cov(X.T,aweights=w); sd=np.sqrt(np.diag(C)); Cr=C/np.outer(sd,sd)
N=a.param_names; ix={n:i for i,n in enumerate(N)}
print('corr matrix\n',np.round(Cr,3))
zb=X[:,ix['zb_km']]; lw=X[:,ix['log10_wOcean_ppt']]
for lo,hi in [(-1,.5),(.5,1.5),(1.5,2.01)]:
    m=(lw>=lo)&(lw<hi); ww=w*m; fm=ww.sum(); ww=ww/fm
    xs=X[m]; ws=ww[m]; Cb=np.cov(xs.T,aweights=ws); sb=np.sqrt(np.diag(Cb))
    print('bin [%g,%g) mass %.3f zb %s  sd_zb %.3f corr(rho_s,C2) %.3f'%(lo,hi,fm,np.round(q(zb[m],ws),2),sb[0],Cb[ix['rho_rock_kgm3'],ix['compensation_C2']]/sb[ix['rho_rock_kgm3']]/sb[ix['compensation_C2']]))
# autocorr proxy: pocoMC iid-ish; report max loglike, best
k=np.argmax(a.log_likelihoods); print('max loglike',a.log_likelihoods[k],dict(zip(N,np.round(X[k],4))))

# --- fine bins + figure ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
zb=X[:,0]; lw=X[:,1]; rs=X[:,ix['rho_rock_kgm3']]; c2=X[:,ix['compensation_C2']]
def q(x,ww):
    i=np.argsort(x); c=np.cumsum(ww[i]); c/=c[-1]; return np.interp([.16,.5,.84],c,x[i])
edges=np.linspace(-1,2,13); rows=[]
for lo,hi in zip(edges[:-1],edges[1:]):
    m=(lw>=lo)&(lw<hi); ww=w[m]
    C=np.cov(np.vstack([rs[m],c2[m]]),aweights=ww); r=C[0,1]/np.sqrt(C[0,0]*C[1,1])
    rows.append((lo,hi,ww.sum(),*q(zb[m],ww),r))
for r in rows: print('lw [%5.2f,%5.2f) mass %.3f zb %.2f/%.2f/%.2f corr(rs,C2) %.3f'%r)
fig,ax=plt.subplots(1,2,figsize=(10,4.2))
h=ax[0].hist2d(lw,zb,bins=[40,40],weights=w,cmap='Greys'); r=np.array(rows); mid=(r[:,0]+r[:,1])/2
ax[0].plot(mid,r[:,4],'C3-o',ms=3,label='median'); ax[0].fill_between(mid,r[:,3],r[:,5],color='C3',alpha=.2,label='16-84%')
for b in (.5,1.5): ax[0].axvline(b,color='k',ls=':',lw=.8)
ax[0].set_xlabel('log10 w_ocean (ppt)'); ax[0].set_ylabel('z_b (km)'); ax[0].legend(); ax[0].set_title('Run (a): z_b jointly with salinity')
ax[1].plot(mid,r[:,6],'C0-o',ms=3); ax[1].set_xlabel('log10 w_ocean (ppt)'); ax[1].set_ylabel('corr(rho_rock, C2) within bin'); ax[1].set_ylim(-1,0)
ax[1].set_title('rho_rock-C2 trade by salinity bin')
fig.tight_layout(); fig.savefig(D+'zb_w_joint.png',dpi=130); print('saved')
