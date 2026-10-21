import sys, numpy as np, json
sys.path.insert(0,'/Users/svance/src/planetprofile-genai/plans/scripts')
sys.path.insert(0,'/Users/svance/src/planetprofile-genai')
import enceladus_zbw_mass_deficit_study as S
import PlanetProfile.Inference.cache_builder as cb
CAP={}
import PlanetProfile.Gravity.Gravity as GG
_orig=GG.SetupGravity
def wrap(P,Pa):
    P,Pa=_orig(P,Pa); CAP['P']=P; CAP['Pa']=Pa; return P,Pa
GG.SetupGravity=wrap
import PlanetProfile.Thermodynamics.LayerPropagators as LP
_SL=LP.SilicateLayers
CALLS=[]
def SLw(P,Pa):
    out=_SL(P,Pa)
    MLs,MAs=out[6],out[7]; r=out[4]
    Mt=(MLs[:,-1]+MAs[:,-1])/P.Bulk.M_kg-1
    iv=out[0]
    CALLS.append(dict(phi=float(P.Sil.fn_phi_frac.multFactor*P.Sil.fn_phi_frac.EOSphiMax_frac), iv=None if __import__('numpy').size(iv)==0 else int(iv), Mt=Mt, rtop=r[:,0]))
    return out
LP.SilicateLayers=SLw
cfg=S._load_config()
zb,w=float(sys.argv[1]),float(sys.argv[2]); dpf=float(sys.argv[3]); nsf=float(sys.argv[4]) if len(sys.argv)>4 else 1.0
s=S.build_node(cfg,zb,w,dp_factor=dpf,nsil_factor=nsf)
P=CAP['P']; M=P.Bulk.M_kg
def outer(r,rho):
    o=np.argsort(r); r=np.asarray(r)[o]; rho=np.asarray(rho)[o]
    ri=np.concatenate(([0],r[:-1])); return np.sum(rho*4/3*np.pi*(r**3-ri**3))
n=len(P.rho_kgm3)
full=outer(P.r_m[:n],P.rho_kgm3)
red=outer(np.array(P.Reduced.r_m),np.array(P.Reduced.rho_kgm3))
ph=np.array(P.phase); 
print('phases unique',np.unique(ph), 'nHydro',P.Steps.nHydro,'nSil',P.Steps.nSil, 'nTot',n)
print('REDUCED sizes',P.Reduced.changeIndices, 'Params RLS',CAP['Pa'].REDUCED_LAYERS_SIZE)
print(json.dumps(dict(PP_Mtot=P.Mtot_kg/M-1, sumMLayer=np.sum(P.MLayer_kg)/M-1, full_outer=full/M-1, reduced_outer=red/M-1, cache=s['mass_residual_frac'],
  CMR2=float(P.CMR2mean), rhoSil=float(P.Sil.rhoMean_kgm3), Rsil_km=float(P.Sil.Rmean_m/1e3), phi=float(P.Sil.phiCalc_frac), Psil0=float(P.P_MPa[P.Steps.nHydro]),
  rho_bounds=[float(x) for x in P.rho_kgm3[P.Steps.nHydro-2:P.Steps.nHydro+3]], r_bounds=[float(x) for x in P.r_m[P.Steps.nHydro-2:P.Steps.nHydro+3]],
  rho_silbot=float(P.rho_kgm3[-1]), Cmeas=P.Bulk.Cmeasured, Cunc=[P.Bulk.CuncertaintyLower,P.Bulk.CuncertaintyUpper], nPoros=P.Steps.nPoros, nSilMax=P.Steps.nSilMax, nRed=len(P.Reduced.r_m)),indent=0))
# segment-by-segment mass comparison
ci=P.Reduced.changeIndices

for c in CALLS:
    iv=c['iv']
    if iv is None: print('none'); continue
    lo=max(iv-2,0)
    print('phi',c['phi'],'iv',iv,'Mt',np.round(c['Mt'][lo:iv+2]*100,4),'dr_km',np.round(-np.diff(c['rtop'][lo:iv+2])/1e3,3))
print('selected phi',P.Sil.phiCalc_frac,'CMR2',P.CMR2mean)
import pickle
pickle.dump(s, open(sys.argv[5],'wb'))
runner=S.make_runner(S.CONFIG_PATH)
def obs(st):
    o=S.observables(runner,st); return {k:o[k][1.0] for k in o}
o0=obs(s)
# pure mass normalization: scale all rho uniformly to target mass
s2=dict(s); s2['rho']=np.asarray(s['rho'])/(1+s['mass_residual_frac'])
o1=obs(s2)
# restore mass by scaling only rock (phase>=50) density
ph=np.asarray(s['phases']); r=np.asarray(s['r_m']); rho=np.asarray(s['rho'],float)
o_=np.argsort(r); ri=np.empty_like(r); rs=r[o_]; ri[o_]=np.concatenate(([0],rs[:-1]))
dm=s['mass_residual_frac']*P.Bulk.M_kg; rock=ph>=50
Mrock=np.sum(rho[rock]*4/3*np.pi*(r[rock]**3-ri[rock]**3))
s3=dict(s); rr=rho.copy(); rr[rock]*= (Mrock-dm)/Mrock; s3['rho']=rr
o3=obs(s3)
sig={'C20':3.699e-5,'C22':1.47e-5,'C30':3.342e-5,'libration_deg':0.003}
print('OBS', {k:o0[k] for k in o0})
print('dsig uniform-rescale', {k:(o1[k]-o0[k])/sig[k] for k in o0})
print('dsig rock-rescale', {k:(o3[k]-o0[k])/sig[k] for k in o0})
