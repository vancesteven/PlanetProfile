function y = L_IceMgSO4(P,T,w,ind)

y = getIcePhaseMgSO4(P,T,w) - ind + .5;