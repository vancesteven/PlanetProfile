function y = L_IceNH3(P,T,w,ind)

y = getIcePhaseNH3(P,T,w) - ind + .5;