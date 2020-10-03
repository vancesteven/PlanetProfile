function [out_Hz,out_nT] = localMax(fr,dat,target_fr)
npks = length(target_fr);
[out_Hz, out_nT] = deal(zeros(1,npks));
for it = 1:npks
    hw = 20; %half-width
    ind = find((abs(fr-target_fr(it))) == min(abs(fr-target_fr(it))));
    inds = ind-hw:ind+hw;
    lmax = find(dat(inds) == max(dat(inds)));
    out_Hz(it) = fr(inds(lmax));
    out_nT(it) = dat(inds(lmax));
end
end %localMax