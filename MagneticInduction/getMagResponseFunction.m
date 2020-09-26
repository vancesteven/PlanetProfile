function [k,Q] = getMagResponseFunction(r,n,nE,w,sig,boundaries,R_m,r0,y0,opts)
    nbound = length(boundaries(:,1));
    nsig = length(sig(:,1));
    lw = length(w);
    Q = zeros(lw,nsig,nbound);
    k = zeros(lw,length(r));
    for lm = 1:nbound
        for jk = 1:nsig
            for ij = 1:lw
                Q(ij,jk,lm) = getQ(n,nE*w(ij),sig(jk,:),boundaries(lm,:),R_m,r0,y0,opts);
            if jk ==1
               k(ij,:) = fnkRad(w(ij),sig(jk,:),r,boundaries(lm,:));
            end
            end
        end
    end
end %getMagResponseFunction
function Q = getQ(n,nEw,sig,boundaries,R_m,r0,y0,opts)
% [~,Qv] = ode15s(@(vr,Q) getdQdr(vr,Q,n,nEw,sig,boundaries),[r0 R_m],y0,opts); % chosen becuase this ode was deemed to be stiff
[~,Qv] = ode45(@(vr,Q) getdQdr(vr,Q,n,nEw,sig,boundaries),[r0 R_m],y0,opts); % this yielded values within 0.1% of those from ode15s. It's the first recommended solution and may work in general. Needs more testing to see if it breaks for some inputs.
% [~,Qv] = ode113(@(vr,Q) getdQdr(vr,Q,n,nEw,sig,boundaries),[r0 R_m],y0,opts); % this yielded values within 0.1% of those from ode15s. It's the first recommended solution and may work in general. Needs more testing to see if it breaks for some inputs.
Q = Qv(end);    
end % get Q
function dQdr = getdQdr(r,Q,n,w,sig,boundaries)
    k = fnkRad(w,sig,r,boundaries);
    dQdr = -n*r.*k.^2/(1+n)/(1+2*n) - Q.*(1+4*n+4*n^2-2*r.^2.*k.^2)/(1+2*n)./r - Q.^2.*(1+n).*r.*k.^2/n/(1+2*n);
%     dQdr = (n+1)/n/(2*n+1)*k.^2.*r.*(n/(n+1)-Q).^2-(2*n+1)./r.*Q; %eckhardt's form is slightly slower
%    dQdr = -r.*k.^2/6 - Q.*(9-2*r.^2.*k.^2)/3./r - Q.^2.*2.*r.*k.^2/3; n =
%    1
%    dipole response
end % fndQrk
function k = sig2k(w,sig)
    mu0 = 4*pi*1e-7;
    k = sqrt(w*mu0*sig*1i);
end % sig2k
function kRad = fnkRad(w,sig,r,boundaries)
    %radial piecewise vector of k
    lr = length(r);
    kRad = ones(1,lr);
    for ij = 1:length(r)
        if r(ij)<=boundaries(1)
            kRad(ij) = sig2k(w,sig(1));
        else
            for jk = 1:length(boundaries)-1
                if r(ij)> boundaries(jk) && r(ij) <= boundaries(jk+1) 
                  kRad(ij) = sig2k(w,sig(jk+1));
                end
            end
        end
    end
end % fnkRad