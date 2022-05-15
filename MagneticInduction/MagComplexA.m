function Ae = MagComplexA(r,sig,w,n, opts)
% r is a three element vector for the base of the ocean, the top of the
% ocean, and the overall radius
% sig is a scalar value for the ocean in S/m
% w is the angular frequency 
prevDigits = digits(300);

mu0=4*pi*1e-7;
% If desired, perform the slower operation of evaluating the k values
% at high precision. When Bessels very nearly cancel, this is an 
% important step to accurately calculate the result.
if exist('opts', 'var')
    if opts.hiprec
        fprintf("Doing high precision k values.\n")
        k = vpa(sqrt(w'.*1i.*mu0.*sig)); % matrix containing the spatially distributed complex wavenumber in frequency, conductivity
    else
        %fprintf("Doing low precision k values.\n")
        k = sqrt(w'.*1i.*mu0.*sig);
    end
else
    %fprintf("Doing (normal) low precision k values.\n")
    k = sqrt(w'.*1i.*mu0.*sig);
end

lr = length(r);
lw = length(w);
Ae = zeros(1,lw);

parfor wi = 1:lw
    %disp("Making Bessels")
    [jn,yn] = sphbess1(k(wi,2:end).*r(1:end-1)); % Non-matching (upper layer) jn,yn
    [jd,yd] = getjydx(n,k(wi,2:end).*r(1:end-1)); % Non-matching (upper layer) jd,yd
    %disp("Halfway done making Bessels")
    [jnm,ynm] = sphbess1(k(wi,:).*r); % Matching (lower layer) jn,yn
    [jdm,ydm] = getjydx(n,k(wi,:).*r); % Matching (lower layer) jd,yd

    T = sym(0);
    for ir = 1:lr-1
        beta = jnm(ir)*yd(ir)-yn(ir)*jdm(ir);
        gamma = ynm(ir)*yd(ir)-yn(ir)*ydm(ir);
        delta = jn(ir)*jdm(ir)-jnm(ir)*jd(ir);
        eps = jn(ir)*ydm(ir)-ynm(ir)*jd(ir);

        T = (delta+T*eps)/(beta+T*gamma);
    end

    % These use matching k,r at outer boundary, hence jnm,ynm, etc.
    beta = jdm(end)-(n+1)*jnm(end);
    gamma = ydm(end)-(n+1)*ynm(end);
    delta = n*jnm(end)+jdm(end);
    eps = n*ynm(end)+ydm(end);

    Ae(wi) = double( (beta+T*gamma)/(delta+T*eps) );
end

digits(prevDigits);
end % MagComplexA

function [jdx,ydx] = getjydxnv(x)
    [j1,y1] = sphbess1nv(x);
    [j2,y2] = sphbess2nv(x);
    jdx = 2*j1-x.*j2;
    ydx = 2*y1-x.*y2;
end % getjydx

function [jdx,ydx] = getjydx(n,x)
    [jn,yn] = sphbessel(n,x);
    [jP,yP] = sphbessel(n+1,x);
    jdx = (n+1)*jn-x.*jP;
    ydx = (n+1)*yn-x.*yP;
end % getjydx

function [jn,yn] = sphbessel(n,x)
    fac = sqrt(pi/2./x);
    jn = fac .* besselj(n+0.5, x);
    yn = -1^(n+1) * fac .* besselj(-n-0.5, x);
end

function [jn,yn] = sphbess1(x)
    [jn,yn] = sphbessel(1,x);
end
