function L_IceNH3_GenLookupTable
%j and ww --> y axis
%i and PP --> x axis
%TT function
ww=linspace(0,22,22);
PP=linspace(0,2000,80);
TT=zeros(15,80);
for j=1:22
    for i=1:80
        display(j/22);
        display(i/80);
        TT(j,i)=fzero(@(T) L_IceNH3(PP(i),T,ww(j),1),[200 400]);
    end
end
save('L_IceNH3_DATA','TT','PP','ww');
end