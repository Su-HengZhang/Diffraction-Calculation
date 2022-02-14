function f_roi=cutcenter(f,Mroi,Nroi)

[M,N]=size(f);

if nargin<2
    Mroi=M;
    Nroi=N;
end

if Mroi>=M
    Mroi=M;
end

if Nroi>=N
   Nroi=N;
end

ml=(0:Mroi-1)-floor(Mroi/2);
nl=(0:Nroi-1)-floor(Nroi/2);

Mc=floor(M/2)+1;
Nc=floor(N/2)+1;

ml_c=ml+Mc;
nl_c=nl+Nc;

f_roi=f(ml_c,nl_c);
