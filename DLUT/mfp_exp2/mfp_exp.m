function [T,rn]=mfp_exp(H,err,non_norm);
%%%%%%%%%  This function evaluate the matrix exponential of H by using the
%%%%%%%%%  algorithm proposed in "High-performance computation of large
%%%%%%%%%  sparse matrix exponential".
%%%%%%%%%  err is the error tolerance
%%%%%%%%%  H is the matrix need to be computed;
%%%%%%%%%  T is the matrix exponential computed using this funciton.
%%%%%%%%%  non_norm==0, H is normal matrix; non_norm==1, H is non-normal
%%%%%%%%%  matrix.
%Wu Feng writed on 2019.8.19 (vonwu@dlut.edu.cn)
%Wu Feng modified on 2021.11.29 (vonwu@dlut.edu.cn)

Nh=size(H,1);
H1=norm_c(H,'fro');   
Ih=speye(Nh);

[M,N]=getMN(H1,err);         
% N=20;
% M=5;
h1=H1/(2^N);
r0=getfi(M,h1);

% [N,M]

if non_norm==0;
    ai=1/(N+1);
elseif non_norm==1;
    ai=1/norm_c(H,'fro');
end
b0=ai*r0/M/exp(2*h1);

%%%%%%Compute T_0;
m=2^N;
H1=H/m;
S=H1;  T=H1;  mt=5;
for i=2:M;
%     i
    S=S*(H1/i);
    %%%fiilter S under the condition that norm(S-S1)<=b0;
%     [S,mt]=flitoutA2(S,b0,mt);
    [S,mt]=filt_norm(S,b0,'fro',mt);
    sc=nnz(S)/Nh^2;
    [i, sc]
    T=sparse(T+S);

end

clear S H1;

rn=[];

mt=1;
for i=1:N
%     i
    r0=2*r0+r0^2
    bi=ai*r0;
    
    T=sparse(2*T+T*T);   
%     T=2*T+T*T;
    
    %%%%%compute the F-norm of Ih+T;
    n1=norm_c(T,'fro');
    dt=diag(T); n2=dt'*dt; 
    dt=dt+1;  n3=dt'*dt;
    nt=(n1^2-n2+n3)^0.5;    
    %%%%%nt is the F-norm of Ih+T;
    [T,mt]=filt_norm(T,bi*nt,'fro',mt);
    sc=nnz(T)/Nh^2 ;   %%%%%%%the sparsity of T;
    rn=[rn,bi*nt];
%     [i, sc]
end
T=Ih+T;

function fi=getfi(M,h)
%%%%%%%%evaluate the relative error of the M-order Taylor series of exp(h),
%%%%%%%%return fi;

%%%%%%s0=h^(M+1)/(M+1)!;
s0=h;
for i=1:M;
    s0=s0*h/(i+1);
end
%%%%%%compute fi;
i=0;  si=s0;  f=s0;  
while si>=1e-30;
    si=si*h*(i+M+1)/(i+1)/(i+M+2);
    if si == inf;
        break;
    end
    f=f+si;
    i=i+1;   
end
fi=f;


function [M,N]=getMN(h,err)

%%%%%%%%%%%In terms of h and the error tolerance err, select M and N.

Nmin=max(floor(log2(h)),0);
Nn=50;
for i=1:50;
    N(i)=Nmin+i;
    hi=h/(2^N(i));
    Mj=max([floor(hi)-1,0]);  rj=1;
    while rj>err;
        Mj=Mj+1;
        r0=getfi(Mj,hi); rj=r0;
        for j=1:N(i);
            rj=2*rj+rj^2;
        end
    end
    rj=1;
    M(i)=Mj;
end
f=M.*(2.^N);
[I,J]=min(f);
M=M(J);
N=N(J);