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

Nh=size(H,1);
H1=norm(H,'fro');   
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
    ai=1/norm(H,'fro');
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
    [S,mt]=flitoutA2(S,b0,mt);
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
    n1=norm(T,'fro');
    dt=diag(T); n2=dt'*dt;
    dt=dt+1;  n3=dt'*dt;
    nt=(n1^2-n2+n3)^0.5;    
    %%%%%nt is the F-norm of Ih+T;
    
    [T,mt]=flitoutA2(T,bi*nt,mt);
    sc=nnz(T)/Nh^2 ;   %%%%%%%the sparsity of T;
    rn=[rn,bi*nt];
%     [i, sc]
end
T=Ih+T;



function [A,m]=flitoutA2(A,eg,m)
%This function filters out the near-zero elements in A to return a sparse
%matrix A_s by using the algorithm 4.1 proposed in "High-performance computing of large
%%%%%%%%%  sparse matrix exponential". 
% If the matrix eliminated after filtering is denoted by B, than b=norm(B,'fro')<=eg
%Wu Feng writed on 2019.8.19

if m==0; m=1; end
eg=eg/2;
b=1;
n=0;
Na=size(A,1);

ef=eg/m;

[Ib,Jb,B]=find(A);     
p=find(abs(B)>ef);

Ia=Ib(p);  Ja=Jb(p);  A=B(p);
B(p)=0;
b=norm(B,2);

n=0;
while b>2*eg;
    m=b/ef;
    ef=eg/m;
    p=find(abs(B)>ef);
    Ia=[Ia;Ib(p)];  
    Ja=[Ja;Jb(p)];  
    A=[A;B(p)];
    B(p)=0;
    b=norm(B,2);
    n=n+1;
end
clear B Ib Jb;
A=sparse(Ia,Ja,A,Na,Na);

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