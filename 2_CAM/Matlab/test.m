%TEST   Simple test of EXPMV_TSPAN.
% A = -gallery('poisson',n);
% CoreNum = 2;
% if isempty(gcp('nocreate'))
%     parpool(CoreNum);
% end

A = load('numbers.txt');
% b = linspace(-1,1,n^2)';
b= load('b.txt');
n = length(A);
t = 450;

tic
Y = expm(t*A)*b;
toc
% 
% 
% tic
% [f,s,m,mv,mvd,unA] = expmv(t,A,b);
% toc



% t0 = 0; tmax = 100;
% q = 9;
% 
% [X,tvals,mv] = expmv_tspan(A,b,t0,tmax,q);

% fprintf('Relative differences between vectors from EXPM and EXPMV_TSPAN.\n')
% fprintf('Should be of order %9.2e.\n', eps/2)
% Y = zeros(size(X));
% for i = 1:length(tvals)
%     Y(:,i) = expm(full(tvals(i)*A))*b;
%     fprintf('%2.0f:  %9.2e\n', i, norm(Y(:,i)-X(:,i),1)/norm(X(:,i),1) )
% end
