function [q_c,Q_c]= findControlSet_nonlinear(x0,x0shMat,xB,Ac,Bc,q0,Q0)
%clear;
%load result
T=4;
% basisMat = [1 zeros(1,9); 0 1 zeros(1,8); 0 0 1 zeros(1,7)]';  % orthogonal basis of (x1, x2) subspace
% psObj_A = rsObj_A.projection(basisMat);  % reach set projection
% psObj_A = psObj_A.cut(T);
% psObj_B = rsObj_B.projection(basisMat);  % reach set projection
% [x0,x0shMat]=x0EllObj_A.double();
%[xB,xB0shMat]=x0EllObj_B.double();
%[ctrMat, ttVec] = rsObj_B.get_center();
%xB=ctrMat(:,end);

% assume lstar. n is number of control variables
l=[0;1;0];
n = 1;
l=l/norm(l);

%%
%Compute int G(t,s) from s=0 to t. G(t,s)=e(t-s)A
% as well as int of sqrt(l'GBQB'G'l)
N=6; % estimated number of product need to take to find Gint 
Gtemp=eye(size(Ac))*T;
for i=1:N
    mAc=mpower(-Ac,i);
    if isequal(mAc,zeros(size(Ac)))
        break
    end
    Gtemp=Gtemp+T^(i+1)/factorial(i+1)*mAc;
end
Gint=expm(Ac*T)*Gtemp;
f= @(s) sqrt(l'*expm(Ac*(T-s))*Bc*Bc'*(expm(Ac*(T-s)))'*l);
GGint=integral(f,0,T,'ArrayValued',1);
        
%%
% define SDP

%f = -2*( xB'*l - l'*expm(Ac*T)*x0 - GGint(1) - sqrt(l'*expm(Ac*T)*x0shMat*...
%    (expm(Ac*T))'*l ))*l'*Gint*Bc;
fc=( xB'*l - l'*expm(Ac*T)*x0 - sqrt(l'*expm(Ac*T)*x0shMat*...
    (expm(Ac*T))'*l ) );
k_v=0.5:0.1:1.2;
for iter=1:length(k_v)
    k=k_v(iter);
cvx_begin SDP
cvx_quiet false
cvx_precision high
    variable q(n)
    variable Q(n,n) semidefinite;
    variable lambda
    lambda > 0
    %fc- GGint(1)*norm(Q,2)-l'*Gint*Bc*q >=0
    [1-lambda zeros(1,n) (q-q0)';
        zeros(n,1) lambda*eye(n) Q;
        q-q0 Q Q0] >=0
    maximize(fc- GGint(1)*norm(Q,2)-l'*Gint*Bc*q + k*log_det(Q))
cvx_end
q_c{iter}=q;
Q_c{iter}=Q'*Q;

%figure(iter)
%plot(ellipsoid(q,Q'*Q))
%hold on
%splot(ellipsoid(q0,Q0),'y')
end
%%
return