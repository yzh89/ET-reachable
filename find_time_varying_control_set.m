function [q,alpha]= find_time_varying_control_set(x0,x0shMat,qB,MB,l,Ac,Bc,u0,U0)
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

% assume lstar
l=[0;1;0;zeros(7,1)];
l=l/norm(l);

%%
% varying factor to obtain good bound on U
%Q0=shMat;
%q0=centVec';

factor=2:0.4:4;
for iter=1:length(factor)
    Q=Q0/factor(iter);
    
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
    f= @(s) sqrt(l'*expm(Ac*(T-s))*Bc*Q*Bc'*(expm(Ac*(T-s)))'*l);
    GGint=integral(f,0,T,'ArrayValued',1);
    
    
    % define QCQP
    n = 3;
    %f = -2*( xB'*l - l'*expm(Ac*T)*x0 - GGint(1) - sqrt(l'*expm(Ac*T)*x0shMat*...
    %    (expm(Ac*T))'*l ))*l'*Gint*Bc;
    fc=( xB'*l - l'*expm(Ac*T)*x0 - GGint(1) - sqrt(l'*expm(Ac*T)*x0shMat*...
        (expm(Ac*T))'*l ) );
    %H = Bc'*Gint'*(l*l')*Gint*Bc;
    
    cvx_begin
    variables q(n)
    (q-q0)'*inv(Q)*(q-q0) <= 1
    maximize(fc-l'*Gint*Bc*q)
    cvx_end
    
    cvx_begin SDP
    variable lambda
    lambda > 0
    [inv(Q0)-lambda*inv(Q), -inv(Q0)*q0+lambda*inv(Q)*q;
        (-inv(Q0)*q0+lambda*inv(Q)*q)', q0'*inv(Q0)*q0-1-lambda*(q'*inv(Q)*q-1)]<=0
    minimize(0)
    cvx_end
    if isnan(lambda)
        %continous
    else
        
        plot(ellipsoid(q,Q))
        hold on
        plot(ellipsoid(q0,Q0))
        break;
    end
end
%%
return