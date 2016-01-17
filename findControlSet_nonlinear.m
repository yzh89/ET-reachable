function [q_c,Q_c]= findControlSet_nonlinear(x0,x0shMat,xB,Ac,Bc,q0,Q0, k_v, method)
%clear;
%load result
T=10;
% basisMat = [1 zeros(1,9); 0 1 zeros(1,8); 0 0 1 zeros(1,7)]';  % orthogonal basis of (x1, x2) subspace
% psObj_A = rsObj_A.projection(basisMat);  % reach set projection
% psObj_A = psObj_A.cut(T);
% psObj_B = rsObj_B.projection(basisMat);  % reach set projection
% [x0,x0shMat]=x0EllObj_A.double();
%[xB,xB0shMat]=x0EllObj_B.double();
%[ctrMat, ttVec] = rsObj_B.get_center();
%xB=ctrMat(:,end);

% assume lstar. n is number of control variables
l=[0;0;0;0;0;1];
n = 2;
l=l/norm(l);

%%
%Compute int G(t,s) from s=0 to t. G(t,s)=e(t-s)A
% as well as int of sqrt(l'GBQB'G'l)
f = @(s) expm(-Ac*s);
Gint = expm(Ac*T)*integral(f,0,T,'ArrayValued',1);
fc=( xB'*l - l'*expm(Ac*T)*x0 - sqrt(l'*expm(Ac*T)*x0shMat*...
    (expm(Ac*T))'*l ) );
q_c = cell(size(k_v));
Q_c = cell(size(k_v));
for iter=1:length(k_v)
    k=k_v(iter);
    % method 1
    
    if method == 1
        f= @(s) sqrt(l'*expm(Ac*(T-s))*(Bc*Bc')*(expm(Ac*(T-s)))'*l);
        GGint=integral(f,0,T,'ArrayValued',1);
        
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
    elseif method ==2
        
        % method 2
        
        f= @(s) sqrt(l'*expm(Ac*(T-s))*Bc*Q0*Bc'*(expm(Ac*(T-s)))'*l);
        GGint=integral(f,0,T,'ArrayValued',1);
        
        cvx_begin SDP
        cvx_quiet false
        cvx_precision high
        variable q(n)
        variable r
        variable lambda
        lambda > 0
        r > 0
        fc- r*GGint(1)-l'*Gint*Bc*q >=0
        [1-lambda zeros(1,n) (q-q0)';
            zeros(n,1) lambda*eye(n) r*sqrt(Q0);
            q-q0 r*sqrt(Q0) Q0] >=0
        maximize(fc- r*GGint(1)-l'*Gint*Bc*q + k*log_det(r*sqrt(Q0)))
        cvx_end
        q_c{iter}=q;
        Q_c{iter}=r^2*Q0;
        
    end
    
    figure
    plot(ellipsoid(q_c{iter},Q_c{iter}))
    hold on
    plot(ellipsoid(q0,Q0),'y')
end

end