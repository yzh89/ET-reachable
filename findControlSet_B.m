clear;
load result_bkup
T=4;
extEllipA=rsObj_A.cut(T).get_ea;
[x0,x0shMat]=x0EllObj_B.double();
Q0=shMat;
q0=centVec';


%seperation ball
d=1;
sepBall=ellipsoid(zeros(10,1),diag([d^2,d^2,d^2 zeros(1,7)]));
dirsMatSimple=[[eye(3);zeros(7,3)],rand(10,7)];

l=[0;1;0;zeros(7,1)];
l=l/norm(l);


for iter=1:length(extEllipA)
    temp=[extEllipA(iter),sepBall];
    extEllipANew(iter,:)=temp.minksum_ea(dirsMatSimple);
end
extEllipANew=reshape(extEllipANew,300,1);
temp=extEllipANew(1);
for iter=1:length(extEllipANew)-1
    temp=temp.intersection_ea(extEllipANew(iter+1));
end
clear extEllipANew;
extEllipANew=temp;
%%
fc=( l'*expm(Ac*T)*x0 -l'*extEllipANew.getCenterVec- sqrt(l'*expm(Ac*T)*x0shMat*...
    (expm(Ac*T))'*l )- sqrt(l'*extEllipANew.getShapeMat*l ));

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
N=i-1;
Gint=expm(Ac*T)*Gtemp;
Gtemp=0;
f= @(s) sqrt(l'*expm(Ac*(T-s))*Bc*Bc'*(expm(Ac*(T-s)))'*l);
GGint=integral(f,0,T,'ArrayValued',1);

%%
n=3;
cvx_begin SDP
cvx_quiet false
cvx_precision high
    variable q(n)
    variable Q(n,n) semidefinite;
    variable lambda
    
    %fc- GGint(1)*norm(Q,2)-l'*Gint*Bc*q >=0
    
    lambda > 0
    [1-lambda zeros(1,n) (q-q0)';
        zeros(n,1) lambda*eye(n) Q;
        q-q0 Q Q0] >=0
    fc- GGint(1)*norm(Q,2)+l'*Gint*Bc*q >0
    maximize(log_det(Q))
cvx_end

%%
% plot the result control set
figure
plot(ellipsoid(q0,Q0),'y')
hold on
plot(ellipsoid(q,Q'*Q));
plot(uBoundsEllObj_cu,'b');

%%
elltool.setconf('accurate')
 uBoundsEllObj_cu = ellipsoid(q, Q'*Q);
lsys_cu = elltool.linsys.LinSysContinuous(Ac, Bc, uBoundsEllObj_cu);  
%dirsMat = [[eye(3);zeros(7,3)],rand(10,7)];
rsObj_B = elltool.reach.ReachContinuous(lsys_cu, x0EllObj_B, dirsMat, timeVec,...
    'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-9); 
%%
basisMat = [1 zeros(1,9); 0 1 zeros(1,8)]';  % orthogonal basis of (x1, x2) subspace
psObj_A = rsObj_A.projection(basisMat);  % reach set projection
psObj_B = rsObj_B.projection(basisMat);  % reach set projection
% plot projection of reach set external approximation:

% to have the use of plObj isn't necessary 
plObj=psObj_A.plotByEa('r');  % external apprx. of reach set 1 (red)
%firstPsObj.plotByEa('r');
hold on
psObj_A.plotByIa('g',plObj);  % internal apprx. of reach set 1 (green)
psObj_B.plotByEa('y',plObj);
psObj_B.plotByIa('b',plObj);
save result_final