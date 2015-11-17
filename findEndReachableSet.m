clear;
load result_bkup
T=4;
extEllipA=rsObj_A.cut(T).get_ea;
extEllipB=rsObj_B.cut(T).get_ea;

%seperation ball
d=0.5;
sepBall=ellipsoid(zeros(3,1),diag([d^2,d^2,d^2]));
dirsMatSimple=[eye(3),rand(3,3)];

l=[0;1;0];
l=l/norm(l);

for iter=1:length(extEllipA)
    temp=[extEllipA(iter).projection([1 zeros(1,9); 0 1 zeros(1,8); 0 0 1 zeros(1,7)]'),sepBall];
    extEllipANew(iter,:)=temp.minksum_ea(dirsMatSimple);
end
extEllipANew=reshape(extEllipANew,180,1);
temp=extEllipANew(1);
for iter=1:length(extEllipANew)-1
    temp=temp.intersection_ea(extEllipANew(iter+1));
end
clear extEllipANew;
extEllipANew=temp;

temp=extEllipB(1);
for iter=1:length(extEllipB)-1
    temp=temp.intersection_ea(extEllipB(iter+1));
end
extEllipBOld=temp;
%%
C=[eye(3),zeros(3,7)];
n=3;
cvx_begin SDP
cvx_quiet false
cvx_precision high
    variable q(n)
    variable Q(n,n) semidefinite;
    variable lambda(length(extEllipB))
    
    %fc- GGint(1)*norm(Q,2)-l'*Gint*Bc*q >=0
    for iter=1:length(extEllipB)
        lambda(iter) > 0
    [1-lambda(iter) zeros(1,n) (q-C*extEllipB(iter).getCenterVec)';
        zeros(n,1) lambda(iter)*eye(n) Q;
        q-C*extEllipB(iter).getCenterVec Q C*extEllipB(iter).getShapeMat*C'] >=0
    end
    l'*q-norm(Q*l)-l'*extEllipANew.getCenterVec-sqrt(l'*extEllipANew.getShapeMat*l)>=0
    maximize(log_det(Q))
cvx_end

qc=q;
Qc=q
%%
n=10;
cvx_begin SDP
cvx_quiet false
cvx_precision high
    variable q(n)
    variable Q(n,n) semidefinite;
    variable lambda(length(extEllipB))
    
    %fc- GGint(1)*norm(Q,2)-l'*Gint*Bc*q >=0
    for iter=1:length(extEllipB)
        lambda(iter) > 0
    [1-lambda(iter) zeros(1,n) (q-extEllipB(iter).getCenterVec)';
        zeros(n,1) lambda(iter)*eye(n) Q;
        q-extEllipB(iter).getCenterVec Q extEllipB(iter).getShapeMat] >=0
    end
    l'*C*q-norm(C*Q*C'*l)-l'*extEllipANew.getCenterVec-sqrt(l'*extEllipANew.getShapeMat*l)>0
    maximize(log_det(C*Q*C'))
cvx_end

%%
basisMat = [1 zeros(1,9); 0 1 zeros(1,8); 0 0 1 zeros(1,7)]';  % orthogonal basis of (x1, x2) subspace
psObj_B = rsObj_B.projection(basisMat);  % reach set projection
psObj_B = psObj_B.cut(4);
% plot projection of reach set external approximation:

% plot the whole reach tube:
psObj_B.plotByEa('y'); % to have the use of plObj isn't necessary 
hold on
plot(extEllipANew)
plot(ellipsoid(qc,Qc'*Qc));
%%
targetSetCenter=[qc;q(4:end)];
targetSetShapeMat=[Qc,Q(1:3,4:10);Q(4:10,1:10)];
%%
timeVecBk = [4 0];
rsBkObj_B = elltool.reach.ReachContinuous(lsys, ellipsoid(targetSetCenter,targetSetShapeMat), dirsMat, timeVecBk,...
    'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-7);  
%%
basisMat = [1 zeros(1,9); 0 1 zeros(1,8)]';  % orthogonal basis of (x1, x2) subspace
psObj_A = rsObj_A.projection(basisMat);  % reach set projection
psObj_B = rsBkObj_B.projection(basisMat);  % reach set projection
% plot projection of reach set external approximation:

% to have the use of plObj isn't necessary 
plObj=psObj_A.plotByEa('r');  % external apprx. of reach set 1 (red)
%firstPsObj.plotByEa('r');
hold on
psObj_A.plotByIa('g',plObj);  % internal apprx. of reach set 1 (green)
psObj_B.plotByEa('y',plObj);
psObj_B.plotByIa('b',plObj);