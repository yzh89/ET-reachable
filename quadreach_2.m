%function quadreach()
clear
elltool.setconf('accurate')
use_old=true;
if use_old
    load result
else
g = 9.81;       % gravity (m/s^2)
m = 1;          % mass (kg)
L = 0.25;       % length of rotor arm (m)
J = 20*m*L^2;   % moment of inertia (kg m^2) (est)

%%
tmp = [0 g; -g 0; 0 0];
Ac = [zeros(3),   eye(3),     zeros(3,2), zeros(3,2); 
      zeros(3),   zeros(3),   tmp,        zeros(3,2); 
      zeros(2,3), zeros(2,3), zeros(2),   eye(2);
      zeros(2,3), zeros(2,3), zeros(2),   zeros(2)];

tmp = [0; 0; 1/m];
Bc = [zeros(3,1) zeros(3,2); 
      tmp        zeros(3,2); 
      zeros(2,1) zeros(2,2); 
      zeros(2,1) L/J*eye(2)];

%Ku = 0.04*[9.935 3.62 3.62; 4.545 3.62 3.62];
Ku = 0.04*[9.935 3.62 3.62;  
           4.545 3.62 3.62];    
       
 centVec = -diff(Ku)/2;
 shMat = diag((Ku(1,:)-centVec).^2,0);
 uBoundsEllObj = ellipsoid(centVec', shMat);
 %uBoundsEllObj=uBoundsEllObj.getShape(1.5);
% plot(uBoundsEllObj)
% linear system
lsys = elltool.linsys.LinSysContinuous(Ac, Bc, uBoundsEllObj);  

timeVec = [0 4];  % time interval% initial conditions:
x0EllObj_A = [zeros(3,1); 0.2 ;zeros(6,1)] + ellipsoid(diag([0.01,0.01,0.01,0,0,0, zeros(1,4)]));
x0EllObj_B= [1.6;0.5;0 ; -0.2 ;zeros(6,1)] + ellipsoid(diag([0.01,0.01,0.01,0,0,0, zeros(1,4)]));

% initial directions (some random vectors in R^4):
% dirsMat = [eye(10,10),rand(10,20)];
%dirsMat=[[eye(3);zeros(7,3)],rand(10,7)];
load('result.mat','dirsMat')
% reach set
rsObj_A = elltool.reach.ReachContinuous(lsys, x0EllObj_A, dirsMat, timeVec,...
    'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-7);  
rsObj_B = elltool.reach.ReachContinuous(lsys, x0EllObj_B, dirsMat, timeVec,...
    'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-7); 

%project to see reachable set at end time
basisMat = [1 zeros(1,9); 0 1 zeros(1,8); 0 0 1 zeros(1,7)]';  % orthogonal basis of (x1, x2) subspace
psObj_A = rsObj_A.projection(basisMat);  % reach set projection
psObj_A = psObj_A.cut(4);
psObj_B = rsObj_B.projection(basisMat);  % reach set projection
psObj_B = psObj_B.cut(4);
% plot projection of reach set external approximation:

% plot the whole reach tube:
plObj=psObj_A.plotByEa('g'); % to have the use of plObj isn't necessary 
hold on
psObj_A.plotByIa('r',plObj);
psObj_B.plotByEa('y',plObj);
psObj_B.plotByIa('b',plObj);

%project to see the whoe reachable set for x y dimension
basisMat = [1 zeros(1,9); 0 1 zeros(1,8)]';  % orthogonal basis of (x1, x2) subspace
psObj_A = rsObj_A.projection(basisMat);  % reach set projection
psObj_B = rsObj_B.projection(basisMat);  % reach set projection
% plot projection of reach set external approximation:

% to have the use of plObj isn't necessary 
plObj=psObj_A.plotByEa('r');  % external apprx. of reach set 1 (red)
hold on
psObj_A.plotByIa('g',plObj);  % internal apprx. of reach set 1 (green)
psObj_B.plotByEa('y',plObj);
psObj_B.plotByIa('b',plObj);

save result
end

[x0,x0shMat]=x0EllObj_A.double();
[ctrMat, ttVec] = rsObj_B.get_center();
xB=ctrMat(:,end);
[qc,Qc]=findControlSet_v2(x0,x0shMat,xB,Ac,Bc,centVec',shMat);


%%
%for iter=1:length(qc)
iter=6;
    q=qc{iter};
    Q=Qc{iter};
 uBoundsEllObj_cu = ellipsoid(q, Q);
 %uBoundsEllObj=uBoundsEllObj.getShape(1.5);
% plot(uBoundsEllObj)
% linear system
lsys_cu = elltool.linsys.LinSysContinuous(Ac, Bc, uBoundsEllObj_cu);  
%dirsMat = [[eye(3);zeros(7,3)],rand(10,7)];
rsObj_A = elltool.reach.ReachContinuous(lsys_cu, x0EllObj_A, dirsMat, timeVec,...
    'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-8); 

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
%end

save result_bkup