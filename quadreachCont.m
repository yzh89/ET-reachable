%function quadreach()
g = 9.81;       % gravity (m/s^2)
m = 1;          % mass (kg)
L = 0.25;       % length of rotor arm (m)
J = 20*m*L^2;   % moment of inertia (kg m^2) (est)

elltool.setconf('accurate')
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
Ku = 0.04*[3 0 1;  
           3 1 1];    
       
 centVec = -diff(Ku)/2;
 shMat = diag((Ku(1,:)-centVec).^2,0);
 uBoundsEllObj = ellipsoid(centVec', shMat);
 %uBoundsEllObj=uBoundsEllObj.getShape(1.5);
% plot(uBoundsEllObj)
% linear system
lsys = elltool.linsys.LinSysContinuous(Ac, Bc, uBoundsEllObj);  

timeVec = [0 4];  % time interval% initial conditions:
x0EllObj = [zeros(3,1); 0.2 ;zeros(6,1)] + ellipsoid(diag([0.01,0.01,0.01,0,0,0, zeros(1,4)]));
% initial directions (some random vectors in R^4):
%dirsMat = [eye(10,10),rand(10,20)];
load('result.mat','dirsMat')
% reach set
rsObj = elltool.reach.ReachContinuous(lsys, x0EllObj, dirsMat, timeVec,...
    'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-7);  
%%
basisMat = [1 zeros(1,9); 0 1 zeros(1,8); 0 0 1 zeros(1,7)]';  % orthogonal basis of (x1, x2) subspace
psObj = rsObj.projection(basisMat);  % reach set projection
psObj = psObj.cut(4);
% plot projection of reach set external approximation:

% plot the whole reach tube:
plObj=psObj.plotByEa('g'); % to have the use of plObj isn't necessary 
hold on
psObj.plotByIa('r',plObj);
%%
basisMat = [1 zeros(1,9); 0 1 zeros(1,8)]';  % orthogonal basis of (x1, x2) subspace
psObj = rsObj.projection(basisMat);  % reach set projection
% plot projection of reach set external approximation:

% to have the use of plObj isn't necessary 
plObj=smartdb.disp.RelationDataPlotter('figureGroupKeySuffFunc', ...
    @(x)sprintf('_forward_reach_set_proj1%d',x));
psObj.plotByEa('r',plObj);  % external apprx. of reach set 1 (red)
%firstPsObj.plotByEa('r');
hold on
psObj.plotByIa('g',plObj);  % internal apprx. of reach set 1 (green)
%firstPsObj.plotByIa('g');
   
   save result
% 
% end
% 
% function [projOrthMatArray,projOrthMatTransArray]=fGetProjMat(projMat,...
%     timeVec,varargin)
%   nTimePoints=length(timeVec);
%   projOrthMatArray=repmat(projMat,[1,1,nTimePoints]);
%   projOrthMatTransArray=repmat(projMat.',[1,1,nTimePoints]);
%  end