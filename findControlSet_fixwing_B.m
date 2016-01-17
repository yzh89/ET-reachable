clear;
load result_fixwing_bkup
T=10;

method =1;
n=2;

extEllipA=rsObj_A.cut(T).get_ea;
[x0,x0shMat]=x0EllObj_A.double(); % B initial set is same as A
Q0=shMat;
q0=centVec';

%seperation ball
% d=10;
% sepBall=ellipsoid(zeros(10,1),diag([d^2,d^2,d^2 zeros(1,7)]));
% dirsMatSimple=[[eye(3);zeros(7,3)],rand(10,7)];

l=[zeros(5,1);1];
l=l/norm(l);

% for iter=1:length(extEllipA)
%     temp=[extEllipA(iter),sepBall];
%     extEllipANew(iter,:)=temp.minksum_ea(dirsMatSimple);
% end
% extEllipANew=reshape(extEllipANew,300,1);

temp=extEllipA(1);
for iter=1:length(extEllipA)-1
    temp=temp.intersection_ea(extEllipA(iter+1));
end
extEllipANew=temp;

fc=( l'*expm(Ac*T)*x0 -l'*extEllipANew.getCenterVec- sqrt(l'*extEllipANew.getShapeMat*l) ...
    -  sqrt(l'*expm(Ac*T)*x0shMat*(expm(Ac*T))'*l ) );

f = @(s) expm(Ac*(T-s));
Gint = integral(f,0,T,'ArrayValued',1);


%% method 1
if method==1
f= @(s) sqrt(l'*expm(Ac*(T-s))*Bc*Bc'*(expm(Ac*(T-s)))'*l);
GGint=integral(f,0,T,'ArrayValued',1);
%fc- GGint(1)*norm(Q,2)+l'*Gint*Bc*q >-1000000
cvx_begin SDP
cvx_quiet false
cvx_precision high
    variable q(n)
    variable Q(n,n);
    variable lambda
    fc- norm(Q,2)*GGint(1)+l'*Gint*Bc*q > 3
    lambda > 0
    [1-lambda zeros(1,n) (q-q0)';
        zeros(n,1) lambda*eye(n) Q;
        q-q0 Q Q0] >=0
    maximize(log_det(Q))
cvx_end
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
fc- r*GGint(1)+l'*Gint*Bc*q > -10
[1-lambda zeros(1,n) (q-q0)';
    zeros(n,1) lambda*eye(n) r*sqrt(Q0);
    q-q0 r*sqrt(Q0) Q0] >=0
maximize(log_det(r*sqrt(Q0)))
cvx_end
Q = r*sqrt(Q0);
end
%%
% plot the result control set
figure
plot(ellipsoid(q0,Q0),'y')
hold on
plot(ellipsoid(q,Q'*Q));
plot(uBoundsEllObj_cu,'b');

%%
elltool.setconf('accurate')
 uBoundsEllObj_cu_B = ellipsoid(q, Q'*Q);
lsys_cu_B = elltool.linsys.LinSysContinuous(Ac, Bc, uBoundsEllObj_cu_B);  
%dirsMat = [[eye(3);zeros(7,3)],rand(10,7)];
rsObj_B = elltool.reach.ReachContinuous(lsys_cu_B, x0EllObj_A, dirsMat, timeVec,...
    'isRegEnabled', true, 'isJustCheck', false, 'regTol', 1e-5); 
%%
basisMat = [1 zeros(1,5); zeros(1,5) 1]';  % orthogonal basis of (x1, x2) subspace
psObj_B = rsObj_B.projection(basisMat);  % reach set projection

s_=size(x_B_star);
psObj_copy_B = psObj_B.getCopyWithCenterModified(x_B_star(1:s_(1)-1,[1,3]));

plotB=psObj_copy_B.cut(T_end);
plObj=plotA.plotByEa('g'); % to have the use of plObj isn't necessary
hold on
plotA.plotByIa('r',plObj);
plotB.plotByEa('y',plObj);
plotB.plotByIa('b',plObj);

% to have the use of plObj isn't necessary
plObj=psObj_copy.plotByEa('g');  % external apprx. of reach set 1 (red)
hold on
psObj_copy.plotByIa('r',plObj);  % internal apprx. of reach set 1 (green)
psObj_copy_B.plotByEa('y',plObj);
psObj_copy_B.plotByIa('b',plObj);
plot3(0:0.1:10,x_star(:,1)',x_star(:,3)');
plot3(0:0.1:10,x_B_star(:,1)',x_B_star(:,3)','k--');
%%
save result_fixwing_final