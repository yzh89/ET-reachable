function [q,a]= find_time_varying_control_set(N,T,Ac,Bc,lstar,c_i,M_i,c_x0,M_x0,c_U,M_U)
%n is number of states 
%m is number of control inputs
%N is number of time samples
%T is Estimated Collision Time
%lstar is n*N matrix 
%mu is weighing vector of length N
%compute Sampling Time
[n,m] = size(Bc);
mu = ones(N,1);
dt = T/N;
%size of Pmatrix
Psize = [n,m];
%Intializing  P(i)
blockP  = cell(N);
lpMat = cell(N,1);

rIntMat = zeros(N,N);
for i=1:1:N
    for j=1:1:N 
        blockP{i,j} = zeros(Psize);
    end
    lpMat{i} = zeros(N,m);
end
%computing P(i) & lP & rInt + matrix equivalent
for k=1:1:N
    for i =1:1:k
        rInt_f = @(s) sqrt((lstar(:,k)')*expm(Ac*(k*dt-s))*Bc*M_U*Bc'*(expm(Ac*(k*dt -s)))'*lstar(:,k));
        rIntMat(k,i) = integral(rInt_f,(i-1)*dt,i*dt,'ArrayValued',1);
        
        blockP_f = @(s) expm(Ac*(k*dt-s))*Bc;
        blockP{k,i}= integral(blockP_f,(i-1)*dt,i*dt,'ArrayValued',1);
        
        lpMat{k}(i,:)=lstar(:,k)'*blockP{k,i};
    end
    
end

% function [h_q] = Calculateh(k_hat,Q)
%     h_q = 0;
%     for l=1:1:k_hat
%         h_q = h_q + lstar(:,k_hat)'*blockP{k_hat,l}*Q(:,l);
%     end 
% end
% function [g_a] = Calculateg(k_hat,alpha_hat)
%     g_a = 0;
%     for l=1:1:k_hat
%         r = @(s) (lstar(:,k_hat)')*expm(Ac*(k_hat*dt-s))*Bc*M_U*Bc'*(expm(Ac*(k_hat*dt -s))')*lstar(:,k_hat);
%         g_a = g_a +alpha_hat(l)*integral(r,(l-1)*dt,l*dt,'ArrayValued',1);
%     end 
% end


cvx_begin SDP
cvx_quiet false
cvx_precision high
variable q(m,N)
variable qdiff(m,N-1)
variable lambda(N,1)
variable a(N,1)
a>0.1
for k =1:1:N
    lambda(k) > 0
    [1-lambda(k),    zeros(1,m),      (q(:,k)-c_U')';
        zeros(m,1),  lambda(k)*eye(m), a(k)*sqrt(M_U);
        q(:,k)-c_U', a(k)*sqrt(M_U),   M_U] >=0
    
    -rIntMat(k,:)*a - sum(diag(lpMat{k}*q))-(lstar(:,k)')*c_i{k} ...
        + sqrt((lstar(:,k)')*M_i{k}*lstar(:,k))-lstar(:,k)'*expm(Ac*k*dt)*c_x0 ...
        - sqrt(lstar(:,k)'*expm(Ac*k*dt)*M_x0*expm(Ac*k*dt)*lstar(:,k)) >0
end
for i=1:1:N-1
    qdiff(:,i) == q(:,i+1) -q(:,i);
end
minimize(-mu'*a + 100*norm(qdiff,'fro'))
cvx_end

end
