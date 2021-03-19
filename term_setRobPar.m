%% Calculating the Terminal Set X_N of Robust MPC
% Monimoy Bujarbaruah

function [Xn,Pinf] = term_setRobPar(Anom, Bnom, delAv, delBv, K, X, U, W, Q, R, nx, nu)

Pinf = dlyap((Anom-Bnom*K)',Q+K'*R*K);                         %  Cost based on dlyap as K is fixed 

%%% extract the set of possible matrices%%%% 
setA = zeros(nx, nx, size(delAv,2)/nx);
setB = zeros(nx, nu, size(delBv,2)/nu);

for i = 1:size(delAv,2)/nx
    setA(:,:,i) = Anom + delAv(:,(i-1)*nx + 1: i*nx);  
end

for i = 1:size(delBv,2)/nu
    setB(:,:,i) = Bnom + delBv(:,(i-1)*nu + 1: i*nu);  
end

%%%%%%%% form all possible pairs %%%%%%
cl_mat = zeros(nx,nx, (size(delAv,2)/nx)*(size(delBv,2)/nu));
count = 0;
for i = 1:(size(delAv,2)/nx)
    for j = 1: size(delBv,2)/nu
        count = count + 1;
        cl_mat(:,:,count) = setA(:,:,i) - setB(:,:,j)*K;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%

%%% check stability
for i = 1: size(cl_mat,3)
   max_abs(i) =  max(abs(eig(cl_mat(:,:,i)))); 
end
%%%%%%

S = Polyhedron('A',X.A,'b',X.b); 
S = S.intersect(Polyhedron('H',[-U.H(:,1:nu)*K U.H(:,nu+1)])); 

Xn = MRPI_Par(cl_mat, S, Polyhedron.emptySet(nx), W);

end
