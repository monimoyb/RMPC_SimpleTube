%% Forming Big Matrices along the horizon for Robust Optimization (see: https://www.sciencedirect.com/science/article/pii/S0005109806000021)
%  Monimoy Bujarbaruah
    
function [capA, capE, capB, capC, capD, Aw_batch, Bu_batch, A_batch] = obtain_matR(A,B,C,D,Xn,nx,nu,N,dim_t) 
                
    capA = eye(nx);

    for k=1:N
        capA = [capA;A^k];
    end

    matE(:,:,1) = [eye(nx),zeros(nx,nx*(N-1))];
    capE = [zeros(nx,nx*N);matE(:,:,1)];

    for k=2:N
        matE(:,:,k) = [A^(k-1),matE(:,1:(end-nx),k-1)];
        capE = [capE;matE(:,:,k)];
    end

    capB = capE*kron(eye(N),B); 
    capC = blkdiag(kron(eye(N),C),Xn.A); 
    capD = [kron(eye(N),D);zeros(dim_t-size(kron(eye(N),D),1),nu*N)];


    %% Batch matrices for simulation check
    Aw_batch = capE(nx+1:end,:); 
    Bu_batch = capB(nx+1:end,:);
    A_batch  = capA(nx+1:end,:);
     
end
