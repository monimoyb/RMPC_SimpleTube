%% Maximal RPI Set Computation for parametric uncertainty 
% Monimoy Bujarbaruah


function P = MRPI_Par(mat_set, X, U,  W)
% computation of a robust invariant set for given LTImodel.  
%
maxIterations = 100;
X0 = X;                     % initial set constraint

for j = 1:maxIterations
    % subtract noise    
    S = X0 - W;
    for i = 1: size(mat_set,3)  % for ALL models 
         model = LTISystem('A',mat_set(:,:,i)); 
         % backward reachable set
         Rm(i) = model.reachableSet('X', S, 'U', U,...
                         'direction', 'backward');
    end
    R = Rm(1); 
    
    for i = 2:size(mat_set,3)
         R = intersect(R,Rm(i)); % take COMMON
    end
    
    % intersect with the state constraints
    P = R.intersect(X0);

    if P==X0
%         j
        break
    else
        X0 = P;
    end
end
