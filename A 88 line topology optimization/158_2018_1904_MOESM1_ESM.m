%%%% An 88 LINE PARAMETERIZED LEVEL SET-BASED TOPOLOGY OPTIMIZATION CODE %%%%
function TOPRBF(nelx,nely,volfrac)
%% LEVEL SET FUNCTION INITIALIZATION
    r = nely*0.1;%RADIUS OF INITIAL HOLES
    hX = nelx*[repmat([1/6,5/6],1,3),repmat([0,1/3,2/3,1],1,2),1/2];
    hY = nely*[kron([0,1/2,1],ones(1,2)),kron([1/4,3/4],ones(1,4)),1/2];
    [X,Y] = meshgrid(0:1:nelx,0:1:nely);
    dX = bsxfun(@minus,repmat(X,[1,1,numel(hX)]),reshape(hX,1,1,numel(hX)));
    dY = bsxfun(@minus,repmat(Y,[1,1,numel(hY)]),reshape(hY,1,1,numel(hY)));
    Phi = max(-3,min(3,min(sqrt(dX.^2+dY.^2)-r,[],3))); 
%% RADIAL BASIS FUNCTION INITIALIZATION
    cRBF = 1e-4;%RBF PARAMETER
    nNode = (nely+1)*(nelx+1);
    Ax = bsxfun(@minus,X(:),X(:)');
    Ay = bsxfun(@minus,Y(:),Y(:)');
    A = sqrt(Ax.^2+Ay.^2+cRBF^2);
    G = [A,ones(nNode,1),X(:),Y(:);[ones(1,nNode);X(:)';Y(:)'],zeros(3,3)];
    pGpX = [Ax./A,repmat([0,1,0],nNode,1);repmat([0;1;0],1,nNode),zeros(3,3)];
    pGpY = [Ay./A,repmat([0,0,1],nNode,1);repmat([0;0;1],1,nNode),zeros(3,3)];
    Alpha = G\[Phi(:);0;0;0];
%% FINITE ELEMENT ANALYSIS PREPARATION
    E0 = 1; Emin = 1e-9; nu = 0.3; %MATERIAL PROPERTIES
    A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
    A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
    B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
    B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
    KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);  
    eleN1 = repmat((1:nely)',1,nelx)+kron(0:nelx-1,(nely+1)*ones(nely,1));
    eleNode = repmat(eleN1(:),1,4)+repmat([0,nely+[1,2],1],nelx*nely,1);
    edofMat = kron(eleNode,[2,2])+repmat([-1,0],nelx*nely,4);    
    iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
    jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);    
%% BOUNDARY CONDITION DEFINITION   
    F = sparse(2*((nely+1)*nelx+ceil(nely/2)+1),1,-100,2*nNode,1);%NODAL LOADS    
    fixeddofs = 1:1:2*(nely+1);%DISPLACEMENT CONSTRAINTS
    freedofs = setdiff(1:2*nNode,fixeddofs);    
    U = zeros(2*nNode,1);
%% ITERATION OPTIMIZATION        
    nLoop = 200; nRelax = 30;
    dt = 0.5; delta = 10; mu = 20; gamma = 0.05;
    comp = zeros(nLoop,1); vol = zeros(nLoop,1);
    for iT = 1:nLoop
       %% FINITE ELEMENT ANALYSIS
        [s,t] = meshgrid(-1:0.1:1,-1:0.1:1);
        tmpPhi = (1-s(:)).*(1-t(:))/4*Phi(eleNode(:,1))'+(1+s(:)).*(1-t(:))/4*...
                  Phi(eleNode(:,2))'+(1+s(:)).*(1+t(:))/4*Phi(eleNode(:,3))'+...
                 (1-s(:)).*(1+t(:))/4*Phi(eleNode(:,4))'; 
        eleVol = sum(tmpPhi>=0,1)'/numel(s);
        vol(iT) = sum(eleVol)/(nelx*nely);
        sK = reshape(KE(:)*(Emin+eleVol'*(E0-Emin)),64*nelx*nely,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        U(freedofs,1) = K(freedofs,freedofs)\F(freedofs,1);
        eleComp = sum((U(edofMat)*KE).*U(edofMat),2).*(Emin+eleVol*(E0-Emin));
        comp(iT) = sum(eleComp);
       %% DISPLAY RESULTS
        fprintf('No.%i,  Obj:%f,  Vol:%f\n',[iT,comp(iT),vol(iT)]);
        figure(1); contourf(Phi,[0,0]);
        colormap([0,0,0]); set(gcf,'color','w'); axis equal; axis off;
        figure(2); surf(Phi); caxis([-12,12]); 
        axis equal; axis([0,nelx,0,nely,-12,12]); view(3);
        figure(3); subplot(2,1,1); plot(comp(1:iT),'-'); title('Compliance');
                   subplot(2,1,2); plot(vol(1:iT),'-'); title('Volume fraction');
       %% CONVERGENCE CHECK
        if iT>nRelax && abs(vol(iT)-volfrac)/volfrac<1e-3 && ...
            all(abs(comp(iT)-comp(iT-9:iT-1))/comp(iT)<1e-3)
            break;
        end
       %% LAGRANGE MULTIPLIER
        if iT<=nRelax
            lag = mu*(vol(iT)-vol(1)+(vol(1)-volfrac)*iT/nRelax);
        else
            lag = lag+gamma*(vol(iT)-volfrac);
            gamma = min(gamma+0.05,5);
        end
       %% LEVEL SET FUNCTION EVOLUTION        
        gradPhi = sqrt((pGpX*Alpha).^2+(pGpY*Alpha).^2); 
        indexDelta = (abs(Phi(:))<=delta);
        DeltaPhi = zeros(size(Phi));
        DeltaPhi(indexDelta) = 0.75/delta*(1-Phi(indexDelta).^2/delta^2);   
        eleComp = reshape(eleComp,nely,nelx);
        eleCompLR = [eleComp(:,1),eleComp]+[eleComp,eleComp(:,end)];
        nodeComp = ([eleCompLR;eleCompLR(end,:)]+[eleCompLR(1,:);eleCompLR])/4;
        B = (nodeComp(:)/median(nodeComp(:))-lag).*DeltaPhi(:)*delta/0.75;
        Alpha = Alpha+dt*(G\[B;0;0;0]);
        Alpha = Alpha/mean(gradPhi(unique(eleNode((eleVol<1 & eleVol>0),:))));
        Phi = reshape(G(1:end-3,:)*Alpha,nely+1,nelx+1);
    end
end

% ======================================================================= %
% TOPRBF, a compact and efficient 88-line MATLAB code for the parameterized
% level set method based topology optimization using radial basis functions
% (RBFs), which is applied to minimize the compliance of a two-dimensional
% linear elastic structure. 
% Developed by: Peng Wei, Zuyu Li, Xueping Li and Michael Yu Wang
% Email: ctpwei@scut.edu.cn
%
% Main references:
% (1) P. Wei, Z. Y. Li, X. P. Li and M. Y. Wang, An 88-line MATLAB code for
% the parameterized level set method based topology optimization using
% radial basis functions, Structural and Multidisciplinary Optimization. 
% DOI: 10.1007/s00158-018-1904-8 
%
% (2) S. Y. Wang, K. M. Lim,  B. C. Khoo, and  M. Y. Wang, An extended
% level set method for shape and topology optimization. Journal of
% Computational Physics, 2007, 221(1), 395¨C421.   
%
% **************************   Disclaimer   ***************************** %
% The authors reserve all rights for the programs. The programs may be
% distributed and used for academic and educational purposes. The authors
% do not guarantee that the code is free from errors, and they shall not be
% liable in any event caused by the use of the program.
% ======================================================================= %
