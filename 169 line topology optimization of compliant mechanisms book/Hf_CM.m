% A 169 LINE MATLAB CODE FOR TOPOLOGY OPTIMIZATION
% OF COMPLIANT MECHANISMS USING A WEIGHTTING SUM
% METHOD WITH SELF - ADJUSTED WEIGHTING FACTORS
function Hf_CM(width , hight , nelx , nely , maxloop ,l_vol)
%% PARAMETERS INITIALIZATION
E1=1; E0 = 1e-3; nu = 0.3;
ew = width/nelx; eh = hight/nely;
 [X,Y] = meshgrid (ew*[-0.5:nelx+0.5],eh*[-0.5:nely +0.5]);
 LSgrid_x = X(:); LSgrid_y = Y(:);
 [M,N] = meshgrid (ew*[0:nelx],eh*[0:nely]);
 FENgrid_x = M(:); FENgrid_y = N(:);
 %% INITIAL CONFIGUATION
phi = init_config(width , hight , nely , nelx , ew , eh , ...
LSgrid_x , LSgrid_y);
 phi_fea = griddata (LSgrid_x , LSgrid_y , phi(:), FENgrid_x , ...
ENgrid_y , 'cubic') ;
 loop = 1;
 alp = 0; bta = 0;
 %% MAIN OPTIMIZATION LOOP
 while loop ≤ maxloop
 %% finite element analysis
 [U,V,C1 ,C2 ,Uout] = FEA(nely , nelx , E1 , E0 , ew , eh , ...
nu , reshape (phi_fea ,nely+1,nelx+1));
 %% construct the velocity
 vel = velocity(E1 , E0 , U, ew , eh , nely , nelx , l_vol , ...
alp , bta , phi , nu);
 %% update the level set equation
 phi = level_set_equation(phi , ew , eh , vel);
 %% re -initialize the level set function to signed ...
distance function
 if mod(loop ,5) ==0; phi = reinit_SDF(phi , ew , eh); end;
 %% calculate the weighting factors of Cin and Cout

 alp = abs(Uout/C1); bta = abs(Uout/C2);
 %% print results
 disp ([ 'It.:' sprintf ('%4i',loop) 'Uout:' ...
sprintf ('%10.2f',full (Uout)) 'V_Ratio:' ...
sprintf ('%10.2f',V/(width*hight)) ])
 %% plot topologies
 phi_fea = griddata (LSgrid_x ,LSgrid_y ,phi(:), ...
FENgrid_x ,FENgrid_y ,'cubic') ;
 colormap ( gray ) ; ...
contourf( reshape (phi_fea ,nely+1,nelx+1) ,[0 0]); ...
axis equal; axis tight; axis off; pause (1e-6);
 loop = loop+1;
 end
 %% SUBFUNCTIONS
 function [phi] = init_config(width , hight , nely , nelx , ...
ew , eh , LSgrid_x , LSgrid_y)
 % Preset an initial configuration with several holes
 cx = width*[ 2/16 4/16 6/16 8/16 10/16 12/16 14/16 ...
 0 2/16 4/16 6/16 8/16 10/16 12/16 14/16 1 ...
 0 2/16 4/16 6/16 8/16 10/16 12/16 14/16 1 ...
 0 2/16 4/16 6/16 8/16 10/16 12/16 14/16 1 ...
 2/16 4/16 6/16 8/16 10/16 12/16 14/16];
 cy = hight*[ 0 0 0 0 0 0 0 ...
 1/4 1/4 1/4 1/4 1/4 1/4 1/4 1/4 1/4 ...
 2/4 2/4 2/4 2/4 2/4 2/4 2/4 2/4 2/4 ...
 3/4 3/4 3/4 3/4 3/4 3/4 3/4 3/4 3/4 ...
 1 1 1 1 1 1 1];
 for i = 1: length (cx)
 tmpPhi(:,i) = -sqrt ((LSgrid_x -cx(i)).^2 + ...
(LSgrid_y -cy(i)) .^2) + hight/11;
 end
 TmpPhi = -( max(tmpPhi.')).';
 phi = reshape (TmpPhi ,nely+2,nelx+2);
 phi = reinit_SDF(phi ,ew ,eh);
 function [velocity] = velocity(E1 , E0 , U, Ew , Eh , nely , ...
nelx , L_Vol , c1 , c2 , phi , nu)
 % Construct the velocity field for updating the
 % level set function phi
 phi = phi(2:nely+1,2:nelx+1);
 Dsub = (1-nu^2)*[1, nu , 0
 nu , 1, 0
 0, 0, (1-nu)/2];
 B = 1/2*[ -1/Ew 0 1/Ew 0 1/Ew 0 -1/Ew 0
 0 -1/Eh 0 -1/Eh 0 1/Eh 0 1/Eh
 -1/Eh -1/Ew -1/Eh 1/Ew 1/Eh 1/Ew 1/Eh -1/Ew];
 for ely = 1:nely
 for elx = 1:nelx
 n1 = (nely+1)*(elx -1)+ely;
 n2 = (nely+1)*elx +ely;
 Ue1 = U([2*n1 -1; 2*n1; 2*n2 -1; 2*n2; 2*n2+1; ...
2*n2+2; 2*n1+1; 2*n1+2],1);
 Ue2 = U([2*n1 -1; 2*n1; 2*n2 -1; 2*n2; 2*n2+1; ...
2*n2+2; 2*n1+1; 2*n1+2],2);
 Stain1 = B*Ue1;
 Stain2 = B*Ue2;
 if phi(ely ,elx) > 0.75*Ew
 ERatio = E1;
 elseif phi(ely ,elx) < -0.75*Ew
 ERatio = E0;
 else
 Value = phi(ely ,elx)/(0.75*Ew);
 ERatio = E1*0.75*(1 - 1e-3)*(Value - ...
Value^3/3) + 0.5*(1+1e-3);
 end
 D = ERatio*Dsub;
 Beta_u1_out(ely ,elx) = Stain2 '*D* Stain1;
 Beta_u1_in(ely ,elx)= Stain1 '*D* Stain1;
 Beta_u2_out(ely ,elx)= Stain2 '*D* Stain2;
 end
 end
 Beta = L_Vol - c1 * Beta_u1_in - c2 * Beta_u2_out + ...
Beta_u1_out;
 velocity(1: nely+2,1:nelx+2) = 0;
 velocity(2: nely+1,2:nelx+1) = Beta/max( max( abs(Beta )));
 velocity = velocity .* exp (1- abs(velocity));
 function [phi]= level_set_equation(phi , dx , dy , Vn)
 % Solve the level set equation using the ENO2 scheme
 Vn_ext = init_normal_ENO2(Vn);
 it = 0;
 t = 0;
 while (it < 5)
 [delt_normal ,H1_abs , H2_abs] = evolve_normal_ENO2(phi , ...
dx , dy , Vn_ext);
 dt = get_dt_normal (0.5,dx,dy ,H1_abs , H2_abs);
 phi = phi -dt*delt_normal;
 it = it+1;
 t = t+dt;
 end
 function [phi] = reinit_SDF(phi , dx , dy)
 % Re -initialize the level set function
 % to be a signed distance function
 s_phi_0 = phi./ sqrt (phi.^2+dx.^2);
 Vn_ext = init_normal_ENO2( s_phi_0);
 it=0;
 while (it < 30)
 [delt_normal ,H1_abs , H2_abs] = evolve_normal_ENO2(phi , ...
dx , dy , Vn_ext);
dt = get_dt_normal (0.5,dx,dy ,H1_abs , H2_abs);
 phi = phi+dt*(s_phi_0 - delt_normal);
 it = it+1;
 end
 function [U,V,u1_in ,u2_out , u1_out] = FEA(nely , nelx , E1 , ...
E0 , ew , eh , nu , phi)
 % Finite element analysis
 K = sparse (2*(nely+1)*(nelx+1), 2*(nely+1)*(nelx+1));
 F = sparse (2*(nely+1)*(nelx+1), 2);
 U = sparse (2*(nely+1)*(nelx+1), 2);
 AreaFactor_V(1: nely ,1: nelx) = 0;
 for ely = 1:nely
 for elx = 1:nelx
 [KE ,AreaFactor] = Real_ES(E1 , E0 , ew , eh , nu , ...
phi(ely:ely+1,elx:elx+1));
 AreaFactor_V(ely ,elx) = AreaFactor;
 n1 = (nely+1)*(elx -1)+ely;
 n2 = (nely+1)*elx +ely;
 edof = [2*n1 -1; 2*n1; 2*n2 -1; 2*n2; 2*n2+1; ...
2*n2+2; 2*n1+1; 2*n1+2];
 K(edof , edof) = K(edof , edof) + KE;
 end
 end
 din = 2*(nely+1) -1;
 dout = 2*(nely+1)*(nelx+1) -1;
 F(din ,1) = 0.1;
 F(dout ,2) = -0.1;
 fixeddofs = ...
union ([2*(nely+1):2*(nely+1):2*(nely+1)*(nelx+1)] ,[1:1:4]);
 alldofs = [1:2*( nely+1)*(nelx+1)];
 freedofs = setdiff( alldofs , fixeddofs);
 U(freedofs ,:) = K(freedofs , freedofs) \ F( freedofs ,:);
 U(fixeddofs ,:)= 0;
 u1_in = U(din , 1);
 u1_out = U(dout , 1);
 u2_out = U(dout , 2);
 V = sum( sum(AreaFactor_V))*ew*eh;
 function [RealKE , AreaFactor ]= Real_ES(E1 , E0 , Ew , Eh , nu , phi)
 % Return the real element stiffness matrix
 if min( min(phi)) > 0
 E = E1; AreaFactor = 1;
 elseif max( max(phi)) < 0
 E = E0; AreaFactor = 0;
 else
 [s,t] = meshgrid ([-1:0.1:1],[-1:0.1:1]);
 tmpphi = (1-s(:)).*(1-t(:))/4*phi(1,1) + ...
(1+s(:)).*(1-t(:))/4*phi(1,2) + ...
(1+s(:)).*(1+t(:))/4*phi(2,1) + ...
(1-s(:)).*(1+t(:))/4*phi(2,2);
 E = length ( find (tmpphi ≥0))/ length (s(:))*E1;
 AreaFactor = length ( find (tmpphi ≥0))/ length (s(:));
 end
 RealKE = Basic_ES(E,nu ,Ew ,Eh);
 function [KE] = Basic_ES(E, nu , a, b)
 % Basic element stiffness matrix
 k = [ -1/6/a/b*(nu*a^2-2*b^2-a^2), 1/8*nu +1/8, ...
-1/12/a/b*(nu*a ^2+4*b^2-a^2), 3/8*nu -1/8, ...
 1/12/a/b*(nu*a^2-2*b^2-a^2), -1/8*nu -1/8, ...
1/6/a/b*(nu*a^2+b^2-a^2), -3/8*nu +1/8];
 KE = E/(1-nu^2)*...
 [ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
 k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
 k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
 k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
 k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
 k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
 k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
 k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];