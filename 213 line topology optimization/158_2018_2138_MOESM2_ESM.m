function TOLD
%initialize---system parameters
dbstop if error
loop=0;wd=cd;
%initialize---material parameters
E0=1e8;nu=0.3;penal=3;c2=5e-4;strain0=0.6;
%initialize---boundary condition and mesh
[nele,nelx,nely,thickness,gobalcoords0,edofMat,fixeddof,loaddof,F,ndof,Kin,Kout,outdof]=initialize2Dinverter;
%initialize---filiter
rmin=1.8;
[H,Hs]=preparefilter(rmin,nele,nelx,nely);
volfrac=0.25;
x = repmat(volfrac,[nele,1]);
xba = (H*x(:))./Hs;  
%iinitialize--MMA parameters
changeobj=1;m=1;a0=1;amma=0;cmma=10000;
d=0;xold1=x;xold2=x;dcdxba(nele,1)=0;
xmax=ones(nele,1);xmin  = 0.001*ones(nele,1);
low=xmin;upp=xmax;
tt=zeros(3,100);
while changeobj>1e-4 && loop < 100   
    tic
    loop = loop+1;   
    %generatecommand 
    generatecommand(F,nele,ndof,thickness,gobalcoords0,E0,nu,xba,edofMat,fixeddof,loaddof,penal,c2,Kin,Kout,outdof);   
    tt(1,loop)=toc;
    %run ANSYS APDL as a subroutine
    cd 'D:\software\ansys14\ANSYS Inc\v140\ansys\bin\winx64'
    !ANSYS140 -b -p ane3fl -i D:\command.txt -o D:\1.out  
    tt(2,loop)=toc;
    cd(wd)    
    %read nodal displacement solution
    U=load('D:\NODALDISPLACEMENT.txt');
    delete('D:\NODALDISPLACEMENT.txt');
    U=U';U=U(:);       
    obj(loop)=U(outdof);
    %read the element logarithmic equivilant strain and update c2
    strain=load('D:\eqvstrain.txt');
    strain=reshape(strain,2,[]);
    maxvms=max(strain');svon=maxvms(2);
    if svon<=strain0
        c2=max((svon/strain0)^0.5,0.8)*c2;
    else
        c2=(svon/strain0)^2*c2;        
    end
    c2=max(c2,1e-6);
    %read the element strain energy 
    strain_energy=load('D:\strainenergy.txt');
    strain_energy99=load('D:\strainenergy99.txt');    
    %sensitivity analysis
    dstrain_energy=reshape(strain_energy99-strain_energy,2,[]);    
    for i=1:nele
    dcdxba(i)=-1/(F(loaddof)/100)*(dstrain_energy(1,i)*penal/xba(i)-...
        dstrain_energy(2,i)*penal*xba(i)^(penal-1)/(1-xba(i)^penal+1e-10));
    end
    dvdxba=ones(nely,nelx);         
    dcdx = H*(dcdxba(:)./Hs);                                                     
    dvdx = H*(dvdxba(:)./Hs);      
    %save data
    lp=num2str(loop);       
    save(lp);         
    %update the design variables by MMA
    ratemma=11000;
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
    mmasub(m,nele,loop,x,xmin,xmax,xold1,xold2,ratemma(1)*obj(loop),ratemma(1)*dcdx,...
    0*dcdx,sum(xba)-nele*volfrac,dvdx',0*dvdx',low,upp,a0,amma,cmma,d);    
    xold2 = xold1;
    xold1 = x;
    x = xmma;  
    xba(:) = (H*x(:))./Hs;       
    if loop>=6
    changeobj=mean((abs(obj(loop-4:loop)-obj(loop-5:loop-1))))/abs(obj(loop));
    end
    fprintf('It.:%5i Obj.:%11.4f Vol.:%2.4f ch.:%2.5f\n',loop,full(obj(loop)),mean(xba),changeobj);          
    colormap(gray); imagesc(1-reshape(xba,nely,nelx)); axis equal; axis tight; axis off;pause(1e-6);  
    tt(3,loop)=toc;
end
plot(1:loop,obj);
end
function [nele,nelx,nely,cc,gobalcoords0,edofMat,fixeddof,loaddof,F,ndof,Kin,Kout,outdof]=initialize2Dinverter
nelx=100;nely=50;
aa=100e-3/nelx;bb=50e-3/nely;cc=1e-3;
[i0,j0] = meshgrid(0:nelx,0:nely);
gobalcoords0x=i0*aa;
gobalcoords0y=bb*nely-j0*bb;
gobalcoords0xy=[gobalcoords0x(:) gobalcoords0y(:)]';
gobalcoords0=gobalcoords0xy(:);
[il,jl] = meshgrid(0,nely);
loadnid = il*(nely+1)+(nely+1-jl); 
loaddof = 2*loadnid(:)-1; 
force=5;
Kin=500;
[ilo,jlo] = meshgrid(nelx,nely);
outnid = ilo*(nely+1)+(nely+1-jlo);  
outdof = 2*outnid(:)-1; 
Kout=100;
[iif,jf] = meshgrid(0,0:nely/5);
fixednid = iif*(nely+1)+(nely+1-jf);
fixeddof = [2*fixednid(:); 2*fixednid(:)-1];
[iif,jf] = meshgrid(0:nelx,nely);
fixednid = iif*(nely+1)+(nely+1-jf);
fixeddof = [2*fixednid(:); fixeddof];
nele = nelx*nely;
ndof = 2*(nelx+1)*(nely+1);
F = sparse(loaddof,1,force,ndof,1);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);                               
edofVec = 2*nodeids(:)+1;                                                              
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*(nely+1) 2*(nely+1)+1 2*(nely+1)-2 2*(nely+1)-1 -2 -1],nele,1);%%%%b=[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],超神节点提取符写法，阐述其相对单元的1节点x方向自由度的位置
end
function [H,Hs]=preparefilter(rmin,nele,nelx,nely)
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
          for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
               for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                   e2 = (i2-1)*nely+j2;
                   k = k+1;
                   iH(k) = e1;
                   jH(k) = e2;
                   sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
                end
           end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
end
function  generatecommand(F,nele,ndof,thickness,gobalcoords0,E0,nu,xba,edofMat,fixeddof,loaddof,penal,c2rate,Kin,Kout,outdof);    
fid = fopen('D:\command.txt','w');
ndot=ndof/2;
fprintf(fid,'/CWD,''D:\\''\n');
fprintf(fid,'/PREP7\n');
%define node
for i=1:ndot
fprintf(fid,'N,%d,%G,%G,%G\n',i,gobalcoords0(2*i-1),gobalcoords0(2*i),0);
end
%define element type and the thickness
fprintf(fid,'et,1,plane182\nKEYOPT,1,1,0\nKEYOPT,1,3,3\nkeyopt,1,6,0\nR,1,%G, \nMPTEMP,1,0 \n',thickness);
edotMat=0.5*edofMat(:,[2,4,6,8]);
%calculate the element parameters
MU0=E0*xba.^penal/(2*(1+nu));
K0=E0*xba.^penal/(2*(1-2*nu));
d=2./K0;
c1=(1-xba.^penal)*E0*1e-9/6;
c2=(1-xba.^penal)*E0*c2rate;
for i=1:nele
        %define original elements
        fprintf(fid,'TB,HYPE,%d,1,2,NEO\nTBTEMP,0\n',2*i-1);
        fprintf(fid,'TBDATA,,%G,%G,,,,\n',MU0(i),d(i));
        fprintf(fid,'type,1\nmat,%d\nreal,1\nesys,0\n',2*i-1);
        fprintf(fid,'e,%d,%d,%d,%d\n',edotMat(i,1),edotMat(i,2),edotMat(i,3),edotMat(i,4));
        %define the additive hyperelastic elements  
        fprintf(fid,'TB,HYPE,%d,1,2,YEOH\nTBTEMP,0\n',2*i);
        fprintf(fid,'TBDATA,,%G,%G,1e-9,1e6,,\n',c1(i),c2(i));
        fprintf(fid,'type,1\nmat,%d\nreal,1\nesys,0\n',2*i);
        fprintf(fid,'e,%d,%d,%d,%d\n',edotMat(i,1),edotMat(i,2),edotMat(i,3),edotMat(i,4));
end
%%%%%define node of the spring
coords0w=reshape(gobalcoords0,2,[]);
lx=max(coords0w(1,:))-min(coords0w(1,:));
fprintf(fid,'N,%d,%f,%f,%f\n',ndot+1,gobalcoords0(loaddof)+2*lx,gobalcoords0(loaddof+1),0);
fprintf(fid,'d,%d,ux,0\n',ndot+1);
fprintf(fid,'d,%d,uy,0\n',ndot+1);
fprintf(fid,'N,%d,%f,%f,%f\n',ndot+2,gobalcoords0(outdof)+2*lx,gobalcoords0(outdof+1),0);
fprintf(fid,'d,%d,ux,0\n',ndot+2);
fprintf(fid,'d,%d,uy,0\n',ndot+2);
%%%%%define the spring
fprintf(fid,'ET,2,LINK180\nKEYOPT,2,2,0\nR,2,1, ,0\n');
fprintf(fid,'MPDATA,ex,%d,,%.15f\nmpdata,prxy,%d,,%f\ntype,2\nmat,%d\nreal,2\nesys,0\ne,%d,%d\n',2*nele+1,Kin*2*lx,2*nele+1,0.3,2*nele+1,(loaddof+1)/2,ndot+1);
fprintf(fid,'MPDATA,ex,%d,,%.15f\nmpdata,prxy,%d,,%f\ntype,2\nmat,%d\nreal,2\nesys,0\ne,%d,%d\n',2*nele+2,Kout*2*lx,2*nele+2,0.3,2*nele+2,(outdof+1)/2,ndot+2);
%apply the displacement
nfix=size(fixeddof,1);
for i=1:nfix
    if mod(fixeddof(i),2)==1
        fprintf(fid,'d,%d,ux,0\n',(fixeddof(i)+1)/2);
    else       
        fprintf(fid,'d,%d,uy,0\n',fixeddof(i)/2);
    end
end
%apply the external load
nload=size(loaddof,1);
for i=1:nload
    if mod(loaddof(i),2)==1
    fprintf(fid,'F,%d,fx,%G\n',(loaddof(i)+1)/2,full(F(loaddof(i))));
    else
    fprintf(fid,'F,%d,fy,%G\n',loaddof(i)/2,full(F(loaddof(i))));
    end
end
%solve 
fprintf(fid,'finish\n/sol\nANTYPE,0\nNLGEOM,1\nNSUBST,1,0,0\n');
fprintf(fid,'CNVTOL,U,-1, \n'); 
fprintf(fid,'CNVTOL,F,%G,0.001,2, ,  \n',sum(abs(full(F))));
fprintf(fid,'OUTRES,ERASE\nOUTRES,ALL,ALL\n/status,solu\nsolve\n');
%apply the addictional external force
fprintf(fid,'F,%d,fx,%f\n',(outdof+1)/2,full(F(loaddof)/100));
fprintf(fid,'NSUBST,1,0,0\n');
%放松收敛准则
fprintf(fid,'CNVTOL,F,%G,0.001,2, ,  \n',1e4*sum(abs(full(F))));
fprintf(fid,'solve\nfinish\n');
%post processing---write the element strain energy 
fprintf(fid,'/POST1\n');
fprintf(fid,'SET, , ,1, ,1, , \n');
fprintf(fid,'esel,all\n');
fprintf(fid,'*dim,STEN,array,%d,1\n',2*nele);
fprintf(fid,'ETABLE,ea,SENE\n');
fprintf(fid,'*vget,STEN,elem,1,etab,ea\n');
fprintf(fid,'*cfopen,strainenergy,txt\n');
fprintf(fid,'*vwrite,STEN(1,1)\n');
fprintf(fid,'(E13.5)\n');
%post processing---the element logarithmic equivilant strain
fprintf(fid,'*dim,strain1,array,%d,1\n',2*nele);
fprintf(fid,'ETABLE,ab,EPEL,EQV\n');
fprintf(fid,'*vget,strain1,elem,1,etab,ab\n');
fprintf(fid,'*cfopen,eqvstrain,txt\n');
fprintf(fid,'*vwrite,strain1(1,1)\n');
fprintf(fid,'(E13.5)\n');
%post processing---the nodal displacement
fprintf(fid,'nsel,all\n');
fprintf(fid,'*dim,UARRAY,array,%d,2\n',ndot);
fprintf(fid,'*vget,UARRAY(1,1),node,1,u,x\n');
fprintf(fid,'*vget,UARRAY(1,2),node,1,u,y\n');
fprintf(fid,'*cfopen,NODALDISPLACEMENT,txt\n');
fprintf(fid,'*vwrite,UARRAY(1,1),UARRAY(1,2)\n');
fprintf(fid,'(E13.5,5X,E13.5)\n');
%post processing---write the element strain energy 
fprintf(fid,'SET, , ,1, ,2, , \n');
fprintf(fid,'*dim,STEN99,array,%d,1\n',2*nele);
fprintf(fid,'ETABLE,ea99,SENE\n');
fprintf(fid,'*vget,STEN99,elem,1,etab,ea99\n');
fprintf(fid,'*cfopen,strainenergy99,txt\n');
fprintf(fid,'*vwrite,STEN99(1,1)\n');
fprintf(fid,'(E13.5)\n');
fprintf(fid,'finish\n');
fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Qi Chen, Xianmin Zhang and Benliang Zhu      %
% Guangdong Key Laboratory of Precision Equipment and Manufacturing Technology,%
% South China University of Technology                                         %
% Please sent your comments to: chenqiscutedu@qq.com                           %
% The code is intended for educational purposes and theoretical details are    %
% discussed in the paper                                                       %
% A 213-line topology optimization code for geometrically nonlinear structures %
% Qi Chen, Xianmin Zhang and Benliang Zhu,                                     %
% Struct Multidisc Optim, 2018                                                 %
% Disclaimer:                                                                  %
% The authors reserves all rights but do not guaranty that the code is free    %
% from errors. Furthermore, we shall not be liable in any event caused by the  %
% use of the program.                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This is the file mmasub.m
%
function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);
%
%    Written in May 1999 by
%    Krister Svanberg <krille@math.kth.se>
%    Department of Mathematics
%    SE-10044 Stockholm, Sweden.
%
%    Modified ("spdiags" instead of "diag") April 2002
%
%
%    This function mmasub performs one MMA-iteration, aimed at
%    solving the nonlinear programming problem:
%         
%      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
%    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
%                xmin_j <= x_j <= xmax_j,    j = 1,...,n
%                z >= 0,   y_i >= 0,         i = 1,...,m
%*** INPUT:
%
%   m    = The number of general constraints.
%   n    = The number of variables x_j.
%  iter  = Current iteration number ( =1 the first time mmasub is called).
%  xval  = Column vector with the current values of the variables x_j.
%  xmin  = Column vector with the lower bounds for the variables x_j.
%  xmax  = Column vector with the upper bounds for the variables x_j.
%  xold1 = xval, one iteration ago (provided that iter>1).
%  xold2 = xval, two iterations ago (provided that iter>2).
%  f0val = The value of the objective function f_0 at xval.
%  df0dx = Column vector with the derivatives of the objective function
%          f_0 with respect to the variables x_j, calculated at xval.
% df0dx2 = Column vector with the non-mixed second derivatives of the
%          objective function f_0 with respect to the variables x_j,
%          calculated at xval. df0dx2(j) = the second derivative
%          of f_0 with respect to x_j (twice).
%          Important note: If second derivatives are not available,
%          simply let df0dx2 = 0*df0dx.
%  fval  = Column vector with the values of the constraint functions f_i,
%          calculated at xval.
%  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
%          f_i with respect to the variables x_j, calculated at xval.
%          dfdx(i,j) = the derivative of f_i with respect to x_j.
%  dfdx2 = (m x n)-matrix with the non-mixed second derivatives of the
%          constraint functions f_i with respect to the variables x_j,
%          calculated at xval. dfdx2(i,j) = the second derivative
%          of f_i with respect to x_j (twice).
%          Important note: If second derivatives are not available,
%          simply let dfdx2 = 0*dfdx.
%  low   = Column vector with the lower asymptotes from the previous
%          iteration (provided that iter>1).
%  upp   = Column vector with the upper asymptotes from the previous
%          iteration (provided that iter>1).
%  a0    = The constants a_0 in the term a_0*z.
%  a     = Column vector with the constants a_i in the terms a_i*z.
%  c     = Column vector with the constants c_i in the terms c_i*y_i.
%  d     = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
%     
%*** OUTPUT:
%
%  xmma  = Column vector with the optimal values of the variables x_j
%          in the current MMA subproblem.
%  ymma  = Column vector with the optimal values of the variables y_i
%          in the current MMA subproblem.
%  zmma  = Scalar with the optimal value of the variable z
%          in the current MMA subproblem.
%  lam   = Lagrange multipliers for the m general MMA constraints.
%  xsi   = Lagrange multipliers for the n constraints alfa_j - x_j <= 0.
%  eta   = Lagrange multipliers for the n constraints x_j - beta_j <= 0.
%   mu   = Lagrange multipliers for the m constraints -y_i <= 0.
%  zet   = Lagrange multiplier for the single constraint -z <= 0.
%   s    = Slack variables for the m general MMA constraints.
%  low   = Column vector with the lower asymptotes, calculated and used
%          in the current MMA subproblem.
%  upp   = Column vector with the upper asymptotes, calculated and used
%          in the current MMA subproblem.
%
epsimin = sqrt(m+n)*10^(-9);
feps = 0.000001;
asyinit = 0.5;
asyincr = 1.2;
asydecr = 0.7;
albefa = 0.1;
een = ones(n,1);
zeron = zeros(n,1);

% Calculation of the asymptotes low and upp :
if iter < 2.5
  low = xval - asyinit*(xmax-xmin);
  upp = xval + asyinit*(xmax-xmin);
else
  zzz = (xval-xold1).*(xold1-xold2);
  factor = een;
  factor(find(zzz > 0)) = asyincr;
  factor(find(zzz < 0)) = asydecr;
  low = xval - factor.*(xold1 - low);
  upp = xval + factor.*(upp - xold1);
end

% Calculation of the bounds alfa and beta :
zzz = low + albefa*(xval-low);
alfa = max(zzz,xmin);
zzz = upp - albefa*(upp-xval);
beta = min(zzz,xmax);

% Calculations of p0, q0, P, Q and b.

ux1 = upp-xval;
ux2 = ux1.*ux1;
ux3 = ux2.*ux1;
xl1 = xval-low;
xl2 = xl1.*xl1;
xl3 = xl2.*xl1;
ul1 = upp-low;
ulinv1 = een./ul1;
uxinv1 = een./ux1;
xlinv1 = een./xl1;
uxinv3 = een./ux3;
xlinv3 = een./xl3;
diap = (ux3.*xl1)./(2*ul1);
diaq = (ux1.*xl3)./(2*ul1);
p0 = zeron;
p0(find(df0dx > 0)) = df0dx(find(df0dx > 0));
p0 = p0 + 0.001*abs(df0dx) + feps*ulinv1;
p0 = p0.*ux2;
q0 = zeron;
q0(find(df0dx < 0)) = -df0dx(find(df0dx < 0));
q0 = q0 + 0.001*abs(df0dx) + feps*ulinv1;
q0 = q0.*xl2;
dg0dx2 = 2*(p0./ux3 + q0./xl3);
del0 = df0dx2 - dg0dx2;
delpos0 = zeron;
delpos0(find(del0 > 0)) = del0(find(del0 > 0));
p0 = p0 + delpos0.*diap;
q0 = q0 + delpos0.*diaq;
P = zeros(m,n);
P(find(dfdx > 0)) = dfdx(find(dfdx > 0));
%%%P = P * diag(ux2);
P = P * spdiags(ux2,0,n,n);
Q = zeros(m,n);
Q(find(dfdx < 0)) = -dfdx(find(dfdx < 0));
%%%Q = Q * diag(xl2);
Q = Q * spdiags(xl2,0,n,n);
%%%dgdx2 = 2*(P*diag(uxinv3) + Q*diag(xlinv3));
dgdx2 = P*spdiags(uxinv3,0,n,n)+Q*spdiags(xlinv3,0,n,n);
dgdx2 = 2*dgdx2;
del = dfdx2 - dgdx2;
delpos = zeros(m,n);
delpos(find(del > 0)) = del(find(del > 0));
%%%P = P + delpos*diag(diap);
P = P + delpos*spdiags(diap,0,n,n);
%%%Q = Q + delpos*diag(diaq);
Q = Q + delpos*spdiags(diaq,0,n,n);
b = P*uxinv1 + Q*xlinv1 - fval ;

%%% Solving the subproblem by a primal-dual Newton method
[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);
end
%    This is the file subsolv.m
%
function [xmma,ymma,zmma,lamma,xsimma,etamma,mumma,zetmma,smma] = ...
subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);
%
%    Written in May 1999 by
%    Krister Svanberg <krille@math.kth.se>
%    Department of Mathematics
%    SE-10044 Stockholm, Sweden.
%
% This function subsolv solves the MMA subproblem:
%         
% minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z +
%          + SUM[ ci*yi + 0.5*di*(yi)^2 ],
%
% subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
%            alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.
%        
% Input:  m, n, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d.
% Output: xmma,ymma,zmma, slack variables and Lagrange multiplers.
%
een = ones(n,1);
eem = ones(m,1);
epsi = 1;
epsvecn = epsi*een;
epsvecm = epsi*eem;
x = 0.5*(alfa+beta);
y = eem;
z = 1;
lam = eem;
xsi = een./(x-alfa);
xsi = max(xsi,een);
eta = een./(beta-x);
eta = max(eta,een);
mu  = max(eem,0.5*c);
zet = 1;
s = eem;
itera = 0;

while epsi > epsimin
  epsvecn = epsi*een;
  epsvecm = epsi*eem;
  ux1 = upp-x;
  xl1 = x-low;
  ux2 = ux1.*ux1;
  xl2 = xl1.*xl1;
  uxinv1 = een./ux1;
  xlinv1 = een./xl1;

  plam = p0 + P'*lam ;
  qlam = q0 + Q'*lam ;
  gvec = P*uxinv1 + Q*xlinv1;
  dpsidx = plam./ux2 - qlam./xl2 ;

  rex = dpsidx - xsi + eta;
  rey = c + d.*y - mu - lam;
  rez = a0 - zet - a'*lam;
  relam = gvec - a*z - y + s - b;
  rexsi = xsi.*(x-alfa) - epsvecn;
  reeta = eta.*(beta-x) - epsvecn;
  remu = mu.*y - epsvecm;
  rezet = zet*z - epsi;
  res = lam.*s - epsvecm;

  residu1 = [rex' rey' rez]';
  residu2 = [relam' rexsi' reeta' remu' rezet res']';
  residu = [residu1' residu2']';
  residunorm = sqrt(residu'*residu);
  residumax = max(abs(residu));

  ittt = 0;
  while residumax > 0.9*epsi & ittt < 100
    ittt=ittt + 1;
    itera=itera + 1;

    ux1 = upp-x;
    xl1 = x-low;
    ux2 = ux1.*ux1;
    xl2 = xl1.*xl1;
    ux3 = ux1.*ux2;
    xl3 = xl1.*xl2;
    uxinv1 = een./ux1;
    xlinv1 = een./xl1;
    uxinv2 = een./ux2;
    xlinv2 = een./xl2;
    plam = p0 + P'*lam ;
    qlam = q0 + Q'*lam ;
    gvec = P*uxinv1 + Q*xlinv1;
    GG = P*spdiags(uxinv2,0,n,n) - Q*spdiags(xlinv2,0,n,n);
    dpsidx = plam./ux2 - qlam./xl2 ;
    delx = dpsidx - epsvecn./(x-alfa) + epsvecn./(beta-x);
    dely = c + d.*y - lam - epsvecm./y;
    delz = a0 - a'*lam - epsi/z;
    dellam = gvec - a*z - y - b + epsvecm./lam;
    diagx = plam./ux3 + qlam./xl3;
    diagx = 2*diagx + xsi./(x-alfa) + eta./(beta-x);
    diagxinv = een./diagx;
    diagy = d + mu./y;
    diagyinv = eem./diagy;
    diaglam = s./lam;
    diaglamyi = diaglam+diagyinv;

    if m < n
      blam = dellam + dely./diagy - GG*(delx./diagx);
      bb = [blam' delz]';
      Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
      AA = [Alam     a
            a'    -zet/z ];
      solut = AA\bb;
      dlam = solut(1:m);
      dz = solut(m+1);
      dx = -delx./diagx - (GG'*dlam)./diagx;
    else
      diaglamyiinv = eem./diaglamyi;
      dellamyi = dellam + dely./diagy;
      Axx = spdiags(diagx,0,n,n) + GG'*spdiags(diaglamyiinv,0,m,m)*GG;
      azz = zet/z + a'*(a./diaglamyi);
      axz = -GG'*(a./diaglamyi);
      bx = delx + GG'*(dellamyi./diaglamyi);
      bz  = delz - a'*(dellamyi./diaglamyi);
      AA = [Axx   axz
            axz'  azz ];
      bb = [-bx' -bz]';
      solut = AA\bb;
      dx  = solut(1:n);
      dz = solut(n+1);
      dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
    end

    dy = -dely./diagy + dlam./diagy;
    dxsi = -xsi + epsvecn./(x-alfa) - (xsi.*dx)./(x-alfa);
    deta = -eta + epsvecn./(beta-x) + (eta.*dx)./(beta-x);
    dmu  = -mu + epsvecm./y - (mu.*dy)./y;
    dzet = -zet + epsi/z - zet*dz/z;
    ds   = -s + epsvecm./lam - (s.*dlam)./lam;
    xx  = [ y'  z  lam'  xsi'  eta'  mu'  zet  s']';
    dxx = [dy' dz dlam' dxsi' deta' dmu' dzet ds']';
    
    stepxx = -1.01*dxx./xx;
    stmxx  = max(stepxx);
    stepalfa = -1.01*dx./(x-alfa);
    stmalfa = max(stepalfa);
    stepbeta = 1.01*dx./(beta-x);
    stmbeta = max(stepbeta);
    stmalbe  = max(stmalfa,stmbeta);
    stmalbexx = max(stmalbe,stmxx);
    stminv = max(stmalbexx,1);
    steg = 1/stminv;

    xold   =   x;
    yold   =   y;
    zold   =   z;
    lamold =  lam;
    xsiold =  xsi;
    etaold =  eta;
    muold  =  mu;
    zetold =  zet;
    sold   =   s;

    itto = 0;
    resinew = 2*residunorm;
    while resinew > residunorm & itto < 50
    itto = itto+1;

    x   =   xold + steg*dx;
    y   =   yold + steg*dy;
    z   =   zold + steg*dz;
    lam = lamold + steg*dlam;
    xsi = xsiold + steg*dxsi;
    eta = etaold + steg*deta;
    mu  = muold  + steg*dmu;
    zet = zetold + steg*dzet;
    s   =   sold + steg*ds;
    ux1 = upp-x;
    xl1 = x-low;
    ux2 = ux1.*ux1;
    xl2 = xl1.*xl1;
    uxinv1 = een./ux1;
    xlinv1 = een./xl1;
    plam = p0 + P'*lam ;
    qlam = q0 + Q'*lam ;
    gvec = P*uxinv1 + Q*xlinv1;
    dpsidx = plam./ux2 - qlam./xl2 ;

    rex = dpsidx - xsi + eta;
    rey = c + d.*y - mu - lam;
    rez = a0 - zet - a'*lam;
    relam = gvec - a*z - y + s - b;
    rexsi = xsi.*(x-alfa) - epsvecn;
    reeta = eta.*(beta-x) - epsvecn;
    remu = mu.*y - epsvecm;
    rezet = zet*z - epsi;
    res = lam.*s - epsvecm;

    residu1 = [rex' rey' rez]';
    residu2 = [relam' rexsi' reeta' remu' rezet res']';
    residu = [residu1' residu2']';
    resinew = sqrt(residu'*residu);
    steg = steg/2;
    end
  residunorm=resinew;
  residumax = max(abs(residu));
  steg = 2*steg;
  end
epsi = 0.1*epsi;
end

xmma   =   x;
ymma   =   y;
zmma   =   z;
lamma =  lam;
xsimma =  xsi;
etamma =  eta;
mumma  =  mu;
zetmma =  zet;
smma   =   s;
end


