%%%%%%%%%% SERA TOPOLOGY OPTIMIZATION CODE %%%%%%%%%%
function sera(nelx,nely,volfrac,rmin)
% INITIALIZE
i=1;Change=1;
x(1:nely,1:nelx)=1.0;
VF(i)=ceil(sum(sum(x)))/(nelx*nely);
PR = 0.03;SR=1.3;B=0.003;
% START ITERATION
 while Change>0.0001
 i=i+1;
 VF(i)=max(VF(i-1)*(1-PR),volfrac);
 [U,Ko,K]=FE(nelx,nely,x);
 c(i)=U'*K*U;
 for ely=1:nely
 for elx=1:nelx
 n1=(nely+1)*(elx-1)+ely;
 n2=(nely+1)*elx+ely;
 Ue=U([2*n1-1;2*n1;2*n2-1;2*n2;2*n2+1;2*n2+2;2*n1+1;2*n1+2],1);
 alfa(ely,elx)=Ue'*Ko*Ue;
 end
end
 % FILTERING TECHNIQUE
 [alfa]=Filter(nelx,nely,rmin,x,alfa);
 % DESIGN UPDATE
 [x,NumElem_Add,NumElem_Rem,PR]=SERA_Update(nelx,nely,alfa,x,VF,
 volfrac,i,SR,B,PR);
 VF(i)=ceil(sum(sum(x)))/(nelx*nely);
% STOPPING CRITERION
if i>20
 Change=abs((sum(c(i-19:i-10))-sum(c(i-9:i)))/sum(c(i-9:i)));
 end
 % PRINT RESULTS
 disp([' Iteration: ' sprintf('%4i', (i-1)) ' Volume fraction: '
 sprintf('%6.4f', VF(i)) ' Compliance: ' sprintf('%6.6f',c(i)) '
Change: ' sprintf('%7.5f', Change) ' Removed: ' sprintf('%5i',
 NumElem_Rem) ' Added: ' sprintf('%5i', NumElem_Add)])
 % PLOT DENSITIES
 colormap(gray); imagesc(-x); axis equal; axis tight;
 axis off; pause(1e-6)
 end
 %%%%%%%%%% FE-ANALYSIS %%%%%%%%%%
 function [U,Ko,K]=FE(nelx,nely,x)
 [Ko] = lk;
 I=zeros(nelx*nely*64,1);
 J=zeros(nelx*nely*64,1);
 X=zeros(nelx*nely*64,1);
 F=sparse(2*(nely+1)*(nelx+1),1);
 U=zeros(2*(nely+1)*(nelx+1),1);
 ntriplets=0;
 for elx = 1:nelx
 for ely = 1:nely
 n1 = (nely+1)*(elx-1)+ely;
 n2 = (nely+1)*e1x + e1y;
 edof = [2*n1-1 2*n1 2*n2-1 2*n2 2*n2+1 2*n2+2 2*n1+1 2*n1+2];
 for krow = 1:8
 for kcol = 1:8
 ntriplets = ntriplets+1;
 I(ntriplets) = edof(krow);
 J(ntriplets) = edof(kcol);
 X(ntriplets) = x(ely,elx)*Ko(krow,kcol);
 end
 end
 end
 end
 K=sparse(I,J,X,2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
 % DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
 F(2,1) = -1;
 fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
 alldofs = [1:2*(nely+1)*(nelx+1)];
 freedofs = setdiff(alldofs,fixeddofs);
 % SOLVING
 U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
 U(fixeddofs,:) = 0;
 %%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%
 function [Ko] = lk
 E = 1; nu = 0.3;
 k=[1/2-nu/6;1/8+nu/8;-1/4-nu/12;-1/8+3*nu/8;
 -1/4+nu/12;-1/8-nu/8;nu/6;1/8-3*nu/8];
 Ko = E/(1-nu^2)*
 [k(1),k(2),k(3),k(4),k(5),k(6),k(7),k(8);
 k(2),k(1),k(8),k(7),k(6),k(5),k(4),k(3);
 k(3),k(8),k(1),k(6),k(7),k(4),k(5),k(2);
 k(4),k(7),k(6),k(1),k(8),k(3),k(2),k(5);
 k(5),k(6),k(7),k(8),k(1),k(2),k(3),k(4);
 k(6),k(5),k(4),k(3),k(2),k(1),k(8),k(7);
 k(7),k(4),k(5),k(2),k(3),k(8),k(1),k(6);
 k(8),k(3),k(2),k(5),k(4),k(7),k(6),k(1)];
 %%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%
 function [alfaNew]=Filter(nelx,nely,rmin,x,alfa)
 alfaNew=zeros(nely,nelx);
 for i = 1:nelx
 for j= 1:nely
 sum=0.0;
 for k = max(i-floor(rmin),1) : min(i+floor(rmin),nelx)
 for l = max(j-floor(rmin),1) : min(j+floor(rmin),nely)
 sum = sum+max(0,rmin-sqrt((i-k)^2+(j-l)^2));
 alfaNew(j,i) = alfaNew(j,i) +
 max(0,rmin-sqrt((i-k)^2+(j-l)^2))*x(l,k)*alfa(l,k);
 end
 end
 alfaNew(j,i) = alfaNew(j,i)/x(j,i)/sum;
 end
 end
 %%%%%%%%%% SERA UPDATE %%%%%%%%%%
 function [x,NumElem_Add,NumElem_Rem,PR]=SERA_Update(nelx,nely,
 alfa,x,VF,volfrac,i,SR,B,PR)
 alfa_min=min(min(alfa));alfa_max=max(max(alfa));
 alfa_V=alfa;alfa_R=alfa;
 alfa_R(x<1.0)=alfa_max; alfa_V(x>1e-9)=alfa_min;
 NumElem_Add=0;NumElem_Rem=0;
if VF(i)>volfrac
    DeltaV(i)=abs(VF(i)-VF(i-1));
 DeltaV_Rem=DeltaV(i)*(SR);
NumElem_Rem=max(1,floor(nelx*nely*DeltaV_Rem));
[x,NumElem_Rem]=Update_R(nelx,nely,x,alfa_R,NumElem_Rem);
NumElem_Add=max(1,floor(NumElem_Rem*(SR-1)));
[x,NumElem_Add]=Update_V(nelx,nely,x,alfa_V,NumElem_Add);
end
 else
if (VF(i-1)-volfrac)>=0.001 PR=PR*0.7;
 else
 DeltaV_Rem=B*volfrac;
 NumElem_Rem=max(1,floor(nelx*nely*DeltaV_Rem));
 [x,NumElem_Rem]=Update_R(nelx,nely,x,alfa_R,NumElem_Rem);
 NumElem_Add=NumElem_Rem;
 [x,NumElem_Add]=Update_V(nelx,nely,x,alfa_V,NumElem_Add);
 end
 end
 %%%%%%%%%% UPDATE_R %%%%%%%%%%
 function[x,NumElem_Rem]=Update_R(nelx,nely,x,alfa_R,NumElem_Rem)
 alfa_R_vec=sort(reshape(alfa_R,(nelx*nely),1),'descend');
 alfa_R_th=alfa_R_vec((nelx*nely)-NumElem_Rem,1);
 Elem_Rem=((alfa_R-alfa_R_th)/abs(alfa_R_th))<1e-6;
 NumElem_Rem=sum(sum(Elem_Rem));
 x(Elem_Rem)=1e-9;
 %%%%%%%%%% UPDATE_V %%%%%%%%%%
 function [x,NumElem_Add]=Update_V(nelx,nely,x,alfa_V,NumElem_Add)
 alfa_V_vec=sort(reshape(alfa_V,(nelx*nely),1), 'descend');
alfa_V_th=alfa_V_vec(NumElem_Add,1);
 Elem_Add=((alfa_V-alfa_V_th)/abs(alfa_V_th))>-1e-6;
 NumElem_Add=sum(sum(Elem_Add));
 x(Elem_Add)=1.0;
    