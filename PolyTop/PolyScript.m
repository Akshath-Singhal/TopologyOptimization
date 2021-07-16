%------------------------------ PolyScript -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyTop: A Matlab %
% implementation of a general topology optimization framework using       %
% unstructured polygonal finite element meshes", Struct Multidisc Optim,  %
% DOI 10.1007/s00158-011-0696-x                                           %
%-------------------------------------------------------------------------%

%% ---------------------------------------------------- CREATE 'fem' STRUCT
[Node,Element,Supp,Load] = PolyMesher(@MbbDomain,5000,30);
fem = struct(...
  'NNode',size(Node,1),...     % Number of nodes
  'NElem',size(Element,1),...  % Number of elements
  'Node',Node,...              % [NNode x 2] array of nodes
  'Element',{Element},...      % [NElement x Var] cell array of elements
  'Supp',Supp,...              % Array of supports
  'Load',Load,...              % Array of loads
  'Nu0',0.3,...                % Poisson's ratio of solid material
  'E0',1.0,...                 % Young's modulus of solid material
  'Reg',0 ...                  % Tag for regular meshes
   );                    
%% ---------------------------------------------------- CREATE 'opt' STRUCT
R = 0.04;
VolFrac = 0.5;
m = @(y)MatIntFnc(y,'SIMP',3);
P = PolyFilter(fem,R);      
zIni = VolFrac*ones(size(P,2),1);
opt = struct(...               
  'zMin',0.0,...               % Lower bound for design variables
  'zMax',1.0,...               % Upper bound for design variables
  'zIni',zIni,...              % Initial design variables
  'MatIntFnc',m,...            % Handle to material interpolation fnc.
  'P',P,...                    % Matrix that maps design to element vars.
  'VolFrac',VolFrac,...        % Specified volume fraction cosntraint
  'Tol',0.01,...               % Convergence tolerance on design vars.
  'MaxIter',150,...            % Max. number of optimization iterations
  'OCMove',0.2,...             % Allowable move step in OC update scheme
  'OCEta',0.5 ...              % Exponent used in OC update scheme
   );              
%% ---------------------------------------------------------- RUN 'PolyTop'
figure;
for penal = 1:0.5:4        %Continuation on the penalty parameter
   disp(['current p: ', num2str(penal)]);
   opt.MatIntFnc = @(y)MatIntFnc(y,'SIMP',penal);
   [opt.zIni,V,fem] = PolyTop(fem,opt);
end
%% ------------------------------------------------------------------------