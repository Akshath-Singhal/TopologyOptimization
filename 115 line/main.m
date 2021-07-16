
function main
 [nx,ny,tol_out,tol_f,iter_max_in,iter_max_out,p,q,e,v,rf] = set_parameters ();
 multitop(nx,ny,tol_out,tol_f,iter_max_in,iter_max_out,p,q,e,v,rf);
end

function [nx,ny,tol,tolf,im_in,im_out,p,q,e,v,rf] = set_parameters ()
  nx = 96; ny = 48; tol = 0.001; tolf = 0.05; im_in = 2; im_out = 200;
  p = 4; q = 3; e = [9 3 1 1e-9]'; v = [0.16 0.08 0.08 0.68]'; rf = 8;
end


 








 
 