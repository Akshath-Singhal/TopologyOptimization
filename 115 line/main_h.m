function main_h
 [nx,ny,tol_out,tol_f,iter_max_in,iter_max_out,p,q,e,v,rf] = set_parameters ();
 multitop_h(nx,ny,tol_out,tol_f,iter_max_in,iter_max_out,p,q,e,v,rf);
end

function [nx,ny,tol,tolf,im_in,im_out,p,q,e,v,rf] = set_parameters ()
 nx = 100; ny = 50; tol = 0.001; tolf = 0.05; im_in = 2; im_out = 400;
 p = 5; q = 3; e = [1000 400 200 100 1]'; v = [0.2 0.1 0.1 0.1 0.5]'; rf = 8;
end


