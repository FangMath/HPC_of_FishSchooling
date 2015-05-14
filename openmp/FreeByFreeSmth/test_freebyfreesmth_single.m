function [pot1, pot2] = test_freebyfreesmth_single(tr, ti, sr, si, delta, distr)
%% test FreeByFreeSmth.c routine with matlab functions
%% Copyright Fang Fang, 04/27/2015, Courant Inst.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mex -I"../Fish/" -I"../eigen/" -v -largeArrayDims ByFreeAccurate_mulSheet.cpp ByFreeAccurate_mulSheet_routine.cpp -lmwblas CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -fopenmp" 
mex -I"../Fish/" -I"../eigen/" -v -largeArrayDims FreeByFreeSmth_single.cpp FreeByFreeSmth_single_routine.cpp -lmwblas CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -fopenmp" 


%%% generate random data for tests %%%%%%%%%
Ns = 1000;
Nt = 900;
%% generate targets
tr = rand(1,Nt); ti = rand(1,Nt);
%% generate sources 
sr = rand(1,Ns); si = rand(1,Ns);
%% generate delta 
delta = rand(1,Ns);
%% generate distribution 
distr = rand(Ns,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tests starting
tc = tic;
[C1,C2] = FreeByFreeSmth_single(tr',ti',sr', si', delta',distr);
pot2 = C2 - C1*1i;
tcc = toc(tc)


%% time matlab calculations %%
tm = tic;
dx = bsxfun(@minus,tr',sr);
dy = bsxfun(@minus,ti',si);
dz = dx.^2 + dy.^2;
ddelta = repmat(delta,Nt,1);
dz = bsxfun(@plus,dz,ddelta.^2);
F = -1/(Ns*2*pi)*(dx-dy*1i)./dz;
pot1 = F*distr/1i;
tmm = toc(tm)


%% check correctness %%
errrr = norm(pot2-pot1)

%% speed up ratio %%
speed_ratio = tcc/tmm

end
