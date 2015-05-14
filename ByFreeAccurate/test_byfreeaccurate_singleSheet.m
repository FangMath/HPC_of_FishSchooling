function [pot2] = test_byfreeaccurate_singleSheet(targetb, zetaf, gammaf)
    close all;
%% test ByFreeAccurate_singleSheet.c routine with matlab functions
%% Copyright Fang Fang, 05/09/2015, Courant Inst.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mex -I"../Fish/" -I"../eigen/" -v -largeArrayDims ByFreeAccurate_singleSheet.cpp ByFreeAccurate_singleSheet_routine.cpp -lmwblas CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -fopenmp" 

%% old compile
%mex -I"../Fish/" -I"../eigen/" -v -largeArrayDims ByFreeAccurate_singleSheet.cpp ByFreeAccurate_singleSheet_routine.cpp -lmwblas CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -fopenmp" 


%%%% generate random data for tests %%%%%%%%%
%nsource = 1000;
%ntarget = 90;
%%% generate targets
%target(1,:) = rand(1,ntarget)-.5; target(2,:) = rand(1,ntarget);
%%% generate sources 
%source(1,:) = (1:nsource)/nsource-.5; 
%source(2,:) = source(1,:).^2;
%%% generate L 
%%% generate distribution 
%gammaf = sin(source(1,:));
%plot(target(1,:), target(2,:), 'ro'); hold on;
%plot(source(1,:), source(2,:), 'b.'); hold on;
%
%% use vortex sheets data  %%%%%%%%%
nsource = length(zetaf);
source(1,:)=real(zetaf);
source(2,:)=imag(zetaf);

ntarget = length(targetb);
target = zeros(2,ntarget);

target(1,:)=real(targetb);
target(2,:)=imag(targetb);

L = ones(1,nsource);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tests starting
tc = tic;
[C1,C2] = ByFreeAccurate_singleSheet(target(1,:)',target(2,:)',source(1,:)', source(2,:)',gammaf.', L');
pot2 = C1+C2*1i;
pot2 = pot2/(2*1i*pi);
tcc = toc(tc)
N=nsource;

%% time matlab calculations %%
tm = tic;
for it = 1:8
if N==1
    X=0;
else
    X(1:N-1) = (gammaf(2:end).' - gammaf(1:end-1).')./(source(1,2:end)'-source(1,1:end-1)'+source(2,2:end)'*1i-source(2,1:end-1)'*1i);
end
X=-X/(2*1i*pi);
distr=X;

N=nsource;
if N > 1
    dx = bsxfun(@minus,target(1,1:ntarget).',source(1,1:nsource));
    dy = bsxfun(@minus,target(2,1:ntarget).',source(2,1:nsource));
    dz = dx+dy*1i;
    F = log(dz(:,1:N-1)./dz(:,2:N));
    pot = F*distr.';
    pot = pot.';
else
    pot = zeros(1,nsource);
end

end
tmm = toc(tm)

speed_ratio = tcc/tmm
%% check correctness %%
dif = (pot2-pot.');
errrr = norm(dif)/sqrt(nsource)

%% speed up ratio %%

pot2 = pot2.';
end

