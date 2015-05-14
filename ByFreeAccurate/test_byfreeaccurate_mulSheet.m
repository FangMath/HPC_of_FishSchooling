function [pot2] = test_byfreeaccurate_mulSheet(targetb, sourcef, gammaf, L, Nw)
    close all;
%% test ByFreeAccurate_singleSheet.c routine with matlab functions
%% Copyright Fang Fang, 05/09/2015, Courant Inst.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mex -I"../Fish/" -I"../eigen/" -v -largeArrayDims ByFreeAccurate_mulSheet.cpp ByFreeAccurate_mulSheet_routine.cpp -lmwblas CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -fopenmp" 

%% old compile
%mex -I"../Fish/" -I"../eigen/" -v -largeArrayDims ByFreeAccurate_singleSheet.cpp ByFreeAccurate_singleSheet_routine.cpp -lmwblas CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -fopenmp" 


%%% generate random data for tests %%%%%%%%%
nsource = 1000;
ntarget = 900;
%% generate targets
target(1,:) = rand(1,ntarget)-.5; target(2,:) = rand(1,ntarget);
%% generate sources 
source1(1,:) = (1:nsource)/nsource-.5; 
source1(2,:) = source1(1,:).^2;
gammaf1 = sin(source1(1,:));
source2(1,:) = (1:nsource)/nsource - .5; 
source2(2,:) = source2(1,:).^2+.5;
gammaf2 = sin(source2(1,:));

source = [source1, source2];
gammaf = [gammaf1, gammaf2];
Nw = 2;
L = size(source1,2);
for ib = 2:Nw
L(ib) = L(ib-1)+size(source2,2); 
end
plot(target(1,:), target(2,:), 'ro'); hold on;
plot(source1(1,:), source1(2,:), 'b.'); hold on;
plot(source2(1,:), source2(2,:), 'b.'); hold on;
%
%%%%% use vortex sheets data  %%%%%%%%%
%nsource = length(sourcef);
%source(1,:)=real(sourcef);
%source(2,:)=imag(sourcef);
%
%ntarget = length(targetb);
%target = zeros(2,ntarget);
%
%target(1,:)=real(targetb);
%target(2,:)=imag(targetb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tests starting
tc = tic;
[C1,C2] = ByFreeAccurate_mulSheet(target(1,:)',target(2,:)',source(1,:)', source(2,:)',gammaf.', L');
pot2 = C1+C2*1i;
pot2 = pot2/(2*1i*pi);
tcc = toc(tc)
N=nsource;

%% time matlab calculations %%
tm = tic;
W(1).zetaf_r = source(1,1:L(1));
W(1).zetaf_i = source(2,1:L(1));
W(1).gammaf = gammaf(1:L(1));
for ib = 2:Nw
W(ib).zetaf_r = source(1,L(ib-1)+1:L(ib));
W(ib).zetaf_i = source(2,L(ib-1)+1:L(ib));
W(ib).gammaf = gammaf(L(ib-1)+1:L(ib));
end
pot = zeros(1,ntarget);
for ib = 1:Nw
pot=pot+bodyByFree(target, [W(ib).zetaf_r;W(ib).zetaf_i], W(ib).gammaf, L);
end
tmm = toc(tm)

%% speed up ratio %%
speed_ratio = tcc/tmm
%% check correctness %%
dif = (pot2-pot.');
errrr = norm(dif)/sqrt(nsource)







function pot=bodyByFree(target, source, gammaf, L)
N=size(source,2);
if N==1
    X=0;
else
    X(1:N-1) = (gammaf(2:end).' - gammaf(1:end-1).')./(source(1,2:end)'-source(1,1:end-1)'+source(2,2:end)'*1i-source(2,1:end-1)'*1i);
end
X=-X/(2*1i*pi);
distr=X;

if N > 1
    dx = bsxfun(@minus,target(1,1:ntarget).',source(1,1:N));
    dy = bsxfun(@minus,target(2,1:ntarget).',source(2,1:N));
    dz = dx+dy*1i;
    F = log(dz(:,1:N-1)./dz(:,2:N));
    pot = F*distr.';
    pot = pot.';
else
    pot = zeros(1,N);
end

end % end of function bodyByFree
pot2 = pot2.';
end % end of test file 

