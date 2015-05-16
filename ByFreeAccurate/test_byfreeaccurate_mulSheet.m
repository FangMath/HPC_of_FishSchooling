function [pot2] = test_byfreeaccurate_mulSheet(targetb, sourcef, gammaf, L, Nw)
    close all;
%% test ByFreeAccurate_mulSheet.cpp routine with matlab functions
%% Copyright Fang Fang, 05/09/2015, Courant Inst.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mex -I"../Fish/" -I"../eigen/" -v -largeArrayDims ByFreeAccurate_mulSheet.cpp ByFreeAccurate_mulSheet_routine.cpp -lmwblas CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -fopenmp" 

%%% generate random data for tests %%%%%%%%%
nsource = 1000;
ntarget = 900;
Nw = 1;
%% generate targets
target(1,:) = rand(1,ntarget)-.5; target(2,:) = rand(1,ntarget);
%% generate sources 
source_r = (1:nsource)/nsource-.5; 
msource_r = source_r;
msource_i = source_r.^2;
mgammaf = sin(source_r);
L = length(msource_r); 
for ib = 2:Nw 
msource_r = [msource_r, source_r];
msource_i = [msource_i, source_r.^2+(ib-1)/Nw];
mgammaf = [mgammaf, sin(source_r)+(ib-1)/Nw];
L(ib) = length(msource_r); 
end

plot(target(1,:), target(2,:), 'ro'); hold on;
plot(msource_r, msource_i, 'b.'); hold on;
source = [msource_r; msource_i];
gammaf = mgammaf;
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
itN = 10;
for it = 1:itN
[C1,C2] = ByFreeAccurate_mulSheet(target(1,:)',target(2,:)',source(1,:)', source(2,:)',gammaf.', L');
pot2 = C1+C2*1i;
pot2 = pot2/(2*1i*pi);
end
tcc = toc(tc)/itN
N=nsource;

%% time matlab calculations %%
tm = tic;
for it = 1:itN
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
end
tmm = toc(tm)/itN

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

