function [pot1, pot2] = test_freebyfreesmth_mul(tr, ti, sr, si, delta, distr)
%% test FreeByFreeSmth.c routine with matlab functions
%% Copyright Fang Fang, 04/27/2015, Courant Inst.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mex -I"../Fish/" -I"../eigen/" -v -largeArrayDims FreeByFreeSmth_mul.cpp FreeByFreeSmth_mul_routine.cpp -lmwblas -lrt CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" CXXOPTIMFLAGS="\$CXXOPTIMFLAGS -fopenmp" 


%%% generate random data for tests %%%%%%%%%
Ns = 1000;
Nt = 900;
Nw = 8;
Lt = zeros(Nw,1); Ls = zeros(Nw,1);
%% generate targets
tr = rand(Nt,1); ti = rand(Nt,1);
%% generate sources 
sr = rand(Ns,1); si = rand(Ns,1);
%% generate delta 
delta = rand(Ns,1);
%% generate distribution 
distr = rand(Ns,1);
Lt(1) = Nt; Ls(1) = Ns;

mtr = tr; mti = ti; msr = sr; msi = si;
mdelta = delta; mdistr = distr;
for ib = 2:Nw
    mtr = [mtr; rand(Nt,1)];
    mti = [mti; rand(Nt,1)];
    msr = [msr; rand(Ns,1)];
    msi = [msi; rand(Ns,1)];
    mdelta = [mdelta; rand(Ns,1)];
    mdistr = [mdistr; rand(Ns,1)];
    Lt(ib) = Lt(ib-1) + Nt;
    Ls(ib) = Ls(ib-1) + Ns;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tests starting
tc = tic;
[C1,C2] = FreeByFreeSmth_mul(mtr,mti, Lt, msr, msi, mdelta, mdistr, Ls);
pot2 = C1 + C2*1i;
pot2 = pot2/1i;
tcc = toc(tc)


%% time matlab calculations %%
tm = tic;
for ib = 1:Nw
if ib == 1
    tr = mtr(1:Lt(1));
    ti = mti(1:Lt(1));
    sr = msr(1:Ls(1));
    si = msi(1:Ls(1));
    delta = mdelta(1:Ls(1));
    distr = mdistr(1:Ls(1));
else
    tr = mtr(Lt(ib-1)+1:Lt(ib));
    ti = mti(Lt(ib-1)+1:Lt(ib));
    sr = msr(Ls(ib-1)+1:Ls(ib));
    si = msi(Ls(ib-1)+1:Ls(ib));
    delta = mdelta(Ls(ib-1)+1:Ls(ib));
    distr = mdistr(Ls(ib-1)+1:Ls(ib));
end
pot1{ib} = smth_matlab(tr, ti, sr, si, delta, distr);
    end
tmm = toc(tm)


%% check correctness %%
for ib = 1:Nw
    if ib == 1
errrr(ib) = norm(pot2(1:Lt(1))-pot1{ib});
else
errrr(ib) = norm(pot2(Lt(ib-1)+1:Lt(ib))-pot1{ib});
    end
end
errrr

%% speed up ratio %%
speed_ratio = tcc/tmm

function pot1 = smth_matlab(tr, ti, sr, si, delta, distr)
Nt = length(tr);
Ns = length(sr);
dx = bsxfun(@minus,tr,sr');
dy = bsxfun(@minus,ti,si');
dz = dx.^2 + dy.^2;
ddelta = repmat(delta',Nt,1);
dz = bsxfun(@plus,dz,ddelta.^2);
F = -1/(Ns*2*pi)*(dx-dy*1i)./dz;
pot1 = F*distr/1i;
end
end
