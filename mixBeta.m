function [] = mixBeta(ID)
% Bayesian mixed beta regression (bounded response)
% use flat prior for beta and Jefferys prior for sigma2
global Y X p n np nonrandom delta2

tot = 1e5; burn = 95e3;
objrate = 0.44; delta2 = 1e4;

ch = str2double(num2str(ID));
ch0 = ch; nChain = 3; nvar = ceil(ch/nChain); ch = ch - (nvar-1)*nChain;
fprintf('nvar = %d, chain = %d:\n', nvar, ch)

rng('default'); rng(ch0*241);

verbose = 0; nonrandom = 0;
load('Ford_NA.mat')
Y = Ytrans(:,nvar); %#ok<*MSNU> <NODEF>
X = [ones(size(X,1),1),X];
p = size(X,2); n = length(unique(np));
ntheta = 2*p+n;
if nonrandom == 1
    ntheta = 2*p;
end
matpara = zeros(tot-burn, ntheta+1);
meantau2 = 0.01; vartau2 = 10^4;
alphatau = 2+meantau2^2/vartau2;
invbetatau = meantau2*(alphatau-1);

accepts = zeros(1,ntheta);
batchLen = 50; batchNum = 0;
batchTot = tot/batchLen;
rates = zeros(batchTot, ntheta);
tunings = -2*ones(1,ntheta);

tau2 = 1e-4;
thetas = [coefs(:,nvar)', normrnd(0, sqrt(tau2), [1,n])]; %#ok<*NODEF>
thetas = thetas.*unifrnd(.9,1.1,[1,length(thetas)]);

tic
for b = 1:tot
    if verbose == 1
        fprintf('%6d', b)
        if(~mod(b,20))
            fprintf('\n')
        end
    end
    for k = 1:ntheta
        re = (k> (2*p));
        loglike = getlik(thetas, re, tau2);
        thetas1 = thetas; thetas1(k) = normrnd(thetas(k), exp(tunings(k)));
        loglike1 = getlik(thetas1, re, tau2);
        MH = exp(loglike1 - loglike);
        u1 = rand(1);
        if u1 <= MH
            thetas = thetas1;
            accepts(k) = accepts(k)+1;
        end
    end
    if nonrandom ~= 1
        tau2 = 1./gamrnd(alphatau+n/2, (invbetatau+sum(thetas((1:n)+2*p).^2)/2)^-1);
    end
    if b > burn
        matpara(b-burn,:) = [thetas, tau2];
    end
    if ~mod(b,batchLen)
        batchNum = batchNum+1;
        accepts = accepts./batchLen;
        rates(batchNum,:) = accepts;
        tunings = tunings + sign((accepts>objrate)-0.5).*min(0.01,1/sqrt(batchNum));
        accepts = zeros(1,ntheta);
    end
end
runtime = toc/3600;
disp(num2str(rates(batchTot-4+(1:4),:)))
fprintf('\n%d iterations are done with elapsed time %.2f hours.\n', tot, runtime)
save(strcat('outBetaCase',num2str(nvar),'_Ch',num2str(ch),'.mat'),'runtime','matpara','rates')
end

function [lik] = getlik(thetas, re, tau2) % likelihood for beta regression
global X p n np Y delta2
us0 = thetas((1:n)+2*p)'; us = us0(np);
mu = F(X*thetas(1:p)' + us); phi = exp(X*thetas(p+(1:p))');
muphi = mu.*phi;
lik = sum(gammaln(phi)) - sum(gammaln(muphi)) - sum(gammaln(phi-muphi))...
    + sum((muphi-1).*log(Y)) + sum((phi-muphi-1).*log(1-Y));
if re == 1
    lik = lik - 0.5*sum(us0.^2)/tau2;
elseif re == 0
    lik = lik - 0.5*sum(thetas(1:(2*p)).^2)/delta2;
end
end

function [y] = F(x)
y = exp(x)./(1+exp(x));
end


