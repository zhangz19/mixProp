function [] = mixProportionalOdds(ID)
% Bayesian Mixed Proportional odds model (ordinal response)
% using Metropolis Hastings algorithm with adaptive MCMC
global Y X C p n np incs nonrandom delta2

simu = 0;
usedat = 1;
reportrate = 0; verbose = 1; nonrandom = 0;
tot = 55e3; burn = 5e4;
delta2 = 1e4; objrate = 0.44;

rng('default')
ch = str2double(num2str(ID));
ch0 = ch; nChain = 5; nvar = ceil(ch/nChain); ch = ch - (nvar-1)*nChain;
fprintf('nvar = %d, chain = %d:\n', nvar, ch)

if usedat == 1; load('Ford_NA.mat'); end
Y = Y(:,nvar);
p = size(X,2); n = length(unique(np)); C = max(Y);
ntheta = C-1+p+n;
tunings = [-2*ones(1,C-1), -3*ones(1,p), -1*ones(1,n)];
if nonrandom == 1; ntheta = C-1+p; tunings = tunings(1:ntheta); end
matpara = zeros(tot-burn, ntheta+1);
alphatau = 0; invbetatau = 0;
if nvar == 1
    thetas = initvec';
    tau2 = inittau2;
else
    thetas = initvec2';
    tau2 = inittau22;
end

if simu == 1 %initial values become true
    alphas = thetas(1:(C-1)); betas = thetas((1:p)+C-1)';
    N = size(X,1);
    rng(24);
    us0 = normrnd(0,sqrt(tau2),[n,1]); us = us0(np);
    if nonrandom == 1
        us0 = us0.*0;
    end
    etas = X*betas + us; Ps = repmat(etas, [1,C-1]) + repmat(alphas,[N,1]);
    Ps = [F(Ps), ones(N,1)]; %CDF
    for c = 0:(C-2)
        Ps(:,C-c) = Ps(:,C-c) - Ps(:,C-c-1); %PDF
    end
    Y = zeros(N,1);
    for i = 1:N; Y(i) = randsample(1:C,1,true,Ps(i,:)); end
    thetas((1:n)+C-1+p) = us0';
    save('simuData.mat','Y','X','thetas')
    load('simuData.mat')
end
incs = cell(1,C);
for c = 1:C
    incs{c} = find(Y==c);
end
if nonrandom == 1
    thetas((1:n)+C-1+p) = 0;
end
thetas1 = thetas;
for i = 2:(C-1)
    thetas1(i) = log(thetas(i) - thetas(i-1));
end
thetas = thetas1;
accepts = zeros(1,ntheta);
batchLen = 50; batchNum = 0;
batchTot = tot/batchLen;
rates = zeros(batchTot, ntheta);

rng(ch0*241);
thetas = thetas.*unifrnd(.9,1.1,[1,length(thetas)]) + 1e-4;
tic
for b = 1:tot
    if verbose == 1
        fprintf('%6d', b)
        if(~mod(b,20))
            fprintf('\n')
        end
    end
    for k = 1:ntheta
        re = (k>C-1+p);
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
    % alphas = cumsum([gammas(1), exp(gammas(2:end))]);
    if nonrandom ~= 1
        tau2 = 1./gamrnd(alphatau+n/2, (invbetatau+sum(thetas((1:n)+C-1+p).^2)/2)^-1);
    end
    if b > burn
        alphas = cumsum([thetas(1), exp(thetas(2:(C-1)))]);
        matpara(b-burn,:) = [alphas, thetas(C:ntheta), tau2];
    end
    if ~mod(b,batchLen)
        if reportrate == 1
            disp(num2str(accepts/batchLen, 2))
        end
        batchNum = batchNum+1;
        accepts = accepts/batchLen;
        rates(batchNum,:) = accepts;
        tunings = tunings + sign((accepts>objrate)-0.5).*min(0.01,1/sqrt(batchNum));
        accepts = zeros(1,ntheta);
    end
end
runtime = toc/3600;
disp(num2str(rates(batchTot-4+(1:4),:)))
fprintf('\n%d iterations are done with elapsed time %.2f hours.\n', tot, runtime)
save(strcat('outCase',num2str(nvar),'_Ch',num2str(ch),'.mat'),'runtime','matpara','rates')
end

function [lik] = getlik(thetas, re, tau2) % likelihood for proportional odds model
global C X p n np incs delta2
lik = 0; betas = thetas((1:p)+C-1)'; us0 = thetas((1:n)+C-1+p)'; us = us0(np);
alphas = cumsum([thetas(1), exp(thetas(2:(C-1)))]);
for c = 1:C
    indc = incs{c};
    if c == 1
        alphac = alphas(c);
        lik = lik + sum(log( F(alphac + X(indc,:)*betas + us(indc)) ));
    elseif c == C
        alphac = alphas(C-1);
        lik = lik + sum(log( 1 - F(alphac + X(indc,:)*betas + us(indc)) ));
    else
        tmp = X(indc,:)*betas + us(indc);
        alphac = alphas(c); alphac1 = alphas(c-1);
        lik = lik + sum(log( F(alphac + tmp) - F(alphac1 + tmp) ));
    end
end
if re == 1
    lik = lik - 0.5*sum(us0.^2)/tau2;
elseif re == 0
    lik = lik - 0.5*sum(thetas(1:(C-1+p)).^2)/delta2;
end
end

function [y] = F(x)
y = exp(x)./(1+exp(x));
end

function [Yhat] = getpred(thetas)
global C X p n np
alphas = cumsum([thetas(1), exp(thetas(2:(C-1)))]); betas = thetas((1:p)+C-1)';
us0 = thetas((1:n)+C-1+p)'; us = us0(np);
N = size(X,1);
etas = X*betas + us; Ps = repmat(etas, [1,C-1]) + repmat(alphas,[N,1]);
Ps = [F(Ps), ones(N,1)]; %CDF
for c = 0:(C-2)
    Ps(:,C-c) = Ps(:,C-c) - Ps(:,C-c-1); %PDF
end
Yhat = zeros(N,1);
for i = 1:N
    Yhat(i) = find(Ps(i,:) == max(Ps(i,:)));
end
end


