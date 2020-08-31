

clear;
% myCluster = parcluster('local');
% myCluster.NumWorkers = 12; % 'Modified' property now TRUE.
% saveProfile(myCluster);
%parpool(12);

%ps = [0.003 0.004 0.006 0.007 0.008 0.009 0.011 0.012 0.013 0.014 0.016 0.017 0.018 0.019 0.02 0.021 0.022 0.023 0.024];
ps = [0.022 0.023 0.024];
idxs = [1 2 3 4];
%idxs = [4];
% CA's counts by 4, Florida's by 3, New York's by 2, and leave Washington's alone
mults = [4 2 3 1];
d = 5;

use_mult = false;
alt_theta = false;
if use_mult
   mult_str = 'mult'; 
else
   mult_str = ''; 
end

if alt_theta
   theta_case = 'alt'; 
else
   theta_case = '';
end

plotting = false;
repeat_boot = true;
states = {'California' 'New York' 'Florida' 'Washington'};
%Ns = [4e7 2e7 2.5e7 1.2e7];
dates = {'15-Mar-2020','18-Mar-2020','22-Mar-2020','17-Mar-2020'};
dates = datetime(dates);
n_states = length(states);

%for r=1:5
for p = ps
    for idx = idxs
        disp(['state: ' states{idx} ' p=' num2str(p)]);
        
        p_str = strrep(num2str(p),'.','_');
        B = 1000;
        
        
        opts = detectImportOptions('~/Documents/GitHub/covid-19-data/us-states.csv');
        dat = readtable('~/Documents/GitHub/covid-19-data/us-states.csv',opts);
        pops = readtable('~/Documents/GitHub/Covid19-undercount/analyze/input/pop_data.csv');
        if alt_theta
            theta = readmatrix('~/Documents/GitHub/Covid19-undercount/analyze/input/theta2.csv');
        else            
            theta = readmatrix('~/Documents/GitHub/Covid19-undercount/analyze/input/theta.csv');
        end
        theta = theta(:,3);
        Ns = zeros(1,n_states);
        for j = 1:n_states
            Ns(j) = pops.Pop(strcmp(states{j},pops.State));
        end        
        
        % weigthing matrix for theta       
        
        id_keep = ismember(dat.state,states);
        
        dat0 = dat(id_keep,:);
        dat0 = sortrows(dat0,[2 1]);
        min_date = datetime('01-Jan-2020');
        %min_date = min(dat0.date);
        %max_date = max(dat0.date);
        max_date = datetime('30-Apr-2020');
        %max_date = datetime('11-May-2020');
        n_days = days(max_date-min_date);
        dat0.dayof = days(dat0.date-min_date)+1;
        dayof_sd = days(dates-min_date)+1;
        %min_state = accumarray(double(categorical(dat0.state)),dat0.dayof,[],@min);
        
        
        D = zeros(n_days+1,n_states);
        
        for state = states
            disp(state);
            D(dat0.dayof(strcmp(state,dat0.state)),strcmp(state,states)) = [0;diff(dat0.deaths(strcmp(state,dat0.state)))];            
        end
        D = D(1:n_days+1,:);
        
        if use_mult
           for j = idxs
              D(:,j) = mults(j).*D(:,j); 
           end
        end
        
        T2 = size(D,1);
        w_theta = zeros(T2,T2);
        n_lag = length(theta);
        for j=1:T2
            w_theta(j,j:min(j+n_lag-1,T2)) = theta(1:min(j+n_lag-1,T2)-j+1);
        end
        
        
        % new york
        %Dn = D(1:end-1,2);
        Dn = D(:,idx);
        %Dn = movmean(Dn,7);
        x = [20.1 6 2.5 0.3 1.5];
        [t,y] = sir_sim(x,dayof_sd(idx),T2,Ns(idx));
        
        nu = -diff(y(:,1));
        nu = [zeros(T2-length(nu),1);nu];
        lam = d_means(w_theta,nu,p,Ns(idx));
        
        
        figure(1); plot(1:T2,Dn,1:T2,lam,'.'); drawnow;
        
        
        zeta = [1 find(Dn>0,1)-1; 1 15; 1 4; 0.2 0.9; 1 4];

        
        [x_hat,H] = mle_sir(x,dayof_sd(idx),T2,p,Ns(idx),Dn,w_theta,zeta);
        disp(num2str(x_hat));
        
        [t,y] = sir_sim(x_hat,dayof_sd(idx),T2,Ns(idx));
        
        nu = -diff(y(:,1));
        nu = [zeros(T2-length(nu),1);nu];
        lam = d_means(w_theta,nu,p,Ns(idx));
        
        figure(2); plot(1:T2,Dn,1:T2,lam,'.'); drawnow;
        
        
        % J = [1 0 0 0; 0 1 0 0; 0 x_hat(3) x_hat(2) 0; 0 0 0 1];
        % C = H\eye(4);
        % S = J'*C*J;
        % save('../output/S.mat','S');
        %
        % save('../output/mles.mat','x_hat','C');
        
        if repeat_boot
            % parametric bootstrap
            Xhat = NaN(B+10,d);
            b = 1;
            parfor b=1:B+10
                if mod(b,100)==0; disp(num2str(b)); end
                D0 = poissrnd(lam');
                try
                    x_hat0 = mle_sir(x_hat,dayof_sd(idx),T2,p,Ns(idx),D0,w_theta,zeta);
                    Xhat(b,:) = x_hat0;
                catch
                    warning(['optimization failed at iterate ' num2str(b) ' repeating this iterate with a new sample']);
                end
            end
            Xhat = Xhat(~isnan(Xhat(:,1)),:);
            Xhat = Xhat(1:B,:);
            
            C = cov(Xhat);
            
            save(strcat('../output/bootC_',states{idx},'_',p_str,'_',mult_str,'_',theta_case,'.mat'),'C');
        else
            load(strcat('../output/bootC_',states{idx},'_',p_str,'_',mult_str,'_',theta_case,'.mat'));
        end
        
        if idx==1
           %c1 = 0.04;
           %c0 = 0.0005;        
           c1 = 0.1;
           c0 = 0.001;
        elseif idx==2
           %c1 = 0.1;
           %c0 = 0.01;
           c1 = 0.3;
           %c0 = 0.03;
           %c1 = 0.05;
           c0 = 0.002;
        elseif idx==3 
           %c1 = 0.1;
           %c0 = 0.01;
           c1 = 0.3;
           c0 = 0.05;
        else
           c1 = 0.2;
           c0 = 0.05;
        end
                   
        
        nmc = 100000;
        adapt = 20000;
        T1 = dayof_sd(idx); N = Ns(idx);
        disp_int = 5000;
        
        
        % mcmc
        % compute log likelihood and prior at the initial state
        % sample the initial state from normal with mean the mle 
        % and covariance the bootstrap estimate,
        % truncated to the prior support
        cont = true;
        while cont
            x = mvnrnd(x_hat,C);
            if ~isinf(logprior(x,zeta))
               cont = false; 
            end
        end
        [l,z,nu,lam] = loglik(x,T1,T2,p,N,Dn,w_theta);
        lp = logprior(x,zeta);
        lt = l+lp;
        
        params = {'$T_0$','$\gamma^{-1}$','$R_0$','$\phi$','$\eta$'};
        ACC = false(nmc,1);
        X = zeros(nmc,d);
        SIR = zeros(T2,2,nmc);
        NU = zeros(nmc,T2);
        LAM = zeros(nmc,T2);
        
        for t = 1:nmc
            if t < adapt
                y = mvnrnd(x,c0.*C);
            elseif t==adapt
                X_bar = mean(X(1:t-1,:),1);
                SXXt = X(1:t-1,:)'*X(1:t-1,:);
                S_X = 1/(t-2).*SXXt - (t-1)/(t-2).*X_bar'*X_bar;
                S_X = c1.*S_X;
                C0 = chol(S_X);
                y = x + normrnd(0,1,[1 d])*C0;
                
                %y = mvnrnd(x,S_X);
            else
                X_bar = ((t-2).* X_bar + x)./(t-1);
                SXXt = SXXt + x'*x;
                S_X = 1/(t-2).*SXXt-(t-1)/(t-2).* X_bar'*X_bar;
                S_X = c1.*S_X;
                C0 = chol(S_X);
                y = x + normrnd(0,1,[1 d])*C0;
                
                %y = mvnrnd(x,S_X);
            end
            
            lp = logprior(y,zeta);
            if ~isinf(lp)
                [l,zp,nup,lamp] = loglik(y,T1,T2,p,N,Dn,w_theta);
                ltprop = l+lp;
                
                acc = exp(ltprop-lt)>rand;
                if acc
                    x = y;
                    lt = ltprop;
                    z = zp;
                    nu = nup;
                    lam = lamp;
                end
            else
                acc = false;
            end
            
            ACC(t) = acc;
            X(t,:) = x;
            SIR(ceil(x(1)):end,:,t) = z(:,1:end-1);
            NU(t,:) = nu;
            LAM(t,:) = lam;
            
            if mod(t,disp_int)==0
                disp(['acceptance rate: ' num2str(mean(ACC(1:t)))]);
                if plotting
                    figure(3);
                    for j=1:d
                        subplot(ceil(sqrt(d)),ceil(sqrt(d)),j); plot(1:t,X(1:t,j),'.'); title(params{j},'interpreter','latex');
                    end
                    drawnow;
                    figure(4); plot((1:T2)',mean(LAM(1:t,:),1),(1:T2)',Dn,(1:T2)',movmean(Dn,7));drawnow;
                end
            end
            
        end
        
        save(strcat('~/Dropbox/Projects/Covid19-undercount/mcmc_',states{idx},'_',p_str,'_',mult_str,'_',theta_case,'.mat'),'X','NU','SIR');
        %save(strcat('~/Dropbox/Projects/Covid19-undercount/mcmc_',states{idx},'_',p_str,'_',mult_str,'_',theta_case,'_',num2str(r),'.mat'),'X','NU','SIR');
        %save(strcat('../output/mcmc_',states{idx},'_',p_str,'_',mult_str,'_',theta_case,'.mat'),'X','NU','SIR');
        
        
    end
end
%end


function [x_hat,H] = mle_sir(x,T1,T2,p,N,D,w_theta,zeta)

    opts = optimoptions('fmincon','display','off');
    [x_hat,~,~,~,~,~,H] = fmincon(@LL,x,[],[],[],[],zeta(:,1),zeta(:,2),[],opts);

    function l = LL(x)
        [~,y] = sir_sim(x,T1,T2,N);
        nu = -diff(y(:,1));
        nu = [zeros(T2-length(nu),1);nu];
        lam = d_means(w_theta,nu,p,N);
        tmp = D.*log(lam')-lam';
        l = -sum(tmp(~isnan(tmp)));   
    end
end



function [t,y] = sir_sim(x,T1,T2,N)
    T0 = x(1); S0 = 1-1/N; gamma = 1/x(2); beta = x(3)*gamma; phi = x(4); eta = x(5);
    tspan = [T0 T1];
    y0 = [S0,1-S0,0];        

    % works   
    sol = ode45(@(t,y) sir(y,beta,gamma), tspan, y0);    
    t1 = (ceil(T0):T1)';
    y1 = deval(sol,t1)';
    
    y0 = y1(end,:);
    tspan = [T1 T2];
    %tspan = [0 T2-T1];
    sol = ode45(@(t,y) sir(y,beta*phi,gamma*eta), tspan, y0);    
    
    t2 = (T1+1:T2)';
    %t2 = (1:(T2-T1))';
    y2 = deval(sol,t2)';
    y = [y1;y2];
    t = [t1;t2];    
end



function dydt = sir(y,beta,gamma)
    dydt = [-beta*y(1)*y(2); beta*y(1)*y(2)-gamma*y(2); gamma*y(2)];
end


function lam = d_means(w_theta,nu,p,N)
    wnu = bsxfun(@times,nu.*p.*N,w_theta);
    lam = sum(wnu,1);
end


function [l,sir,nu,lam] = loglik(x,T1,T2,p,N,D,w_theta)
    [~,sir] = sir_sim(x,T1,T2,N);
    nu = -diff(sir(:,1));
    nu = [zeros(T2-length(nu),1);nu];
    lam = d_means(w_theta,nu,p,N);
    tmp = D.*log(lam')-lam';
    l = sum(tmp(~isnan(tmp)));
end


function l = logprior(x,zeta)
    is_gtr = x(:)>zeta(:,1);
    is_less = x(:)<zeta(:,2);
    if all(is_gtr & is_less)
       l = 0; 
    else
       l = -Inf;
    end
end











