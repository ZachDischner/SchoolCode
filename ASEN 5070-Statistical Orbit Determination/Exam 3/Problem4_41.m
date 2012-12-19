% Problem 4-41 in StatoD book

clc;clear all;close all
set(0,'defaulttextinterpreter','latex')

obs     =load('Exam3_Problem5_data.txt');
time    =obs(:,1);
Y       =obs(:,2);

% Beta    = 0.02;
% sigma   = 0.67;
Beta = 1.45;
sigma = 10.49;

Eta(1) = Y(1);

Truth = sin(2*pi*time/10);

for ii = 1:length(Y)
    
    Htilde(ii) = 1;
    R(ii) = 1;
    Q(ii) = 1;
    
    if ii > 1
        m(ii) = exp(-Beta*(time(ii)-time(ii-1)));
        Gamma(ii) = sqrt(sigma^2/(2*Beta)*(1-m(ii)^2));
        PhiStep(ii) = m(ii);
        EtaBar(ii)  = PhiStep(ii)*EtaHat(ii-1);
        Pbar(ii)    = PhiStep(ii)*P(ii-1)*PhiStep(ii)' + Gamma(ii)*Q(ii-1)*Gamma(ii)';
    else
        m(ii) = exp(-Beta*(time(ii)));
        Gamma(ii) = sqrt(sigma^2/(2*Beta)*(1-m(ii)^2));
        Eta(ii) = m;
        Pbar(ii) = 1;
        EtaBar(ii) = 1;
    end

    K(ii) = Pbar(ii)/(Pbar(ii)+1);
    EtaHat(ii) = EtaBar(ii) + K(ii)*(Y(ii)-EtaBar(ii));
    P(ii) = (eye(size(K(ii)*Htilde(ii))) - K(ii)*Htilde(ii))*Pbar(ii);
    
    
end



figure
plot(Truth)
hold on
plot(Y,'r.')
plot(EtaHat,'k')


legend('Truth','Observations','Best Estimate, \eta^{\wedge}')
xlabel('Time')
ylabel('Force ($\eta$)')

figure
hist((Y-EtaHat'),50)
xlabel('Bins')
ylabel('Y - $\eta$')
title('Gaussian Noise of Residules')

% Fit Result to Sin 
[fitresult, gof] = createFit(time, EtaHat);
% Eta_fcn = feval(fitresult,time);
%Convert cfit to a function handle
cfit2mfile(fitresult,'eta_force');
% Evaluate 'backed out' force function
T2 = eta_force(time);

figure
plot(time,Truth)
hold on
plot(time,T2,'k')
xlabel('Time')
ylabel('Amplitude')
legend('True Force','Sin Fit Function')


%% Now do it again, correcting for the modeled force
clear y Htilde R Q m Gamma PhiStep Pbar EtaBar Eta Hat K P
Y       =obs(:,2) - eta_force(time);

Beta    = 0.02;
sigma   = 0.67;

Eta(1) = Y(1);

Truth = zeros(length(time),1);%sin(2*pi*time/10);

for ii = 1:length(Y)
    
    Htilde(ii) = 1;
    R(ii) = 1;
    Q(ii) = 1;
    
    if ii > 1
        m(ii) = exp(-Beta*(time(ii)-time(ii-1)));
        Gamma(ii) = sqrt(sigma^2/(2*Beta)*(1-m(ii)^2));
        PhiStep(ii) = m(ii);
        EtaBar(ii)  = PhiStep(ii)*EtaHat(ii-1);
        Pbar(ii)    = PhiStep(ii)*P(ii-1)*PhiStep(ii)' + Gamma(ii)*Q(ii-1)*Gamma(ii)';
    else
        m(ii) = exp(-Beta*(time(ii)));
        Gamma(ii) = sqrt(sigma^2/(2*Beta)*(1-m(ii)^2));
        Eta(ii) = m;
        Pbar(ii) = 1;
        EtaBar(ii) = 1;
    end

    K(ii) = Pbar(ii)/(Pbar(ii)+1);
    EtaHat(ii) = EtaBar(ii) + K(ii)*(Y(ii)-EtaBar(ii));
    P(ii) = (eye(size(K(ii)*Htilde(ii))) - K(ii)*Htilde(ii))*Pbar(ii);
    
    
end

figure
plot(Truth)
hold on
plot(Y,'r.')
plot(EtaHat,'k')


legend('Truth','Observations','Best Estimate, \eta^{\wedge}')
xlabel('Time')
ylabel('Force ($\eta$)')

figure
hist((Y-EtaHat'),50)
xlabel('Bins')
ylabel('Y - $\eta$')
title('Gaussian Noise of Residules')

figure_awesome('save')


