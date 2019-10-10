% Bauer Alexander
% #ok<*ASGLU>
% #ok<*NASGU>

clear all java          % clears workspace
close all               % closes figures


% Script for evaluating the experiment "Enzymkinetik"
% Input: spreadsheet with:
% 1) calibration curve
% 2) enzyme kinetics
% 3) inhibitor kinetics

% load all data from spreadsheetfile:
% Lamber-Beer calibration curve data:

temp = input('Spreadsheet to process?');         % name spreadsheet to process
str_load = strcat(temp,'.xlsx');                 % format string for loading   
cal_con = xlsread(str_load, 1, 'A2:A7');         % load concentrations from spreadsheet
cal_abs = xlsread(str_load, 1, 'B2:B7');         % load absorption from speadsheet 




% 1) calibration curve:
% This section evaluates the p-nitrophenol calibration curve and calculates
% the extinction coefficient

% load data from spreadsheet file

d_cuv = 1;                                       % cuvette thickness [cm]

% perform least square fit on linear equation:
xx0 = [1 0];                                     % starting conditions
F = @(xx0,cal_con)(xx0(1)*cal_con+xx0(2));       % setup linear equation to be fitted: y = mx + c
[ahat]=lsqcurvefit(F,xx0,cal_con,cal_abs);       % perform least square fit of data on equation y = mx + c

% create idealized calibration curve:
x = linspace(min(cal_con),max(cal_con),100);     % create evenenly spaced concentration vector
y = x.*ahat(1)+ahat(2);                          % absorption vector calculated from fitted values

% calculate extinction coefficient
NP_e_true = 18.3*10^3;                           % [1/(M*cm)] 
NP_e = ahat(1)/d_cuv;                            % [1/(yM*cm)]
NP_e = NP_e/10^(-6);                             % [1/(M*cm)]
NP_error = abs((NP_e-NP_e_true)/NP_e_true*100);  % calculate percentage error



% plot measured data and calibration curve:
figure;                                          % opens new figure
hold on                                          % holds current figure
plot(x,y,'linewidth',2,'color','r')              % plots calculated calibration curve 
plot(cal_con,cal_abs,'*','markersize',8)         % plots measured calibration values
xlabel('concentration [yM]','fontsize',14)       % x-axis label
ylabel('absorption [a.u.]','fontsize',14)        % y-axis label

% clear all unnecessary values
clear exitflag jacobian lambda output residual resnorm temp xx0 F x y





% 2) Michaelis-Menten:
% This section evaluates the enzyme activity from the absorption change

% enzyme parameters:
ezy_mw = 80000*1.66*10^(-21);                    % enzyme weight [mg]
U_true = 3992;                                   % true specific activity [yM/min]
stock_con = 17.77;                               % enzyme concentration [mg/ml]
avogadro = 6.02*10^(23);                         % Avogadro's number
ezy_con = stock_con/ezy_mw;                      % enzyme concentration [n/ml]
ezy_con = ezy_con*10^3;                          % enzyme concentration [n/l]
ezy_con = ezy_con/avogadro;                      % enzyme concentration [M]
ezy_con = ezy_con*10^(6);                        % enzyme concentration [yM]
ezy_con = ezy_con/3;                             % enzyme concentration after dilution [yM]
ezy_con = ezy_con*2/1002;                        % enzyme concentration in measurement cuvette [yM]
disp_str = strcat('The enzyme concentration is: ',int2str(ezy_con),' [yM]');
disp(disp_str)

% Michaelis-Menten enzyme kinetics data
NPP_con = xlsread(str_load, 1, 'D2:J2');         % load concentration from speadsheet [yM]
temp = xlsread(str_load, 1, 'D3:J7');            % load absorption change from spreadsheet [1/min]
taqtaq = mean(temp);
NPP_speed = mean(temp)./(NP_e*d_cuv);            % change for each concentration [M/min]
NPP_speed = NPP_speed/60;                        % change for each concentration [M/s]
NPP_speed = NPP_speed*10^6;                      % convert speed [ym/s]

% Lineweaver-Burke evaluation:
i_NPP_con = 1./NPP_con;                                   % inverse substrate concentration [1/yM]
i_NPP_speed = 1./NPP_speed;                               % inverse reaction speed [s/yM]
xx0 = [1 1];                                              % initialize starting conditions
F = @(xx0,i_NPP_con)(xx0(1).*i_NPP_con+xx0(2));           % set up equation to be fitted: y = mx+c
[ahat] = lsqcurvefit(F,xx0,i_NPP_con,i_NPP_speed);        % perform least square fit of data
LB_v_max = 1/ahat(2);                                     % maximum reaction speed
LB_K_M = ahat(1)/LB_v_max;                                % Michaelis-Menten constant

% Michaelis-Menten evaluation:
xx0 = [100 10];                                           % initialize starting conditions
G = @(xx0,NPP_con)((xx0(1).*NPP_con)./(xx0(2)+NPP_con));  % setup equation to be fitted: y = (mx)/(c+x)
[para] = lsqcurvefit(G,xx0,NPP_con,NPP_speed);            % perform least square fit of data
MM_v_max = para(1);                                       % m = maximum reaction speed [ym/s]
MM_K_M = para(2);                                         % c = Michaelis-Menten constant [ym]
A_m = MM_v_max/ezy_con;                                   % molar activity [1/s]
A_s = A_m/ezy_mw;                                         % specific activity [1/(s*mg)]
A_s = A_s*60;                                             % specific activity [1/(min*mg)]
U_cal = A_s/avogadro;                                     % calculated activity [mol/(min*mg)]
U_cal = U_cal*10^6;                                       % calculated activity [ymol/(min*mg)] = [U/mg]


% calculate and plot enzyme kinetics plot:
xx = linspace(NPP_con(1),NPP_con(end));                             % assign vector of substrate concentration
MM_yy = (MM_v_max*xx)./(MM_K_M+xx);                                 % enzyme kinetics plot with Michaelis-Menten estimations
LB_yy = (LB_v_max*xx)./(LB_K_M+xx);                                 % enzyme kinetics plot with Linewaever-Burke estimations
leg = ['Michaelis-Menten';'Lineweaver-Burke';'   measured data'];   % creates charcer array for legend

figure;                                          % opens new figure
hold on                                          % hold current figure
plot(xx,MM_yy,'linewidth',2,'color','r')         % plots MM calculated values in red
plot(xx,LB_yy,'linewidth',2,'color','c')         % plots LB calculated values in cyan
plot(NPP_con,NPP_speed,'*','markersize',8)       % plots measured values as blue stars
xlabel('concentration [yM]','fontsize',14)       % x-axis label
ylabel('speed [yM/s]','fontsize',14)             % y-axis label
legend(leg)                                      % creates legend
hold off                                         % stops holding current plot

% Reaktionshemmung
% This section evaluates the enzyme reaction with constant substrate and
% varying inhibitor concentration

% load data from spreadsheet file
sub_con = 2500;                                  % substrate concentration [yM]
KHPO_con = xlsread(str_load, 1, 'L2:P2');        % load inhibitor concentration from spreadsheet [mM]
KHPO_con = 1000.*KHPO_con;                       % convert inhibitor concantration [yM]
temp = xlsread(str_load, 1, 'L3:P7');            % load absorption change from spreadsheet [1/min]
KHPO_speed = mean(temp);                         % average absorption change [1/min]
KHPO_speed = KHPO_speed./(NP_e*d_cuv);           % change of concentration [M/min]
KHPO_speed = KHPO_speed/60;                      % change of concentration [M/sec]
KHPO_speed = KHPO_speed*10^(6);                  % change of concentration [yM/sec]
i_KHPO = 1./KHPO_speed;                          % inverse change of concentration [s/yM]


% perform least square fit on concentration vs inverse reaction speed
xxx0 = [0.001 0];                                % starting conditions [m c]
H = @(xxx0,KHPO_con)(xxx0(1)*KHPO_con+xxx0(2));  % setup equation to be fitted: y = mx + c
[ahat_i] = lsqcurvefit(H,xxx0,KHPO_con,i_KHPO);  % least square fit algorithm

% calculate inhibitor constant:
K_i = MM_K_M/(ahat_i(1)*MM_v_max*sub_con);       % Inhibitor constant [yM]

base = linspace(KHPO_con(1),KHPO_con(end));      % create interpolated inhibitor concentration vector
taq1 = base*ahat_i(1)+ahat_i(2);                 % calculate inverse reaction speed graph 
taq2 = 1./taq1;                                  % convert to normal reaction speed graph

figure;                                          % opens new figure
hold on                                          % holds current figure
plot(KHPO_con,i_KHPO,'*')                        % plot inhibitor concentration vs inverse reaction speed
plot(base,taq1,'linewidth',2,'color','r')        % plot inhibitor concentration vs calculated inverse reaction speed
xlabel('inhibitor concentration [yM]')           % x-axis label
ylabel('inverse reaction speed [s/yM]')          % y-axis label

figure;                                          % opens new figure
hold on                                          % holds current figure
plot(KHPO_con,KHPO_speed,'*')                    % plot inhibitor concentration vs measured reaction speed
plot(base,taq2,'color','r')                      % plot inhibitor concentration vs caclulated reaction speed
xlabel('inhibitor concentration [yM]')           % x-axis label
ylabel('reaction speed [yM/s]')                  % y-axis label


% display results on screen
disp_str = strcat('The measured extinktion coefficient is: ',num2str(NP_e),' 1/(M*cm)');       % create string for displaying extinction coefficient 
disp(disp_str);                                                                                % displays extinction coefficient [1/cm  1/M]
temp = NP_e*10^(-6);                                                                           % convert extinction coefficient to 1/(yM*cm)
disp_str = strcat('The measured extinktion coefficient is: ',num2str(temp),' 1/(yM*cm)');      % create string for displaying extinction coefficient 
disp(disp_str);                                                                                % displays extinction coefficient [1/cm  1/M]
disp_str = strcat('The error is: ',num2str(NP_error),' [%]');                                  % create string for displaying error
disp(disp_str);                                                                                % displays percentage error [%]


disp_str = strcat('MM_v_max is: ',num2str(MM_v_max),' yM/s');              % create string for displaying maximum reaction speed
disp(disp_str);                                                            % displays MM-maximum reaction speed
disp_str = strcat('MM_v_max is: ',num2str(MM_v_max/1000*60),' mM/min');    % create string for displaying maximum reaction speed
disp(disp_str);                                                            % displays MM-maximum reaction speed
disp_str = strcat('LB_v_max is: ',num2str(LB_v_max),' yM/s');              % create string for displaying maximum reaction speed
disp(disp_str);                                                            % displays LB-maximum reaction speed
disp_str = strcat('LB_v_max is: ',num2str(LB_v_max/1000*60),' mM/min');    % create string for displaying maximum reaction speed
disp(disp_str);                                                            % displays LB-maximum reaction speed


disp_str = strcat('MM_K_M is: ',num2str(MM_K_M),' yM');                    % create string for displaying Michaelis-Menten constant
disp(disp_str);                                                            % displays MM constant
disp_str = strcat('LB_K_M is: ',num2str(LB_K_M),' yM');                    % create string for displaying Michaelis-Menten constant
disp(disp_str);                                                            % displays MM constant






