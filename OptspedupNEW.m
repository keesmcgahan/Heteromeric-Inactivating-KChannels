
%[ACTdata]= ActCurveKv14;
 [locmins]=kvchannelopt;


function [locmins]=kvchannelopt

data=[2.53037669594676e-05,0.000306955085983109,0.00318267445758426,0.0234391649515618,0.106428124809781,0.290572936654734,0.525175499431972,0.722675172818586,0.850818511396516,0.922745879274974,0.960272475330649,0.979234472674357,0.988734287155502,0.993540899092112,0.996026444204975];
dataIN=[0.995179240676961,0.990716294531252,0.967944920421038,0.890890472070149,0.749604733812209,0.614277138879068,0.533447723324825,0.495647131351706,0.479488388991891,0.472601465525254,0.469569551475566,0.468193673639994,0.467556335198404,0.467283050872215,0.467261588400973];
dataTaus=[0.100000000000000,0.177999846575680,1.86604216362275,2.79595967407082,3.46019335091557,3.69083640913037,3.52230921494050,3.12680744773204,2.66443203472987,2.22170669696517,1.83284580742685,1.50404164150028,1.23083692094569,1.00594468106846,0.821646096258852];
dataRec=[0.461386870435747,0.461386870435747,0.461386870435747,0.461896118243928,0.462686464304096,0.464263305723347,0.468965623230917,0.476710678673144,0.534704053615712,0.598159984081556,0.700034381852640,0.873855863466406,0.928047813567986,0.967469631892407,0.993515524263076];
lbb = 0;
ubb = 1000;


[parameters,locmins] = neuron_parameterization(data,dataTaus,dataIN,dataRec,lbb,ubb);


[ACTdata,Vsteps,Taus]=ActCurveKv14(parameters);
[INACTdata,VstepsIN]=InActCurve(parameters);
[RECdata,VstepsRec]=InacRecov(parameters);

disp('Parameters')
disp(parameters)

figure(1)
hold on
plot(Vsteps,data,'bo')
plot(Vsteps,ACTdata,'r*')
hold off

figure(2)
hold on
plot(VstepsIN,dataIN,'bo')
plot(VstepsIN,INACTdata,'r*')
hold off

figure(3)
hold on
plot(Vsteps,Taus,'r*')
plot(Vsteps,dataTaus,'bo')
hold off

figure(4)
hold on
plot(VstepsRec,dataRec,'bo')
plot(VstepsRec,RECdata,'r*')
hold off

end

function [ACTdata,Vsteps,Taus]=ActCurveKv14(parms)
    prot=1;
    t=[0 50];
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
    output=zeros(1,length(Vsteps));
    ACTdata=zeros(1,length(Vsteps));
    Taus=zeros(1,length(Vsteps));
    Ninfs=zeros(1,length(Vsteps));  
    
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    [time,S]=voltageclamp(t,Vstep,prot,parms);
    output(i)=S(end,5);
    k = find(time>=20,1);
    j = find(time>=50,1)-1;
    ACTdata(i)=max(S(k:j,5));
    
    fo = fitoptions('Method','NonlinearLeastSquares','Lower',[ACTdata(i).^(1/4),0],'Upper',[ACTdata(i).^(1/4),10],...
               'StartPoint',[ACTdata(i).^(1/4),.1]);
    ft = fittype('(a*(1-exp(-(x-20)/b)))^4','options',fo);
    [curve4,fof1]=fit(time(k:j),(S(k:j,5)),ft);
    Taus(i)=curve4.b;
    Ninfs(i)=curve4.a;
    
    
 end

 
end  

function [INACTdata,Vsteps]=InActCurve(parms)
    prot=2;
    
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
    %output=zeros(1,length(Vsteps));
    INACTdata=zeros(1,length(Vsteps));
    
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    t=[0 5000];
        
    [sol]=voltageclamp2(t,Vstep,prot,parms);
    endtime=6000;
    [solext]=vclampextend(t,Vstep,prot,sol,endtime,parms);
    
    time=solext.x;
    yvals=solext.y;
     k = find(time>=5000,1);
     INACTdata(i)=max(yvals(5,k:end));
   
 end
 
 
end  

function [Recover,Vsteps]=InacRecov(parms)
    prot=3;
    %t=[0:0.1:500];
    Vsteps=[.1,.5,1,5,10,20,50,100,500,1000,2000,5000,7000,10000,20000];
    %Vsteps=[.1];
     Recover=zeros(1,length(Vsteps));
     t=[0 5200];
     Vstep=0;
     [sol]=voltageclamp2(t,Vstep,prot,parms);
    
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
        
        
    
    endtime=6000+Vstep;
    [solext]=vclampextend(t,Vstep,prot,sol,endtime,parms);
    
    time=solext.x;
    yvals=solext.y;
     k = find(time>=5200,1);
     Recover(i)=max(yvals(5,k:end));
   
 end

 
end 

function [parameters,locmins] = neuron_parameterization(data,dataTaus,dataIN,dataRec,lbb,ubb)


 % % Establish initial conditions for the ODE
parms=[700*rand(1),5*rand(1),500*rand(1),5*rand(1),rand(1),rand(1)];% % Initial guess for parameters
%parms=[561.9305,    2.1720,   15.6349  ,  4.9635,    0.4396,    0.5006];
lb = [0,0,0,0,0,0];   % % lower bound vector
ub = [700,10,500,10,1,1];   % % upper bound vector


% % % % fmincon is the optimization function used to produce MLE estimates
% % % % in our MLE procedure. The error function called is based on the NLL

tic
%parms = fmincon(@(x) error(data,dataTaus,dataIN,dataRec,...
%   [x(1) x(2) x(3) x(4) x(5) x(6)]),parms,[],[],[],[],lb,ub);

% gs=GlobalSearch;
% problem=createOptimProblem('fmincon','x0',parms,...
%      'objective',@(x) error(data,dataTaus,dataIN,dataRec,...
%    [x(1) x(2) x(3) x(4) x(5) x(6)]),'lb',lb,'ub',ub);
% 
% [parms,fg,exitflag,output,locmins]=run(gs,problem);

[xga,fval,efga,outga,population,scores] = ga(@(x)error(data,dataTaus,dataIN,dataRec,...
    [x(1) x(2) x(3) x(4) x(5) x(6)]),6,[],[],[],[],lb,ub);

 locmins=horzcat(population,scores);

toc

err=error(data,dataTaus,dataIN,dataRec,parms)
% % % % To calculate BIC an estimate of varience (stderror) is needed
% % % % as described in the BIC Estimation Section.
% [~,hz]=ode45(@(t,y)hzode(t,y,parms),t,IC); % % 
% N = length(h); % % number of data points
% v = length(parms); % % number of fit parameters
% errs2 = (hz(:,1)-h).^2+(hz(:,3)-z).^2; 
% stderror = 1/(2*N)*sum(errs2);%Variance of the data from the model estimate
% NLL=N*log(2*pi*stderror)+1/(2*stderror)*sum(errs2);
% BIC=2*NLL+v*log(N); % % BIC for parameterized model given the data

parameters=parms; % % parameter estimates for normalized data


end

%%%%%I think i need to fee parms and IC into ACTCURVEKV14 and the rest of
%%%%%the funcs?
function err=error(data,dataTaus,dataIN,dataRec,parms)
    
% % % % For a given parameter set, the ODE is solved numerically using
% % % % ode45
[ACTdata,Vsteps,Taus]=ActCurveKv14(parms);
[INACTdata,Vsteps]=InActCurve(parms);
[Recover,Vsteps]=InacRecov(parms);
% % % % the error between the model estimates and the data is calculated,
% % % % and subsequently minimized by fmincon

err=sum((INACTdata-dataIN).^2)+sum((ACTdata-data).^2)+sum((Recover-dataRec).^2)+sum(((dataTaus-Taus).^2)/3.7);
%err=sum((ACTdata-data).^2);
%err=sum((INACTdata-dataIN).^2);
%err=sum((dataTaus-Taus).^2);
%+sum(((dataTaus-Taus).^2)/3.7)+sum((INACTdata-dataIN).^2)
%sum((ACTdata-data).^2)+sum(((dataTaus-Taus).^2)/3.7)+
end


function [time,S]=voltageclamp(t,Vstep,prot,parms)

S0 = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]'; 
options = odeset('RelTol',1e-6,'AbsTol',1e-6); 
[time,S]=ode15s(@(t,S)deRHS(t,S,Vstep,prot,parms),t,S0,options);
S(end,:);

end

function [sol]=voltageclamp2(t,Vstep,prot,parms)

S0 = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]'; 
options = odeset('RelTol',1e-6,'AbsTol',1e-6); 
[sol]=ode15s(@(t,S)deRHS(t,S,Vstep,prot,parms),t,S0,options);


end

function [solext]=vclampextend(t,Vstep,prot,sol,endtime,parms)

solext=odextend(sol,@(t,S)deRHS(t,S,Vstep,prot,parms),endtime);

end


function [V]=SSActkv14paper(t,Vstep)

    if t<20
        V=-90;
    elseif t>20 && t<50
        V=Vstep;
    elseif t>50 && t<80
        V=-50;
    else 
        V=-90;
    end
end

function [V]=SSinact(t,Vstep)

if t<5000
    V=Vstep;
else
    V=50;    
end

end

function [V]=RecovProt(t,Vstep)

if t<20
    V=-90;
elseif t>=20 && t<5200
    V=50;
elseif t>=5200 && t<(Vstep+5200)
    V=-90;
else
    V=50;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%ODE%%%%%%%%%%%%%%%%%%%5

function s_prime=deRHS(t,S,Vstep,prot,parms) 

    if prot==1
        V=SSActkv14paper(t,Vstep);
    elseif prot==2
        V=SSinact(t,Vstep);
    elseif prot==3    
        V=RecovProt(t,Vstep);
    end

%Variables
% p0=S(1);
% p1=S(2);
% p2=S(3);
% p3=S(4);
% O=S(5);
%I=S(6);

p00=S(1);
p10=S(2);
p20=S(3);
p30=S(4);
%P40 is open state
p40=S(5);
p11=S(6);
p21=S(7);
p31=S(8);
p41=S(9);
p22=S(10);
p32=S(11);
p42=S(12);
p33=S(13);
p43=S(14);
p44=S(15);


%Parameters
% m1=.8;
% m2=.01;
% n1=.0001;
% n2=.0001;
% c1=.025;
% c2=10.35;
% d1=.01;
% d2=.01;

% k1=parms(1);
% k2=parms(2);
% k3=parms(3);
% k4=parms(4);
% k5=parms(5);
% k6=parms(6);
% k7=parms(7);


m1=parms(1);
n1=parms(2);
c1=parms(3);
d1=parms(4);
m2=parms(5);
n2=parms(6);

 a1=m1*exp(n1*V*.01)/1000;
 %a2=m2*exp(-n2*V*.01)/1000;
 b1=c1*exp(-d1*V*.01)/1000;
 %b2=c2*exp(d2*V*.01)/1000;
 a2=m2/1000;
 b2=n2/1000;
% a1=k1*exp((V-(100*k3))/(100*k4));
% b1=k2*exp((V-(100*k3))/(100*k4)).*exp(-(V+(100*k5))/(100*k6))./(1+(k7/100)*exp(-(V+(100*k5))/(100*k6)));

%think about spitting out rates to look at time constants?


s_prime=[-4*a1*p00+b1*p10,...
    4*a1*p00-b1*p10-3*a1*p10+2*b1*p20-a2*p10+b2*p11,...
    3*a1*p10-2*b1*p20-2*a1*p20+3*b1*p30-2*a2*p20+b2*p21,...
    2*a1*p20-3*b1*p30-1*a1*p30+4*b1*p40-3*a2*p30+b2*p31,...
    a1*p30-4*b1*p40-4*a2*p40+b2*p41,...
    a2*p10-b2*p11-3*a1*p11+b1*p21,...
    2*a2*p20-b2*p21+3*a1*p11-b1*p21-2*a1*p21+2*b1*p31-a2*p21+2*b2*p22,...
    3*a2*p30-b2*p31+2*a1*p21-2*b1*p31-a1*p31+3*b1*p41-2*a2*p31+2*b2*p32,...
    4*a2*p40-b2*p41+a1*p31-3*b1*p41-3*a2*p41+2*b2*p42,...
    a2*p21-2*b2*p22-2*a1*p22+b1*p32,...
    2*a2*p31-2*b2*p32+2*a1*p22-b1*p32-a2*p32+3*b2*p33-a1*p32+2*b1*p42,...
    3*a2*p41-2*b2*p42+a1*p32-2*b1*p42-2*a2*p42+3*b2*p43,...
    a2*p32-3*b2*p33-a1*p33+b1*p43,...
    2*a2*p42-3*b2*p43+a1*p33-b1*p43-a2*p43+4*b2*p44,...
    a2*p43-4*b2*p44]';
end