clear all
clc
%[time,S]=voltageclamp(t);
%[Actend]=ActCurve;
%[Inactend]=InActCurve;
 %InacRecovCurve
 %[output,ACTdata]= ActCurveKv14;
 [INACTdata]=InActCurve;
%[Recover]=InacRecovCurve;

 
function [output]=ActCurve
    prot=1;
    t=[0:0.1:500];
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50,60,70,90];
    output=zeros(1,length(Vsteps));
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    [time,S]=voltageclamp(t,Vstep,prot);
    output(i)=S(end,5);
    figure(1)
    hold on
    plot(time,S(:,5))
  end
end  

function [INACTdata,output]=InActCurve
    prot=2;
    t=[0:0.1:6000];
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
    output=zeros(1,length(Vsteps));
    INACTdata=zeros(1,length(Vsteps));
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    [time,S]=voltageclamp(t,Vstep,prot);
    output(i)=S(end,5);
    INACTdata(i)=max(S(50000:end,5));
    figure(3)
    hold on
    plot(time,S(:,5))
 end
%   figure(4)
%  plot(Vsteps,INACTdata,'o')
end  

function [Recover]=InacRecovCurve
    prot=3;
    %t=[0:0.1:500];
    Vsteps=[.1,.5,1,5,10,20,50,100,500,1000,2000,5000,7000,10000,20000];
    %Vsteps=[.1];
     Recover=zeros(1,length(Vsteps));
     
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
        t=[0:.1:6000+Vstep];
    [time,S]=voltageclamp(t,Vstep,prot);
    Recover(i)=max(S(52000:end,5));
    figure(5)
    hold on
    plot(time,S(:,5))
 end
%   figure(6)
%    plot(Vsteps,Recover,'o')
end  

function [output,ACTdata]=ActCurveKv14
    prot=4;
    t=[0:0.1:120];
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
    output=zeros(1,length(Vsteps));
    ACTdata=zeros(1,length(Vsteps));
    Taus=zeros(1,length(Vsteps));
    Ninfs=zeros(1,length(Vsteps));
    fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0],'Upper',[1,10],...
               'StartPoint',[.1,.1]);
    ft = fittype('(a*(1-exp(-(x-20)/b)))^4','options',fo);
    
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    [time,S]=voltageclamp(t,Vstep,prot);
    output(i)=S(end,5);
    ACTdata(i)=max(S(200:500,5));
    
    [curve4,fof1]=fit(time(200:500),(S(200:500,5)),ft);
    Taus(i)=curve4.b;
    Ninfs(i)=curve4.a;
    figure(1)
    hold on
    %plot(time(200:500),(S(200:500,5)),'o')
    plot(time,S(:,5))
    %plot(time(200:500),curve4)
    
    
 end
    Taus
    Ninfs
   figure(2)
 plot(Vsteps,ACTdata,'o')
end  



function [time,S]=voltageclamp(t,Vstep,prot)

S0 = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]'; 
[time,S]=ode45(@(t,S)deRHS(t,S,Vstep,prot),t,S0);
S(end,:);

end



function [V]=SSact(t,Vstep)

if t<300
    V=-100;
else
    V=Vstep;    
end

end



function [V]=SSinact(t,Vstep)

if t<5000
    V=Vstep;
else
    V=50;    
end

end

function [V]=InacRecov(t,Vstep)

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

function [V]=SSActkv14paper(t,Vstep)

    if t<20
        V=-90;
    elseif t>=20 && t<=50
        V=Vstep;
    elseif t>50 && t<100
        V=-50;
    else 
        V=-90;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%ODE%%%%%%%%%%%%%%%%%%%5

function s_prime=deRHS(t,S,Vstep,prot) 

    if prot==1
        V=SSact(t,Vstep);
    elseif prot==2
        V=SSinact(t,Vstep);
    elseif prot==3
        V=InacRecov(t,Vstep);
    elseif prot==4
        V=SSActkv14paper(t,Vstep);
    end

%Variables
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


parms=[112.0378    2.1154    1.3872    4.9630    0.0881    0.9845    7.7309    1.4008];
m1=parms(1);
n1=parms(2);
c1=parms(3);
d1=parms(4);
m2=parms(5);
n2=parms(6);
c2=parms(7);
d2=parms(8);


%%%%%%%rates of transitions

 a1=m1*exp(n1*V*.01);
 a2=m2*exp(-n2*V*.01);
 b1=c1*exp(-d1*V*.01);
 b2=c2*exp(d2*V*.01);
 %a2=m2;
 %b2=n2;

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