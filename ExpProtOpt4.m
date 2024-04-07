% clear all
% clc
%[time,S]=voltageclamp(t);
%[Actend,time]=ActCurve;
%[Inactend]=InActCurve;
 %InacRecovCurve
   [output,ACTdata,time]= ActCurveKv14;
%   [INACTdata,timeIN]=InActCurve;
%[Recover,time]=InacRecovCurve;
 [Recoverext,timeext,yvals]=InacRecovExt;
  [INACTdataext,timeext]=InActCurveExt;

function [Recover,time,yvals]=InacRecovExt
    prot=3;
    %t=[0:0.1:500];
    Vsteps=[.1,.5,1,5,10,20,50,100,500,1000,2000,5000,7000,10000,20000];
    %Vsteps=[.1];
     Recover=zeros(1,length(Vsteps));
     t=[0 5200];
     Vstep=0;
     [sol]=voltageclamp2(t,Vstep,prot);
 tic    
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
        
        
    
    endtime=6000+Vstep;
    [solext]=vclampextend(t,Vstep,prot,sol,endtime);
    
    time=solext.x;
    yvals=solext.y;
     k = find(time>=5200,1);
     Recover(i)=max(yvals(5,k:end));
     figure(5)
     hold on
     plot(time,yvals(5,:))
 end
toc
  figure(6)
  hold on
   plot(Vsteps,Recover,'o')
  xlabel('\Delta t (ms)')
  ylabel('Ratio P2 to P1')
end   


function [output,time]=ActCurve
    prot=1;
    t=[0 500];
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

function [INACTdata,time]=InActCurve
    prot=2;
    t=[0 6000];
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
    output=zeros(1,length(Vsteps));
    INACTdata=zeros(1,length(Vsteps));
    tic
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    [time,S]=voltageclamp(t,Vstep,prot);
    output(i)=S(end,5);
    k = find(time>=5000,1);
    INACTdata(i)=max(S(k:end,5));
    figure(3)
    hold on
    plot(time,S(:,5))
 end
 toc
  figure(4)
   hold on
 plot(Vsteps,INACTdata,'o')
end  

function [INACTdata,time]=InActCurveExt
    prot=2;
    
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
    %output=zeros(1,length(Vsteps));
    INACTdata=zeros(1,length(Vsteps));
    tic
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    t=[0 5000];
        
    [sol]=voltageclamp2(t,Vstep,prot);
    endtime=6000;
    [solext]=vclampextend(t,Vstep,prot,sol,endtime);
    
    time=solext.x;
    yvals=solext.y;
     k = find(time>=5000,1);
     INACTdata(i)=max(yvals(5,k:end));
    figure(3)
    hold on
    plot(time,yvals(5,:))
 end
 toc
  figure(4)
   hold on
 plot(Vsteps,INACTdata,'o')
  xlabel('Voltage (mV)')
  ylabel('Ratio P2 to P1')
end  

function [Recover,time]=InacRecovCurve
    prot=3;
    %t=[0:0.1:500];
    Vsteps=[.1,.5,1,5,10,20,50,100,500,1000,2000,5000,7000,10000,20000];
    %Vsteps=[.1];
     Recover=zeros(1,length(Vsteps));
     tic
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
        t=[0:.1: 6000+Vstep];
    [time,S]=voltageclamp(t,Vstep,prot);
    k = find(time>=5200,1);
    Recover(i)=max(S(k:end,5));
    figure(5)
    hold on
    plot(time,S(:,5))
 end
 toc
  figure(6)
   plot(Vsteps,Recover,'o')
end  

function [output,ACTdata,time]=ActCurveKv14
    prot=4;
    t=[0 50];
    Vsteps=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
    output=zeros(1,length(Vsteps));
    ACTdata=zeros(1,length(Vsteps));
    Taus=zeros(1,length(Vsteps));
    Ninfs=zeros(1,length(Vsteps));  
    tic
 for i=1:length(Vsteps)
        Vstep=Vsteps(i);
    [time,S]=voltageclamp(t,Vstep,prot);
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
 toc
 figure(1)
   hold on
    plot(Vsteps,Taus,'o')
    xlabel('Voltage (mV)')
    ylabel('Tau Activation (ms)')
    
   figure(2)
    hold on
 plot(Vsteps,ACTdata,'o')
     xlabel('Voltage (mV)')
    ylabel('Open Probability')
end  



function [time,S]=voltageclamp(t,Vstep,prot)

S0 = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]'; 
options = odeset('RelTol',1e-6,'AbsTol',1e-6); 
[time,S]=ode15s(@(t,S)deRHS(t,S,Vstep,prot),t,S0,options);
S(end,:);

end

function [sol]=voltageclamp2(t,Vstep,prot)

S0 = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
options = odeset('RelTol',1e-6,'AbsTol',1e-6); 
[sol]=ode15s(@(t,S)deRHS(t,S,Vstep,prot),t,S0,options);


end

function [solext]=vclampextend(t,Vstep,prot,sol,endtime)

solext=odextend(sol,@(t,S)deRHS(t,S,Vstep,prot),endtime);

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
    elseif t>20 && t<50
        V=Vstep;
    elseif t>50 && t<80
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


%parms=[452.4278    1.3168    3.2838    7.5995    0.1221    0.5359];
%parms=[424.9,  2.12,  8.516,   5.36,    .099,  .453];
parms=[429.42,  2.15,  8.56,   5.36,    .102,  .46];
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