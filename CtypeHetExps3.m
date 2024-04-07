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
    xlabel('Voltage (mV)')
   ylabel('Ratio P2 to P1')
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
 
   figure(2)
    hold on
 plot(Vsteps,ACTdata,'o')
 
end  



function [time,S]=voltageclamp(t,Vstep,prot)

S0 = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]'; 
options = odeset('RelTol',1e-6,'AbsTol',1e-6); 
[time,S]=ode15s(@(t,S)deRHS(t,S,Vstep,prot),t,S0,options);
S(end,:);

end

function [sol]=voltageclamp2(t,Vstep,prot)

S0 = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
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



%%%for a 1 kv11: 3kv14Ctype het channel
c000=S(1);
c001=S(2);
c101=S(3);
c201=S(4);
c301=S(5);
c111=S(6);
c211=S(7);
c311=S(8);
c221=S(9);
c321=S(10);
c331=S(11);
c100=S(12);
c200=S(13);
c300=S(14);
c110=S(15);
c210=S(16);
c310=S(17);
c220=S(18);
c320=S(19);
c330=S(20);


%Best Fit Parms so far
%parms=[452.4278    1.3168    3.2838    7.5995    0.1221    0.5359];
%parms=[424.9,  2.12,  8.516,   5.36,    .099,  .453];
parms=[429.42,  2.15,  8.56,   5.36,    .102,  .46];
m1=parms(1);
n1=parms(2);
c1=parms(3);
d1=parms(4);
m2=parms(5);
n2=parms(6);

%Rates Kv14
 a1=m1*exp(n1*V*.01)/1000;
 b1=c1*exp(-d1*V*.01)/1000;
 aI=m2/1000;
 bI=n2/1000;
 %Rates Kv11
 a2=.12889*exp((V+45)/33.90877);
 b2=.12889*exp(-(V+45)/12.42101);
 
%%1Kv1.1 and 3Kv1.4
 s_prime=[-3*a1*c000-a2*c000+b1*c100+b2*c001,...
          a2*c000-3*a1*c001-b2*c001+b1*c101,...
          3*a1*c001-2*a1*c101+a2*c100-aI*c101+bI*c111+2*b1*c201-b1*c101-b2*c101,...
          2*a1*c101-a1*c201-2*aI*c201+a2*c200-2*b1*c201+3*b1*c301+bI*c211-b2*c201,...
          a1*c201-3*aI*c301+a2*c300-b2*c301-3*b1*c301+bI*c311,...
          -bI*c111+aI*c101-2*a1*c111+b1*c211-b2*c111+a2*c110,...
          2*a1*c111-a1*c211-aI*c211+2*aI*c201+a2*c210-b1*c211+2*b1*c311-bI*c211+2*bI*c221-b2*c211,...
          a1*c211+3*aI*c301-2*aI*c311+a2*c310-2*b1*c311-bI*c311+2*bI*c321-b2*c311,...
          aI*c211-a1*c221+a2*c220-2*bI*c221+b1*c321-b2*c221,...
          a1*c221+2*aI*c311-aI*c321+a2*c320-b1*c321-2*bI*c321+3*bI*c331-b2*c321,...
          aI*c321+a2*c330-3*bI*c331-b2*c331,...
          3*a1*c000-2*a1*c100-aI*c100-a2*c100-b1*c100+2*b1*c200+bI*c110+b2*c101,...
          2*a1*c100-a1*c200-2*aI*c200-a2*c200-2*b1*c200+3*b1*c300+bI*c210+b2*c201,...
          a1*c200-3*aI*c300-a2*c300-3*b1*c300+b2*c301+bI*c310,...
          -2*a1*c110+aI*c100-a2*c110+b1*c210-bI*c110+b2*c111,...
          2*a1*c110-a1*c210+2*aI*c200-aI*c210-a2*c210-b1*c210+2*b1*c310-bI*c210+2*bI*c220+b2*c211,...
          a1*c210+3*aI*c300-2*aI*c310-a2*c310-2*b1*c310-bI*c310+2*bI*c320+b2*c311,...
          aI*c210-a1*c220-a2*c220+b1*c320-2*bI*c220+b2*c221,...
          a1*c220+2*aI*c310-aI*c320-a2*c320-b1*c320-2*bI*c320+3*bI*c330+b2*c321,...
          aI*c320-a2*c330-3*bI*c330+b2*c331]';


end