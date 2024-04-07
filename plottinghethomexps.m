Voltage=[-90,-80,-70,-60,-50,-40,-30,-20,-10,-0,10,20,30,40,50];
 Timedelta=[.1,.5,1,5,10,20,50,100,500,1000,2000,5000,7000,10000,20000];
 
 Mat1114allvals=Kv11Kv14HetHomallEXPsSheet1;
 Mat1114QSS=Kv11Kv14HetHomallEXPsQSS1;
% for i=1:15
%    figure(1)
%    hold on
%    plot(Voltage(i),Mat1114allvals(2,i),"bo");
%    plot(Voltage(i),Mat1114allvals(6,i),"ro");
%    plot(Voltage(i),Mat1114allvals(10,i),"go");
%    plot(Voltage(i),Mat1114allvals(14,i),"ko");
%    plot(Voltage(i),Mat1114allvals(18,i),"yo");
%    %plot(Voltage(i),Mat1114allvals(33,i),"ro");
%    hold off
%    legend('kv11hom','3:1','2:2','1:3','hv14Hom')
%    xlabel('Voltage (mV)')
%    ylabel('Open Probability')
%    %legend('14full','14QSS')
% end
% for i=1:15
%    figure(2)
%    hold on
%    plot(Voltage(i),Mat1114allvals(8,i),"ro");
%    plot(Voltage(i),Mat1114allvals(12,i),"go");
%    plot(Voltage(i),Mat1114allvals(16,i),"ko");
%    plot(Voltage(i),Mat1114allvals(20,i),"yo");
%    hold off
%    legend('3:1','2:2','1:3','hv14Hom')
% end
% for i=1:15
%    figure(3)
%    hold on
%    plot(Timedelta(i),Mat1114allvals(7,i),"ro");
%    plot(Timedelta(i),Mat1114allvals(11,i),"go");
%    plot(Timedelta(i),Mat1114allvals(15,i),"ko");
%    plot(Timedelta(i),Mat1114allvals(19,i),"yo");
%    hold off
%    legend('3:1','2:2','1:3','hv14Hom')
%      
%   
% end
% 
% for i=1:15
%    figure(1)
%    hold on
%    plot(Timedelta(i),Mat1114allvals(23,i),"ro");
%    plot(Timedelta(i),Mat1114allvals(31,i),"go");
%    plot(Timedelta(i),Mat1114allvals(27,i),"ko");
%    plot(Timedelta(i),Mat1114allvals(19,i),"yo");
%    
%    plot(Timedelta(i),Mat1114QSS(8,i),"rv");
%    plot(Timedelta(i),Mat1114QSS(6,i),"gv");
%    plot(Timedelta(i),Mat1114QSS(4,i),"kv");
%    plot(Timedelta(i),Mat1114QSS(2,i),"yv");
%    hold off
%    legend('3:1','2:2','1:3','hv14Hom','3:1QSS','2:2QSS','1:3QSS','0:4QSS')
%    %legend('3:1','2:2','1:3','kv14')
%    xlabel('\Delta t (ms)')
%    ylabel('Ratio P2 to P1')
%   
% end
% % % % 
% for i=1:15
%    figure(2)
%    hold on
%    plot(Voltage(i),Mat1114allvals(24,i),"ro");
%    plot(Voltage(i),Mat1114allvals(32,i),"go");
%    plot(Voltage(i),Mat1114allvals(28,i),"ko");
%    plot(Voltage(i),Mat1114allvals(20,i),"yo");
%    
%    plot(Voltage(i),Mat1114QSS(7,i),"rv");
%    plot(Voltage(i),Mat1114QSS(5,i),"gv");
%    plot(Voltage(i),Mat1114QSS(3,i),"kv");
%    plot(Voltage(i),Mat1114QSS(1,i),"yv");
%    hold off
%    legend('3:1','2:2','1:3','kv14','3:1QSS','2:2QSS','1:3QSS','0:4QSS')
%    %legend('3:1','2:2','1:3','kv14')
%    xlabel('Voltage (mV)')
%    ylabel('Ratio P2 to P1')
% end
% % 
% for i=1:15
%    figure(5)
%    hold on
%    plot(Voltage(i),Mat1114QSS(7,i),"ro");
%    plot(Voltage(i),Mat1114QSS(5,i),"go");
%    plot(Voltage(i),Mat1114QSS(3,i),"ko");
%    plot(Voltage(i),Mat1114QSS(1,i),"yo");
%    hold off
%    legend('3:1','2:2','1:3','hv14Hom')
% end
% % 
% % 
% % for i=1:15
% %    figure(4)
% %    hold on
% %    plot(Timedelta(i),Mat1114QSS(8,i),"ro");
% %    plot(Timedelta(i),Mat1114QSS(6,i),"go");
% %    plot(Timedelta(i),Mat1114QSS(4,i),"ko");
% %    plot(Timedelta(i),Mat1114QSS(2,i),"yo");
% %    hold off
% %    legend('3:1','2:2','1:3','hv14Hom')
% %      
% %   
% % end