
clc;
coff_X=[-10.0000;-15.0000;-20.0000;-25.0000;-30.0000;-35.0000;-40.0000;-45.0000;-50.0000];
plot(coff_X,Chebyshev_coff_2Y,'k-*');
hold on;
% grid on;
plot(coff_X,Chebyshev_coff_4Y,'k:*');
hold on;
% grid on;
plot(coff_X,SOCP_coff_2y,'k-o');
hold on;
% grid on;
plot(coff_X,SOCP_coff_4y,'k:o');
hold on;
% grid on;
h = gca;
set(h,'FontSize',10,'FontName','Times New Roman');
set(h,'FontName','Times New Roman');
xlabel('Sidelobe Level(dB)');
ylabel('Mainlobe Width(Degree)');
title('Relationship Between Mainlobe Width And Sidelobe Level');
legend('Dolph-chebyshev  (1/2 lamda)','Dolph-chebyshev  (1/4 lamda)','SOCP  (1/2 lamda)','SOCP  (1/4 lamda)');
