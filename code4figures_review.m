%% this is the script to make the additional figure for the first review of the paper

load('all_results_11aug.mat')
addpath(genpath('/home/np451/Desktop/dataset alignment'));

%%


%% polarisation vs cilia density for small shear stress

figure(10)
ms=8 ; lw=2
fa=[0.052,0.085]
x_15= shear_tap;
y_15= POL(ind_tap); ey_15= emPOL(ind_tap);
 x= shear_tap_div;
y= POL(ind_tap_div); ey= emPOL(ind_tap_div);
shear_ind=[1,4,6]
for ii=1:3
jj=shear_ind(ii)    
errorbar(fa,[y(jj),y_15(jj)],[ey(jj),ey_15(jj)],'d','MarkerSize',ms,'LineWidth',lw);hold on;
end

set(gca,'XScale', 'lin', 'YScale', 'lin');
xlabel('Ciliated fraction area','Interpreter','latex') ;
ylabel('Alignment $\Phi$ ','Interpreter','latex');
 ylim([-0.5,1.4]); xlim([fa(1)-fa(1)*0.2,fa(2)+fa(2)*0.2])
 x0=0;y0=0;width=600;height=300;
set(gcf,'position',[x0,y0,width,height])
set(gca,'FontSize',13);
legend({strcat('shear=',num2str(x(shear_ind(1)),1),' dyne/cm$^2$'),strcat('shear=',num2str(x(shear_ind(2)),1),' dyne/cm$^2$'),strcat('shear=',num2str(x(shear_ind(3)),2),' dyne/cm$^2$')},'Interpreter','latex','Location','sw')


