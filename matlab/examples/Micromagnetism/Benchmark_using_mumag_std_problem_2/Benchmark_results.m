clearvars
close all

i_arr = 1:9; 

figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
fig1 = axes('Parent',figure1,'Layer','top','FontSize',10);
hold all
grid on
box on
set(fig1,'Yscale','log')

figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
fig2 = axes('Parent',figure2,'Layer','top','FontSize',10);
hold all
grid on
box on

figure3= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
fig3 = axes('Parent',figure3,'Layer','top','FontSize',10);
hold all
grid on
box on

figure4= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
fig4 = axes('Parent',figure4,'Layer','top','FontSize',10);
hold all
grid on
box on

%% Plot results of the MagTense benchmarking
linestr = {'-','-.','--',':'};
clear elapsedTime; 
k = 1;

Color_arr = parula(length(i_arr));

for m = 2%[1 2]
    for i = 1:length(i_arr)
        resolution = [i*20,i*5,1];
        switch m
            case 1
                load(['CPU\Solution_' num2str(resolution(1)) '_' num2str(resolution(2)) '_' num2str(resolution(3)) '_no_CUDA.mat']); 
                plot_sym = 'o';
                work_str = 'CUDA';
            case 2
                load(['CUDA\Solution_' num2str(resolution(1)) '_' num2str(resolution(2)) '_' num2str(resolution(3)) '_with_CUDA.mat']); 
                plot_sym = 'd';
                work_str = 'CPU';
        end

%         if ((m == 2) && ((i == 3) || (i == 5) || (i == 7) || (i == 9)))
%             plot(fig6,problem_dym.t,mean(solution_dym.M(:,:,1),2),'k','Linestyle',linestr{k}); 
%             plot(fig6,problem_dym.t,mean(solution_dym.M(:,:,2),2),'k','Linestyle',linestr{k}); 
%             plot(fig6,problem_dym.t,mean(solution_dym.M(:,:,3),2),'k','Linestyle',linestr{k});
%             xlim(fig6,[min(problem_dym.t) max(problem_dym.t)])
%             k = k+1;
%         end
% 
%         %--- Calculate the difference from the NIST published results
%         for j = 1:3
%             NIST_avg_interp = interp1(1e-9*t,mean_avg(1,:,j),problem_dym.t);
%             abs_res_diff(i,j) = trapz(problem_dym.t,abs(mean(solution_dym.M(:,:,j),2)-NIST_avg_interp'));
%         end

        elapsedTime_arr(i) = [elapsedTime]; 
        
        plot(fig2,results.dlex,results.Mxr,'Marker',plot_sym,'Color',Color_arr(i,:));
        plot(fig3,results.dlex,results.Myr,'Marker',plot_sym,'Color',Color_arr(i,:));
        plot(fig4,results.dlex,abs(results.Hc),'Marker',plot_sym,'Color',Color_arr(i,:)); 
        legend_arr{k} = [work_str '-' num2str(resolution(1)) 'x' num2str(resolution(2)) 'x' num2str(resolution(3))];
        k = k+1;
    end
  
    switch m
        case 1
            elapsedTime_CPU = elapsedTime_arr;
            plot(fig1,i_arr,elapsedTime_arr,'r-s');
        case 2
            elapsedTime_CUDA = elapsedTime_arr;
            plot(fig1,i_arr,elapsedTime_arr,'b-s');  
    end
    

end



%% Legend, titles and labels
legend(fig1,'CPU','CUDA','Location','NorthWest')
title(fig1,['Standard problem 2'])
xlabel(fig1,'n with resolution=[n*20,n*5,1]')
ylabel(fig1,'Simulation time [s]')
% figure(figure1)
% print('-dpng',['Field_' num2str(NIST_field) '_Benchmark_problem_4_time_res.png'])

xlabel(fig2,'d/l_{ex}'); 
ylabel(fig2,'M_{xr}/M_s');
legend(fig2,legend_arr)

xlabel(fig3,'d/l_{ex}'); 
ylabel(fig3,'M_{yr}/M_s');
legend(fig3,legend_arr)
  
xlabel(fig4,'d/l_{ex}'); 
ylabel(fig4,'|H_c|/M_s');
legend(fig4,legend_arr)