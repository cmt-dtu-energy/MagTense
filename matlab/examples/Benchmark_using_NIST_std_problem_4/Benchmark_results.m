clearvars
close all

NIST_field = 2;
i_arr = 1:9; 

addpath('../MagTense/matlab/MEX_files');
addpath('../MagTense/matlab/util');
        
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

figure6= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
fig6 = axes('Parent',figure6,'Layer','top','FontSize',10);
hold all
grid on
box on

%--- Fake legend for figure 6
plot(fig6,-1,0,'k-')
plot(fig6,-1,0,'k-.')
plot(fig6,-1,0,'k--')
plot(fig6,-1,0,'k:')

%% Compare with published solutions available from NIST webpage
t=linspace(0,1,1000);
data = load(['Published_solutions_field' num2str(NIST_field)]);
mean_avg=data.mean_avg; 
err=data.err;
colours = [[1 0 0];[0 1 0];[0 0 1]];
weak_colours = colours + ~colours*0.75;
fill_ts=[t,fliplr(t)];  
for j=1:3
    std_errors{j}(1:2,:)=[mean_avg(1,:,j)+err(1,:,j);mean_avg(1,:,j)-err(1,:,j)];
    interval = [std_errors{j}(1,:),fliplr(std_errors{j}(2,:))];
    fill(fig6,1e-9*fill_ts,interval,weak_colours(j,:),'linestyle','none')
%     plot(fig6,1e-9*t,mean_avg(1,:,j),'color',colours(j,:))
end

%% Plot results of the MagTense benchmarking
linestr = {'-','-.','--',':'};
clear elapsedTime; 
k = 1;

for m = [1 2]
    for i = 1:length(i_arr)
        switch m
            case 1
                load(['Field_' num2str(NIST_field) '\CPU\Solution_' num2str(i) '_x_20_5_1_no_CUDA.mat']); 
            case 2
                load(['Field_' num2str(NIST_field) '\CUDA\Solution_' num2str(i) '_x_20_5_1_with_CUDA.mat']); 
         end

        if ((m == 2) && ((i == 3) || (i == 5) || (i == 7) || (i == 9)))
            plot(fig6,problem_t.t,mean(solution_t.M(:,:,1),2),'k','Linestyle',linestr{k}); 
            plot(fig6,problem_t.t,mean(solution_t.M(:,:,2),2),'k','Linestyle',linestr{k}); 
            plot(fig6,problem_t.t,mean(solution_t.M(:,:,3),2),'k','Linestyle',linestr{k});
            xlim(fig6,[min(problem_t.t) max(problem_t.t)])
            k = k+1;
        end

        %--- Calculate the difference from the NIST published results
        for j = 1:3
            NIST_avg_interp = interp1(1e-9*t,mean_avg(1,:,j),problem_t.t);
            abs_res_diff(i,j) = trapz(problem_t.t,abs(mean(solution_t.M(:,:,j),2)-NIST_avg_interp'));
        end

        elapsedTime(i,:) = [elapsedTime_part1,elapsedTime_part2]; 
    end
  
    switch m
        case 1
            elapsedTime_CPU = elapsedTime;
            plot(fig1,i_arr,elapsedTime(:,2),'r-s'); 
        case 2
            elapsedTime_CUDA = elapsedTime;
 
            plot(fig1,i_arr,elapsedTime(:,2),'b-s'); 
            plot(fig2,i_arr,sum(abs_res_diff(:,:),2),'k-x');
    end
end

plot(fig3, i_arr,elapsedTime_CPU(:,2)./elapsedTime_CUDA(:,2),'k-d');

%% Load and plot the Matlab simulation data
clear abs_res_diff
for m = 1:2
    for i = 1:length(i_arr)
        switch m 
            case 1
                load(['Field_' num2str(NIST_field) '\Matlab_simulations_CPU\Matlab_resolution_' num2str(i) '_x_20_5_1\PhysParams_TestDynamicsStdProb4_1st_fieldsolution.mat'])
                load(['Field_' num2str(NIST_field) '\Matlab_simulations_CPU\Matlab_resolution_' num2str(i) '_x_20_5_1\PhysParams_TestDynamicsStdProb4_1st_field.mat'])
            case 2
                load(['Field_' num2str(NIST_field) '\Matlab_simulations_GPU\Matlab_resolution_' num2str(i) '_x_20_5_1\PhysParams_TestDynamicsStdProb4_1st_fieldsolution.mat'])
                load(['Field_' num2str(NIST_field) '\Matlab_simulations_GPU\Matlab_resolution_' num2str(i) '_x_20_5_1\PhysParams_TestDynamicsStdProb4_1st_field.mat'])
        end
        
        for k=1:size(SigmaSol,1) 
            Sigma = SigmaSol(k,:).' ;
            NN = round(numel(Sigma)/3) ;

            SigmaX = Sigma(0*NN+[1:NN]) ;
            SigmaY = Sigma(1*NN+[1:NN]) ;
            SigmaZ = Sigma(2*NN+[1:NN]) ;
            SigmaN = sqrt(SigmaX.^2+SigmaY.^2+SigmaZ.^2) ;
            Mx(k) = mean(SigmaX./SigmaN) ;
            My(k) = mean(SigmaY./SigmaN) ;
            Mz(k) = mean(SigmaZ./SigmaN) ;
            Mk(k) = Mx(k)*ProblemSetupStruct.HystDir(1) + My(k)*ProblemSetupStruct.HystDir(2) + Mz(k)*ProblemSetupStruct.HystDir(3) ;
        end
        M_matlab = [Mx; My; Mz];

        for j = 1:3
            NIST_avg_interp = interp1(1e-9*t,mean_avg(1,:,j),ProblemSetupStruct.t);
            abs_res_diff(i,j) = trapz(ProblemSetupStruct.t,abs(M_matlab(j,:)-NIST_avg_interp));
        end    

        %--- This is the elapsed time of part 2, i.e. the dynamic part.
        elapsedTime_Matlab(i) = elapsedTime;

    %     plot(fig6, ProblemSetupStruct.t,Mx,'ro')
    %     plot(fig6, ProblemSetupStruct.t,My,'go')
    %     plot(fig6, ProblemSetupStruct.t,Mz,'bo')
    end

    switch m
        case 1
            plot(fig2,i_arr,sum(abs_res_diff(:,:),2),'k-d');
            plot(fig1, i_arr,elapsedTime_Matlab,'g-s');
        case 2
            plot(fig1, i_arr,elapsedTime_Matlab,'m-s');
    end
end

%% Legend, titles and labels
legend(fig1,'CPU','CUDA','Matlab','Matlab GPU','Location','NorthWest')
title(fig1,['Standard problem 4, field ' num2str(NIST_field)])
xlabel(fig1,'n with resolution=[n*20,n*5,1]')
ylabel(fig1,'Simulation time dynamic part [s]')
figure(figure1)
print('-dpng',['Field_' num2str(NIST_field) '_Benchmark_problem_4_time_res.png'])

legend(fig2,'Fortran','Matlab','Location','NorthEast')
title(fig2,['Standard problem 4, field ' num2str(NIST_field)])
xlabel(fig2,'n with resolution=[n*20,n*5,1]')
ylabel(fig2,'Error, i.e. \Sigma_i \int{}|m_{i}-m_{i,NIST}|dt')
figure(figure2)
print('-dpng',['Field_' num2str(NIST_field) '_Benchmark_problem_4_error_res.png'])

legend(fig3,'t_{CPU}/t_{CUDA}','t_{CPU}/t_{CUDA, thres frac 0.2}','t_{CUDA}/t_{CUDA, thres frac 0.2}','Location','NorthWest')
title(fig3,['Standard problem 4, field ' num2str(NIST_field)])
xlabel(fig3,'n with resolution=[n*20,n*5,1]')
ylabel(fig3,'Time ratio [-]')
% figure(figure3)
% print('-dpng',['Field_' num2str(NIST_field) '_Benchmark_problem_4_time_rel_res.png'])


title(fig6,['Standard problem 4, field ' num2str(NIST_field)])
xlabel(fig6,'Time [s]')
ylabel(fig6,'m_i [-]')
legend(fig6,'MagTense [60 15 1]','MagTense [100 25 1]','MagTense [140 35 1]','MagTense [180 45 1]', 'NIST \sigma{}(Mx)','NIST \sigma{}(My)','NIST \sigma{}(Mz)');
figure(figure6)
print('-dpng',['Field_' num2str(NIST_field) '_Benchmark_problem_4_t_specific_solutions.png'])