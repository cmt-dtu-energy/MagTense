clearvars
close all

%--- Use either field 1 or field 2 from the NIST example
NIST_field = 1;
resolution = [36,9,1];

figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on; box on
mu0 = 4*pi*1e-7;

addpath('../../MEX_files');
addpath('../util');

%% Setup the problem for the initial configuration
%takes the size of the grid as arguments (nx,ny,nz) and a function handle
%that produces the desired field (if not present zero applied field is inferred)
tic
problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));

problem.setUseCuda( true );

problem.alpha = 4.42e3;
problem.gamma = 0;

%initial magnetization
problem.m0(:) = 1/sqrt(3);

%time grid on which to solve the problem
problem = problem.setTime( linspace(0,100e-9,200) );
problem.setTimeDis = int32(10);
HystDir = 1/mu0*[1,1,1] ;

%time-dependent applied field
HextFct = @(t) (1e-9-t)' .* HystDir .* (t<1e-9)';

problem = problem.setHext( HextFct );

solution = struct();
%convert the class obj to a struct so it can be loaded into fortran
prob_struct = struct(problem);


solution = MagTenseLandauLifshitzSolver_mex( prob_struct, solution );
figure; M_end = squeeze(solution.M(end,:,:)); quiver(solution.pts(:,1),solution.pts(:,2),M_end(:,1),M_end(:,2)); axis equal; title('Starting state - Fortran')


%% Setup problem for the time-dependent solver
problem = DefaultMicroMagProblem(resolution(1),resolution(2),resolution(3));
% problem.gamma = -2.21e5;
% problem.alpha = -4.42e3/problem.Ms;
problem.alpha = 4.42e3 ;
problem.gamma = 2.21e5 ;
problem.dem_thres = 0;%1e-6;
problem = problem.setTime( linspace(0,1e-9,200) ); %
problem.setTimeDis = int32(10);

if (NIST_field == 1)
    %field 1
    HystDir = 1/mu0*[-24.6,4.3,0]/1000 ;
end
if (NIST_field == 2)
    %field 2
    HystDir = 1/mu0*[-35.5,-6.3,0]/1000 ;
end

HextFct = @(t) (t>-1)' .*HystDir;
problem = problem.setHext( HextFct );

problem.m0(:) = solution.M(end,:,:);

%convert the class obj to a struct so it can be loaded into fortran
solution_t = struct();
prob_struct = struct(problem);

solution_t = MagTenseLandauLifshitzSolver_mex( prob_struct, solution_t );
toc

plot(fig1,problem.t,mean(solution_t.M(:,:,1),2),'rx'); 
plot(fig1,problem.t,mean(solution_t.M(:,:,2),2),'gx'); 
plot(fig1,problem.t,mean(solution_t.M(:,:,3),2),'bx'); 
% figure; hold all; for i=2:4; plot(problem.Hext(:,1),problem.Hext(:,i),'.'); end;

%% Run the Matlab version of the micromagnetism code
addpath('..\..\micromagnetism')
tic
Matlab_model_params.nGrid = resolution';
Matlab_model_params.Field_dir = HystDir;
Script_3D_Std_Problem_4(fig1,Matlab_model_params);
toc

%% Compare with published solutions available from NIST webpage
t=linspace(0,1,1000);
data = load(['Published_solutions_field' num2str(NIST_field)]);
mean_avg=data.mean_avg; err=data.err;
colours = [[1 0 0];[0 1 0];[0 0 1]];
weak_colours = colours + ~colours*0.75;
fill_ts=[t,fliplr(t)];  
for j=1:3
    std_errors{j}(1:2,:)=[mean_avg(1,:,j)+err(1,:,j);mean_avg(1,:,j)-err(1,:,j)];
    interval = [std_errors{j}(1,:),fliplr(std_errors{j}(2,:))];
    fill(fig1,1e-9*fill_ts,interval,weak_colours(j,:),'linestyle','none')
    plot(fig1,1e-9*t,mean_avg(1,:,j),'color',colours(j,:))
end

legend(fig1,'Fortran Mx','Fortran My','Fortran Mz','Matlab Mx','Matlab My','Matlab Mz','NIST \sigma{}(Mx)','NIST <Mx>','NIST \sigma{}(My)','NIST <My>','NIST \sigma{}(Mz)','NIST <Mz>');
ylabel(fig1,'<M_i>/M_s')
xlabel(fig1,'Time [ns]')

xlim(fig1,[0 1e-9])
figure(figure1)


