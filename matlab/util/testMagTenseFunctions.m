util_dir = pwd;

%% Test magnetostatic part by running validation scripts and comparing with FEM simulations
cd(util_dir)
cd('..\examples\Magnetostatics\Validation_field_cylindrical_slice')
[rel_int_error_cylindrical_slice_example_1] = MagTense_Validation_cylindrical_slice_example_1();
assert(all(rel_int_error_cylindrical_slice_example_1 < 5))

cd(util_dir)
cd('..\examples\Magnetostatics\Validation_field_cylindrical_slice')
[rel_int_error_cylindrical_slice_example_2] = MagTense_Validation_cylindrical_slice_example_2();
assert(all(rel_int_error_cylindrical_slice_example_2 < 5))

cd(util_dir)
cd('..\examples\Magnetostatics\Validation_field_prism')
[rel_int_error_prism] = MagTense_Validation_prism();
assert(all(rel_int_error_prism < 5))

cd(util_dir)
cd('..\examples\Magnetostatics\Validation_field_sphere')
[rel_int_error_sphere] = MagTense_Validation_sphere();
assert(all(rel_int_error_sphere < 5))

cd(util_dir)
cd('..\examples\Magnetostatics\Validation_field_spheroid')
[rel_int_error_spheroid] = MagTense_Validation_spheroid();
assert(all(rel_int_error_spheroid < 5))

cd(util_dir)
cd('..\examples\Magnetostatics\Validation_field_tetrahedron')
[rel_int_error_tetrahedron] = MagTense_Validation_tetrahedron();
assert(all(rel_int_error_tetrahedron < 5))


%% Test micromagnetic part by running the mumag standard problems and comparing with published solutions

%% Standard problem 4
cd(util_dir)
cd('..\examples\Micromagnetism\mumag_micromag_Std_problem_4')
[~,~,~,~,~,~,rel_int_error_Standard_problem_4] = Standard_problem_4( 1 );
assert(all(rel_int_error_Standard_problem_4 < 100))


%% Standard problem 6
cd(util_dir)
cd('..\examples\Micromagnetism\mumag_micromag_Std_problem_6')
parameter_variation = {'akj','ak','aj','a','kj','k'};
HPs=containers.Map(parameter_variation,[1.568,1.089,1.206,0.838,1.005,0.565]); % Theoretical pinning fields for the different cases.

for i = 1:length(parameter_variation)
    %% AKJ parameters
    HP = Standard_problem_6(parameter_variation{i}, 80, 201, 'x','ShowTheResult',false); 
    HP_theory = HPs(parameter_variation{i});
    rel_int_error_Standard_problem_6_AKJ(i) = abs(HP-HP_theory)/HP_theory*100;
    disp(rel_int_error_Standard_problem_6_AKJ(i))
    assert(rel_int_error_Standard_problem_6_AKJ(i) < 5)
end

directions = {'x','y','z'};
for i = 1:length(directions)
    %% Physical directions
    HP = Standard_problem_6('akj', 80, 201, directions{i}); 
    HP_theory = HPs('akj');
    rel_int_error_Standard_problem_6_dir(i) = abs(HP-HP_theory)/HP_theory*100;
    disp(rel_int_error_Standard_problem_6_dir(i))
    assert(rel_int_error_Standard_problem_6_dir(i) < 5)
end

%% Unstructured mesh (only works in the x-direction)
Standard_problem_6('akj', 80, 201, 'x', 'use_uniform_mesh', false); 
HP_theory = HPs('akj');
rel_int_error_Standard_problem_6_mesh = abs(HP-HP_theory)/HP_theory*100;
disp(rel_int_error_Standard_problem_6_mesh)
assert(rel_int_error_Standard_problem_6_mesh < 5)


%% Display the errors
disp("Relative errors for all tests [in %]")
disp(rel_int_error_cylindrical_slice_example_1)
disp(rel_int_error_cylindrical_slice_example_2)
disp(rel_int_error_prism)
disp(rel_int_error_sphere)
disp(rel_int_error_spheroid)
disp(rel_int_error_tetrahedron)
disp(rel_int_error_Standard_problem_4)
disp(rel_int_error_Standard_problem_6_AKJ)
disp(rel_int_error_Standard_problem_6_dir)
disp(rel_int_error_Standard_problem_6_mesh)