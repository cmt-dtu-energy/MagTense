function model = getGeneralComsolCompare( model, tile )

box_size = 10;

model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('pos', [tile.offset(1) tile.offset(2) tile.offset(3)]);
model.component('comp1').geom('geom1').feature('blk1').set('base', 'center');
largest_abc = max(tile.abc);
model.component('comp1').geom('geom1').feature('blk1').set('size', [largest_abc largest_abc largest_abc]*box_size);

model.component('comp1').geom('geom1').run;

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});

model.component('comp1').physics.create('mfnc', 'MagnetostaticsNoCurrents', 'geom1');
model.component('comp1').physics('mfnc').create('mfc2', 'MagneticFluxConservation', 3);
model.component('comp1').physics('mfnc').feature('mfc2').selection.set([2]);
model.component('comp1').physics('mfnc').feature('mfc2').set('ConstitutiveRelationBH', 'RemanentFluxDensity');
model.component('comp1').physics('mfnc').feature('mfc2').set('e_crel_BH_RemanentFluxDensity', [tile.u_ea(1); tile.u_ea(2); tile.u_ea(3)]);
model.component('comp1').physics('mfnc').feature('mfc2').set('normBr_crel_BH_RemanentFluxDensity_mat', 'userdef');
model.component('comp1').physics('mfnc').feature('mfc2').set('normBr_crel_BH_RemanentFluxDensity', tile.Mrem*4*pi*1e-7);

model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
% model.component('comp1').mesh('mesh1').create('ref1', 'Refine');
% model.component('comp1').mesh('mesh1').feature('ref1').selection.geom('geom1', 3);
% model.component('comp1').mesh('mesh1').feature('ref1').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('size').set('hauto', 3);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'cg');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'amg');
model.sol('sol1').runAll;

model.result.dataset.create('cln3', 'CutLine3D');
model.result.dataset.create('cln2', 'CutLine3D');
model.result.dataset.create('cln1', 'CutLine3D');
model.result.dataset.create('surf1', 'Surface');
model.result.dataset('cln1').set('genpoints', [-largest_abc*box_size tile.offset(2) tile.offset(3); largest_abc*box_size tile.offset(2) tile.offset(3)]);
model.result.dataset('cln2').set('genpoints', [tile.offset(1) -largest_abc*box_size tile.offset(3); tile.offset(1) largest_abc*box_size tile.offset(3)]);
model.result.dataset('cln3').set('genpoints', [tile.offset(1) tile.offset(2) -largest_abc*box_size; tile.offset(1) tile.offset(2) largest_abc*box_size]);

model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('expr', '1');
model.result.export.create('plot10', 'Plot');
model.result.export('plot10').set('filename', 'Comsol_spheroid_surface.txt');

model.result.create('pg2', 'PlotGroup1D');
model.result('pg2').create('lngr1', 'LineGraph');
model.result('pg2').create('lngr2', 'LineGraph');
model.result('pg2').create('lngr3', 'LineGraph');
model.result('pg2').create('lngr4', 'LineGraph');
model.result('pg2').create('lngr5', 'LineGraph');
model.result('pg2').create('lngr6', 'LineGraph');
model.result('pg2').create('lngr7', 'LineGraph');
model.result('pg2').create('lngr8', 'LineGraph');
model.result('pg2').create('lngr9', 'LineGraph');
model.result('pg2').feature('lngr1').set('xdata', 'expr');
model.result('pg2').feature('lngr1').set('expr', 'mfnc.Hx');
model.result('pg2').feature('lngr2').set('xdata', 'expr');
model.result('pg2').feature('lngr2').set('expr', 'mfnc.Hy');
model.result('pg2').feature('lngr3').set('xdata', 'expr');
model.result('pg2').feature('lngr3').set('expr', 'mfnc.Hz');
model.result('pg2').feature('lngr4').set('xdata', 'expr');
model.result('pg2').feature('lngr4').set('expr', 'mfnc.Hx');
model.result('pg2').feature('lngr5').set('xdata', 'expr');
model.result('pg2').feature('lngr5').set('expr', 'mfnc.Hy');
model.result('pg2').feature('lngr6').set('xdata', 'expr');
model.result('pg2').feature('lngr6').set('expr', 'mfnc.Hz');
model.result('pg2').feature('lngr7').set('xdata', 'expr');
model.result('pg2').feature('lngr7').set('expr', 'mfnc.Hx');
model.result('pg2').feature('lngr8').set('xdata', 'expr');
model.result('pg2').feature('lngr8').set('expr', 'mfnc.Hy');
model.result('pg2').feature('lngr9').set('xdata', 'expr');
model.result('pg2').feature('lngr9').set('expr', 'mfnc.Hz');
model.result.export.create('plot1', 'Plot');
model.result.export.create('plot2', 'Plot');
model.result.export.create('plot3', 'Plot');
model.result.export.create('plot4', 'Plot');
model.result.export.create('plot5', 'Plot');
model.result.export.create('plot6', 'Plot');
model.result.export.create('plot7', 'Plot');
model.result.export.create('plot8', 'Plot');
model.result.export.create('plot9', 'Plot');


model.result('pg2').set('xlabel', 'z-coordinate (m)');
model.result('pg2').set('ylabel', 'Magnetic field, z component (A/m)');
model.result('pg2').set('xlabelactive', false);
model.result('pg2').set('ylabelactive', false);
model.result('pg2').feature('lngr1').set('data', 'cln1');
model.result('pg2').feature('lngr1').set('xdataexpr', 'x');
model.result('pg2').feature('lngr2').set('data', 'cln1');
model.result('pg2').feature('lngr2').set('xdataexpr', 'x');
model.result('pg2').feature('lngr3').set('data', 'cln1');
model.result('pg2').feature('lngr3').set('xdataexpr', 'x');
model.result('pg2').feature('lngr4').set('data', 'cln2');
model.result('pg2').feature('lngr4').set('xdataexpr', 'y');
model.result('pg2').feature('lngr5').set('data', 'cln2');
model.result('pg2').feature('lngr5').set('xdataexpr', 'y');
model.result('pg2').feature('lngr6').set('data', 'cln2');
model.result('pg2').feature('lngr6').set('xdataexpr', 'y');
model.result('pg2').feature('lngr7').set('data', 'cln3');
model.result('pg2').feature('lngr7').set('xdataexpr', 'z');
model.result('pg2').feature('lngr8').set('data', 'cln3');
model.result('pg2').feature('lngr8').set('xdataexpr', 'z');
model.result('pg2').feature('lngr9').set('data', 'cln3');
model.result('pg2').feature('lngr9').set('xdataexpr', 'z');

model.result.export('plot1').set('plotgroup', 'pg2');
model.result.export('plot1').set('plot', 'lngr1');
model.result.export('plot1').set('filename', 'Comsol_spheroid_Hx_x.txt');
model.result.export('plot1').run;

model.result.export('plot2').set('plotgroup', 'pg2');
model.result.export('plot2').set('plot', 'lngr2');
model.result.export('plot2').set('filename', 'Comsol_spheroid_Hy_x.txt');
model.result.export('plot2').run;

model.result.export('plot3').set('plotgroup', 'pg2');
model.result.export('plot3').set('plot', 'lngr3');
model.result.export('plot3').set('filename', 'Comsol_spheroid_Hz_x.txt');
model.result.export('plot3').run;

model.result.export('plot4').set('plotgroup', 'pg2');
model.result.export('plot4').set('plot', 'lngr4');
model.result.export('plot4').set('filename', 'Comsol_spheroid_Hx_y.txt');
model.result.export('plot4').run;

model.result.export('plot5').set('plotgroup', 'pg2');
model.result.export('plot5').set('plot', 'lngr5');
model.result.export('plot5').set('filename', 'Comsol_spheroid_Hy_y.txt');
model.result.export('plot5').run;

model.result.export('plot6').set('plotgroup', 'pg2');
model.result.export('plot6').set('plot', 'lngr6');
model.result.export('plot6').set('filename', 'Comsol_spheroid_Hz_y.txt');
model.result.export('plot6').run;

model.result.export('plot7').set('plotgroup', 'pg2');
model.result.export('plot7').set('plot', 'lngr7');
model.result.export('plot7').set('filename', 'Comsol_spheroid_Hx_z.txt');
model.result.export('plot7').run;

model.result.export('plot8').set('plotgroup', 'pg2');
model.result.export('plot8').set('plot', 'lngr8');
model.result.export('plot8').set('filename', 'Comsol_spheroid_Hy_z.txt');
model.result.export('plot8').run;

model.result.export('plot9').set('plotgroup', 'pg2');
model.result.export('plot9').set('plot', 'lngr9');
model.result.export('plot9').set('filename', 'Comsol_spheroid_Hz_z.txt');
model.result.export('plot9').run;

model.component('comp1').view('view1').hideObjects.create('hide1');
model.component('comp1').view('view1').hideObjects('hide1').init;
model.component('comp1').view('view1').hideObjects('hide1').add({'blk1'});
model.result('pg1').run;
model.result.export('plot10').run;

end
