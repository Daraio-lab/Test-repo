%3D LINEAR-HOMOGENIZATION
%
%This code implements a static homogenization for a 3D Unit Cell, RVE using
%COMSOL Multiphysics
%Running this file with the MatlabLink of Comsol will save an .mph file
%which shows the used geometry.
%The model fits to an orthotropic material, therefore 6 load cases have to
%be considered. The implementation is similiar to the following paper: Steven, G.P. (1997), ?Homogenization of multicomponent orthotropic materials using FEA?,
%Comm. Num. Meth. Engng, Vol. 13, pp. 517-31.
%However it has been extended to 3D.
%Please let me know if you have problems running the file or want to adapt
%it to your geometry.

clear all;
close all;
clc;

filename = 'HomogenizationLinearElasticIsotropicTest.mph';

a = 0.1;
dx = a;
dy = a;
dz = a;

Ax = dx*dz;
Ay = dy*dz;
Az = dx*dy;

%Material properties of linear elastic material

E=210e9;
nu=0.3;

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/Users/skroedel/Desktop/Projects/SiSuperlatticePhononAnalogy/07-05-15SL-Layervariation-simulation/temporary');

model.param.set('a', num2str(a));
model.param.set('deltax', 'a');
model.param.set('deltay', 'a');
model.param.set('deltaz', 'a');

model.modelNode.create('comp1');

model.geom.create('geom1', 3);
model.geom('geom1').feature.create('blk1', 'Block');
model.geom('geom1').feature('blk1').set('size', {'a' 'a' 'a'});
model.geom('geom1').run;

model.material.create('mat1');

model.physics.create('solid', 'SolidMechanics', 'geom1');
model.physics('solid').feature.create('disp1', 'Displacement2', 2);
model.physics('solid').feature('disp1').selection.set([1]);
model.physics('solid').feature.create('disp2', 'Displacement2', 2);
model.physics('solid').feature('disp2').selection.set([2]);
model.physics('solid').feature.create('disp3', 'Displacement2', 2);
model.physics('solid').feature('disp3').selection.set([6]);
model.physics('solid').feature.create('disp4', 'Displacement2', 2);
model.physics('solid').feature('disp4').selection.set([5]);
model.physics('solid').feature.create('disp5', 'Displacement2', 2);
model.physics('solid').feature('disp5').selection.set([4]);
model.physics('solid').feature.create('disp6', 'Displacement2', 2);
model.physics('solid').feature('disp6').selection.set([3]);

model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').feature.create('ftet1', 'FreeTet');
model.mesh('mesh1').feature('size').set('hauto', 5);

model.material('mat1').propertyGroup('def').set('youngsmodulus', num2str(E));
model.material('mat1').propertyGroup('def').set('poissonsratio', num2str(nu));
model.material('mat1').propertyGroup('def').set('density', '7800');

model.mesh('mesh1').run;


model.study.create('std1');
model.study('std1').feature.create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature.create('s1', 'Stationary');
model.sol('sol1').feature('s1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.study('std1').feature('stat').set('initstudyhide', 'on');
model.study('std1').feature('stat').set('initsolhide', 'on');
model.study('std1').feature('stat').set('notstudyhide', 'on');
model.study('std1').feature('stat').set('notsolhide', 'on');

model.batch.create('p1', 'Parametric');
model.batch.create('p2', 'Parametric');
model.batch('p1').feature.create('so1', 'Solutionseq');
model.batch('p2').feature.create('so1', 'Solutionseq');
model.batch('p1').study('std1');
model.batch('p2').study('std1');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').name('Compile Equations: Stationary');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').feature('s1').set('control', 'stat');
model.sol('sol1').runAll;

disp('Runnning 6 load cases.');
tic;

C = zeros(6,6);
mphsave(model,filename);
for i = 1:6 %%%%Define all 6 load cases for the static homogenization
    %%%%%%
    if i == 1 % eps 1
        
        % Set BC
        model.physics('solid').feature('disp1').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp1').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp2').set('Direction', {'0'; '1'; '1'});
        model.physics('solid').feature('disp2').set('U0', {'0'; '0';'0'});
        model.physics('solid').feature('disp3').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp3').set('U0', {'deltax'; '0'; '0'});
        model.physics('solid').feature('disp4').set('Direction', {'0'; '1'; '1'});
        model.physics('solid').feature('disp4').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp5').set('Direction', {'0'; '1'; '1'});
        model.physics('solid').feature('disp5').set('U0', {'0'; '0';'0'});
        model.physics('solid').feature('disp6').set('Direction', {'0'; '1'; '1'});
        model.physics('solid').feature('disp6').set('U0', {'0'; '0';'0'});
        
        % Run Model
        model.study('std1').run;
        
        %Fetch By integration and calculate stress related to the dimension
        %of the UC
        Sxx = (1/Ax)*mphint2(model,'solid.sx','surface','selection',[6]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Sxy = (1/Ax)*mphint2(model,'solid.sxy','surface','selection',[6]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Sxz = (1/Ax)*mphint2(model,'solid.sxz','surface','selection',[6]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Syy = (1/Ay)*mphint2(model,'solid.sy','surface','selection',[5]);%*mphint2(model,'1','surface','selection',[50 52 171 336 338]);
        Syz = (1/Az)*mphint2(model,'solid.syz','surface','selection',[3]);%*mphint2(model,'1','surface','selection',[10 45 160 326 329]);
        Szz = (1/Az)*mphint2(model,'solid.sz','surface','selection',[3]);%*mphint2(model,'1','surface','selection',[10 45 160 326 329]);
        
        C(1,1) = Sxx;
        C(2,1) = Syy;
        C(3,1) = Szz;
        C(4,1) = Sxy;
        C(5,1) = Syz;
        C(6,1) = Sxz;
        disp(C);
    end
    
    if i == 2 % eps 2
        
         % Set BC
        
        model.physics('solid').feature('disp1').set('Direction', {'1'; '0'; '1'});
        model.physics('solid').feature('disp1').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp2').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp2').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp3').set('Direction', {'1'; '0'; '1'});
        model.physics('solid').feature('disp3').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp4').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp4').set('U0', {'0'; 'deltay'; '0'});
        model.physics('solid').feature('disp5').set('Direction', {'1'; '0'; '1'});
        model.physics('solid').feature('disp5').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp6').set('Direction', {'1'; '0'; '1'});
        model.physics('solid').feature('disp6').set('U0', {'0'; '0'; '0'});
        
        % Run Model
        model.study('std1').run;
        
        %Fetch By integration and calculate stress related to the dimension
        %of the UC
        Sxx = (1/Ax)*mphint2(model,'solid.sx','surface','selection',[6]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Sxy = (1/Ax)*mphint2(model,'solid.sxy','surface','selection',[2]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Sxz = (1/Ax)*mphint2(model,'solid.sxz','surface','selection',[3]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Syy = (1/Ay)*mphint2(model,'solid.sy','surface','selection',[5]);%*mphint2(model,'1','surface','selection',[50 52 171 336 338]);
        Syz = (1/Az)*mphint2(model,'solid.syz','surface','selection',[2]);%*mphint2(model,'1','surface','selection',[10 45 160 326 329]);
        Szz = (1/Az)*mphint2(model,'solid.sz','surface','selection',[3]);%*mphint2(model,'1','surface','selection',[10 45 160 326 329]);
        
        C(1,2) = Sxx;
        C(2,2) = Syy;
        C(3,2) = Szz;
        C(4,2) = Sxy;
        C(5,2) = Syz;
        C(6,2) = Sxz;
         disp(C)
    end
        
    if i == 3 % eps 3
        
        model.physics('solid').feature('disp1').set('Direction', {'1'; '1'; '0'});
        model.physics('solid').feature('disp1').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp2').set('Direction', {'1'; '1'; '0'});
        model.physics('solid').feature('disp2').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp3').set('Direction', {'1'; '1'; '0'});
        model.physics('solid').feature('disp3').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp4').set('Direction', {'1'; '1'; '0'});
        model.physics('solid').feature('disp4').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp5').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp5').set('U0', {'0'; '0'; 'deltaz'});
        model.physics('solid').feature('disp6').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp6').set('U0', {'0'; '0'; '0'});
        
        % Run Model
        model.study('std1').run;
        
        %Fetch By integration and calculate stress related to the dimension
        %of the UC
        Sxx = (1/Ax)*mphint2(model,'solid.sx','surface','selection',[1]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Sxy = (1/Ax)*mphint2(model,'solid.sxy','surface','selection',[5]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Sxz = (1/Ax)*mphint2(model,'solid.sxz','surface','selection',[4]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Syy = (1/Ay)*mphint2(model,'solid.sy','surface','selection',[5]);%*mphint2(model,'1','surface','selection',[50 52 171 336 338]);
        Syz = (1/Az)*mphint2(model,'solid.syz','surface','selection',[4]);%*mphint2(model,'1','surface','selection',[10 45 160 326 329]);
        Szz = (1/Az)*mphint2(model,'solid.sz','surface','selection',[4]);%*mphint2(model,'1','surface','selection',[10 45 160 326 329]);
        
        C(1,3) = Sxx;
        C(2,3) = Syy;
        C(3,3) = Szz;
        C(4,3) = Sxy;
        C(5,3) = Syz;
        C(6,3) = Sxz;
         disp(C)
    end
    
    if i == 4 % eps 12
        
        
        model.physics('solid').feature('disp1').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp1').set('U0', {'deltax*0.5*y/(a)'; '0'; '0'});
        model.physics('solid').feature('disp2').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp2').set('U0', {'0'; 'deltax*0.5*x/(a)'; '0'});
        model.physics('solid').feature('disp3').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp3').set('U0', {'deltax*0.5*y/(a)'; '0.5*deltax'; '0'});
        model.physics('solid').feature('disp4').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp4').set('U0', {'0.5*deltax'; 'deltax*0.5*x/(a)'; '0'});
        model.physics('solid').feature('disp5').set('Direction', {'0'; '0'; '1'});
        model.physics('solid').feature('disp5').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp6').set('Direction', {'0'; '0'; '1'});
        model.physics('solid').feature('disp6').set('U0', {'0'; '0'; '0'});
        mphsave(model,'temporary/sheartest.mph');
              % Run Model
        model.study('std1').run;
        
        %Fetch By integration and calculate stress related to the dimension
        %of the UC
        Sxx = (1/Ax)*mphint2(model,'solid.sx','surface','selection',[1]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Sxy = (1/Ax)*mphint2(model,'solid.sxy','surface','selection',[1]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Sxz = (1/Ax)*mphint2(model,'solid.sxz','surface','selection',[1]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Syy = (1/Ay)*mphint2(model,'solid.sy','surface','selection',[5]);%*mphint2(model,'1','surface','selection',[50 52 171 336 338]);
        Syz = (1/Az)*mphint2(model,'solid.syz','surface','selection',[3]);%*mphint2(model,'1','surface','selection',[10 45 160 326 329]);
        Szz = (1/Az)*mphint2(model,'solid.sz','surface','selection',[3]);%*mphint2(model,'1','surface','selection',[10 45 160 326 329]);
        
        C(1,4) = Sxx;
        C(2,4) = Syy;
        C(3,4) = Szz;
        C(4,4) = Sxy;
        C(5,4) = Syz;
        C(6,4) = Sxz;
         disp(C);
    end
    
    if i == 5 % eps 13
        
        model.physics('solid').feature('disp1').set('Direction', {'1'; '0'; '0'});
        model.physics('solid').feature('disp1').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp2').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp2').set('U0', {'0'; 'deltax*0.5*z/(a)'; '0'});
        model.physics('solid').feature('disp3').set('Direction', {'1'; '0'; '0'});
        model.physics('solid').feature('disp3').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp4').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp4').set('U0', {'0'; 'deltax*0.5*z/(a)'; '0.5*deltaz'});
        model.physics('solid').feature('disp5').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp5').set('U0', {'0'; '0.5*deltax'; '0.5*y/(a)*deltaz'});
        model.physics('solid').feature('disp6').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp6').set('U0', {'0'; '0'; '0.5*y/(a)*deltaz'});
        
              % Run Model
        model.study('std1').run;
        
        %Fetch By integration and calculate stress related to the dimension
        %of the UC
        Sxx = (1/Ax)*mphint2(model,'solid.sx','surface','selection',[1 ]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Sxy = (1/Ax)*mphint2(model,'solid.sxy','surface','selection',[1]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Sxz = (1/Ax)*mphint2(model,'solid.sxz','surface','selection',[1 ]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Syy = (1/Ay)*mphint2(model,'solid.sy','surface','selection',[5]);%*mphint2(model,'1','surface','selection',[50 52 171 336 338]);
        Syz = (1/Az)*mphint2(model,'solid.syz','surface','selection',[3]);%*mphint2(model,'1','surface','selection',[10 45 160 326 329]);
        Szz = (1/Az)*mphint2(model,'solid.sz','surface','selection',[3]);%*mphint2(model,'1','surface','selection',[10 45 160 326 329]);
        
        C(1,5) = Sxx;
        C(2,5) = Syy;
        C(3,5) = Szz;
        C(4,5) = Sxy;
        C(5,5) = Syz;
        C(6,5) = Sxz;
         disp(C);
    end
    
    if i == 6 % eps 23
        
        model.physics('solid').feature('disp1').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp1').set('U0', {'0.5*z/(a)*deltax'; '0'; '0'});
        model.physics('solid').feature('disp2').set('Direction', {'0'; '1'; '0'});
        model.physics('solid').feature('disp2').set('U0', {'0'; '0'; '0'});
        model.physics('solid').feature('disp3').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp3').set('U0', {'0.5*z/(a)*deltax'; '0'; '0.5*deltaz'});
        model.physics('solid').feature('disp4').set('Direction', {'0'; '1'; '0'});
        model.physics('solid').feature('disp4').set('U0', {'0'; '0'; '0.5*deltaz'});
        model.physics('solid').feature('disp5').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp5').set('U0', {'0.5*deltax'; '0'; 'deltax*0.5*x/(a)'});
        model.physics('solid').feature('disp6').set('Direction', {'1'; '1'; '1'});
        model.physics('solid').feature('disp6').set('U0', {'0'; '0'; 'deltaz*0.5*x/(a)'});
        
              % Run Model
        model.study('std1').run;
        
        %Fetch By integration and calculate stress related to the dimension
        %of the UC
        Sxx = (1/Ax)*mphint2(model,'solid.sx','surface','selection',[1]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Sxy = (1/Ax)*mphint2(model,'solid.sxy','surface','selection',[1]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Sxz = (1/Ax)*mphint2(model,'solid.sxz','surface','selection',[1 ]);%*mphint2(model,'1','surface','selection',[1 6 18 37 40]);
        Syy = (1/Ay)*mphint2(model,'solid.sy','surface','selection',[5]);%*mphint2(model,'1','surface','selection',[50 52 171 336 338]);
        Syz = (1/Az)*mphint2(model,'solid.syz','surface','selection',[3]);%*mphint2(model,'1','surface','selection',[10 45 160 326 329]);
        Szz = (1/Az)*mphint2(model,'solid.sz','surface','selection',[3]);%*mphint2(model,'1','surface','selection',[10 45 160 326 329]);
        
        C(1,6) = Sxx;
        C(2,6) = Syy;
        C(3,6) = Szz;
        C(4,6) = Sxy;
        C(5,6) = Syz;
        C(6,6) = Sxz;
        disp(C);
    end
    
end
toc;

disp('Compare to linear elastic material:')

c = E * (1-nu)/((1+nu)*(1-2*nu));

Cmat = c.* [1 nu/(1-nu) nu/(1-nu) 0 0 0;
    nu/(1-nu) 1 nu/(1-nu) 0 0 0;
    nu/(1-nu) nu/(1-nu) 1 0 0 0;
    0 0 0 (1-2*nu)/(2*(1-nu)) 0 0;
    0 0 0 0 (1-2*nu)/(2*(1-nu)) 0;
    0 0 0 0 0 (1-2*nu)/(2*(1-nu))]
