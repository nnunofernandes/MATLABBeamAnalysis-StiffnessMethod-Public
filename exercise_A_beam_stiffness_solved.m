%% EXERCISE A: CONTINUOUS BEAM BY THE DIRECT STIFFNESS METHOD
% Standalone teaching script. Run this file; no toolbox or extra file needed.
% Use MATLAB sections (%%) to walk through the solution during the lecture.
%
% BEAM: fixed support 1 at x=0; roller 2 at x=6 m; fixed support 3 at x=10 m.
% LOADS: 25 kN downward at x=3 m; 15 kN/m downward over span 2 (6 <= x <= 10).
% REQUIRED: moments at fixed supports 1 and 3. EI is constant but unspecified.
% MODEL: small-displacement, prismatic Euler-Bernoulli beam; no shear/axial DOF.
% UNITS: kN, m, radians; consequently EI is in kN*m^2.
%
% EXPECTED ANSWER:
% Reaction moments ON the beam: support 1 = +18.500 kN*m (CCW),
%                              support 3 = -20.375 kN*m (CW).
% Internal support moments (sagging positive): [-18.500; -19.250; -20.375].

clear; clc; close all;
format short g;

%% 1. Define the physical model
xNode = [0 6 10];                  % Node locations, m.
connectivity = [1 2; 2 3];          % One beam element per loaded span.
q = [0; -15];                      % Uniform loads, kN/m; downward negative.
P = [-25; 0];                      % Interior point loads, kN.
a = [3; 0];                        % Point-load distance from element left, m.
settlement = 0;                    % All support displacements are zero.

% Because every element has the SAME EI and there are no support movements,
% the common EI cancels from the force solution. A reference EI=1 is valid
% for finding reactions and moments. It is NOT a measured beam stiffness.
% The reference rotations/deflections are not physical numerical predictions.
% Enter an actual common EI here if physical displacements are required.
EI = 1;                            % Reference value, kN*m^2.
referenceMoments = [-18.5; -19.25; -20.375];
exerciseName = 'Exercise A';
supportLabel = 'Support 2: roller; v2 = 0; theta2 free';
deflectionFactor = EI;              % Plot EI*v: independent of reference EI.
deflectionLabel = 'EI v (kN m^3)';
deflectionTitle = 'Stiffness-scaled deflection EI*v: divide by actual EI to obtain v';

%% 2. Degrees of freedom and boundary conditions
% Each Euler-Bernoulli beam node has TWO bending degrees of freedom:
% vertical displacement v, and cross-section rotation theta = dv/dx.
% Positive v and forces: upward. Positive theta and nodal moments: CCW.
% The global order is d = [v1; theta1; v2; theta2; v3; theta3].
% Element 1 uses DOFs [1 2 3 4]; element 2 uses DOFs [3 4 5 6].
% Both end supports are fixed. The roller restrains v2, but NOT theta2.
% A roller under a continuous beam is NOT an internal hinge: M2 need not be 0.

nNode = numel(xNode);
nElement = size(connectivity,1);
nDOF = 2*nNode;
elementDOF = zeros(nElement,4);
for e = 1:nElement
    i = connectivity(e,1);
    j = connectivity(e,2);
    elementDOF(e,:) = [2*i-1, 2*i, 2*j-1, 2*j];
end

prescribedDOF = [1 2 3 5 6];
freeDOF = 4;                       % Only theta2 is unknown.
d = zeros(nDOF,1);
d(3) = settlement;                 % Downward settlement is NEGATIVE.

%% 3. Element stiffness matrices and consistent load vectors
% For an element of length L, ordered as [vi; thetai; vj; thetaj]:
%
% ke = EI/L^3 * [ 12    6L   -12    6L
%                6L   4L^2  -6L   2L^2
%               -12   -6L    12   -6L
%                6L   2L^2  -6L   4L^2 ].
%
% A uniform transverse load q (upward positive) has equivalent nodal load
% fe,UDL = q * [L/2; L^2/12; L/2; -L^2/12].
% These loads act on the NODES in K*d = F. Fixed-end member actions have
% the OPPOSITE sign: s_fixed = -fe when all element displacements are zero.
%
% An interior point load P is distributed using cubic Hermite functions:
% fe,P = P*[N1(a); N2(a); N3(a); N4(a)].
% Equivalent nodal MOMENTS must be included as well as nodal forces.

K = zeros(nDOF);
F = zeros(nDOF,1);
elementK = zeros(4,4,nElement);
elementLoad = zeros(4,nElement);
L = zeros(nElement,1);

for e = 1:nElement
    nodes = connectivity(e,:);
    L(e) = xNode(nodes(2)) - xNode(nodes(1));
    Le = L(e);

    ke = EI/Le^3 * ...
        [12,       6*Le,      -12,       6*Le;
         6*Le,     4*Le^2,    -6*Le,     2*Le^2;
         -12,      -6*Le,     12,        -6*Le;
         6*Le,     2*Le^2,    -6*Le,     4*Le^2];

    fe = q(e)*[Le/2; Le^2/12; Le/2; -Le^2/12];

    if P(e) ~= 0
        assert(a(e) > 0 && a(e) < Le, 'Point load must lie inside its span.');
        xi = a(e)/Le;
        N = [1 - 3*xi^2 + 2*xi^3;
             Le*(xi - 2*xi^2 + xi^3);
             3*xi^2 - 2*xi^3;
             Le*(-xi^2 + xi^3)];
        fe = fe + P(e)*N;
    end

    elementK(:,:,e) = ke;
    elementLoad(:,e) = fe;

    % ASSEMBLY: add contributions at the corresponding global DOFs.
    dof = elementDOF(e,:);
    K(dof,dof) = K(dof,dof) + ke;
    F(dof) = F(dof) + fe;
end

disp('Global DOF order: [v1, theta1, v2, theta2, v3, theta3]');
disp('Global stiffness matrix K (mixed force/displacement units):');
disp(K);
disp('Consistent global load vector F (kN and kN*m):');
disp(F);

%% 4. Partition and solve the stiffness equations
% Reordering into free (f) and prescribed (c) DOFs gives:
%
% [Kff Kfc] [df] = [Ff       ]
% [Kcf Kcc] [dc]   [Fc + Rc  ].
%
% Therefore: Kff*df = Ff - Kfc*dc.
% The term -Kfc*dc is essential for a nonzero support settlement.
% Use MATLAB backslash to solve the system; do not explicitly invert K.

Kff = K(freeDOF,freeDOF);
Kfc = K(freeDOF,prescribedDOF);
effectiveLoad = F(freeDOF) - Kfc*d(prescribedDOF);
d(freeDOF) = Kff \ effectiveLoad;

%% 5. Recover support reactions and member end actions
% R is the force/moment exerted by supports ON the beam.
% At the free rotation DOF, R(4) must be zero.
% Each column of endAction is [Vi; Mi; Vj; Mj], with forces upward
% and moments CCW positive at BOTH ends (external end-action convention).

R = K*d - F;
endAction = zeros(4,nElement);
for e = 1:nElement
    de = d(elementDOF(e,:));
    endAction(:,e) = elementK(:,:,e)*de - elementLoad(:,e);
end

% Internal bending moment uses positive SAGGING throughout the beam:
% M(0) = -Mi and M(L) = +Mj.
M_support = [-endAction(2,1); endAction(4,1); endAction(4,2)];

fprintf('\nSUPPORT REACTIONS (upward / CCW positive)\n');
supportReactions = table((1:3)', R(1:2:5), R(2:2:6), ...
    'VariableNames', {'Support','Vertical_kN','ReactionMoment_kNm'});
disp(supportReactions);

fprintf('INTERNAL SUPPORT MOMENTS (sagging positive)\n');
disp(table((1:3)', M_support, ...
    'VariableNames', {'Support','InternalMoment_kNm'}));

fprintf('MEMBER END ACTIONS (upward / CCW positive)\n');
disp(table((1:nElement)', endAction(1,:)', endAction(2,:)', ...
    endAction(3,:)', endAction(4,:)', ...
    'VariableNames', {'Element','Vi_kN','Mi_kNm','Vj_kN','Mj_kNm'}));

% EI*theta2 is known even though the physical EI is unspecified.
fprintf('theta2 = %.6f / EI radians, with EI entered in kN*m^2.\n', EI*d(4));
fprintf('EI is a reference value here; raw d values are not physical predictions.\n');
fprintf('Requested: support 1 = %.3f kN*m CCW; support 3 = %.3f kN*m CW.\n', ...
    R(2), -R(6));

%% 6. Recover continuous shear, moment, and deflection fields
% Let x be measured from the left end of an element, and <x-a> = max(x-a,0).
% Equilibrium gives the EXACT fields for these piecewise constant loads:
% V(x) = Vi + q*x + P*H(x-a)
% M(x) = -Mi + Vi*x + q*x^2/2 + P*<x-a>.
% For this sign convention, EI*v'' = M. Integrating twice gives:
% v(x) = vi + thetai*x + ...
%        (-Mi*x^2/2 + Vi*x^3/6 + q*x^4/24 + P*<x-a>^3/6)/EI.
%
% This load-aware recovery includes bending BETWEEN the nodes.
% A plot made only by joining nodal displacements would miss that bending.
% With constant EI and the exact consistent loads, the nodal solution and
% these recovered fields are exact for the ideal Euler-Bernoulli model.

xPlot = cell(nElement,1);
vPlot = cell(nElement,1);
VPlot = cell(nElement,1);
MPlot = cell(nElement,1);

for e = 1:nElement
    x = linspace(0,L(e),301);
    if P(e) ~= 0
        % Repeated x at the load produces a vertical shear-force jump.
        x = sort([x, a(e), a(e)]);
    end
    ramp = max(x-a(e),0);
    H = double(x > a(e));
    if P(e) ~= 0
        atLoad = find(abs(x-a(e)) < 1e-12);
        H(atLoad) = 1;
        H(atLoad(1)) = 0;
    end

    s = endAction(:,e);
    de = d(elementDOF(e,:));
    V = s(1) + q(e)*x + P(e)*H;
    M = -s(2) + s(1)*x + q(e)*x.^2/2 + P(e)*ramp;
    v = de(1) + de(2)*x + ...
        (-s(2)*x.^2/2 + s(1)*x.^3/6 + q(e)*x.^4/24 + ...
         P(e)*ramp.^3/6)/EI;
    theta = de(2) + ...
        (-s(2)*x + s(1)*x.^2/2 + q(e)*x.^3/6 + P(e)*ramp.^2/2)/EI;

    % Check that the recovered field meets both nodal end conditions.
    assert(abs(v(end)-de(3)) < 1e-8*max(1,norm(d,inf)), ...
        'Recovered right-end displacement does not match the nodal value.');
    assert(abs(theta(end)-de(4)) < 1e-8*max(1,norm(d,inf)), ...
        'Recovered right-end rotation does not match the nodal value.');

    xPlot{e} = xNode(connectivity(e,1)) + x;
    vPlot{e} = v;
    VPlot{e} = V;
    MPlot{e} = M;
end

%% 7. Verify equilibrium, continuity, and the reference solution
% Global equilibrium must use the ORIGINAL physical loads, not a mixture
% of physical loads and their equivalent nodal loads (which double counts).
externalVerticalLoad = sum(q(:).*L + P(:));
loadMomentAboutLeft = 0;
for e = 1:nElement
    xLeft = xNode(connectivity(e,1));
    loadMomentAboutLeft = loadMomentAboutLeft + ...
        q(e)*L(e)*(xLeft + L(e)/2) + P(e)*(xLeft + a(e));
end
forceResidual = sum(R(1:2:5)) + externalVerticalLoad;
momentResidual = sum(R(1:2:5).*xNode(:)) + ...
    sum(R(2:2:6)) + loadMomentAboutLeft;

assert(norm(K-K','fro') < 1e-10*max(1,norm(K,'fro')), ...
    'The assembled stiffness matrix must be symmetric.');
assert(abs(forceResidual) < 1e-8*max(1,abs(externalVerticalLoad)), ...
    'Global vertical equilibrium failed.');
assert(abs(momentResidual) < 1e-8*max(1,abs(loadMomentAboutLeft)), ...
    'Global moment equilibrium failed.');
assert(norm(R(freeDOF),inf) < 1e-8, ...
    'An unrestrained rotation cannot develop a support reaction moment.');
assert(abs(endAction(4,1)+endAction(2,2)) < 1e-8, ...
    'Member end moments must balance at the unloaded roller joint.');
% This last benchmark applies to the original exercise input values.
% Disable it when deliberately changing the geometry, loading, EI, or settlement.
assert(norm(M_support-referenceMoments,inf) < 1e-7, ...
    'Support moments differ from the independent slope-deflection solution.');
fprintf('Checks passed: force/moment equilibrium, joint balance, and reference moments.\n');
fprintf('Vertical residual = %.3e kN; moment residual = %.3e kN*m.\n', ...
    forceResidual, momentResidual);

%% 8. Draw the beam and the response diagrams
% Positive shear follows V = dM/dx; positive bending moment is sagging.
% Both shear and moment are plotted positive ABOVE their zero axes.
figure('Color','w','Name',exerciseName,'NumberTitle','off', ...
       'Position',[100 50 1100 900]);
spanLength = xNode(end)-xNode(1);
supportWidth = 0.018*spanLength;

subplot(4,1,1); hold on;
plot(xNode,[0 0 0],'k-','LineWidth',4);
for xWall = [xNode(1), xNode(end)]
    plot([xWall xWall],[-0.5 0.5],'k-','LineWidth',5);
end
patch(xNode(2)+[-supportWidth 0 supportWidth], ...
    [-0.45 -0.07 -0.45],[0.85 0.85 0.85],'EdgeColor','k');
for e = 1:nElement
    xLeft = xNode(connectivity(e,1));
    if q(e) ~= 0
        xArrow = linspace(xLeft+0.06*L(e),xLeft+0.94*L(e),9);
        quiver(xArrow,1.05*ones(size(xArrow)),zeros(size(xArrow)), ...
            -0.92*ones(size(xArrow)),0,'k','LineWidth',1);
        plot([xLeft xLeft+L(e)],[1.05 1.05],'k-');
        text(xLeft+L(e)/2,1.4,sprintf('%.0f kN/m',abs(q(e))), ...
            'HorizontalAlignment','center');
    end
    if P(e) ~= 0
        quiver(xLeft+a(e),1.05,0,-0.92,0,'k','LineWidth',1.8);
        text(xLeft+a(e),1.4,sprintf('%.0f kN',abs(P(e))), ...
            'HorizontalAlignment','center');
    end
    text(xLeft+L(e)/2,-0.5,sprintf('Element %d: L = %g m',e,L(e)), ...
        'HorizontalAlignment','center');
end
text(xNode(2),-0.9,supportLabel,'HorizontalAlignment','center');
text(xNode(1),-0.75,'Support 1','HorizontalAlignment','center');
text(xNode(3),-0.75,'Support 3','HorizontalAlignment','center');
xlim([xNode(1)-0.04*spanLength, xNode(end)+0.04*spanLength]);
ylim([-1.15 1.8]); axis off;
title([exerciseName, ' | fixed - roller - fixed']);

subplot(4,1,2); hold on; grid on;
plot([xNode(1) xNode(end)],[0 0],'k:');
for e = 1:nElement
    plot(xPlot{e},deflectionFactor*vPlot{e},'b-','LineWidth',1.8);
end
plot(xNode,deflectionFactor*d(1:2:5),'ko','MarkerFaceColor','k');
xlim(xNode([1 end])); ylabel(deflectionLabel);
title(deflectionTitle);

subplot(4,1,3); hold on; grid on;
plot([xNode(1) xNode(end)],[0 0],'k:');
for e = 1:nElement
    plot(xPlot{e},VPlot{e},'Color',[0 0.5 0.3],'LineWidth',1.8);
end
plot([xNode(2) xNode(2)],[VPlot{1}(end) VPlot{2}(1)], ...
    'Color',[0 0.5 0.3],'LineWidth',1.8);
xlim(xNode([1 end])); ylabel('V (kN)');
title('Shear-force diagram');

subplot(4,1,4); hold on; grid on;
plot([xNode(1) xNode(end)],[0 0],'k:');
for e = 1:nElement
    plot(xPlot{e},MPlot{e},'r-','LineWidth',1.8);
end
plot(xNode,M_support,'ko','MarkerFaceColor','k');
xlim(xNode([1 end])); xlabel('Position along beam, x (m)'); ylabel('M (kN m)');
title('Bending-moment diagram: sagging positive, hogging negative');

%% 9. Suggested lecture questions
% 1. Why does the roller leave theta2 free but still transmit bending moment?
% 2. Why must the distributed/point load vector contain nodal moments?
% 3. Why are a left reaction moment and the left internal moment opposite?
% 4. Which quantities change if the common EI is multiplied by two?
% 5. What changes in Kff*df = Ff-Kfc*dc if support 2 moves downward?
% Workspace variables K, F, d, R, endAction, and M_support remain available.
