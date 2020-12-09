clc; clear


b = 40;
h = 40;
t_f = 3;
t_w = 1;
E = 70000;
G = E/(2*(1+0.3));

S = 100;


A = 2*b*t_f + h*t_w;
x_cg = b^2*t_f/A;
y_cg = b/2;

I_xx = h^3*t_w/12 + b*h^2*t_f/2;
I_yy = h*t_w*x_cg^2 + 2* ( b^3*t_f/12 + b*t_f * (b/2 - x_cg)^2 );

J = h*t_w^3/3 + 2*b*t_f^3/3;

x_sc = (h^2*b^2*t_f)/(4*I_xx);
y_sc = h/2;

T = x_sc * S;

w_lc = h*x_sc/2 * T/(G*J);
w_rc = -(h*x_sc/2 - h*b/2) * T/(G*J);

I_os = I_xx + I_yy + A*((x_sc + x_cg)^2+(y_sc-y_cg)^2);

%% Gamma_R

% Taking 0 as the edge of the cross section

% gamma_R1 = (2/3)*h^2*b^2 * t_f;
% gamma_R2 = -h^3/(3*x_sc) * ( (b/2 - x_sc)^3 - b^3 / 8 ) * t_w;
% gamma_R3 = 2*h^2/3 * ( (b-x_sc)^3 - (b/2 - x_sc)^3 ) * t_f;
% 
% gamma_R = gamma_R1 + gamma_R2 + gamma_R3;


% Taking 0 as the center of symmetry

gamma_R1 = (x_sc^2*h^3)/24 * t_w;
gamma_R2 = h^2/12 * ( (b-x_sc)^3 + x_sc^3 ) * t_f;

gamma_R = 2* (gamma_R1 + gamma_R2);

%% Euler buckling

L = [500, 1000, 2000];

buckling_horizontal = (pi^2 * E * I_yy)./ (4*L.^2);
buckling_vertical = (pi^2 * E * I_xx)./ (4*L.^2);


%% Local buckling

sigma_crit_flange = 0.385 * E * (t_f/b)^2;
sigma_crit_web = 3.6 * E * (t_w/h)^2;

buckling_local = A * min(sigma_crit_flange, sigma_crit_web);

buckling_local_alt = (sigma_crit_flange * 2*b*t_f  +  sigma_crit_web * t_w*h);


%% Torsion buckling

buckling_torsion = (A/I_os) * (G*J + (pi^2*E*gamma_R)./L.^2);


%% Combined buckling

buckling_combined = zeros(3,1);

x_s = x_sc + x_cg;
y_s = 0;

for i=1:3

    
    P_cr_x = buckling_vertical(i);
    P_cr_y = buckling_horizontal(i);
    P_cr_theta = buckling_torsion(i);

    coefficients = [x_s^2+y_s^2 - I_os/A, ...
                   -(P_cr_x*y_s^2 + P_cr_y*x_s^2) + (P_cr_x + P_cr_y + P_cr_theta) * I_os/A, ...
                   - I_os/A * (P_cr_x*P_cr_y + P_cr_x*P_cr_theta + P_cr_y*P_cr_theta), ...
                   I_os/A * (P_cr_x + P_cr_y + P_cr_theta)];

    provisional = roots(coefficients);
    buckling_combined(i) = provisional(2);
    
end

%% Print results to file

fIDoutput = fopen("numerics.txt", "w");

fIDVars = fopen("variables.txt", "r");

i = 0;
names = {};
vars = {};
units = {};

while ~feof(fIDVars)
    i = i+1;
    line  = fgetl(fIDVars);
    [splitStr, perdoned] = regexp(line, "([^\t\s]*)\s+([^\t\s]*)\s+(.*)" ,'tokens','match');
    names{i} = splitStr{1,1}{1,3};
    units{i} = splitStr{1,1}{1,2};
    vars{i} = splitStr{1,1}{1,1};
end

space_name = 38;
space_number = 19;

fprintf(fIDoutput, "\n\tNumerical values for the problem (all in mm, N):\n\n\n");

for i=1:length(vars)
    n_tabs1 = (space_name - length(names{i}))/4;
    n_tabs2 = fix((space_number - strlength(sprintf("%.3f",eval(vars{i}))))/4); 
    
    fprintf(fIDoutput, "%s:", names{i});
    
    for j=1:n_tabs1
        fprintf(fIDoutput,"\t");
    end
    fprintf(fIDoutput, "%.3f", eval(vars{i}));
    
    for j=1:n_tabs2
        fprintf(fIDoutput,"\t");
    end
    fprintf(fIDoutput, "%s\n", units{i});
    
end

fclose(fIDVars);


fprintf(fIDoutput, "\n\n\n\n\tBuckling loads:\n\n\n");
fprintf(fIDoutput, "Direction\t|\t%d mm\t\t\t%d mm\t\t\t%d mm\n", L(1), L(2), L(3));
fprintf(fIDoutput, "------------+---------------------------------------------\n");
fprintf(fIDoutput, "Euler horz.\t|\t%.1f N\t\t%.1f N\t\t%.1f N\n", buckling_horizontal(1), ...
        buckling_horizontal(2), buckling_horizontal(3));
fprintf(fIDoutput, "Euler vert.\t|\t%.1f N\t\t%.1f N\t\t%.1f N\n", buckling_vertical(1), ...
        buckling_vertical(2), buckling_vertical(3));
fprintf(fIDoutput, "Torsional\t|\t%.1f N\t\t%.1f N\t\t%.1f N\n", buckling_torsion(1), ...
        buckling_torsion(2), buckling_torsion(3)); 
fprintf(fIDoutput, "Combined\t|\t%.1f N\t\t%.1f N\t\t%.1f N\n", buckling_combined(1), ...
        buckling_combined(2), buckling_combined(3));
fprintf(fIDoutput, "Local\t\t|\t%.1f N\t\t%.1f N\t\t%.1f N\n", buckling_local, ...
        buckling_local, buckling_local);    
    
    
    
fclose(fIDoutput);


