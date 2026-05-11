clc;
clear all;
close all;
script_dir = fileparts(mfilename('fullpath'));  % directory where this script lives
%% %%%%%%%%%% Geometric parameters %%%%%%%%%%%%
%% VALIDATION: Kim et al. (2014) tripod centrifuge test T2
%% Paper: J. Geotech. Geoenviron. Eng. 2014.140:04014008
%% Centrifuge at 70g, prototype scale values used below
%% Expected: My = 93 MN.m at 0.6 deg rotation
do_b = 6.5; %bucket outer diameter [m] - Kim: Dt=6.5m
thickness = 0.025; %bucket thickness [m] - Kim: 25mm prototype
l_b = 8.0; %length of the bucket skirt [m] - Kim: Lt=8.0m (L/D=1.23)
Z_lid =0.0001; % thickness of the lid of bucket
n = 40; %number of strips in a single circumference
n_l =25; % number of rings in the bottom
n_q =10;
eccentricity = 33.39; % [m] - Kim: load at 33m above seabed (M/H=33m)
th1= 60;
th2 = 180;
th3 = 300;
Deadload=4650; %kN - Kim: tripod weight 4.65 MN
l = 15.5; % [m] - Kim: CtC=26.85m, circumradius=26.85/sqrt(3)=15.5m
%% Group interaction p-multiplier (Barari et al. 2023 parametric calibration)
SD_ratio = l / do_b;  % spacing-to-diameter ratio (S/D)
delta_ref = 0.05;  % reference displacement (m) for full f_g reduction
bucket_angles = [th1, th2, th3];  % [60, 180, 300]
%% %%%%%%%%%% Soil parameters %%%%%%%%%%%%
%% Yellow Sea silty sand SM layer (0 to -11m), Dr~70-83%
gamma = 9.65; % buoyant unit weight [kN/m3] - estimated from gd=1.55, Gs=2.65
S_gamma = 0.5;
e_max =0;
e_min =0;
D_r = 0.827; % relative density - Kim T2: Dr=82.7% for SM layer
m_c = 161.42*D_r^2+199.8*D_r+36.877; % 80-150 soft sand % 150-250 medium sand % 250-400 stiff sand  (rewrite based on dr)
A_c = 0.3474*D_r^2+0.4222*D_r+0.328;
fi_crit = 33; % friction angle [deg] - Kim: phi'=33 from direct shear test
m =0; % for finding the FI_preak %m is either 0 or 3
delta = 22; % interface friction [deg] - Kim: 2/3*phi'=2/3*33=22
%%
ro_b = do_b/2; %outer radious of the bucket
di_b =do_b-2*thickness; %bucket inner diameter
ri_b = di_b/2; %outer radious of the bucket
r_pile = di_b/2 + thickness/2; %from the center of buckt to the center of the thickness
r_b = abs(di_b-do_b)/4; %radius of the pile
penetration =l_b; %fully penetrated
circle_split = linspace(0,360,n+1);
R_split = linspace(0,r_pile,n_q);
penetration1 =l_b; %fully penetrated
penetration2 =l_b; %fully penetrated
penetration3 =l_b; %fully penetrated
%%
ro_b = do_b/2; %outer radious of the bucket
di_b =do_b-2*thickness; %bucket inner diameter
ri_b = di_b/2; %outer radious of the bucket
r_pile = ri_b/2 + thickness/2; %from the center of buckt to the center of the thickness
r_b = abs(di_b-do_b)/4; %radius of the pile

circle_split = linspace(0,360,n+1);
R_split = linspace(0,r_pile,n_q);
R_sec = R_split(2)-R_split(1);
%% Code initiation
QQ = zeros(1,10000);
kkk =0;
abc=0;
cc1=0;
% full penetration action
    abc=abc+1;
 cc1=cc1+1;
%    if abs(l_b - penetration)<=0.01
%        FP =1;
%    else
%         FP=0;
%    end
count =0;

%%  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ dealload effect with no tilt angle $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  %%
dz=0;
SUMz=0;
while true
    if SUMz>Deadload
    dz=dz-0.00001;
    else
     dz=dz+0.00001;
    end
for dd=0
   q=[0 -dd];
phi = q(1);     % Rotation around x-axis
theta = q(2);   % Rotation around y-axis
psi = 0;     % Rotation around z-axis
j=0;
for k=0:l_b/n_l:l_b-(l_b/n_l)
for i=2:n+1 
    j=j+1;
    first_s_end(:,j) = [l*cosd(th1)+ro_b*cosd(circle_split(i)) l*sind(th1)+ro_b*sind(circle_split(i)) -k-(l_b/n_l)/2]'; %  bucket points on the skirt  w. r. t. the body frame
    second_s_end(:,j) = [l*cosd(th2)+ro_b*cosd(circle_split(i)) l*sind(th2)+ro_b*sind(circle_split(i)) -k-(l_b/n_l)/2]'; %  bucket points on the skirt  w. r. t. the body frame
    third_s_end(:,j) = [l*cosd(th3)+ro_b*cosd(circle_split(i)) l*sind(th3)+ro_b*sind(circle_split(i)) -k-(l_b/n_l)/2]'; %  bucket points on the skirt  w. r. t. the body frame

    % notice: the stars are in the middle of layers.
end
end
% w. r. t.  fixed frame
A1 = zeros(size(first_s_end));
A2 = zeros(size(second_s_end));
A3 = zeros(size(third_s_end));
A1(3,:) = -dz;
A2(3,:) = -dz;
A3(3,:) = -dz;
first_s_end_fixed       = A1 + first_s_end;
second_s_end_fixed = A2+ second_s_end;
third_s_end_fixed     = A3 + third_s_end;

j=0;
for i=2:n+1 
    j=j+1;
    Pile_first_s_end(:,j) = [l*cosd(th1)+r_pile*cosd(circle_split(i)) l*sind(th1)+r_pile*sind(circle_split(i)) -l_b]'; % tip of the bucket points  w. r. t the body frame
        Pile_second_s_end(:,j) = [l*cosd(th2)+r_pile*cosd(circle_split(i)) l*sind(th2)+r_pile*sind(circle_split(i)) -l_b]'; % tip of the bucket points  w. r. t the body frame
            Pile_third_s_end(:,j) = [l*cosd(th3)+r_pile*cosd(circle_split(i)) l*sind(th3)+r_pile*sind(circle_split(i)) -l_b]'; % tip of the bucket points  w. r. t the body frame
end
B1 = zeros(size(Pile_first_s_end));
B2 = zeros(size(Pile_second_s_end));
B3 = zeros(size(Pile_third_s_end));

B1(3,:) = -dz;
B2(3,:) = -dz;
B3(3,:) = -dz;

Pile_first_s_end_fixed = Pile_first_s_end + B1;
Pile_second_s_end_fixed = Pile_second_s_end + B2;
Pile_third_s_end_fixed = Pile_third_s_end + B3;

j=0;
jj = 0;
co =0;
R_sec = R_split(2)-R_split(1);
for jj =2:size(R_split,2)
j=0;
% w.r.t body frame
for i=2:n+1 
    j=j+1;
r_Pile_first_s_end(:,j+co*n) = [l*cosd(th1)+(R_split(jj)-R_sec/2)*cosd(circle_split(i)) l*sind(th1)+(R_split(jj)-R_sec/2)*sind(circle_split(i)) -Z_lid]';
r_Pile_second_s_end(:,j+co*n) = [l*cosd(th2)+(R_split(jj)-R_sec/2)*cosd(circle_split(i)) l*sind(th2)+(R_split(jj)-R_sec/2)*sind(circle_split(i)) -Z_lid]';
r_Pile_third_s_end(:,j+co*n) = [l*cosd(th3)+(R_split(jj)-R_sec/2)*cosd(circle_split(i)) l*sind(th3)+(R_split(jj)-R_sec/2)*sind(circle_split(i)) -Z_lid]';

end
co = co+1;
end
j=0;
C1 = zeros(size(r_Pile_first_s_end));
C2 = zeros(size(r_Pile_second_s_end));
C3 = zeros(size(r_Pile_third_s_end));

C1(3,:) = -dz;
C2(3,:) = -dz;
C3(3,:) = -dz;

r_Pile_first_s_end_fixed = r_Pile_first_s_end+C1;
r_Pile_second_s_end_fixed = r_Pile_second_s_end+C2;
r_Pile_third_s_end_fixed = r_Pile_third_s_end+C3;

active_z_values1 = Zvalues(l_b,n,n_l,penetration1);
active_z_values2 = Zvalues(l_b,n,n_l,penetration2);
active_z_values3 = Zvalues(l_b,n,n_l,penetration3);

FI_peak1 = fi_finder(gamma,active_z_values1,D_r,fi_crit,m);
FI_peak2 = fi_finder(gamma,active_z_values2,D_r,fi_crit,m);
FI_peak3 = fi_finder(gamma,active_z_values3,D_r,fi_crit,m);

fi1 = FI_peak1;
fi2= FI_peak2;
fi3 = FI_peak3;

fi_end1 = fi1(end);
fi_end2 = fi2(end);
fi_end3 = fi3(end);

displacement1 = (first_s_end(3, :) - first_s_end_fixed(3, :)) * 1e3;  
displacement2 = (second_s_end(3, :) - second_s_end_fixed(3, :)) * 1e3;  
displacement3 = (third_s_end(3, :) - third_s_end_fixed(3, :)) * 1e3;  

%% py and tz forces self-weight
[~,~,P_i_out1] = coeffinder_deadload(first_s_end,first_s_end_fixed,active_z_values1,l_b,n_l,fi1,n,do_b,gamma); %KN
[~,~,P_i_out2] = coeffinder_deadload(second_s_end,second_s_end_fixed,active_z_values2,l_b,n_l,fi2,n,do_b,gamma); %KN
[~,~,P_i_out3] = coeffinder_deadload(third_s_end,third_s_end_fixed,active_z_values3,l_b,n_l,fi3,n,do_b,gamma); %KN

[tz_in1,tz_out1]=tzcurve_noOCR(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in2,tz_out2]=tzcurve_noOCR(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in3,tz_out3]=tzcurve_noOCR(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);

KKz1 = tz_in1+tz_out1;
KKz2 = tz_in2+tz_out2;
KKz3 = tz_in3+tz_out3;

KKz_in = tz_in1+tz_in2+tz_in3;
KKz_out = tz_out1+tz_out2+tz_out3;
KKz =KKz_in+KKz_out;
%% qz forces self-weight
[qb1_1] = qzcurve2(fi_end1,Pile_first_s_end,Pile_first_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_2] = qzcurve2(fi_end2,Pile_second_s_end,Pile_second_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_3] = qzcurve2(fi_end3,Pile_third_s_end,Pile_third_s_end_fixed,gamma,n,thickness,di_b,do_b);

[qb2_1] = qzcurve4(fi1(1),Z_lid,r_Pile_first_s_end, r_Pile_first_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_2] = qzcurve4(fi2(1),Z_lid,r_Pile_second_s_end, r_Pile_second_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_3] = qzcurve4(fi3(1),Z_lid,r_Pile_third_s_end, r_Pile_third_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);

SUMz1= sum(KKz1)+sum(qb1_1)+sum(qb2_1);
SUMz2= sum(KKz2)+sum(qb1_2)+sum(qb2_2);
SUMz3= sum(KKz3)+sum(qb1_3)+sum(qb2_3);
SUMz = SUMz1+SUMz2+SUMz3;
end
(Deadload - SUMz)
 if abs(Deadload - SUMz)<Deadload
    break
 end

end 
%% defining the local frams and relative
 %dz=0
mask1 = (active_z_values1 == 0);
mask2 = (active_z_values2 == 0);
mask3 = (active_z_values3 == 0);
%% Pre-compute bucket adjacency (passive wedge overlap)
passive_extent = do_b/2 + l_b*tand(60);
overlap_threshold = 2 * passive_extent;
n_buckets = length(bucket_angles);
bucket_centers = l * [cosd(bucket_angles); sind(bucket_angles)];
adj_matrix = false(n_buckets);
for ii = 1:n_buckets
    for jj = ii+1:n_buckets
        d_ij = norm(bucket_centers(:,ii) - bucket_centers(:,jj));
        adj_matrix(ii,jj) = d_ij < overlap_threshold;
        adj_matrix(jj,ii) = adj_matrix(ii,jj);
    end
end
%%  $$$$$$$$$$$$$$$$$$$$$$$$$   bucket rotation angle change step by step $$$$$$$$$$$$$$$$$$$$$$$$$ %%
%% Storage for Fig 19, 20, 21 equivalents (per-bucket)
dd_vec = 0.005*pi/180:0.025*pi/180:1.5*pi/180;  % <-- keep in sync with the for-loop below
n_steps_total = length(dd_vec);
M_py_pct = zeros(n_steps_total, 3);
M_tz_pct = zeros(n_steps_total, 3);
M_qztip_pct = zeros(n_steps_total, 3);
M_qzlid_pct = zeros(n_steps_total, 3);
tz_right_out_store = zeros(n_steps_total, 3);
tz_right_in_store  = zeros(n_steps_total, 3);
tz_left_out_store  = zeros(n_steps_total, 3);
tz_left_in_store   = zeros(n_steps_total, 3);
KKx_store1 = zeros(n_steps_total, n*n_l);
KKx_store2 = zeros(n_steps_total, n*n_l);
KKx_store3 = zeros(n_steps_total, n*n_l);
% Precompute right/left strip masks per bucket
strip_angles = (1:n) * (360/n);
th_vec_all = [th1, th2, th3];
right_mask_cell = cell(1,3);
left_mask_cell  = cell(1,3);
for b_idx = 1:3
    rel = cosd(strip_angles - th_vec_all(b_idx));
    right_mask_cell{b_idx} = repmat(rel > 0, 1, n_l);
    left_mask_cell{b_idx}  = repmat(rel < 0, 1, n_l);
end
Q = [0,0,0]';     %initial acenter of rotation
for dd = dd_vec
R_tilt = rodrigues_rotation(dd, 180);
first_s_lo = [l*cosd(th1), l*sind(th1), 0];
second_s_lo = [l*cosd(th2), l*sind(th2), 0];
third_s_lo = [l*cosd(th3), l*sind(th3), 0];
translation_vector = [l,0,0];

% Bucket1_local_frame = LocalFrameFunction(translation_vector,th1);
% Bucket2_local_frame = LocalFrameFunction(translation_vector,th2);
% Bucket3_local_frame = LocalFrameFunction(translation_vector,th3);

% relative_first_s_lo = first_s_lo - Bucket1_local_frame';
% relative_second_s_lo = second_s_lo - Bucket2_local_frame';
% relative_third_s_lo = third_s_lo - Bucket3_local_frame';
%%%%%%%%%%%Suction arrangement angle and geometrica vectors%%%%%%%%%%%
tower_tip = [0,0,eccentricity]';
tower_tip_fixed = transformation(tower_tip,Q,R_tilt);
angle_indicator = atan2d(tower_tip_fixed(2),tower_tip_fixed(1));
intensity_indicator = sqrt(tower_tip_fixed(1)^2+tower_tip_fixed(2)^2);
SUMz =0; 
count = count + 1;
ddd(count) = dd*180/pi;    
%%%%%%%%%%%

active_z_values1 = Zvalues(l_b,n,n_l,penetration1);
active_z_values2 = Zvalues(l_b,n,n_l,penetration2);
active_z_values3 = Zvalues(l_b,n,n_l,penetration3);
ee=0;
ee_tar = eccentricity;
nnn=0;
cc=0;
mmm=0;
QQQ = zeros(1,10000);
ccc =0;
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  Outer While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %% 
while  abs((ee)-ee_tar)>0.001*ee_tar %while condition to find the z-axis of the center of rotation
        tower_tip_fixed = transformation(tower_tip,Q,R_tilt);
angle_indicator = atan2d(tower_tip_fixed(2),tower_tip_fixed(1));
intensity_indicator = sqrt(tower_tip_fixed(1)^2+tower_tip_fixed(2)^2);
ccc = ccc+1;

 if ((ee)-ee_tar)>0.0
 mmm=mmm-0.015;%*abs(((ee)-ee_tar));
 else
 mmm=mmm+0.015;%*abs(((ee)-ee_tar));
 end
if cc==1
 Q = [nnn*cosd(angle_indicator), nnn*sind(angle_indicator), -(1.0*l_b)+mmm]';
else
 Q = [(nnn)*cosd(angle_indicator), nnn*sind(angle_indicator), -(1.0*l_b)+mmm]';
end
if rem(ccc,5)==0
fprintf('\rdegree: %0.3f  | ((ee)-ee_tar): %0.2f |(SUMz-Deadload): %0.3f | Qx: %0.3f| Qy: %0.2f | Qz: %0.3f', dd*180/pi, CON2,CON1,Q(1),Q(2),Q(3));  
end
QQQ(ccc+2)=Q(3);
if abs(abs(QQQ(ccc+2))-abs(QQQ(ccc)))==0
     break
end
count2 =0;
j=0;
fi1=FI_peak1;
fi2=FI_peak2;
fi3=FI_peak3;
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  inner While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %% 
 while abs(SUMz-Deadload)>(Deadload/1000) 
CON1 = abs(SUMz-Deadload); %track the condition
tower_tip = [0,0,eccentricity]';
tower_tip_fixed = transformation(tower_tip,Q,R_tilt);
angle_indicator = atan2d(tower_tip_fixed(2),tower_tip_fixed(1));
intensity_indicator = sqrt(tower_tip_fixed(1)^2+tower_tip_fixed(2)^2);
cc=cc+1;
if SUMz-Deadload>0
nnn=nnn+0.01;
else
nnn=nnn-0.01;
end
Q = [(0+nnn)*cosd(angle_indicator), (0+nnn)*sind(angle_indicator), Q(3)]';% iterative increasing or decreasing the x-axis of the rotation center
fprintf('\rdegree: %0.3f  | ((ee)-ee_tar): %0.3f |(Deadload - SUMz): %0.3f | Qx: %0.3f| Qy: %0.2f | Qz: %0.3f', dd*180/pi,((ee)-ee_tar) ,CON1,Q(1),Q(2),Q(3));

QQ(cc+2)=Q(1);
if abs(abs(QQ(cc+2))-abs(QQ(cc)))==0
     break
 end

[first_s_end,first_s_end_in,first_s_end_fixed,first_s_end_fixed_in,Pile_first_s_end,Pile_first_s_end_fixed,r_Pile_first_s_end,r_Pile_first_s_end_fixed]=points(R_split,th1,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[second_s_end,second_s_end_in,second_s_end_fixed,second_s_end_fixed_in,Pile_second_s_end,Pile_second_s_end_fixed,r_Pile_second_s_end,r_Pile_second_s_end_fixed]=points(R_split,th2,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[third_s_end,third_s_end_in,third_s_end_fixed,third_s_end_fixed_in,Pile_third_s_end,Pile_third_s_end_fixed,r_Pile_third_s_end, r_Pile_third_s_end_fixed]=points(R_split,th3,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
%% relative points calculation outer while
[relative_first_s_end,relative_first_s_end_fixed,relative_first_s_end_in,relative_first_s_end_fixed_in]...
    =relative_points(first_s_end, first_s_end_in, first_s_end_fixed, first_s_end_fixed_in, Pile_first_s_end, Pile_first_s_end_fixed, r_Pile_first_s_end, r_Pile_first_s_end_fixed, th1,l,dz,Q,R_tilt);

[relative_second_s_end,relative_second_s_end_fixed,relative_second_s_end_in,relative_second_s_end_fixed_in]...
    =relative_points(second_s_end, second_s_end_in, second_s_end_fixed, second_s_end_fixed_in, Pile_second_s_end, Pile_second_s_end_fixed, r_Pile_second_s_end, r_Pile_second_s_end_fixed, th2,l,dz,Q,R_tilt);

[relative_third_s_end,relative_third_s_end_fixed,relative_third_s_end_in,relative_third_s_end_fixed_in]...
    =relative_points(third_s_end, third_s_end_in, third_s_end_fixed, third_s_end_fixed_in, Pile_third_s_end, Pile_third_s_end_fixed, r_Pile_third_s_end, r_Pile_third_s_end_fixed, th3,l,dz,Q,R_tilt);


%% py and tz forces inner while
displacement1 = (first_s_end(3, :) - first_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement2 = (second_s_end(3, :) - second_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement3 = (third_s_end(3, :) - third_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm

displacement_in1 = (first_s_end_in(3, :) - first_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in2 = (second_s_end_in(3, :) - second_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in3 = (third_s_end_in(3, :) - third_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm

[~,~,P_i_out1] = coeffinder_LD123(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1); %KN
[~,~,P_i_out2] = coeffinder_LD123(second_s_end,second_s_end_fixed,abs(second_s_end_fixed(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2); %KN
[~,~,P_i_out3] = coeffinder_LD123(third_s_end,third_s_end_fixed,abs(third_s_end_fixed(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3); %KN
P_i_out1(mask1) = 0;
P_i_out2(mask2) = 0;
P_i_out3(mask3) = 0;

[~,tz_out1]=tzcurve_noOCR(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out2]=tzcurve_noOCR(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out3]=tzcurve_noOCR(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);

[tz_in1,~]=tzcurve_noOCR(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement_in1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in2,~]=tzcurve_noOCR(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement_in2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in3,~]=tzcurve_noOCR(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement_in3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);

tz_out1(mask1) = 0;
tz_out2(mask2) = 0;
tz_out3(mask3) = 0;

tz_in1(mask1) = 0;
tz_in2(mask2) = 0;
tz_in3(mask3) = 0;

KKz_in1 = tz_in1;
KKz_in2 = tz_in2;
KKz_in3 = tz_in3;

KKz_out1 = tz_out1;
KKz_out2 = tz_out2;
KKz_out3 = tz_out3;

KKz1 =KKz_in1+KKz_out1;
KKz2 =KKz_in2+KKz_out2;
KKz3 =KKz_in3+KKz_out3;

KKZ1(count,:) = tz_out1+tz_in1;
KKZ2(count,:) = tz_out2+tz_in2;
KKZ3(count,:) = tz_out3+tz_in3;
%% qz forces inner while
[qb1_1] = qzcurve2(fi_end1,Pile_first_s_end,Pile_first_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_2] = qzcurve2(fi_end2,Pile_second_s_end,Pile_second_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_3] = qzcurve2(fi_end3,Pile_third_s_end,Pile_third_s_end_fixed,gamma,n,thickness,di_b,do_b);

[qb2_1] = qzcurve4(fi_end1,Z_lid,r_Pile_first_s_end_fixed, r_Pile_first_s_end_fixed,gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_2] = qzcurve4(fi_end2,Z_lid,r_Pile_second_s_end,r_Pile_second_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_3] = qzcurve4(fi_end3,Z_lid,r_Pile_third_s_end,r_Pile_third_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);

SUMz= sum(KKz1)+sum(qb1_1)+sum(qb2_1)+...
             sum(KKz2)+sum(qb1_2)+sum(qb2_2)+...
             sum(KKz3)+sum(qb1_3)+sum(qb2_3);
         

 end
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  End of inner While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %% 
%% points  calculation outer while
[first_s_end,first_s_end_in,first_s_end_fixed,first_s_end_fixed_in,Pile_first_s_end,Pile_first_s_end_fixed,r_Pile_first_s_end,r_Pile_first_s_end_fixed]=points(R_split,th1,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[second_s_end,second_s_end_in,second_s_end_fixed,second_s_end_fixed_in,Pile_second_s_end,Pile_second_s_end_fixed,r_Pile_second_s_end,r_Pile_second_s_end_fixed]=points(R_split,th2,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[third_s_end,third_s_end_in,third_s_end_fixed,third_s_end_fixed_in,Pile_third_s_end,Pile_third_s_end_fixed,r_Pile_third_s_end,r_Pile_third_s_end_fixed]=points(R_split,th3,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
%% relative points calculation outer while
% [relative_first_s_end,relative_first_s_end_fixed,relative_first_s_end_in,relative_first_s_end_fixed_in]...
%     =relative_points(first_s_end, first_s_end_in, first_s_end_fixed, first_s_end_fixed_in, Pile_first_s_end, Pile_first_s_end_fixed, r_Pile_first_s_end, r_Pile_first_s_end_fixed, th1,l,dz,Q,R_tilt);
% 
% [relative_second_s_end,relative_second_s_end_fixed,relative_second_s_end_in,relative_second_s_end_fixed_in]...
%     =relative_points(second_s_end, second_s_end_in, second_s_end_fixed, second_s_end_fixed_in, Pile_second_s_end, Pile_second_s_end_fixed, r_Pile_second_s_end, r_Pile_second_s_end_fixed, th2,l,dz,Q,R_tilt);
% 
% [relative_third_s_end,relative_third_s_end_fixed,relative_third_s_end_in,relative_third_s_end_fixed_in]...
%     =relative_points(third_s_end, third_s_end_in, third_s_end_fixed, third_s_end_fixed_in, Pile_third_s_end, Pile_third_s_end_fixed, r_Pile_third_s_end, r_Pile_third_s_end_fixed, th3,l,dz,Q,R_tilt);

%% py and tz forces outer while
% Group effect deactivated (f_g defaults to 1.0 inside coeffinder_LD123)
% fg_vec = detect_fg_per_bucket( ...
%     {first_s_end,second_s_end,third_s_end}, ...
%     {first_s_end_fixed,second_s_end_fixed,third_s_end_fixed}, ...
%     adj_matrix, SD_ratio, delta_ref, l_b/do_b, 'tri');
[KKx1,KKy1,~] = coeffinder_LD123(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1); %KN
[KKx2,KKy2,~] = coeffinder_LD123(second_s_end,second_s_end_fixed,abs(second_s_end_fixed(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2); %KN
[KKx3,KKy3,~] = coeffinder_LD123(third_s_end,third_s_end_fixed,abs(third_s_end_fixed(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3); %KN

KKx1(mask1) = 0;
KKy1(mask1) = 0;
KKx2(mask2) = 0;
KKy2(mask2) = 0;
KKx3(mask3) = 0;
KKy3(mask3) = 0;

SUMx = sum(KKx1)+sum(KKx2)+sum(KKx3);
SUMy = sum(KKy1)+sum(KKy2)+sum(KKy3);

displacement1 = (first_s_end(3, :) - first_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement2 = (second_s_end(3, :) - second_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement3 = (third_s_end(3, :) - third_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm

displacement_in1 = (first_s_end_in(3, :) - first_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in2 = (second_s_end_in(3, :) - second_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in3 = (third_s_end_in(3, :) - third_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm

[~,~,P_i_out1] = coeffinder_LD123(first_s_end,first_s_end_fixed,abs(first_s_end(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1); %KN
[~,~,P_i_out2] = coeffinder_LD123(second_s_end,second_s_end_fixed,abs(second_s_end(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2); %KN
[~,~,P_i_out3] = coeffinder_LD123(third_s_end,third_s_end_fixed,abs(third_s_end(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3); %KN
P_i_out1(mask1) = 0;
P_i_out2(mask2) = 0;
P_i_out3(mask3) = 0;

[~,tz_out1]=tzcurve_noOCR(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out2]=tzcurve_noOCR(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out3]=tzcurve_noOCR(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);

[tz_in1,~]=tzcurve_noOCR(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in2,~]=tzcurve_noOCR(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in3,~]=tzcurve_noOCR(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);


tz_out1(mask1) = 0;
tz_out2(mask2) = 0;
tz_out3(mask3) = 0;

tz_in1(mask1) = 0;
tz_in2(mask2) = 0;
tz_in3(mask3) = 0;

KKz_in1 = tz_in1;
KKz_in2 = tz_in2;
KKz_in3 = tz_in3;

KKz_out1 = tz_out1;
KKz_out2 = tz_out2;
KKz_out3 = tz_out3;

KKz1 =KKz_in1+KKz_out1;
KKz2 =KKz_in2+KKz_out2;
KKz3 =KKz_in3+KKz_out3;

KKZ1(count,:) = tz_out1+tz_in1;
KKZ2(count,:) = tz_out2+tz_in2;
KKZ3(count,:) = tz_out3+tz_in3;
%% qz forces outer while
[qb1_1] = qzcurve2(fi_end1,Pile_first_s_end,Pile_first_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_2] = qzcurve2(fi_end2,Pile_second_s_end,Pile_second_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_3] = qzcurve2(fi_end3,Pile_third_s_end,Pile_third_s_end_fixed,gamma,n,thickness,di_b,do_b);

[qb2_1] = qzcurve4(fi_end1,Z_lid,r_Pile_first_s_end,r_Pile_first_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_2] = qzcurve4(fi_end2,Z_lid,r_Pile_second_s_end,r_Pile_second_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_3] = qzcurve4(fi_end3,Z_lid,r_Pile_third_s_end,r_Pile_third_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);

SUMz= sum(KKz1)+sum(qb1_1)+sum(qb2_1)+...
             sum(KKz2)+sum(qb1_2)+sum(qb2_2)+...
             sum(KKz3)+sum(qb1_3)+sum(qb2_3);
%% py moments outer while
forces_py1 = [KKx1;KKy1;0*KKz1]';
forces_py2 = [KKx2;KKy2;0*KKz2]';
forces_py3 = [KKx3;KKy3;0*KKz3]';

positions1 = first_s_end_fixed';
positions2 = second_s_end_fixed';
positions3 = third_s_end_fixed';

moments_py1 = cross(positions1, forces_py1);
moments_py2 = cross(positions2, forces_py2);
moments_py3 = cross(positions3, forces_py3);

total_moment_py1 = sum(moments_py1);
total_moment_py2 = sum(moments_py2);
total_moment_py3 = sum(moments_py3);

total_moment_py_vec = total_moment_py1+total_moment_py2+total_moment_py3;
total_moment_py = sqrt(total_moment_py_vec(1)^2+total_moment_py_vec(2)^2+total_moment_py_vec(3)^2);

for i=1:n_l
    pc(i)= n*i;
end

%% tz moments outer while
forces_tz_out1 = [0*KKx1;0*KKy1;tz_out1]';
forces_tz_out2 = [0*KKx2;0*KKy2;tz_out2]';
forces_tz_out3 = [0*KKx3;0*KKy3;tz_out3]';

forces_tz_in1 = [0*KKx1;0*KKy1;tz_in1]';
forces_tz_in2 = [0*KKx2;0*KKy2;tz_in2]';
forces_tz_in3 = [0*KKx3;0*KKy3;tz_in3]';

positions_out1 = first_s_end_fixed';
positions_out2 = second_s_end_fixed';
positions_out3 = third_s_end_fixed';

positions_in1 = first_s_end_fixed_in';
positions_in2 = second_s_end_fixed_in';
positions_in3 = third_s_end_fixed_in';

moments_tz_out1 = cross(positions_out1, forces_tz_out1);
moments_tz_out2 = cross(positions_out2, forces_tz_out2);
moments_tz_out3 = cross(positions_out3, forces_tz_out3);

moments_tz_in1 = cross(positions_in1, forces_tz_in1);
moments_tz_in2= cross(positions_in2, forces_tz_in2);
moments_tz_in3 = cross(positions_in3, forces_tz_in3);

total_moment_tz1 = sum(moments_tz_in1)+sum(moments_tz_out1);
total_moment_tz2 = sum(moments_tz_in2)+sum(moments_tz_out2);
total_moment_tz3 = sum(moments_tz_in3)+sum(moments_tz_out3);

total_moment_tz_vec = total_moment_tz1+total_moment_tz2+total_moment_tz3;
total_moment_tz = sqrt(total_moment_tz_vec(1)^2+total_moment_tz_vec(2)^2+total_moment_tz_vec(3)^2);
%% qz lid moments outer while

positions_pile1 = Pile_first_s_end_fixed';
positions_pile2 = Pile_second_s_end_fixed';
positions_pile3 = Pile_third_s_end_fixed';

XX = zeros(size(qb1_1));
YY = zeros(size(qb1_1));

force_vectors1 = [XX;YY;qb1_1]';
force_vectors2 = [XX;YY;qb1_2]';
force_vectors3 = [XX;YY;qb1_3]';

moments_qz_1 = cross(positions_pile1, force_vectors1);
moments_qz_2 = cross(positions_pile2, force_vectors2);
moments_qz_3 = cross(positions_pile3, force_vectors3);

total_moment_qz_1 = sum(moments_qz_1);
total_moment_qz_2 = sum(moments_qz_2);
total_moment_qz_3 = sum(moments_qz_3);

total_moment_qz_tip_vec = total_moment_qz_1+total_moment_qz_2+total_moment_qz_3;
total_moment_tip = sqrt(total_moment_qz_tip_vec(1)^2+total_moment_qz_tip_vec(2)^2+total_moment_qz_tip_vec(3)^2);

%% qz tip moments outer while
r_positions_pile1 = r_Pile_first_s_end_fixed';
r_positions_pile2 = r_Pile_second_s_end_fixed';
r_positions_pile3 = r_Pile_third_s_end_fixed';

XXX = zeros(size(qb2_1));
YYY = zeros(size(qb2_1));

force_vectors_2_1 = [XXX;YYY;qb2_1]';
force_vectors_2_2 = [XXX;YYY;qb2_2]';
force_vectors_2_3 = [XXX;YYY;qb2_3]';

moments_qz_2_1 = cross(r_positions_pile1, force_vectors_2_1);
moments_qz_2_2 = cross(r_positions_pile2, force_vectors_2_2);
moments_qz_2_3 = cross(r_positions_pile3, force_vectors_2_3);

total_moment_qz_2_1 = sum(moments_qz_2_1);
total_moment_qz_2_2 = sum(moments_qz_2_2);
total_moment_qz_2_3 = sum(moments_qz_2_3);

total_moment_qz_lid_vec = total_moment_qz_2_1+total_moment_qz_2_2+total_moment_qz_2_3;
total_moment_lid = sqrt(total_moment_qz_lid_vec(1)^2+total_moment_qz_lid_vec(2)^2+total_moment_qz_lid_vec(3)^2);
%% outer While condition calculation
MOMENT =total_moment_lid + total_moment_tip + total_moment_tz + total_moment_py;

F_Hx = SUMx;
F_Hy = SUMy;
F_HH  = sqrt (F_Hx^2+F_Hy^2);
ee = abs(MOMENT/F_HH);

CON2=abs((ee)-ee_tar); 

end 
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  End of outer While  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ %% 

j=0;
%   Q=[0,0,0]';
fi1=FI_peak1;
fi2=FI_peak2;
fi3=FI_peak3;
%% points  calculation outer while
[first_s_end,first_s_end_in,first_s_end_fixed,first_s_end_fixed_in,Pile_first_s_end,Pile_first_s_end_fixed,r_Pile_first_s_end,r_Piles_first_s_end_fixed]=points(R_split,th1,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[second_s_end,second_s_end_in,second_s_end_fixed,second_s_end_fixed_in,Pile_second_s_end,Pile_second_s_end_fixed,r_Pile_second_s_end,r_Piles_second_s_end_fixed]=points(R_split,th2,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
[third_s_end,third_s_end_in,third_s_end_fixed,third_s_end_fixed_in,Pile_third_s_end,Pile_third_s_end_fixed,r_Pile_third_s_end,r_Piles_third_s_end_fixed]=points(R_split,th3,dz,n,n_l,Q,R_tilt,Z_lid,l_b,l,ro_b,ri_b,r_pile);
%% relative points calculation final
% [relative_first_s_end,relative_first_s_end_fixed,relative_first_s_end_in,relative_first_s_end_fixed_in]...
%     =relative_points(first_s_end, first_s_end_in, first_s_end_fixed, first_s_end_fixed_in, Pile_first_s_end, Pile_first_s_end_fixed, r_Pile_first_s_end, r_Pile_first_s_end_fixed, th1,l,dz,Q,R_tilt);
% 
% [relative_second_s_end,relative_second_s_end_fixed,relative_second_s_end_in,relative_second_s_end_fixed_in]...
%     =relative_points(second_s_end, second_s_end_in, second_s_end_fixed, second_s_end_fixed_in, Pile_second_s_end, Pile_second_s_end_fixed, r_Pile_second_s_end, r_Pile_second_s_end_fixed, th2,l,dz,Q,R_tilt);
% 
% [relative_third_s_end,relative_third_s_end_fixed,relative_third_s_end_in,relative_third_s_end_fixed_in]...
%     =relative_points(third_s_end, third_s_end_in, third_s_end_fixed, third_s_end_fixed_in, Pile_third_s_end, Pile_third_s_end_fixed, r_Pile_third_s_end, r_Pile_third_s_end_fixed, th3,l,dz,Q,R_tilt);

%% py and tz forces final
% Group effect deactivated (f_g defaults to 1.0 inside coeffinder_LD123)
% fg_vec = detect_fg_per_bucket( ...
%     {first_s_end,second_s_end,third_s_end}, ...
%     {first_s_end_fixed,second_s_end_fixed,third_s_end_fixed}, ...
%     adj_matrix, SD_ratio, delta_ref, l_b/do_b, 'tri');
[KKx1,KKy1,~] = coeffinder_LD123(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1); %KN
[KKx2,KKy2,~] = coeffinder_LD123(second_s_end,second_s_end_fixed,abs(second_s_end_fixed(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2); %KN
[KKx3,KKy3,~] = coeffinder_LD123(third_s_end,third_s_end_fixed,abs(third_s_end_fixed(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3); %KN

KKx1(mask1) = 0;
KKy1(mask1) = 0;
KKx2(mask2) = 0;
KKy2(mask2) = 0;
KKx3(mask3) = 0;
KKy3(mask3) = 0;

% KKX1(count,:) = KKx1;
% KKY1(count,:) = KKy1;
SUMx1= sum(KKx1);
SUMx2 = sum(KKx2);
SUMx3 = sum(KKx3);
SUMx = SUMx1+SUMx2+SUMx3;

SUMy1= sum(KKy1);
SUMy2 = sum(KKy2);
SUMy3 = sum(KKy3);
SUMy = SUMy1+SUMy2+SUMy3;

displacement1 = (first_s_end(3, :) - first_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement2 = (second_s_end(3, :) - second_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm
displacement3 = (third_s_end(3, :) - third_s_end_fixed(3, :)) * 1e3; %1e3 is to make it to mm

Reshaped_first_s_end = reshape(first_s_end(3, :),[n,n_l]);
Reshaped_second_s_end = reshape(second_s_end(3, :),[n,n_l]);
Reshaped_third_s_end = reshape(third_s_end(3, :),[n,n_l]);

Reshaped_first_s_end_fixed = reshape(first_s_end_fixed(3, :) ,[n,n_l]);
Reshaped_second_s_end_fixed  = reshape(second_s_end_fixed(3, :) ,[n,n_l]);
Reshaped_third_s_end_fixed  = reshape(third_s_end_fixed(3, :) ,[n,n_l]);

Reshaped_displacement1 = reshape(displacement1,[n,n_l]);
Reshaped_displacement2 = reshape(displacement2,[n,n_l]);
Reshaped_displacement3 = reshape(displacement3,[n,n_l]);

displacement_in1 = (first_s_end_in(3, :) - first_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in2 = (second_s_end_in(3, :) - second_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm
displacement_in3 = (third_s_end_in(3, :) - third_s_end_fixed_in(3, :)) * 1e3; %1e3 is to make it to mm

[~,~,P_i_out1] = coeffinder_LD123(first_s_end,first_s_end_fixed,abs(first_s_end_fixed(3, :)),l_b,n_l,fi1,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th1); %KN
[~,~,P_i_out2] = coeffinder_LD123(second_s_end,second_s_end_fixed,abs(second_s_end_fixed(3, :)),l_b,n_l,fi2,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th2); %KN
[~,~,P_i_out3] = coeffinder_LD123(third_s_end,third_s_end_fixed,abs(third_s_end_fixed(3, :)),l_b,n_l,fi3,n,do_b,gamma,abs(Q(3)),A_c,m_c,translation_vector,th3); %KN
P_i_out1(mask1) = 0;
P_i_out2(mask2) = 0;
P_i_out3(mask3) = 0;

[~,tz_out1]=tzcurve_noOCR(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out2]=tzcurve_noOCR(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[~,tz_out3]=tzcurve_noOCR(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);

[tz_in1,~]=tzcurve_noOCR(P_i_out1, fi1, gamma, do_b, di_b, (first_s_end_fixed(3, :)), displacement_in1, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in2,~]=tzcurve_noOCR(P_i_out2, fi2, gamma, do_b, di_b, (second_s_end_fixed(3, :)), displacement_in2, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);
[tz_in3,~]=tzcurve_noOCR(P_i_out3, fi3, gamma, do_b, di_b, (third_s_end_fixed(3, :)), displacement_in3, l_b, n_l, n, fi_crit,delta,D_r,e_max,e_min);

tz_out1(mask1) = 0;
tz_out2(mask2) = 0;
tz_out3(mask3) = 0;

tz_in1(mask1) = 0;
tz_in2(mask2) = 0;
tz_in3(mask3) = 0;

Reshaped_tz_out1 = reshape(tz_out1,[n,n_l]);
Reshaped_tz_out2 = reshape(tz_out2,[n,n_l]);
Reshaped_tz_out3 = reshape(tz_out3,[n,n_l]);
Reshaped_tz_in1 = reshape(tz_in1,[n,n_l]);
Reshaped_tz_in2 = reshape(tz_in2,[n,n_l]);
Reshaped_tz_in3 = reshape(tz_in3,[n,n_l]);

KKz_in1 = tz_in1;
KKz_in2 = tz_in2;
KKz_in3 = tz_in3;

KKz_out1 = tz_out1;
KKz_out2 = tz_out2;
KKz_out3 = tz_out3;



KKz1 =KKz_in1+KKz_out1;
KKz2 =KKz_in2+KKz_out2;
KKz3 =KKz_in3+KKz_out3; 

KKZ1(count,:) = tz_out1+tz_in1;
KKZ2(count,:) = tz_out2+tz_in2;
KKZ3(count,:) = tz_out3+tz_in3;
%% qz forces final
[qb1_1] = qzcurve2(fi_end1,Pile_first_s_end,Pile_first_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_2] = qzcurve2(fi_end2,Pile_second_s_end,Pile_second_s_end_fixed,gamma,n,thickness,di_b,do_b);
[qb1_3] = qzcurve2(fi_end3,Pile_third_s_end,Pile_third_s_end_fixed,gamma,n,thickness,di_b,do_b);

[qb2_1] = qzcurve4(fi_end1,Z_lid,r_Pile_first_s_end, r_Pile_first_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_2] = qzcurve4(fi_end2,Z_lid,r_Pile_second_s_end, r_Pile_second_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);
[qb2_3] = qzcurve4(fi_end3,Z_lid,r_Pile_third_s_end, r_Pile_third_s_end_fixed, gamma, n, n_q, R_split,di_b,D_r,S_gamma);

SUMz= sum(KKz1)+sum(qb1_1)+sum(qb2_1)+...
             sum(KKz2)+sum(qb1_2)+sum(qb2_2)+...
             sum(KKz3)+sum(qb1_3)+sum(qb2_3);
%% py moments final
forces_py1 = [KKx1;KKy1;0*KKz1]';
forces_py2 = [KKx2;KKy2;0*KKz2]';
forces_py3 = [KKx3;KKy3;0*KKz3]';

positions1 = first_s_end_fixed';
positions2 = second_s_end_fixed';
positions3 = third_s_end_fixed';

moments_py1 = cross(positions1, forces_py1);
moments_py2 = cross(positions2, forces_py2);
moments_py3 = cross(positions3, forces_py3);

total_moment_py1 = sum(moments_py1);
total_moment_py2 = sum(moments_py2);
total_moment_py3 = sum(moments_py3);

total_moment_py_vec = total_moment_py1+total_moment_py2+total_moment_py3;
total_moment_py = sqrt(total_moment_py_vec(1)^2+total_moment_py_vec(2)^2+total_moment_py_vec(3)^2);

for i=1:n_l
    pc(i)= n*i;
end

%% tz moments final
forces_tz_out1 = [0*KKx1;0*KKy1;tz_out1]';
forces_tz_out2 = [0*KKx2;0*KKy2;tz_out2]';
forces_tz_out3 = [0*KKx3;0*KKy3;tz_out3]';

forces_tz_in1 = [0*KKx1;0*KKy1;tz_in1]';
forces_tz_in2 = [0*KKx2;0*KKy2;tz_in2]';
forces_tz_in3 = [0*KKx3;0*KKy3;tz_in3]';

positions_out1 = first_s_end_fixed';
positions_out2 = second_s_end_fixed';
positions_out3 = third_s_end_fixed';

positions_in1 = first_s_end_fixed_in';
positions_in2 = second_s_end_fixed_in';
positions_in3 = third_s_end_fixed_in';

moments_tz_out1 = cross(positions_out1, forces_tz_out1);
moments_tz_out2 = cross(positions_out2, forces_tz_out2);
moments_tz_out3 = cross(positions_out3, forces_tz_out3);

moments_tz_in1 = cross(positions_in1, forces_tz_in1);
moments_tz_in2= cross(positions_in2, forces_tz_in2);
moments_tz_in3 = cross(positions_in3, forces_tz_in3);

total_moment_tz1 = sum(moments_tz_in1)+sum(moments_tz_out1);
total_moment_tz2 = sum(moments_tz_in2)+sum(moments_tz_out2);
total_moment_tz3 = sum(moments_tz_in3)+sum(moments_tz_out3);

total_moment_tz_vec = total_moment_tz1+total_moment_tz2+total_moment_tz3;
total_moment_tz = sqrt(total_moment_tz_vec(1)^2+total_moment_tz_vec(2)^2+total_moment_tz_vec(3)^2);
%% qz lid moments final

positions_pile1 = Pile_first_s_end_fixed';
positions_pile2 = Pile_second_s_end_fixed';
positions_pile3 = Pile_third_s_end_fixed';

XX = zeros(size(qb1_1));
YY = zeros(size(qb1_1));

force_vectors1 = [XX;YY;qb1_1]';
force_vectors2 = [XX;YY;qb1_2]';
force_vectors3 = [XX;YY;qb1_3]';

moments_qz_1 = cross(positions_pile1, force_vectors1);
moments_qz_2 = cross(positions_pile2, force_vectors2);
moments_qz_3 = cross(positions_pile3, force_vectors3);

total_moment_qz_1 = sum(moments_qz_1);
total_moment_qz_2 = sum(moments_qz_2);
total_moment_qz_3 = sum(moments_qz_3);

total_moment_qz_tip_vec = total_moment_qz_1+total_moment_qz_2+total_moment_qz_3;
total_moment_qz_tip = sqrt(total_moment_qz_tip_vec(1)^2+total_moment_qz_tip_vec(2)^2+total_moment_qz_tip_vec(3)^2);
%% qz tip moments final
r_positions_pile1 = r_Pile_first_s_end_fixed';
r_positions_pile2 = r_Pile_second_s_end_fixed';
r_positions_pile3 = r_Pile_third_s_end_fixed';

XXX = zeros(size(qb2_1));
YYY = zeros(size(qb2_1));

force_vectors_2_1 = [XXX;YYY;qb2_1]';
force_vectors_2_2 = [XXX;YYY;qb2_2]';
force_vectors_2_3 = [XXX;YYY;qb2_3]';

moments_qz_2_1 = cross(r_positions_pile1, force_vectors_2_1);
moments_qz_2_2 = cross(r_positions_pile2, force_vectors_2_2);
moments_qz_2_3 = cross(r_positions_pile3, force_vectors_2_3);

total_moment_qz_2_1 = sum(moments_qz_2_1);
total_moment_qz_2_2 = sum(moments_qz_2_2);
total_moment_qz_2_3 = sum(moments_qz_2_3);

total_moment_qz_lid_vec = total_moment_qz_2_1+total_moment_qz_2_2+total_moment_qz_2_3;
total_moment_qz_lid = sqrt(total_moment_qz_lid_vec(1)^2+total_moment_qz_lid_vec(2)^2+total_moment_qz_lid_vec(3)^2);
%% Moment Calculation final
MOMENT =total_moment_qz_tip+ total_moment_qz_lid + total_moment_tz + total_moment_py;

MOMENTT(count) =total_moment_qz_tip+ total_moment_qz_lid + total_moment_tz + total_moment_py;

F_Hx = sum(KKx1)+sum(KKx2)+sum(KKx3);
F_Hy = sum(KKy1)+sum(KKy2)+sum(KKy3);
F_HH  = sqrt (F_Hx^2+F_Hy^2);
M_1(count)  = F_HH  *eccentricity;
M_2(count)  = MOMENT;

%% --- Store data for Fig 19, 20, 21 per bucket ---
% Fig 19: per-bucket moment contribution percentages
M_py_b  = [norm(total_moment_py1), norm(total_moment_py2), norm(total_moment_py3)];
M_tz_b  = [norm(total_moment_tz1), norm(total_moment_tz2), norm(total_moment_tz3)];
M_tip_b = [norm(total_moment_qz_1), norm(total_moment_qz_2), norm(total_moment_qz_3)];
M_lid_b = [norm(total_moment_qz_2_1), norm(total_moment_qz_2_2), norm(total_moment_qz_2_3)];
for b_idx = 1:3
    Mtot_b = M_py_b(b_idx) + M_tz_b(b_idx) + M_tip_b(b_idx) + M_lid_b(b_idx);
    if Mtot_b > 0
        M_py_pct(count, b_idx)    = M_py_b(b_idx)  / Mtot_b * 100;
        M_tz_pct(count, b_idx)    = M_tz_b(b_idx)  / Mtot_b * 100;
        M_qztip_pct(count, b_idx) = M_tip_b(b_idx) / Mtot_b * 100;
        M_qzlid_pct(count, b_idx) = M_lid_b(b_idx) / Mtot_b * 100;
    end
end
% Fig 20: per-bucket t-z right/left forces (normalized)
tz_out_cell = {tz_out1, tz_out2, tz_out3};
tz_in_cell  = {tz_in1,  tz_in2,  tz_in3};
for b_idx = 1:3
    rm = right_mask_cell{b_idx}; lm = left_mask_cell{b_idx};
    tz_right_out_store(count,b_idx) = sum(tz_out_cell{b_idx}(rm)) / (gamma * do_b^3);
    tz_right_in_store(count,b_idx)  = sum(tz_in_cell{b_idx}(rm))  / (gamma * do_b^3);
    tz_left_out_store(count,b_idx)  = sum(tz_out_cell{b_idx}(lm)) / (gamma * do_b^3);
    tz_left_in_store(count,b_idx)   = sum(tz_in_cell{b_idx}(lm))  / (gamma * do_b^3);
end
% Fig 21: store p-y forces per step
KKx_store1(count,:) = KKx1;
KKx_store2(count,:) = KKx2;
KKx_store3(count,:) = KKx3;
% Store per-step data for animation (used by animate_3d_heatmaps.m)
anim_s_fixed{count}    = {first_s_end_fixed, second_s_end_fixed, third_s_end_fixed};
anim_s_fixed_in{count} = {first_s_end_fixed_in, second_s_end_fixed_in, third_s_end_fixed_in};
anim_KKx{count} = {KKx1, KKx2, KKx3};
anim_KKy{count} = {KKy1, KKy2, KKy3};
anim_KKz{count} = {KKz1, KKz2, KKz3};
anim_KKz_out{count} = {KKz_out1, KKz_out2, KKz_out3};
anim_KKz_in{count}  = {KKz_in1, KKz_in2, KKz_in3};

colors2 = [
    0.111, 0.111, 0.111;    % Dark Grey
    0.8500, 0.3250, 0.0980; % Orange
    0.4940, 0.1840, 0.5560; % Purple
    0.4660, 0.6740, 0.1880; % Green
    0.3010, 0.7450, 0.9330; % Light Blue
    0.9290, 0.6940, 0.1250; % Yellow
    0.6350, 0.0780, 0.1840  % Dark Red
];

figure(4)
plot(dd*180/pi,M_1(count),'*r','LineWidth',5)
hold on
plot(dd*180/pi,M_2(count),'*b','LineWidth',5)

grid on
hold on
set(gca,'TickLabelInterpreter','latex');
set(gca,'fontweight','bold','fontsize',22)
xlabel('Rotation (degree)', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
ylabel('Moment Load (MNm)', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
figure(107)
scatter3(Q(1),Q(2),Q(3),'o','LineWidth',8)
xlabel('X', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
ylabel('Y', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
ylabel('Z', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')

hold on
% inner_right(count) = (l_b/n_l)*(pi*do_b/n)*sum(KKz1_in(pc)+KKz1_in(pc-11)+KKz1_in(pc-10)+KKz1_in(pc-1)+KKz1_in(pc-2))/gamma/do_b^3;
% inner_left(count) = (l_b/n_l)*(pi*do_b/n)*sum(KKz1_in(pc-8)+KKz1_in(pc-7)+KKz1_in(pc-6)+KKz1_in(pc-5)+KKz1_in(pc-4))/gamma/do_b^3;
% outer_left(count) = (l_b/n_l)*(pi*do_b/n)*sum(KKz1_out(pc-8)+KKz1_out(pc-7)+KKz1_out(pc-6)+KKz1_out(pc-5)+KKz1_out(pc-4))/gamma/do_b^3;
% outer_right(count) = (l_b/n_l)*(pi*do_b/n)*sum(KKz1_out(pc)+KKz1_out(pc-11)+KKz1_out(pc-10)+KKz1_out(pc-1)+KKz1_out(pc-2))/gamma/do_b^3;
% figure(1)
% plot(dd*180/pi,inner_right(count), '>b', 'LineWidth', 6); %inner right
hold on
% plot(dd*180/pi,inner_left(count), '>r', 'LineWidth', 6); %inner left
% hold on
% %% 
% plot(dd*180/pi,outer_left(count), 'or', 'LineWidth', 6); %outer left
% hold on
% %%
% plot(dd*180/pi,outer_right(count), 'ob', 'LineWidth', 6); %outer right
% hold on
% set(gca, 'XScale', 'log');  % Set the x-axis to logarithmic scale
% hold on;
grid on;
xlabel('Angle (degrees)','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
ylabel('Skirt wall friction (F_t_z/\gamma^\prime/D^3) ','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontweight','bold','fontsize',22)
% set(gca, 'YAxisLocation', 'left');
legend({'Inner right side friction','Inner left side friction','Outer left side friction','Outer right side friction'},'fontsize',16,'interpreter','latex','location','northeast')

 CENTE_OF_ROTATION(:,count) = Q;

% figure(3)
%  plot(dd*180/pi,total_moment_tz/total_moment(count)*100,'-*b','LineWidth',5)
%   plot(dd*180/pi,total_moment_py/total_moment(count)*100,'-*r','LineWidth',5)
%     plot(dd*180/pi,2*(total_moment_qz_1+total_moment_qz_2)/total_moment(count)*100,'-*k','LineWidth',5)
% %       plot(dd*180/pi,total_moments_qz_2(2)/total_moment(count)*100,'-*m','LineWidth',5)   
%  hold on
% % grid on
% set(gca,'TickLabelInterpreter','latex');
% % set(gca,'fontweight','bold','fontsize',22)
% % set(gca, 'XScale', 'log');  % Set the x-axis to logarithmic scale
% hold on;
% xlabel('Rotation Angel (degree)', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
% ylabel('share moment percentage', 'FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman')
% % xlim([0.1, 0.5]);
% figure(11)
% HH(count) = sum((KKx1));
% VV(count) = sum(qb1_1)+sum(qb2_1)+sum(sqrt(KKz1.^2));
% MM(count) = total_moment(count);
% e = eccentricity;
% SUMx = sum(KKx1);
% SUMy = sum(KKy1);
% plot(dd*do_b,u(count),'-ob','LineWidth',5)
% hold on
% cc=0;
% figure(99)
% for e = [5,10,20,40,70,100]
%     cc=cc+1;
% colors1 = [0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.2050, 0.1111, 0.7880; 0.3010, 0.6740, 0.2880; 0.1010, 0.9740, 0.8888];
% plot(dd*180/pi,-total_moment(count)/(abs(Q(3))+e)*e/1000,'-o','LineWidth',4,'Color', colors1(cc,:))
% hold on
% grid on
% 
% end
% xlabel('Rotation(\theta)','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
% ylabel('Moment (MN.m)','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontweight','bold','fontsize',22)
% legend({'h=5','h=10','h=20','h=40','h=70','h=100'},'fontsize',16,'interpreter','latex','location','northeast')
% figure(100)
% cc=0;
% for e = eccentricity
%     cc=cc+1;
% colors1 = [0.011, 0.115, 0.172; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.2050, 0.1111, 0.7880; 0.3010, 0.6740, 0.2880; 0.1010, 0.9740, 0.8888];
% plot(dd*180/pi,F_H(count)/1000,'-o','LineWidth',4,'Color', colors1(cc,:))
% hold on
% grid on
% end
% legend({'h=0','h=5','h=10','h=20','h=40','h=70','h=100'},'fontsize',16,'interpreter','latex','location','northeast')
% xlabel('Rotation(\theta)','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
% ylabel('Horizontal Load(H)','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontweight','bold','fontsize',22)
% 
% figure(101)
% cc=0;
% for e = eccentricity
%     cc=cc+1;
% colors1 = [0.011, 0.115, 0.172; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.2050, 0.1111, 0.7880; 0.3010, 0.6740, 0.2880; 0.1010, 0.9740, 0.8888];
% plot(dd*180/pi,total_moment(count)/(abs(Q(3))+e)/1000,'-o','LineWidth',4,'Color', colors1(cc,:))
% hold on
% grid on
% end
% legend({'h=0','h=5','h=10','h=20','h=40','h=70','h=100'},'fontsize',16,'interpreter','latex','location','northeast')
% xlabel('Rotation(\theta)','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
% ylabel('Horizontal Load(H)','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontweight','bold','fontsize',22)
% colors2 = [
%     0.111, 0.111, 0.111;    % Dark Grey
%     0.8500, 0.3250, 0.0980; % Orange
%     0.4940, 0.1840, 0.5560; % Purple
%     0.4660, 0.6740, 0.1880; % Green
%     0.3010, 0.7450, 0.9330; % Light Blue
%     0.9290, 0.6940, 0.1250; % Yellow
%     0.6350, 0.0780, 0.1840  % Dark Red
% ];
% 
% total_moment_py_count(count) = sum(total_moment_py)
% total_moment_tz_count(count) = sum(total_moment_tz)
end
M_2_f =M_2'/1000;
% 
% %% ========== Fig 19 equivalent: Spring Contribution % per Bucket ==========
% figure('Name','Fig 19 - Spring Contribution %','Position',[50 50 1600 450]);
% bucket_names = {'Bucket 1 (60\circ)', 'Bucket 2 (180\circ)', 'Bucket 3 (300\circ)'};
% for b_idx = 1:3
%     subplot(1,3,b_idx);
%     semilogx(ddd(1:count), M_tz_pct(1:count,b_idx), '-o', 'LineWidth', 2, 'Color', [0.494 0.184 0.556]); hold on;
%     semilogx(ddd(1:count), M_py_pct(1:count,b_idx), '-s', 'LineWidth', 2, 'Color', [0.850 0.325 0.098]);
%     semilogx(ddd(1:count), M_qztip_pct(1:count,b_idx), '-^', 'LineWidth', 2, 'Color', [0.466 0.674 0.188]);
%     semilogx(ddd(1:count), M_qzlid_pct(1:count,b_idx), '-d', 'LineWidth', 2, 'Color', [0.301 0.745 0.933]);
%     grid on; ylim([0 100]);
%     xlabel('Rotation Angle (degree)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     ylabel('Sharing percentage of springs contribution (%)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     title(bucket_names{b_idx}, 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'Interpreter', 'tex');
%     legend({'t-z Springs','p-y Springs','q-z Springs (tip)','q-z Springs (lid)'}, 'FontSize', 10, 'Location', 'best');
%     set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% end
% exportgraphics(gcf, 'Fig19_tripod.png', 'Resolution', 600)
% 
% %% ========== Fig 20 equivalent: t-z Force Variation per Bucket ==========
% figure('Name','Fig 20 - t-z Force Variation','Position',[50 50 1600 450]);
% for b_idx = 1:3
%     subplot(1,3,b_idx);
%     semilogx(ddd(1:count), tz_right_out_store(1:count,b_idx), '-o', 'LineWidth', 2, 'Color', [0.111 0.111 0.111]); hold on;
%     semilogx(ddd(1:count), tz_right_in_store(1:count,b_idx), '-s', 'LineWidth', 2, 'Color', [0.494 0.184 0.556]);
%     semilogx(ddd(1:count), tz_left_in_store(1:count,b_idx), '-^', 'LineWidth', 2, 'Color', [0.850 0.325 0.098]);
%     semilogx(ddd(1:count), tz_left_out_store(1:count,b_idx), '-d', 'LineWidth', 2, 'Color', [0.466 0.674 0.188]);
%     grid on;
%     xlabel('Rotation angle (deg.)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     ylabel('Skirt Wall Friction, (F_{tz}/\gamma^{\prime}/D_o^3)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     title(bucket_names{b_idx}, 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'Interpreter', 'tex');
%     legend({'Outer right side','Inner right side','Inner left side','Outer left side'}, 'FontSize', 10, 'Location', 'best');
%     set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% end
% exportgraphics(gcf, 'Fig20_tripod.png', 'Resolution', 600)
% 
% %% ========== Fig 21 equivalent: Lateral Soil Resistance per Bucket ==========
% figure('Name','Fig 21 - Lateral Soil Resistance','Position',[50 50 1600 450]);
% z_layer_centers = ((1:n_l) - 0.5) * (l_b/n_l);
% z_norm_fig21 = -z_layer_centers / l_b;  % 0 at surface, -1 at tip
% plot_steps_fig21 = 1:count;  % plot all rotation steps
% KKx_stores_all = {KKx_store1, KKx_store2, KKx_store3};
% fig21_colors = lines(length(plot_steps_fig21));
% for b_idx = 1:3
%     subplot(1,3,b_idx);
%     KKx_s = KKx_stores_all{b_idx};
%     leg_entries = cell(1, length(plot_steps_fig21));
%     for si = 1:length(plot_steps_fig21)
%         step_idx = plot_steps_fig21(si);
%         kkx_mat = reshape(KKx_s(step_idx,:), n, n_l);
%         kkx_per_layer = sum(kkx_mat, 1);
%         kkx_normalized = kkx_per_layer / (gamma * do_b^3);
%         plot(kkx_normalized, z_norm_fig21, '-o', 'LineWidth', 2, 'Color', fig21_colors(si,:)); hold on;
%         leg_entries{si} = sprintf('\\theta = %g\\circ', round(ddd(step_idx),2));
%     end
%     grid on;
%     xlabel('F_{py}/\gamma^{\prime}/D^3', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     ylabel('z/L', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     title(bucket_names{b_idx}, 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman', 'Interpreter', 'tex');
%     legend(leg_entries, 'FontSize', 9, 'Location', 'best', 'Interpreter', 'tex');
%     set(gca, 'FontSize', 12, 'FontWeight', 'bold');
% end
% exportgraphics(gcf, 'Fig21_tripod.png', 'Resolution', 600)
% 
% % legend({'q-z','q-z sliding','t-z','p-y','total'},'fontsize',12,'interpreter','latex','location','northeast')
% %% ========== Figure 15: 3D p-y Surface (left) + 2D Heatmaps (right column) ==========
% figure(15);
% set(gcf, 'Position', [50 50 1800 1000], 'Color', 'w');
% 
% % --- (a) 3D p-y Surface (left) ---
% ax_3d = axes('Position', [0.03 0.08 0.50 0.86]);
% hold(ax_3d, 'on');
% s_fixed_fig15 = {first_s_end_fixed, second_s_end_fixed, third_s_end_fixed};
% KKx_fig15 = {KKx1, KKx2, KKx3};
% KKy_fig15 = {KKy1, KKy2, KKy3};
% for b_f = 1:3
%     Xg = reshape(s_fixed_fig15{b_f}(1,:), n, n_l);
%     Yg = reshape(s_fixed_fig15{b_f}(2,:), n, n_l);
%     Zg = reshape(s_fixed_fig15{b_f}(3,:), n, n_l);
%     Cg = reshape(sqrt(KKx_fig15{b_f}.^2 + KKy_fig15{b_f}.^2), n, n_l);
%     Xg = [Xg; Xg(1,:)]; Yg = [Yg; Yg(1,:)]; Zg = [Zg; Zg(1,:)]; Cg = [Cg; Cg(1,:)];
%     surf(ax_3d, Xg, Yg, Zg, Cg, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3, 'FaceAlpha', 0.95);
% end
% scatter3(ax_3d, Q(1), Q(2), Q(3), 200, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.2 0.9 0.2], 'LineWidth', 2);
% plot3(ax_3d, CENTE_OF_ROTATION(1,1:count), CENTE_OF_ROTATION(2,1:count), CENTE_OF_ROTATION(3,1:count), ...
%     '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.5, 'MarkerSize', 7, 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'k');
% text(ax_3d, l*cosd(th1), l*sind(th1), 4, 'B1 (60°)',   'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(ax_3d, l*cosd(th2), l*sind(th2), 4, 'B2 (180°)', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(ax_3d, l*cosd(th3), l*sind(th3), 4, 'B3 (300°)', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% colormap(ax_3d, jet); cb15 = colorbar(ax_3d);
% cb15.Label.String = '|F_{py}| (kN)'; cb15.Label.FontSize = 16; cb15.Label.FontWeight = 'bold'; cb15.FontSize = 13;
% t15 = title(ax_3d, sprintf('(a) Tripod p-y Force Distribution at \\theta = %.2f°', dd*180/pi), ...
%     'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(t15, 'Units', 'normalized', 'Position', [0.5, 1.08, 0]);
% xlabel(ax_3d, 'X (m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% ylabel(ax_3d, 'Y (m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% zlabel(ax_3d, 'Z (m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(ax_3d, 'FontSize', 14, 'FontWeight', 'bold');
% grid(ax_3d, 'on'); axis(ax_3d, 'tight'); axis(ax_3d, 'equal'); view(ax_3d, 135, 25);
% hold(ax_3d, 'off');
% 
% % --- (b) 2D p-y Heatmap panels (right column, stacked vertically) ---
% hm_strip_angles = (1:n) * (360/n);
% hm_z_norm = -((1:n_l) - 0.5) * (l_b/n_l) / l_b;
% KKx_hm = {KKx1, KKx2, KKx3};
% KKy_hm = {KKy1, KKy2, KKy3};
% bucket_names_hm = {'Bucket 1 (60°)', 'Bucket 2 (180°)', 'Bucket 3 (300°)'};
% cmax_py = max(sqrt([KKx1.^2+KKy1.^2, KKx2.^2+KKy2.^2, KKx3.^2+KKy3.^2]));
% hm_positions = {[0.57 0.70 0.20 0.22], [0.57 0.40 0.20 0.22], [0.57 0.10 0.20 0.22]};
% hm_labels = {'(b1)', '(b2)', '(b3)'};
% for b_f = 1:3
%     ax_hm = axes('Position', hm_positions{b_f});
%     Cpy = reshape(sqrt(KKx_hm{b_f}.^2 + KKy_hm{b_f}.^2), n, n_l)';
%     pcolor(ax_hm, hm_strip_angles, hm_z_norm, Cpy); shading(ax_hm, 'interp');
%     colormap(ax_hm, jet); clim(ax_hm, [0 cmax_py]);
%     if b_f == 3
%         cb30 = colorbar(ax_hm);
%         cb30.Label.String = '|F_{py}| (kN)'; cb30.Label.FontSize = 12; cb30.Label.FontWeight = 'bold'; cb30.FontSize = 10;
%     end
%     xlabel(ax_hm, 'Azimuthal angle (°)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     ylabel(ax_hm, 'z/L', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     title(ax_hm, [hm_labels{b_f} ' ' bucket_names_hm{b_f}], 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     set(ax_hm, 'FontSize', 10, 'FontWeight', 'bold');
% end
% 
% %% ========== Figure 16: 3D t-z Surface (left) + 2D Heatmaps (right column) ==========
% figure(16);
% set(gcf, 'Position', [50 50 1800 1000], 'Color', 'w');
% 
% % Diverging colormap for signed t-z
% n_cmap = 256;
% blue_white_red = [linspace(0,1,n_cmap/2)', linspace(0,1,n_cmap/2)', ones(n_cmap/2,1); ...
%                   ones(n_cmap/2,1), linspace(1,0,n_cmap/2)', linspace(1,0,n_cmap/2)'];
% 
% % --- (a) 3D t-z Surface (left) ---
% ax_3d16 = axes('Position', [0.03 0.08 0.50 0.86]);
% hold(ax_3d16, 'on');
% KKz_fig16 = {KKz1, KKz2, KKz3};
% for b_f = 1:3
%     Xg = reshape(s_fixed_fig15{b_f}(1,:), n, n_l);
%     Yg = reshape(s_fixed_fig15{b_f}(2,:), n, n_l);
%     Zg = reshape(s_fixed_fig15{b_f}(3,:), n, n_l);
%     Cg = reshape(KKz_fig16{b_f}, n, n_l);
%     Xg = [Xg; Xg(1,:)]; Yg = [Yg; Yg(1,:)]; Zg = [Zg; Zg(1,:)]; Cg = [Cg; Cg(1,:)];
%     surf(ax_3d16, Xg, Yg, Zg, Cg, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3, 'FaceAlpha', 0.95);
% end
% scatter3(ax_3d16, Q(1), Q(2), Q(3), 200, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.2 0.9 0.2], 'LineWidth', 2);
% plot3(ax_3d16, CENTE_OF_ROTATION(1,1:count), CENTE_OF_ROTATION(2,1:count), CENTE_OF_ROTATION(3,1:count), ...
%     '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.5, 'MarkerSize', 7, 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'k');
% text(ax_3d16, l*cosd(th1), l*sind(th1), 4, 'B1 (60°)',   'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(ax_3d16, l*cosd(th2), l*sind(th2), 4, 'B2 (180°)', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(ax_3d16, l*cosd(th3), l*sind(th3), 4, 'B3 (300°)', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% cmax_tz = max(abs([KKz1, KKz2, KKz3]));
% clim(ax_3d16, [-cmax_tz, cmax_tz]);
% colormap(ax_3d16, blue_white_red);
% cb16 = colorbar(ax_3d16);
% cb16.Label.String = 'F_{tz} (kN)'; cb16.Label.FontSize = 16; cb16.Label.FontWeight = 'bold'; cb16.FontSize = 13;
% t16 = title(ax_3d16, sprintf('(a) Tripod t-z Force Distribution at \\theta = %.2f°', dd*180/pi), ...
%     'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(t16, 'Units', 'normalized', 'Position', [0.5, 1.08, 0]);
% xlabel(ax_3d16, 'X (m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% ylabel(ax_3d16, 'Y (m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% zlabel(ax_3d16, 'Z (m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(ax_3d16, 'FontSize', 14, 'FontWeight', 'bold');
% grid(ax_3d16, 'on'); axis(ax_3d16, 'tight'); axis(ax_3d16, 'equal'); view(ax_3d16, 135, 25);
% hold(ax_3d16, 'off');
% 
% % --- (b) 2D t-z Heatmap panels (right column, stacked vertically) ---
% cmax_tz_hm = max(abs([KKz1, KKz2, KKz3]));
% hm_positions16 = {[0.57 0.70 0.20 0.22], [0.57 0.40 0.20 0.22], [0.57 0.10 0.20 0.22]};
% hm_labels16 = {'(b1)', '(b2)', '(b3)'};
% for b_f = 1:3
%     ax_hm16 = axes('Position', hm_positions16{b_f});
%     Ctz = reshape(KKz_fig16{b_f}, n, n_l)';
%     pcolor(ax_hm16, hm_strip_angles, hm_z_norm, Ctz); shading(ax_hm16, 'interp');
%     clim(ax_hm16, [-cmax_tz_hm, cmax_tz_hm]);
%     colormap(ax_hm16, blue_white_red);
%     if b_f == 3
%         cb31 = colorbar(ax_hm16);
%         cb31.Label.String = 'F_{tz} (kN)'; cb31.Label.FontSize = 12; cb31.Label.FontWeight = 'bold'; cb31.FontSize = 10;
%     end
%     xlabel(ax_hm16, 'Azimuthal angle (°)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     ylabel(ax_hm16, 'z/L', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     title(ax_hm16, [hm_labels16{b_f} ' ' bucket_names_hm{b_f}], 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
%     set(ax_hm16, 'FontSize', 10, 'FontWeight', 'bold');
% end
% 
% %% ========== Figure 18: Outer t-z (left) + Inner t-z (right) Side by Side ==========
% figure(18);
% set(gcf, 'Position', [50 50 1600 800], 'Color', 'w');
% 
% % --- (a) Outer t-z Surface (left) ---
% ax_outer = axes('Position', [0.03 0.08 0.43 0.86]);
% hold(ax_outer, 'on');
% KKz_out_fig = {KKz_out1, KKz_out2, KKz_out3};
% cmax_tz_out = max(abs([KKz_out1, KKz_out2, KKz_out3]));
% for b_f = 1:3
%     Xg = reshape(s_fixed_fig15{b_f}(1,:), n, n_l);
%     Yg = reshape(s_fixed_fig15{b_f}(2,:), n, n_l);
%     Zg = reshape(s_fixed_fig15{b_f}(3,:), n, n_l);
%     Cg = reshape(KKz_out_fig{b_f}, n, n_l);
%     Xg = [Xg; Xg(1,:)]; Yg = [Yg; Yg(1,:)]; Zg = [Zg; Zg(1,:)]; Cg = [Cg; Cg(1,:)];
%     surf(ax_outer, Xg, Yg, Zg, Cg, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3, 'FaceAlpha', 0.95);
% end
% scatter3(ax_outer, Q(1), Q(2), Q(3), 200, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.2 0.9 0.2], 'LineWidth', 2);
% plot3(ax_outer, CENTE_OF_ROTATION(1,1:count), CENTE_OF_ROTATION(2,1:count), CENTE_OF_ROTATION(3,1:count), ...
%     '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.5, 'MarkerSize', 7, 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'k');
% text(ax_outer, l*cosd(th1), l*sind(th1), 4, 'B1 (60°)',   'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(ax_outer, l*cosd(th2), l*sind(th2), 4, 'B2 (180°)', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(ax_outer, l*cosd(th3), l*sind(th3), 4, 'B3 (300°)', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% clim(ax_outer, [-cmax_tz_out, cmax_tz_out]);
% colormap(ax_outer, blue_white_red);
% cb18 = colorbar(ax_outer);
% cb18.Label.String = 'F_{tz,outer} (kN)'; cb18.Label.FontSize = 16; cb18.Label.FontWeight = 'bold'; cb18.FontSize = 13;
% t18 = title(ax_outer, sprintf('(a) Outer t-z at \\theta = %.2f°', dd*180/pi), ...
%     'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(t18, 'Units', 'normalized', 'Position', [0.5, 1.08, 0]);
% xlabel(ax_outer, 'X (m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% ylabel(ax_outer, 'Y (m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% zlabel(ax_outer, 'Z (m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(ax_outer, 'FontSize', 14, 'FontWeight', 'bold');
% grid(ax_outer, 'on'); axis(ax_outer, 'tight'); axis(ax_outer, 'equal'); view(ax_outer, 135, 25);
% hold(ax_outer, 'off');
% 
% % --- (b) Inner t-z Surface (right) ---
% ax_inner = axes('Position', [0.53 0.08 0.43 0.86]);
% hold(ax_inner, 'on');
% s_fixed_in_fig = {first_s_end_fixed_in, second_s_end_fixed_in, third_s_end_fixed_in};
% KKz_in_fig = {KKz_in1, KKz_in2, KKz_in3};
% cmax_tz_in = max(abs([KKz_in1, KKz_in2, KKz_in3]));
% for b_f = 1:3
%     Xg = reshape(s_fixed_in_fig{b_f}(1,:), n, n_l);
%     Yg = reshape(s_fixed_in_fig{b_f}(2,:), n, n_l);
%     Zg = reshape(s_fixed_in_fig{b_f}(3,:), n, n_l);
%     Cg = reshape(KKz_in_fig{b_f}, n, n_l);
%     Xg = [Xg; Xg(1,:)]; Yg = [Yg; Yg(1,:)]; Zg = [Zg; Zg(1,:)]; Cg = [Cg; Cg(1,:)];
%     surf(ax_inner, Xg, Yg, Zg, Cg, 'EdgeColor', [0.3 0.3 0.3], 'EdgeAlpha', 0.3, 'FaceAlpha', 0.95);
% end
% scatter3(ax_inner, Q(1), Q(2), Q(3), 200, 'p', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.2 0.9 0.2], 'LineWidth', 2);
% plot3(ax_inner, CENTE_OF_ROTATION(1,1:count), CENTE_OF_ROTATION(2,1:count), CENTE_OF_ROTATION(3,1:count), ...
%     '-s', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.5, 'MarkerSize', 7, 'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'k');
% text(ax_inner, l*cosd(th1), l*sind(th1), 4, 'B1 (60°)',   'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(ax_inner, l*cosd(th2), l*sind(th2), 4, 'B2 (180°)', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% text(ax_inner, l*cosd(th3), l*sind(th3), 4, 'B3 (300°)', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
% clim(ax_inner, [-cmax_tz_in, cmax_tz_in]);
% colormap(ax_inner, blue_white_red);
% cb19 = colorbar(ax_inner);
% cb19.Label.String = 'F_{tz,inner} (kN)'; cb19.Label.FontSize = 16; cb19.Label.FontWeight = 'bold'; cb19.FontSize = 13;
% t19 = title(ax_inner, sprintf('(b) Inner t-z at \\theta = %.2f°', dd*180/pi), ...
%     'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(t19, 'Units', 'normalized', 'Position', [0.5, 1.08, 0]);
% xlabel(ax_inner, 'X (m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% ylabel(ax_inner, 'Y (m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% zlabel(ax_inner, 'Z (m)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% set(ax_inner, 'FontSize', 14, 'FontWeight', 'bold');
% grid(ax_inner, 'on'); axis(ax_inner, 'tight'); axis(ax_inner, 'equal'); view(ax_inner, 135, 25);
% hold(ax_inner, 'off');
% 
% figure(17);
% scatter3(r_Pile_first_s_end(1,:), r_Pile_first_s_end(2,:), r_Pile_first_s_end(3,:), 'ko', 'LineWidth', 1);  % 'ko' is for black circles
% hold on
% scatter3(r_Pile_first_s_end_fixed(1,:), r_Pile_first_s_end_fixed(2,:), r_Pile_first_s_end_fixed(3,:), 'bo', 'LineWidth', 1);
% hold on
% % force_directions = [zeros(1, size(qb2_1,2)); zeros(1, size(qb2_1,2)); qb2_1*100]
% % hold on
% % scatter3(first_s_end_fixed(1,:),first_s_end_fixed(2,:),first_s_end_fixed(3,:),'*blue', 'LineWidth', 5)
% % % Plotting the force vectors
% % hold on
% % quiver3(r_Pile_first_s_end_fixed(1,:), r_Pile_first_s_end_fixed(2,:), r_Pile_first_s_end_fixed(3,:)+0.998*l_b, ... 
% %          zeros(1, size(qb2_1,2)), zeros(1, size(qb2_1,2)), qb2_1*1, 'r', 'LineWidth', 3, 'AutoScale', 'off');  % 'r' is for red color
% 
% title('3D Plot of Nodes and Force Vectors');
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');


% figure(17)
% % Assume nodes and forces are already defined
% x = r_Pile_first_s_end_fixed(1,:)
% y = r_Pile_first_s_end_fixed(2,:)
% z = r_Pile_first_s_end_fixed(3,:)
% force_magnitude = qb2_1;  % This should be the same length as x, y, z
% 
% % Create a grid over the node space
% [xq, yq] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
% 
% % Interpolate the z-force values on this grid
% % Assuming the force is directed along z and is proportional to the magnitude
% zq = griddata(x, y, force_magnitude .* z, xq, yq, 'natural');
% 
% % Plot the surface
% figure;
% surf(xq, yq, zq);
% shading interp;  % Interpolate colors across the surface
% colorbar;  % Adds a color bar to indicate the scale of force magnitudes
% title('Surface Plot of Force Magnitudes');
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Force Magnitude');
%% KKX1 plot (disabled - KKX1 not stored)
% i=0;
% for i=1:n_l
%     pc(i)= n*i;
% end
% K = [KKX1; first_s_end_fixed(3,:)];
% KK = zeros(size(K, 1), n_l);
% for i = 1:n_l
%     for j = 0:n-1
%         if pc(i) - j > 0
%             KK(:, i) = KK(:, i) + K(:, pc(i) - j);
%         end
%     end
% end
% KK = KK';
% KKfinal = [-KK(:, 1:end-1), KK(:, end) / n];
% figure(19)
% for i =1:size(ddd,2)
% plot(KKfinal(:,i)/(do_b*tan(pi/n))/(l_b/n_l)/gamma/do_b^3,KKfinal(:,end)/(l_b),'-*','LineWidth',3)
% hold on
% grid on
% end
% xlabel('(F_p_y/\gamma^\prime/D_o^3)','FontSize', 22, 'FontWeight', 'bold','FontName','Times New Roman');
% ylabel('(z/l_b)','FontSize', 22, 'FontWeight', 'bold','FontName','Times New Roman');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontweight','bold','fontsize',22)
% legend({'$0^\circ$','$0.5^\circ$','$1^\circ$','$1.5^\circ$','$2.0^\circ$','$2.5^\circ$','$3^\circ$'},'fontsize',20,'interpreter','latex','location','northeast')

%%
% colors1 = [0.001, 0.15, 0.175; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.1010, 0.6740, 0.1880;0.001, 0.15, 0.175; 0.75, 0, 0.75; 0.6350, 0.0080, 0.1840; 0.2510, 0.2510, 0.5519; 0.3010, 0.7550, 0.5530; 0.4411, 0.4155, 0.0350];
% colors2 = [0.001, 0.15, 0.175; 0.75, 0, 0.75; 0.6350, 0.0080, 0.1840; 0.2510, 0.2510, 0.5519; 0.3010, 0.7550, 0.5530; 0.4411, 0.4155, 0.0350];
% 
% Q1 = [KKZ1;first_s_end_fixed(3,:)];
% QQ1 = Q1(:,pc)'+Q1(:,pc-1)'+Q1(:,pc-7)';
% QQfinal1 = [QQ1(:,1:end-1),QQ1(:,end)/3 ];
% 
% Q2 = [KKZ1;first_s_end_fixed(3,:)];
% QQ2 = Q2(:,pc-3)'+Q2(:,pc-4)'+Q2(:,pc-5)';
% QQfinal2 = [QQ2(:,1:end-1),QQ2(:,end)/3];

% figure(20)
% % for i =1:20
% % plot(QQfinal1(:,i)/(pi*do_b/n)/(l_b/n_l)/gamma/do_b^3/2,QQfinal1(:,end)/l_b,'-o','LineWidth',4)
% % hold on
% % plot(QQfinal2(:,i)/(pi*do_b/n)/(l_b/n_l)/gamma/do_b^3/2,QQfinal2(:,end)/l_b,'-o','LineWidth',4)
% %  hold on
% % grid on
% % end

% xlabel('(F_t_z/\gamma^\prime D_o^2)','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
% ylabel('Depth,(z/L)','FontSize', 24, 'FontWeight', 'bold','FontName','Times New Roman');
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'fontweight','bold','fontsize',22)
% set(gcf, 'Position', [100, 100, 1200, 600]); % Wider figure to accommodate legend
% set(gca, 'YAxisLocation', 'right');
% 
% legend({'$0^\circ$ left side','$0.01^\circ$ left side','$0.05^\circ$ left side','$0.15^\circ$ left side',...
%     '$0.25^\circ$ left side','$0.50^\circ$ left side'},'fontsize',16,'interpreter','latex','location','northeast')
% % legend({'$0^\circ$ right side','$0.01^\circ$ right side','$0.05^\circ$ right side','$0.15^\circ$ right side',...
% %     '$0.25^\circ$ right side','$0.50^\circ$ right side'},'fontsize',16,'interpreter','latex','location','northeast')
% figure(222)
% scatter3(first_s_end(1,:),first_s_end(2,:),first_s_end(3,:),'*k')
% hold on
% scatter3(second_s_end(1,:),second_s_end(2,:),second_s_end(3,:),'*k')
% hold on
% scatter3(third_s_end(1,:),third_s_end(2,:),third_s_end(3,:),'*k')
% 
% % grid on
%  hold on
% scatter3(first_s_end_fixed(1,:),first_s_end_fixed(2,:),first_s_end_fixed(3,:),'*r')
% hold on
% scatter3(second_s_end_fixed(1,:),second_s_end_fixed(2,:),second_s_end_fixed(3,:),'*r')
% hold on
% scatter3(third_s_end_fixed(1,:),third_s_end_fixed(2,:),third_s_end_fixed(3,:),'*r')
% 
% 
% 
% 
% figure(223)
% scatter3(relative_first_s_end(1,:),relative_first_s_end(2,:),relative_first_s_end(3,:),'*k')
% hold on
% scatter3(relative_second_s_end(1,:),relative_second_s_end(2,:),relative_second_s_end(3,:),'*k')
% hold on
% scatter3(relative_third_s_end(1,:),relative_third_s_end(2,:),relative_third_s_end(3,:),'*k')
% 
% % grid on
%  hold on
% scatter3(relative_first_s_end_fixed(1,:),relative_first_s_end_fixed(2,:),relative_first_s_end_fixed(3,:),'*r')
% hold on
% scatter3(relative_second_s_end_fixed(1,:),relative_second_s_end_fixed(2,:),relative_second_s_end_fixed(3,:),'*r')
% hold on
% scatter3(relative_third_s_end_fixed(1,:),relative_third_s_end_fixed(2,:),relative_third_s_end_fixed(3,:),'*r')

