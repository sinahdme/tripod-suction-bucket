clc; clear all; close all;
%% ====== TRIPOD BUCKET FOUNDATION — 3D VISUALIZATION ======
% Standalone script: shows untilted (gray) vs tilted (colored) buckets
% Two buckets in compression (red/orange), one in tension (blue)

%% Geometric parameters (Kim 2014 tripod)
do_b = 6.5;         % bucket outer diameter [m]
ro_b = do_b/2;
thickness = 0.025;
di_b = do_b - 2*thickness;
ri_b = di_b/2;
l_b = 8.0;          % skirt length [m]
l = 15.5;           % circumradius [m]
n = 32;             % circumferential divisions (finer for smooth look)
n_l = 20;           % vertical layers
eccentricity = 33.39;

% Pod angles — tripod at 0, 120, 240
th = [0, 120, 240];

% Tilt angle (degrees) — large enough to see clearly
tilt_deg = 3;
dd = tilt_deg * pi/180;

%% Generate bucket geometry (body frame)
circle_split = linspace(0, 360, n+1);
ca = cosd(circle_split(2:end));
sa = sind(circle_split(2:end));
z_centers = -((0:n_l-1) + 0.5) * (l_b/n_l);

% Store per-bucket: skirt nodes
s_end = cell(3,1);      % body frame
for b = 1:3
    cx = l*cosd(th(b));
    cy = l*sind(th(b));
    nodes = zeros(3, n*n_l);
    idx = 1;
    for j = 1:n_l
        zc = z_centers(j);
        nodes(:, idx:idx+n-1) = [cx + ro_b*ca; cy + ro_b*sa; zc*ones(1,n)];
        idx = idx + n;
    end
    s_end{b} = nodes;
end

%% Generate lid circle nodes (for drawing lid disks)
lid_top = cell(3,1);
lid_bot = cell(3,1);
n_lid_r = 8; % radial divisions for lid
for b = 1:3
    cx = l*cosd(th(b));
    cy = l*sind(th(b));
    lid_top{b} = [cx + ro_b*ca; cy + ro_b*sa; zeros(1,n)];
    lid_bot{b} = [cx + ro_b*ca; cy + ro_b*sa; -l_b*ones(1,n)];
end

%% Set up rotation
% Rotation center Q (approximate: at skirt tip depth, offset from center)
Q = [5.5; 0; -l_b];  % approximate center of rotation

% Rodrigues rotation: tilt of magnitude dd in the alpha=0 direction
% (alpha=0 means tower tilts toward -x, so -x side goes down = compression)
R_tilt = rodrigues_rotation(dd, 0);

%% Transform points
s_end_fixed = cell(3,1);
lid_top_fixed = cell(3,1);
lid_bot_fixed = cell(3,1);

for b = 1:3
    nodes = s_end{b};
    nodes_f = zeros(size(nodes));
    for i = 1:size(nodes,2)
        nodes_f(:,i) = transformation(nodes(:,i), Q, R_tilt);
    end
    s_end_fixed{b} = nodes_f;

    lt = lid_top{b};
    lt_f = zeros(size(lt));
    for i = 1:size(lt,2)
        lt_f(:,i) = transformation(lt(:,i), Q, R_tilt);
    end
    lid_top_fixed{b} = lt_f;

    lb_ = lid_bot{b};
    lb_f = zeros(size(lb_));
    for i = 1:size(lb_,2)
        lb_f(:,i) = transformation(lb_(:,i), Q, R_tilt);
    end
    lid_bot_fixed{b} = lb_f;
end

%% Determine compression/tension per bucket
% Based on average z-displacement of skirt tip (bottom ring)
bucket_type = zeros(3,1); % +1 = compression (down), -1 = tension (up)
for b = 1:3
    % bottom ring: last n nodes
    z_orig = mean(s_end{b}(3, end-n+1:end));
    z_tilt = mean(s_end_fixed{b}(3, end-n+1:end));
    dz_avg = z_tilt - z_orig;
    if dz_avg < 0
        bucket_type(b) = 1;   % pushed further down = compression
    else
        bucket_type(b) = -1;  % pulled up = tension
    end
end

fprintf('Bucket states: ');
for b = 1:3
    if bucket_type(b) == 1
        fprintf('Pod%d=COMPRESSION  ', b);
    else
        fprintf('Pod%d=TENSION  ', b);
    end
end
fprintf('\n');

%% ====== PLOTTING ======
figure('Position', [50, 50, 1400, 900], 'Color', 'w');
hold on;

% Colors
color_compression = [0.85, 0.15, 0.15];   % red
color_tension     = [0.15, 0.40, 0.85];   % blue
color_untilted    = [0.6, 0.6, 0.6];      % gray
alpha_untilted    = 0.15;
alpha_tilted      = 0.45;

%% --- Draw ground plane (mudline at z=0) ---
ground_size = 25;
gx = [-ground_size ground_size ground_size -ground_size];
gy = [-ground_size -ground_size ground_size ground_size];
gz = [0 0 0 0];
patch(gx, gy, gz, [0.82, 0.76, 0.62], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

%% --- Draw untilted buckets (gray wireframe + transparent surface) ---
for b = 1:3
    nodes = s_end{b};

    % Draw rings
    for j = 1:n_l
        idx = (j-1)*n + (1:n);
        ring = nodes(:, idx);
        ring = [ring, ring(:,1)]; % close the ring
        plot3(ring(1,:), ring(2,:), ring(3,:), '-', 'Color', [color_untilted, 0.4], 'LineWidth', 0.8);
    end

    % Draw vertical lines (every 4th strip for clarity)
    for i = 1:4:n
        vline = zeros(3, n_l);
        for j = 1:n_l
            vline(:,j) = nodes(:, (j-1)*n + i);
        end
        plot3(vline(1,:), vline(2,:), vline(3,:), '-', 'Color', [color_untilted, 0.3], 'LineWidth', 0.5);
    end

    % Draw top lid (untilted) as transparent disk
    lt = lid_top{b};
    lt_closed = [lt, lt(:,1)];
    patch(lt_closed(1,:), lt_closed(2,:), lt_closed(3,:), color_untilted, ...
          'FaceAlpha', alpha_untilted, 'EdgeColor', color_untilted, 'EdgeAlpha', 0.3);

    % Draw bottom ring (untilted)
    lb_ = lid_bot{b};
    lb_closed = [lb_, lb_(:,1)];
    plot3(lb_closed(1,:), lb_closed(2,:), lb_closed(3,:), '-', 'Color', [color_untilted, 0.4], 'LineWidth', 0.8);
end

%% --- Draw tilted buckets (colored by compression/tension) ---
for b = 1:3
    nodes = s_end_fixed{b};

    % Choose color
    if bucket_type(b) == 1
        col = color_compression;
    else
        col = color_tension;
    end

    % Draw filled surface rings (cylinder wall)
    for j = 1:n_l-1
        idx1 = (j-1)*n + (1:n);
        idx2 = j*n + (1:n);
        for i = 1:n
            i2 = mod(i, n) + 1;
            verts = [nodes(:, idx1(i))'; nodes(:, idx1(i2))'; ...
                     nodes(:, idx2(i2))'; nodes(:, idx2(i))'];
            patch('Vertices', verts, 'Faces', [1 2 3 4], ...
                  'FaceColor', col, 'FaceAlpha', alpha_tilted, ...
                  'EdgeColor', 'none');
        end
    end

    % Draw rings on tilted bucket (solid colored lines)
    for j = 1:n_l
        idx = (j-1)*n + (1:n);
        ring = nodes(:, idx);
        ring = [ring, ring(:,1)];
        plot3(ring(1,:), ring(2,:), ring(3,:), '-', 'Color', col, 'LineWidth', 1.2);
    end

    % Draw top lid (tilted)
    lt = lid_top_fixed{b};
    lt_closed = [lt, lt(:,1)];
    patch(lt_closed(1,:), lt_closed(2,:), lt_closed(3,:), col, ...
          'FaceAlpha', alpha_tilted+0.1, 'EdgeColor', col, 'LineWidth', 1.5);

    % Draw bottom ring (tilted)
    lb_ = lid_bot_fixed{b};
    lb_closed = [lb_, lb_(:,1)];
    patch(lb_closed(1,:), lb_closed(2,:), lb_closed(3,:), col, ...
          'FaceAlpha', alpha_tilted, 'EdgeColor', col, 'LineWidth', 1.5);
end

%% --- Draw connection frame (lines from center to each bucket top) ---
center_body = [0; 0; 0];
center_fixed = transformation(center_body, Q, R_tilt);

for b = 1:3
    cx = l*cosd(th(b));
    cy = l*sind(th(b));
    pod_top = transformation([cx; cy; 0], Q, R_tilt);
    plot3([center_fixed(1), pod_top(1)], ...
          [center_fixed(2), pod_top(2)], ...
          [center_fixed(3), pod_top(3)], ...
          '-', 'Color', [0.3, 0.3, 0.3], 'LineWidth', 2.5);
end

% Draw untilted frame (dashed gray)
for b = 1:3
    cx = l*cosd(th(b));
    cy = l*sind(th(b));
    plot3([0, cx], [0, cy], [0, 0], '--', 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1.5);
end

%% --- Draw center of rotation ---
plot3(Q(1), Q(2), Q(3), 'p', 'MarkerSize', 15, 'MarkerFaceColor', [0.2, 0.8, 0.2], ...
      'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

%% --- Draw load arrow (tower eccentricity direction) ---
tower_base = center_fixed;
tower_top = transformation([0;0;eccentricity], Q, R_tilt);
% Scale the arrow for visibility (show partial tower)
arrow_top = tower_base + 0.4*(tower_top - tower_base);
quiver3(arrow_top(1), arrow_top(2), arrow_top(3), ...
        3*sind(tilt_deg), 0, 0, 0, ...
        'LineWidth', 3, 'Color', [0.9, 0.5, 0.0], 'MaxHeadSize', 2);
text(arrow_top(1)+2, arrow_top(2), arrow_top(3)+1, 'H (load)', ...
     'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.9, 0.5, 0.0]);

%% --- Labels for each bucket ---
for b = 1:3
    cx = l*cosd(th(b));
    cy = l*sind(th(b));
    pod_top = transformation([cx; cy; 0], Q, R_tilt);

    if bucket_type(b) == 1
        label = sprintf('Pod %d\n(Compression)', b);
        tcol = color_compression;
    else
        label = sprintf('Pod %d\n(Tension)', b);
        tcol = color_tension;
    end
    text(pod_top(1), pod_top(2), pod_top(3)+1.5, label, ...
         'FontSize', 13, 'FontWeight', 'bold', 'Color', tcol, ...
         'HorizontalAlignment', 'center');
end

%% --- Axis formatting ---
axis equal;
grid on;
box on;
xlabel('X (m)', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylabel('Y (m)', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
zlabel('Z (m)', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
title(sprintf('Tripod Suction Bucket Foundation — Tilt = %.1f%c', tilt_deg, char(176)), ...
      'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
set(gca, 'FontSize', 13, 'FontName', 'Times New Roman');
view(135, 25);
xlim([-22 22]);
ylim([-22 22]);
zlim([-12 8]);

% Legend
h1 = patch(NaN, NaN, NaN, color_untilted, 'FaceAlpha', 0.3);
h2 = patch(NaN, NaN, NaN, color_compression, 'FaceAlpha', 0.5);
h3 = patch(NaN, NaN, NaN, color_tension, 'FaceAlpha', 0.5);
h4 = plot3(NaN, NaN, NaN, 'p', 'MarkerSize', 12, 'MarkerFaceColor', [0.2,0.8,0.2], 'MarkerEdgeColor', 'k');
legend([h1, h2, h3, h4], {'Untilted position', 'Compression bucket', 'Tension bucket', 'Center of rotation'}, ...
       'FontSize', 13, 'Location', 'northeast', 'FontName', 'Times New Roman');

lighting gouraud;
camlight('headlight');

fprintf('Tilt angle: %.1f degrees\n', tilt_deg);
fprintf('Visualization complete.\n');
