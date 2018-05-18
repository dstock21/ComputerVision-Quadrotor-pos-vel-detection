function [vel, omg] = estimate_vel(sensor, varargin)
%ESTIMATE_VEL 6DOF velocity estimator
%   sensor - struct stored in provided dataset, fields include
%          - is_ready: logical, indicates whether sensor data is valid
%          - t: timestamp
%          - rpy, omg, acc: imu readings, you should not use these in this phase
%          - img: uint8, 240x376 grayscale image
%          - id: 1xn ids of detected tags
%          - p0, p1, p2, p3, p4: 2xn pixel position of center and
%                                four corners of detected tags
%            Y
%            ^ P3 == P2
%            | || P0 ||
%            | P4 == P1
%            o---------> X
%   varargin - any variables you wish to pass into the function, could be
%              a data structure to represent the map or camera parameters,
%              your decision. But for the purpose of testing, since we don't
%              know what inputs you will use, you have to specify them in
%              init_script by doing
%              estimate_vel_handle = ...
%                  @(sensor) estimate_vel(sensor, your personal input arguments);
%   vel - 3x1 velocity of the quadrotor in world frame
%   omg - 3x1 angular velocity of the quadrotor

persistent PT cnsp tp dt
Kinv = varargin{1};
estimate_H_handle = varargin{2};

%% initialize
fs = 17;
Nmin = 100;
nRANSAC = 50;
thresh = 5;

if isempty(PT)
    tp = sensor.t;
    dt = tp;
    
    PT = vision.PointTracker();
    cp = detectHarrisFeatures(sensor.img, 'FilterSize', fs);
    cnsp = cp.Location;
    initialize(PT, cnsp, sensor.img);
    
    vel = zeros(3,1);
    omg = zeros(3,1);
    return
end

dt = 0.9*dt+0.1*(sensor.t-tp);
tp = sensor.t;

%%
[cns, val] = PT(sensor.img);

cns_c = cns(val,:);
cnsp_c = cnsp(val,:);

%% convert uv to xydot
xy = zeros(size(cns_c,1), 3);
xyp = zeros(size(xy));
xydot = zeros(size(xy));
for i = 1:size(xy,1)
    xy(i,:) = [cns_c(i,:) 1]*Kinv';
    xyp(i,:) = [cnsp_c(i,:) 1]*Kinv';
    xydot(i,:) = (xy(i,:) - xyp(i,:)) / dt;
end

%%
H = estimate_H_handle(sensor);
Hinv = inv(H);
Nin_opt = 0;
Vopt = [];

%% run RANSAC
for k = 1:nRANSAC
    sam = randi(size(xydot,1), [3,1]);
    
    xydot_sam = xydot(sam,1:2)';
    A = zeros(6,6);
    for i=1:3
        Zi = 1/(Hinv(3,1:3)*[cnsp_c(sam(i),:) 1]');
        xi = xy(sam(i),1);
        yi = xy(sam(i),2);
        A(2*i-1:2*i, :) = [-1/Zi, 0, xi/Zi, xi*yi, -(1+xi^2), yi; 0, -1/Zi, yi/Zi, 1+yi^2, -xi*yi, -xi];
    end
    
    V = A\xydot_sam(:);
    
    %%
    cmp = zeros(size(cns_c,1));
    for i=1:size(cmp,1)
        Zi = 1/(Hinv(3,1:3)*[cns_c(i,:) 1]');
        xi = xy(i,1);
        yi = xy(i,2);
        Ai = [-1/Zi, 0, xi/Zi, xi*yi, -(1+xi^2), yi; 0, -1/Zi, yi/Zi, 1+yi^2, -xi*yi, -xi];
        cmp(i) = norm(Ai*V / Zi);
    end
    
    inliers = cmp < thresh;
    Nin = sum(inliers);
    
    if Nin > Nin_opt
        Nin_opt = Nin;
        Vopt = V;
    end
end

%%
vel = Vopt(1:3);
omg = Vopt(4:6);

if size(cns_c, 1) < Nmin
    release(PT);
    PT = vision.PointTracker();
    cp = detectHarrisFeatures(sensor.img, 'FilterSize', fs);
    cnsp = cp.Location;
    initialize(PT, cnsp, sensor.img);
else
    cnsp = cns;
end

end
