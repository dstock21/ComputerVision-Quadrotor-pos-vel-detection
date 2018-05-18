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

persistent PT cnsp tp dt cnt
Kinv = varargin{1};
estimate_H_handle = varargin{2};
Tr = varargin{3};

%% initialize
fs = 21;
Nmin = 70;
nRANSAC = 100;
thresh = (1/311.052);

if isempty(sensor.id)
    vel = [];
    omg = [];
    return
end

if isempty(PT)
    tp = sensor.t;
    dt = tp;
    cnt=0;
    
    PT = vision.PointTracker();
    cp = detectHarrisFeatures(sensor.img, 'FilterSize', fs);
    cnsp = cp.Location;
    initialize(PT, cnsp, sensor.img);
    
    vel = zeros(3,1);
    omg = zeros(3,1);
    return
end

dt = 0.85*dt+0.15*(sensor.t-tp);
tp = sensor.t;

%%
[cns, val] = PT(sensor.img);

cns_c = cns(val,:);
cnsp_c = cnsp(val,:);
N = size(cns_c, 1);

%% convert uv to xydot
xy = Kinv*[cns_c'; ones(1,N)];
xyp = Kinv*[cnsp_c'; ones(1,N)];
xydot = (xy(1:2,:) - xyp(1:2,:)) / dt;

%% precompute
[H, R] = estimate_H_handle(sensor);
Hinv = inv(H);
Rinv = inv(R);

A = zeros(2*N, 6);
for i=1:N
    Zi = 1/(Hinv(3,1:3)*[cnsp_c(i,:) 1]');
    xi = xy(1,i);
    yi = xy(2,i);
    A(2*i-1:2*i, :) = [-1/Zi, 0, xi/Zi, xi*yi, -(1+xi^2), yi; 
        0, -1/Zi, yi/Zi, 1+yi^2, -xi*yi, -xi];
end

%%
Nin_opt = 0;
in_opt = [];

%% run RANSAC
for k = 1:nRANSAC
    sam = randperm(N, 3);
    sami = [2*sam-1; 2*sam];
    sami = sami(:);
    
    V = A(sami,:)\xydot(sami);
    
    %%
    cmp = sum(reshape(A*V - xydot(:), [2, N]).^2, 1);
    
    inliers = find(cmp < thresh);
    Nin = length(inliers);
    
    if Nin > Nin_opt
        Nin_opt = Nin;
        in_opt = inliers;
        
        if Nin_opt > 0.95*N
            break;
        end
    end
end

%% recompute V
in_opti = [2*in_opt-1; 2*in_opt];
in_opti = in_opti(:);

Vopt = A(in_opti,:)\xydot(in_opti);

%%
omg = Rinv*Vopt(4:6);
vel = Rinv*Vopt(1:3) + cross(omg, Rinv*Tr);

cnt = cnt+1;
if N < Nmin || mod(cnt, 5) == 0
    release(PT);
    cp = detectHarrisFeatures(sensor.img, 'FilterSize', fs);
    cnsp = cp.Location;
    initialize(PT, cnsp, sensor.img);
else
    cnsp = cns;
end

end
