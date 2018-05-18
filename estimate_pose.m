function [pos, q] = estimate_pose(sensor, varargin)
%ESTIMATE_POSE 6DOF pose estimator based on apriltags
%   sensor - struct stored in provided dataset, fields include
%          - is_ready: logical, indicates whether sensor data is valid
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
%              estimate_pose_handle = ...
%                  @(sensor) estimate_pose(sensor, your personal input arguments);
%   pos - 3x1 position of the quadrotor in world frame
%   q   - 4x1 quaternion of the quadrotor [w, x, y, z] where q = w + x*i + y*j + z*k

%%
id = sensor.id+1;

N = size(id, 2);
%% compute tag map
if N >=1
    
    %% compute tag map
    Kinv = varargin{1};
    Hr = varargin{2};
    Xtag1 = varargin{3};
    Xtag2 = varargin{4};
    Xtag3 = varargin{5};
    Xtag4 = varargin{6};
    Xtagc = varargin{7};

    %%
    Xi1 = Xtag1(:,id);
    Xi2 = Xtag2(:,id);
    Xi3 = Xtag3(:,id);
    Xi4 = Xtag4(:,id);
    Xic = Xtagc(:,id);

    %%

    pcC = zeros(3,N);
    p1C = zeros(3,N);
    p2C = zeros(3,N);
    p3C = zeros(3,N);
    p4C = zeros(3,N);

    for i = 1:N
        pcC(:,i) = [sensor.p0(1,i); sensor.p0(2,i); 1];
        p1C(:,i) = [sensor.p1(1,i); sensor.p1(2,i); 1];
        p2C(:,i) = [sensor.p2(1,i); sensor.p2(2,i); 1];
        p3C(:,i) = [sensor.p3(1,i); sensor.p3(2,i); 1];
        p4C(:,i) = [sensor.p4(1,i); sensor.p4(2,i); 1];
    end

    %%
    A = zeros(N*10, 9);
    for i = 1:N
        A(10*i-9:10*i-8, :) = [-Xi1(1,i), -Xi1(2,i), -1, 0, 0, 0, ...
            p1C(1,i)*Xi1(1,i), p1C(1,i)*Xi1(2,i), p1C(1,i);
            0, 0, 0, -Xi1(1,i), -Xi1(2,i), -1, ...
            p1C(2,i)*Xi1(1,i), p1C(2,i)*Xi1(2,i), p1C(2,i)];
        A(10*i-7:10*i-6, :) = [-Xi2(1,i), -Xi2(2,i), -1, 0, 0, 0, ...
            p2C(1,i)*Xi2(1,i), p2C(1,i)*Xi2(2,i), p2C(1,i);
            0, 0, 0, -Xi2(1,i), -Xi2(2,i), -1, ...
            p2C(2,i)*Xi2(1,i), p2C(2,i)*Xi2(2,i), p2C(2,i)];
        A(10*i-5:10*i-4, :) = [-Xi3(1,i), -Xi3(2,i), -1, 0, 0, 0, ...
            p3C(1,i)*Xi3(1,i), p3C(1,i)*Xi3(2,i), p3C(1,i);
            0, 0, 0, -Xi3(1,i), -Xi3(2,i), -1, ...
            p3C(2,i)*Xi3(1,i), p3C(2,i)*Xi3(2,i), p3C(2,i)];
        A(10*i-3:10*i-2, :) = [-Xi4(1,i), -Xi4(2,i), -1, 0, 0, 0, ...
            p4C(1,i)*Xi4(1,i), p4C(1,i)*Xi4(2,i), p4C(1,i);
            0, 0, 0, -Xi4(1,i), -Xi4(2,i), -1, ...
            p4C(2,i)*Xi4(1,i), p4C(2,i)*Xi4(2,i), p4C(2,i)];
        A(10*i-1:10*i, :) = [-Xic(1,i), -Xic(2,i), -1, 0, 0, 0, ...
            pcC(1,i)*Xic(1,i), pcC(1,i)*Xic(2,i), pcC(1,i);
            0, 0, 0, -Xic(1,i), -Xic(2,i), -1, ...
            pcC(2,i)*Xic(1,i), pcC(2,i)*Xic(2,i), pcC(2,i)];
    end

    [~, ~, V] = svd(A);
    H1 = reshape(V(:,end),3,3);
    H1 = H1/V(9,9);
    H1 = H1';
    H1 = Kinv*H1;

    [U, ~, V] = svd([H1(:,1), H1(:,2), cross(H1(:,1), H1(:,2))]);

    Rc = U*[1, 0, 0; 0, 1, 0; 0, 0, det(U*V')]*V';
    Tc = H1(:,3)/norm(H1(:,1));
    Hc = [Rc, Tc; 0,0,0,1];

    H = Hr*Hc;
    H = inv(H);
    R = H(1:3,1:3);
    T = H(1:3,4);

    %%
    pos = T;

    q = RotToQuat(R);
else
    pos = [];
    q = [];
end
