% Add additional inputs after sensor if you want to
% Example:
% your_input = 1;
% estimate_vel_handle = @(sensor) estimate_vel(sensor, your_input);
%
% We will only call estimate_vel_handle in the test function.
% Note that thise will only create a function handle, but not run the function

clear estimate_vel;

K = [311.0520 0        201.8724;
     0         311.3885 113.6210;
     0         0        1];

Kinv = inv(K);

Rr = [1/sqrt(2), -1/sqrt(2), 0;
    -1/sqrt(2), -1/sqrt(2), 0;
    0, 0, -1];
Tr = [-0.04; 0.0; -0.03];
Hr = [Rr, -Rr*Tr; 0,0,0,1];

Xtag4 = ones(3, 108);
for i = 1:9
    if i <= 3
        Xtag4(2, 12*i-11:12*i) = (i-1)*2*0.152;
    elseif i <= 6
        Xtag4(2, 12*i-11:12*i) = ((i-1)*2-1)*0.152 + 0.178;
    else
        Xtag4(2, 12*i-11:12*i) = ((i-1)*2-2)*0.152 + 2*0.178;
    end
end


for i = 1:12
    Xtag4(1, i:12:108) = (i-1)*2*0.152;
end

Xtag1 = Xtag4;
Xtag1(1, :) = Xtag1(1, :) + 0.152;

Xtag2 = Xtag4;
Xtag2(1:2,:) = Xtag4(1:2,:) + 0.152;

Xtag3 = Xtag4;
Xtag3(2, :) = Xtag3(2, :) + 0.152;

Xtagc = Xtag4;
Xtagc(1:2,:) = Xtag4(1:2,:) + 0.152/2;

estimate_H_handle = @(sensor) estimate_H(sensor, Kinv, Hr, Xtag1, Xtag2, Xtag3, Xtag4, Xtagc);

estimate_vel_handle = @(sensor) estimate_vel(sensor, Kinv, estimate_H_handle, Tr);
