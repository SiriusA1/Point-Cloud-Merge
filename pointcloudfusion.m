%% clear workspace
clc
clear all
close all

%% Part 1: applying filters
depth1 = double(imread('depth1.png'));
depth2 = double(imread('depth2.png'));
depth3 = double(imread('depth3.png'));

rgb1 = imread('rgb1.png');
rgb2 = imread('rgb2.png');
rgb3 = imread('rgb3.png');

gray_1 = double(rgb2gray(rgb1));
gray_2 = double(rgb2gray(rgb2));
gray_3 = double(rgb2gray(rgb3));

x_deriv = [-1 0 1; -1 0 1; -1 0 1];
y_deriv = x_deriv';

Ix_1 = imfilter(gray_1, x_deriv, "replicate");
Ix_2 = imfilter(gray_2, x_deriv, "replicate");
Ix_3 = imfilter(gray_3, x_deriv, "replicate");

Iy_1 = imfilter(gray_1, y_deriv, "replicate");
Iy_2 = imfilter(gray_2, y_deriv, "replicate");
Iy_3 = imfilter(gray_3, y_deriv, "replicate");

Ixx_1 = rescale(abs(imfilter(Ix_1, x_deriv, "replicate")), 0, 255);
Ixx_2 = rescale(abs(imfilter(Ix_2, x_deriv, "replicate")), 0, 255);
Ixx_3 = rescale(abs(imfilter(Ix_3, x_deriv, "replicate")), 0, 255);

Iyy_1 = rescale(abs(imfilter(Iy_1, y_deriv, "replicate")), 0, 255);
Iyy_2 = rescale(abs(imfilter(Iy_2, y_deriv, "replicate")), 0, 255);
Iyy_3 = rescale(abs(imfilter(Iy_3, y_deriv, "replicate")), 0, 255);

Ixy_1 = rescale(abs(imfilter(Ix_1, y_deriv, "replicate")), 0, 255);
Ixy_2 = rescale(abs(imfilter(Ix_2, y_deriv, "replicate")), 0, 255);
Ixy_3 = rescale(abs(imfilter(Ix_3, y_deriv, "replicate")), 0, 255);

gaussian = [1 4 7 4 1; 4 16 26 16 4; 7 26 41 26 7; 4 16 26 16 4; 1 4 7 4 1] * (1/273);

Ixx_1_gauss = imfilter(Ixx_1, gaussian, "replicate");
Ixx_2_gauss = imfilter(Ixx_2, gaussian, "replicate");
Ixx_3_gauss = imfilter(Ixx_3, gaussian, "replicate");

Iyy_1_gauss = imfilter(Iyy_1, gaussian, "replicate");
Iyy_2_gauss = imfilter(Iyy_2, gaussian, "replicate");
Iyy_3_gauss = imfilter(Iyy_3, gaussian, "replicate");

Ixy_1_gauss = imfilter(Ixy_1, gaussian, "replicate");
Ixy_2_gauss = imfilter(Ixy_2, gaussian, "replicate");
Ixy_3_gauss = imfilter(Ixy_3, gaussian, "replicate");


%% Harris response
k = 0.05;
[rows, cols] = size(Ixx_1_gauss);
det_1 = zeros(rows, cols);
det_2 = zeros(rows, cols);
det_3 = zeros(rows, cols);
trace_1 = zeros(rows, cols);
trace_2 = zeros(rows, cols);
trace_3 = zeros(rows, cols);
for x = 1:rows
    for y = 1:cols
        det_1(x,y) = Ixx_1_gauss(x,y) * Iyy_1_gauss(x,y) - Ixy_1_gauss(x,y) * Ixy_1_gauss(x,y);
        det_2(x,y) = Ixx_2_gauss(x,y) * Iyy_2_gauss(x,y) - Ixy_2_gauss(x,y) * Ixy_2_gauss(x,y);
        det_3(x,y) = Ixx_3_gauss(x,y) * Iyy_3_gauss(x,y) - Ixy_3_gauss(x,y) * Ixy_3_gauss(x,y);
        trace_1(x,y) = k * (Ixx_1_gauss(x,y) + Iyy_1_gauss(x,y))^2;
        trace_2(x,y) = k * (Ixx_2_gauss(x,y) + Iyy_2_gauss(x,y))^2;
        trace_3(x,y) = k * (Ixx_3_gauss(x,y) + Iyy_3_gauss(x,y))^2;
    end
end

    
harris_response_1 = rescale(abs(det_1 - trace_1), 0, 255);
harris_response_2 = rescale(abs(det_2 - trace_2), 0, 255);
harris_response_3 = rescale(abs(det_3 - trace_3), 0, 255);

%imwrite(harris_response_1, "harris_response_1.png")

%% non-maximum suppression
pad_1 = padarray(harris_response_1, [1 1], 0,'both');
pad_2 = padarray(harris_response_2, [1 1], 0,'both');
pad_3 = padarray(harris_response_3, [1 1], 0,'both');
[rows, cols] = size(pad_1);

for x = 2:rows-1
    for y = 2:cols-1
        window_1 = pad_1( (x - 1 : x + 1) , (y - 1 : y + 1) );
        window_2 = pad_2( (x - 1 : x + 1) , (y - 1 : y + 1) );
        window_3 = pad_3( (x - 1 : x + 1) , (y - 1 : y + 1) );
        [M_1,I_1] = max(window_1, [], "all");
        [M_2,I_2] = max(window_2, [], "all");
        [M_3,I_3] = max(window_3, [], "all");
        if M_1 ~= pad_1(x, y) || I_1 ~= 5
            pad_1(x, y) = 0;
        end
        if M_2 ~= pad_2(x, y) || I_2 ~= 5
            pad_2(x, y) = 0;
        end
        if M_3 ~= pad_3(x, y) || I_3 ~= 5
            pad_3(x, y) = 0;
        end
    end
end

suppressed_1 = pad_1( (2 : rows-1), (2 : cols-1) );
suppressed_2 = pad_2( (2 : rows-1), (2 : cols-1) );
suppressed_3 = pad_3( (2 : rows-1), (2 : cols-1) );
%imwrite(suppressed_1, "suppressed_1.png");

%% keep top 100
[rows, cols] = size(suppressed_1);
[sorted_1, ids_1] = sort(suppressed_1(:), "descend");
[sorted_2, ids_2] = sort(suppressed_2(:), "descend");
[sorted_3, ids_3] = sort(suppressed_3(:), "descend");
top_100_1 = ids_1(1:100);
top_100_2 = ids_2(1:100);
top_100_3 = ids_3(1:100);

for i = 101:(rows*cols)
    suppressed_1(ids_1(i)) = 0;
    suppressed_2(ids_2(i)) = 0;
    suppressed_3(ids_3(i)) = 0;
end

%imwrite(suppressed_1, "top_100_1.png");

%% Part 2: Corners to 3D points
K = [525.0 0.0 329.5; 0.0 525.0 239.5; 0.0 0.0 1.0];

corners_1 = zeros(3, 100);
corners_2 = zeros(3, 100);
corners_3 = zeros(3, 100);
for i = 1:100
    x_1 = ceil(ids_1(i) / cols);
    y_1 = mod(ids_1(i), cols);
    x_2 = ceil(ids_2(i) / cols);
    y_2 = mod(ids_2(i), cols);
    x_3 = ceil(ids_3(i) / cols);
    y_3 = mod(ids_3(i), cols);
    if depth1(ids_1(i)) ~= 0
        corners_1(:, i) = (1/5000) * depth1(ids_1(i)) * K\[x, y, 1]';
    end
    if depth2(ids_2(i)) ~= 0
        corners_2(:, i) = (1/5000) * depth2(ids_2(i)) * K\[x, y, 1]';
    end
    if depth3(ids_3(i)) ~= 0
        corners_3(:, i) = (1/5000) * depth3(ids_3(i)) * K\[x, y, 1]';
    end
end

%% Part 3: Corner Matching
rt_image_1 = ranktransform(gray_1, 5);
rt_image_2 = ranktransform(gray_2, 5);
rt_image_3 = ranktransform(gray_3, 5);

%imwrite(rt_image_1, "rt_image_1.png");

[matches_21, match_ids_1] = SAD(rt_image_2, rt_image_1, top_100_2, corners_2, top_100_1, corners_1, 11);
[matches_23, match_ids_3] = SAD(rt_image_2, rt_image_3, top_100_2, corners_2, top_100_3, corners_3, 11);
[sorted_matches_21, match_ids_21] = sort(matches_21, "descend");
[sorted_matches_23, match_ids_23] = sort(matches_23, "descend");
top_10_21 = match_ids_21(1:10);
top_10_23 = match_ids_23(1:10);

%% Part 4: Pose Estimation
% gets really messy here so I'll give an example:
% corners_1(match_ids_1(top_10_21(1))) is the 3d points of the best
% matching corner from image 1 that matches to image 2
v_11 = corners_1(:, match_ids_1(top_10_21(1))) - corners_1(:, match_ids_1(top_10_21(5)));
v_12 = corners_1(:, match_ids_1(top_10_21(5))) - corners_1(:, match_ids_1(match_ids_21(14)));
v_21 = corners_2(:, top_10_21(1)) - corners_2(:, top_10_21(5));
v_22 = corners_2(:, top_10_21(5)) - corners_2(:, match_ids_21(14));
R1_prime = [v_21, v_22, (v_21.*v_22)]/[v_11, v_12, (v_11.*v_12)];
[U,S,V] = svd(R1_prime);
R1 = U*V';
t1 = corners_2(:, top_10_21(1)) - R1*corners_1(:, match_ids_1(top_10_21(1)));

v_31 = corners_3(:, match_ids_3(top_10_23(6))) - corners_3(:, match_ids_3(top_10_23(7)));
v_32 = corners_3(:, match_ids_3(top_10_23(7))) - corners_3(:, match_ids_3(match_ids_23(12)));
v_21 = corners_2(:, top_10_23(6)) - corners_2(:, top_10_23(7));
v_22 = corners_2(:, top_10_23(7)) - corners_2(:, match_ids_23(12));
R3_prime = [v_21 v_22 (v_21.*v_22)]/[v_31 v_32 (v_31.*v_32)];
[U,S,V] = svd(R3_prime);
R3 = U*V';
t3 = corners_2(:, top_10_23(6)) - R3*corners_3(:, match_ids_3(top_10_23(6)));

%% Part 5: Merge images
coords_3d = zeros(921600, 3);
colors = zeros(921600, 3);
[rows, cols] = size(gray_1);
i = 1;
for x = 1:rows
    for y = 1:cols
        if depth2(x, y) ~= 0
            coords_3d(i, :) = (1/5000) * depth2(x, y) * K\[x, y, 1]';
            colors(i, :) = rgb2(x, y, :);
            i = i + 1;
        end
        if depth1(x, y) ~= 0
            temp_1 = (1/5000) * depth1(x, y) * K\[x, y, 1]';
            coords_3d(i, :) = R1*temp_1 + t1;
            colors(i, :) = rgb1(x, y, :);
            i = i + 1;
        end
        if depth3(x, y) ~= 0
            temp_3 = (1/5000) * depth3(x, y) * K\[x, y, 1]';
            coords_3d(i, :) = R3*temp_3 + t3;
            colors(i, :) = rgb3(x, y, :);
            i = i + 1;
        end
    end
end

coloredPtCloud = pointCloud(coords_3d);
coloredPtCloud.Color = uint8(colors);
pcwrite(coloredPtCloud,'coloredPtCloud','PLYFormat','ascii');

%% Problem 2
[rows, cols] = size(gray_1);
image1_3d = zeros(rows, cols, 3);
image2_3d = zeros(rows, cols, 3);
image3_3d = zeros(rows, cols, 3);
for x = 1:rows
    for y = 1:cols
        if depth1(x, y) ~= 0
            image1_3d(x, y, :) = (1/5000) * depth1(x, y) * K\[x, y, 1]';
        end
        if depth2(x, y) ~= 0
            image2_3d(x, y, :) = (1/5000) * depth2(x, y) * K\[x, y, 1]';
        end
        if depth3(x, y) ~= 0
            image3_3d(x, y, :) = (1/5000) * depth3(x, y) * K\[x, y, 1]';
        end
    end
end

%% Find normals
image1_norms = zeros(rows, cols, 3);
image2_norms = zeros(rows, cols, 3);
image3_norms = zeros(rows, cols, 3);
c = [0.5; 0.5; 0.5];
for x = 4:rows-3
    for y = 4:cols-3
        % image 1
        i = 1;
        points = zeros(3, 3);
        for u = -3:3
            for v = -3:3
                if i == 4
                    normal = cross(points(1, :)-points(2, :), points(1, :)-points(3, :));
                    d = -points(1,:)*normal';
                    if d < 0
                        normal = -normal;
                    end
                    color_vec = ((normal' / 2*norm(normal')) + c)*255;
                    image1_norms(x, y, :) = color_vec;
                    i = i + 1;
                elseif image1_3d(x+u, y+v, 1) ~= 0 && i < 4
                    points(i, :) = image1_3d(x+u, y+v, :);
                    i = i + 1;
                end
            end
        end

        % image 2
        i = 1;
        points = zeros(3, 3);
        for u = -3:3
            for v = -3:3
                if i == 4
                    normal = cross(points(1, :)-points(2, :), points(1, :)-points(3, :));
                    d = -points(1,:)*normal';
                    if d < 0
                        normal = -normal;
                    end
                    color_vec = ((normal' / 2*norm(normal')) + c)*255;
                    image2_norms(x, y, :) = color_vec;
                    i = i + 1;
                elseif image2_3d(x+u, y+v, 1) ~= 0 && i < 4
                    points(i, :) = image2_3d(x+u, y+v, :);
                    i = i + 1;
                end
            end
        end

        % image 3
        i = 1;
        points = zeros(3, 3);
        for u = -3:3
            for v = -3:3
                if i == 4
                    normal = cross(points(1, :)-points(2, :), points(1, :)-points(3, :));
                    d = -points(1,:)*normal';
                    if d < 0
                        normal = -normal;
                    end
                    color_vec = ((normal' / 2*norm(normal')) + c)*255;
                    image3_norms(x, y, :) = color_vec;
                    i = i + 1;
                elseif image3_3d(x+u, y+v, 1) ~= 0 && i < 4
                    points(i, :) = image3_3d(x+u, y+v, :);
                    i = i + 1;
                end
            end
        end
    end
end

%% make images
imwrite(image1_norms, "image1_norms.png");
imwrite(image2_norms, "image2_norms.png");
imwrite(image3_norms, "image3_norms.png");