function [matches, row_ids] = SAD(image_2, other_image, corners_2, corners_2_3d, corners_other, corners_other_3d, filter_width)
%SAD Summary of this function goes here
%   Detailed explanation goes here
  
  [rows, cols] = size(image_2);
  matches = zeros(100, 100);
  filter_radius = ((filter_width-1)/2);
  image_2 = padarray(image_2, [5 5], 0,'both');
  other_image = padarray(other_image, [5 5], 0,'both');

  % Loop over corners of image 2
  for i = 1: size(corners_2)
    if corners_2_3d(1, i) == 0
      continue
    end
    x_2 = ceil(corners_2(i) / cols);
    y_2 = mod(corners_2(i), cols);
      % Loop over corners of other image
      for j = 1: size(corners_other)
        if corners_other_3d(1, j) == 0
          continue
        end
        x_other = ceil(corners_other(j) / cols);
        y_other = mod(corners_other(j), cols);
        % Loop over window
        for u = -filter_radius: filter_radius
          for v = -filter_radius: filter_radius
            matches(i,j) = matches(i,j) + abs(image_2(x_2+u, y_2+v) - other_image(x_other+u, y_other+v));
          end
        end
      end
  end
  
  [matches, row_ids] = max(matches);

end

