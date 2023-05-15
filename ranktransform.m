function rt_image = ranktransform(image, window_width)
%RANKTRANSFORM Summary of this function goes here
%   Detailed explanation goes here
    %% pad image
    window_radius = ((window_width-1)/2);
    rt_image = zeros(size(image));
    [H,W] = size(image);

    %% Loop over internal image in steps of 5
    for y = window_radius+1: 5: H - window_radius
        for x = window_radius+1: 5: W - window_radius
            % Loop over window
            pixel_vals = zeros(window_width*window_width, 1);

            % make list of all pixel values in window
            i = 1;
            for v = -window_radius: window_radius
                for u = -window_radius: window_radius
                    pixel_vals(i) = image(y+v, x+u);
                    i = i+1;
                end
            end
      
            % create a sorted list from pixel values retrieved above
            ranked_pixels = sort(unique(pixel_vals));

            % step through window again, replacing pixel values with index ("rank") of pixel val from above list
            for v = -window_radius: window_radius
                for u = -window_radius: window_radius
                    rt_image(y+v, x+u) = find(ranked_pixels == image(y+v, x+u));
                end
            end
        
        end
    end

end

