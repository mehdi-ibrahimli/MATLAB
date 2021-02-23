clear all
clc
image = imread('woodring.png');
image = im2gray(image);

%%
orientation = [45 30];
wavelength = [5,20,40,60];
aspect_ratio = [1, 2, 5];
gb = {};
for aspct = 1:length(aspect_ratio)
    for w = 1:length(wavelength)
        for or = 1:length(orientation)
            
            gamma = aspect_ratio(aspct);
            lambda = wavelength(w);
            theta = orientation(or);
            
            BW = 2;
            sigma_x = lambda/pi*sqrt(log(2)/2)*(2^BW+1)/(2^BW-1);
            sigma_y = sigma_x / gamma;

            nstds = 3;
            xmax = max(abs(nstds * sigma_x * cos(theta)), abs(nstds * sigma_y * sin(theta)));
            xmax = ceil(max(1, xmax));
            ymax = max(abs(nstds * sigma_x * sin(theta)), abs(nstds * sigma_y * cos(theta)));
            ymax = ceil(max(1, ymax));
            xmin = -xmax; ymin = -ymax;
            [x,y] = meshgrid(xmin:xmax, ymin:ymax);
 
            x_theta = x * cos(theta) + y * sin(theta);
            y_theta = -x * sin(theta) + y * cos(theta);
            gab = exp(-(x_theta.^2 + y_theta.^2 * gamma^2)/(2 * sigma_x^2)).*cos((2 * pi * x_theta)/lambda);
            gb{or,w,aspct} = gab;
        end
    end
    %figure, montage(gb(:,:,1),'DisplayRange',[])
end
figure, montage(gb(:,:,1),'DisplayRange',[])
%%
[u,v,z] = size(gb);
filtered = {};
for i = 1:u
    for j = 1:v
        for k = 1:z
        filtered{i,j,k} = imfilter(image, gb{i,j,k});
        end
    end
end
%%
for m = 1:z
    hold on
    figure, montage(filtered(:,:,m),'DisplayRange',[])
end
