
% % Numerical Mathematics
% % Matlab Sheet 4
% % RSGI
% % WS20/21
% % Mehdi Ibrahimli


% this version fixes the thresholding bugs and does not follow the
% matlab sheet steps explicitly for speed

clear all;
clc;
image = imread('rug1.jpg');
image = rgb2gray(image);
image = im2double(image);

%% creating a scale space
% 0. creating sigma dimension as described in the paper (Lowe, 2004)
tic
s = 2;
sigma = 0.707;
k = 2^(1/s);
stack = s + 3;
octave_count = 3;
scale_space = {octave_count,stack};
DoG = {stack - 1,octave_count};
for i = 1:octave_count
    sigma_1 = sigma * 2^(i - 1);
    for j = 1:stack
        k_layer = k^(j-1);
        sigma_layer = sigma_1 * k_layer;
        scale_space{i,j} = imgaussfilt(image,sigma_layer);
        if j >= 2
            DoG{j-1,i} = scale_space{i,j} - scale_space{i,j - 1};
        end
    end
end

%% finding interest points
% 1. finding local min and max
keypoint = [];
counter = 1;

for si_oct = 1:octave_count
    octave_stacks = cat(3,DoG{:,si_oct});
    [a , b] = size(DoG{1,si_oct});
    for si_s = 2:stack - 2
        for i = 2:a-1
            for j = 2:b-1
                cp = octave_stacks(i,j,si_s);
                if (cp<octave_stacks(i-1,j-1,si_s-1) && cp<octave_stacks(i-1,j,si_s-1)  && cp<octave_stacks(i-1,j+1,si_s-1) && cp<octave_stacks(i,j-1,si_s-1) && cp<octave_stacks(i,j,si_s-1) && cp<octave_stacks(i,j+1,si_s-1) && cp<octave_stacks(i+1,j-1,si_s-1)...
                        && cp<octave_stacks(i+1,j,si_s-1) && cp<octave_stacks(i+1,j+1,si_s-1) && cp<octave_stacks(i-1,j-1,si_s) && cp<octave_stacks(i-1,j,si_s) && cp<octave_stacks(i-1,j+1,si_s) && cp<octave_stacks(i,j-1,si_s) && cp<octave_stacks(i,j+1,si_s)...
                        && cp<octave_stacks(i+1,j-1,si_s) && cp<octave_stacks(i+1,j,si_s) && cp<octave_stacks(i+1,j+1,si_s) && cp<octave_stacks(i-1,j-1,si_s+1) && cp<octave_stacks(i-1,j,si_s+1) && cp<octave_stacks(i-1,j+1,si_s+1) && cp<octave_stacks(i,j-1,si_s+1)...
                        && cp<octave_stacks(i,j,si_s+1) && cp<octave_stacks(i,j+1,si_s+1) && cp<octave_stacks(i+1,j-1,si_s+1) && cp<octave_stacks(i+1,j,si_s+1) && cp<octave_stacks(i+1,j+1,si_s+1))|| ...
                        (cp>octave_stacks(i-1,j-1,si_s-1) && cp>octave_stacks(i-1,j,si_s-1)  && cp>octave_stacks(i-1,j+1,si_s-1) && cp>octave_stacks(i,j-1,si_s-1) && cp>octave_stacks(i,j,si_s-1) && cp>octave_stacks(i,j+1,si_s-1) && cp>octave_stacks(i+1,j-1,si_s-1)...
                        && cp>octave_stacks(i+1,j,si_s-1) && cp>octave_stacks(i+1,j+1,si_s-1) && cp>octave_stacks(i-1,j-1,si_s) && cp>octave_stacks(i-1,j,si_s) && cp>octave_stacks(i-1,j+1,si_s) && cp>octave_stacks(i,j-1,si_s) && cp>octave_stacks(i,j+1,si_s)...
                        && cp>octave_stacks(i+1,j-1,si_s) && cp>octave_stacks(i+1,j,si_s) && cp>octave_stacks(i+1,j+1,si_s) && cp>octave_stacks(i-1,j-1,si_s+1) && cp>octave_stacks(i-1,j,si_s+1) && cp>octave_stacks(i-1,j+1,si_s+1) && cp>octave_stacks(i,j-1,si_s+1)...
                        && cp>octave_stacks(i,j,si_s+1) && cp>octave_stacks(i,j+1,si_s+1) && cp>octave_stacks(i+1,j-1,si_s+1) && cp>octave_stacks(i+1,j,si_s+1) && cp>octave_stacks(i+1,j+1,si_s+1))
                    keypoint(counter,:)=[i,j,si_s,si_oct];
                    counter = counter+1;
                 end
            end
        end
    end
end
%% taylor interpolation, accurate position and edge elimination
big_DOG = cat(3,DoG{:,:});
X = [];
r = 10;
cntr = 1;
for kp = 1:length(keypoint)
    x = keypoint(kp,1);
    y = keypoint(kp,2);
    sig = (keypoint(kp,4)-1) * (stack -1) + keypoint(kp,3);
  % increment
    x1 = (big_DOG(x+1,y,sig) - big_DOG(x,y,sig)) / (big_DOG(x+1,y,sig) - 2*big_DOG(x,y,sig) + big_DOG(x-1,y,sig));
    y1 = (big_DOG(x,y+1,sig) - big_DOG(x,y,sig)) / (big_DOG(x,y+1,sig) - 2*big_DOG(x,y,sig) + big_DOG(x,y-1,sig));
  sig1 = (big_DOG(x,y,sig+1) - big_DOG(x,y,sig)) / (big_DOG(x,y,sig+1) - 2*big_DOG(x,y,sig) + big_DOG(x,y,sig-1));
  % hesse
    H = zeros(2);
    H(1,1) =  (big_DOG(x+1,y,sig) - 2*big_DOG(x,y,sig) + big_DOG(x-1,y,sig));
    H(1,2) =  (big_DOG(x+1,y+1,sig) - big_DOG(x+1,y,sig) - big_DOG(x,y+1,sig) + 2*big_DOG(x,y,sig) - big_DOG(x-1,y,sig) - big_DOG(x,y-1,sig) + big_DOG(x-1,y-1,sig));
    H(2,1) =  H(1,2);
    H(2,2) =  (big_DOG(x,y+1,sig) - 2*big_DOG(x,y,sig) + big_DOG(x,y-1,sig));
    % plug the accurate location into taylor expansion
    res = big_DOG(x,y,sig) + (big_DOG(x+1,y,sig) - big_DOG(x,y,sig)) * x1 + ((big_DOG(x,y+1,sig) - big_DOG(x,y,sig)) * y1)...
        + ((big_DOG(x,y,sig+1) - big_DOG(x,y,sig)) * sig1) + 1/2*(big_DOG(x+1,y,sig) - 2*big_DOG(x,y,sig) + big_DOG(x-1,y,sig))*x1^2 ...
        + 1/2*(big_DOG(x,y+1,sig) - 2*big_DOG(x,y,sig) + big_DOG(x,y-1,sig))*y1^2 + 1/2*(big_DOG(x,y,sig+1) - 2*big_DOG(x,y,sig) + big_DOG(x,y,sig-1))*sig1^2 ...
        + (big_DOG(x+1,y+1,sig) - big_DOG(x+1,y,sig) - big_DOG(x,y+1,sig) + 2*big_DOG(x,y,sig) - big_DOG(x-1,y,sig) - big_DOG(x,y-1,sig) + big_DOG(x-1,y-1,sig))*x1*y1...
        + (big_DOG(x+1,y,sig+1) - big_DOG(x+1,y,sig) - big_DOG(x,y,sig+1) + 2*big_DOG(x,y,sig) - big_DOG(x-1,y,sig) - big_DOG(x,y,sig-1) + big_DOG(x-1,y,sig-1))*x1*sig1...
        + (big_DOG(x,y+1,sig+1) - big_DOG(x,y+1,sig) - big_DOG(x,y,sig+1) + 2*big_DOG(x,y,sig) - big_DOG(x,y-1,sig) - big_DOG(x,y,sig-1) + big_DOG(x,y-1,sig-1))*y1*sig1;
    
    % eliminating points that don't satisfy 3. and 4.
    
    if (abs(res)> 0.003) && (((H(1,1)+H(2,2))^2/det(H)) < ((r+1)^2)/r)
         X(cntr,1) = keypoint(kp,1) + x1;
         X(cntr,2) = keypoint(kp,2) + y1;
         X(cntr,3) = keypoint(kp,3) + sig1;
         cntr = cntr + 1;
    end
end
toc
%% results
figure, imshow(image)
hold on
axis on
plot(X(:,2), X(:,1),'r.');



