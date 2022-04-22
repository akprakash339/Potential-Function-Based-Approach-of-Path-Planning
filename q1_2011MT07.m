clc;
clear all;
close all;

%% read image 
I = rgb2gray(imread('cs.png'));
imshow(I);

%% Parameters
q_goal = [20,50];
zeta = 0.1;
dstar = 8;
eta = 10;
qstar = 8;
N=25;
%neighbors = [1 0;1 1;0 1;-1 1;-1 0;-1 -1; 0 -1;1 -1];  %define pixel connectivity
neighbors = [1 0;0 1;-1 0; 0 -1];
%% Attractive Potential

figure;
[X,Y] = meshgrid(1:100,1:100);

U_att = zeros(size(X));
gradU_att = zeros([size(X),2]);

for i = 1:100
    for j = 1:100
        q = [j,i];
        if norm(q-q_goal)<= dstar
            U_att(i,j) = 0.5*zeta*norm(q-q_goal)^2;
            gradU_att(i,j,:) = zeta*(q-q_goal);
        else
            U_att(i,j) = dstar*zeta*norm(q-q_goal) - 0.5*zeta*dstar^2;
            gradU_att(i,j,:) = dstar*zeta*(q-q_goal)/norm(q-q_goal);
        end       
    end
end
[gradU_att(:,:,1),gradU_att(:,:,2)] = gradient(U_att);
subplot(1,2,1);
title('Attracive potential');
surf(X,Y,U_att);
subplot(1,2,2);
title('Gradient of Attractive Potential');
hold on;
imshow(I);
quiver(X,Y,gradU_att(:,:,1),gradU_att(:,:,2));
hold off;
axis([1 100 1 100]);
daspect([1 1 1]);

%% Brushfire 

figure;
subplot(1,3,1);
title('Brushfire');
Inorm = (I/255);
world = zeros(N,N);
hold on
for i=1:(100/N):100
    for j=1:(100/N):100
        if min(min(Inorm(j:j+(100/N)-1, i:i+(100/N)-1)))==0
            hold on
            fill([i i+(100/N)-1 i+(100/N)-1 i],[j j j+(100/N)-1 j+(100/N)-1],'r', 'facealpha', 0.7);
            world(1+floor(i/(100/N)), 1+floor(j/(100/N)))=1;     
            text(i-1+(100/N)/2,j-1+(100/N)/2, '1');
            hold off
        else
            hold on
            fill([i i+(100/N)-1 i+(100/N)-1 i],[j j j+(100/N)-1 j+(100/N)-1],'g', 'facealpha', 0.1);
            world(1+floor(i/(100/N)), 1+floor(j/(100/N)))=0;
            hold off
        end        
    end
end
hold off

counter = 1;
while min(min(world))==0 
    for i=1:N
        for j=1:N            
            if world(i,j)==counter
               %scan each neighbor
               for n=1:size(neighbors,1)
                   if i+neighbors(n,1)>0 && j+neighbors(n,2)>0 &&...
                      i+neighbors(n,1)<=N && j+neighbors(n,2)<=N    
                       if world(i+neighbors(n,1), j+neighbors(n,2))==0
                           world(i+neighbors(n,1), j+neighbors(n,2)) = ...
                               world(i, j)+1;
                             % write the values of pixels
                            hold on
                            text( (100/N)*(i+neighbors(n,1))-(100/N)/2,...
                                (100/N)*(j+neighbors(n,2))-...
                                (100/N)/2, num2str(world((i+neighbors(n,1)), (j+neighbors(n,2)))));
                            hold off
                       end
                   end
               end

            end            
        end
    end
    
    counter = counter + 1;
end
daspect([1 1 1]);

%% Repulsive Potential

world = world';
[dworldx,dworldy] = gradient(world);
[X,Y] = meshgrid(1:100,1:100);
U_rep = zeros(size(X));
gradU_rep = zeros([size(X),2]);

for i=1:99
    for j=1:99
        Dq = world(1+floor(i/(100/N)), 1+floor(j/(100/N)));
        gradDq = [dworldx(1+floor(i/(100/N)), 1+floor(j/(100/N)));
            dworldy(1+floor(i/(100/N)), 1+floor(j/(100/N)))];
        if Dq <= qstar
            U_rep(i,j) = 0.5*eta*((1/Dq)-(1/qstar))^2;
            gradU_rep(i,j,:) = eta*((1/qstar)-(1/Dq))*(1/Dq^2)*gradDq;
        else
            U_rep(i,j) = 0;
            gradU_rep(i,j,:) = [0;0];
        end
    end
end
[gradU_rep(:,:,1),gradU_rep(:,:,2)] = gradient(U_rep);
subplot(1,3,2);
title('Repulsive potential');
surf(X,Y,U_rep);
subplot(1,3,3);
title('Gradient of Repulsive Potential');
hold on;
imshow(I);
quiver(X,Y,gradU_rep(:,:,1),gradU_rep(:,:,2));
hold off;
axis([1 100 1 100]);
daspect([1 1 1]);

%% Total Potential
 
figure;
U = U_att+(5*U_rep);
[gradU(:,:,1),gradU(:,:,2)] = gradient(U);

subplot(1,2,1);
title('Total Potential');
surf(X,Y,U);
subplot(1,2,2);
title('Gradient of Total Potential');
hold on
imshow(I);
quiver(X,Y,gradU(:,:,1),gradU(:,:,2));
hold off

