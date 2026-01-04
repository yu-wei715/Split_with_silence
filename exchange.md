%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test
edge_pts= [];
edge_pts_3d= [];
edge_pts_proj=[];
[x,y] = meshgrid(0:288:2879, 0:288:2879);

%        4
%   |<---------
%   |         ^
% 1 V         | 3
%   --------->|
%        2
for j = 1:size(x,1)-1
    points = [x(j,1), y(j,1)];
    edge_pts = [edge_pts;points];
end
for i = 1:size(x,2)-1
    points = [x(size(x,1),i), y(size(x,1),i)];
    edge_pts = [edge_pts;points];
end
for j = size(x,1):-1:2
    points = [x(j, size(x,2)), y(j, size(x,2))];
    edge_pts = [edge_pts;points];
end
for i = size(x,2):-1:2
    points = [x(1,i), y(1,i)];
    edge_pts = [edge_pts;points];
end

for i = 1:size(edge_pts,1)
    points = edge_pts(i,:);
    uv = 1\(points-[1454, 1439])';
    rho = vecnorm(uv,2,1);
    D = [863.3761 -4.921e-4 3.1235e-7 -2.0134e-10];
    Zc = D(1)+D(2)*rho.^2+D(3)*rho.^3+D(4)*rho.^4;
    Xc = uv(1,:);Yc = uv(2,:);
    worldPoints = [Xc;Yc;Zc];
    nw = vecnorm(worldPoints,2,1); % lambda
    nw(nw == 0) = eps;
    uvz = worldPoints ./ [nw;nw;nw];
    plot3(uvz(1), uvz(2), uvz(3), 'b*'); hold on;
    edge_pts_3d = [edge_pts_3d;uvz'];
end

for i=1:size(edge_pts_3d,1)
    theta = acos(edge_pts_3d(i,2));      % 與北極夾角
    phi   = atan2(-edge_pts_3d(i,1), -edge_pts_3d(i,3));  % 方位角
    new_x = theta * cos(phi);
    new_y = theta * sin(phi);
    edge_pts_proj =[[new_x new_y]; edge_pts_proj]
end

% %make the sphere
for j = 0:0.2:pi
    for i = 0:0.2:2*pi
        x= -sin(j)*sin(i);
        y = -cos(j);
        z = -sin(j)*cos(i);
        point = [x y z];
        is_front = in_edge(point, edge_pts_3d)
        if is_front
            plot3(x, y, z, 'r*'); hold on;
        else
            plot3(x, y, z, 'g*'); hold on;
        end

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% in edge
function inside = in_edge(point, edge)
    sum = 0;
    polar = [0 0 1];
    % for i=1:size(edge,1)
        % next_idx = mod((i),size(edge,1))+1;
        % dot_val = dot(edge(i, :),point);
        % da = dot_val*point
        % a =  edge(i,:) - da;
        % dot_val = dot(edge(next_idx, :),point);
        % db = dot_val*point
        % b =  edge(next_idx,:) - db;
        % na = norm(a);
        % nb = norm(b);
        % if na < eps || nb < eps
        %     continue;
        % end
        % a = a / na;
        % b = b / nb;
        % theta = atan2(dot(point, cross(a,b)), dot(a,b));
        % sum = sum + theta;
    % end
    polar = [0 0 1];
    for i=1:size(edge,1)
        theta = acos(edge(i,3));      % 與北極夾角
        phi   = atan2(edge(i,2), edge(i,1));  % 方位角
    
        % 投影到平面
        new_x = theta * cos(phi);
        new_y = theta * sin(phi);
        plot(new_x, new_y, 'b*'); hold on;
    end


    if sum > 1e-16
        inside = 0;
    else
        inside = 1;
    end
end
