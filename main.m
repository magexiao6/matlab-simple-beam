%% 定义基本信息
E = 210e9; % 弹性模量（pa）
h = 2; % 梁高（m）
l = 10; % 梁长（m）
b = 1; % 梁宽（m）
I = 1/12*b*h^3; % 惯性矩（m^4） 矩形截面

%% 有限元模型参数
n_elements = 10; % 单元数量
n_nodes = n_elements + 1; % 节点数量
l_elem = l/ n_elements; % 单元长度

% 在txt文件中写入节点坐标
file_nodes = fopen("location of nodes.txt", "w");
no_node = 1;
for i = 1: n_nodes
    x = l_elem* (i-1);
    y = 0;
    fprintf(file_nodes, "%d, %f, %f\n", no_node, x, y);
    no_node = no_node+ 1;
end
fclose(file_nodes);

% 在txt文件中写入单元坐标
file_elems = fopen("location of elems.txt", "w");
no_elem = 1;
for i = 1: n_elements
    fprintf(file_elems, "%d, %d\n", i, i+1);
    no_elem = no_elem+ 1;
end
fclose(file_elems);

%% 刚度矩阵
% 单元刚度矩阵
k = E*I/ (l_elem^3)* [12 6*l_elem -12 6*l_elem;
                      6*l_elem 4*l_elem^2 -6*l_elem 2*l_elem^2;
                      -12 -6*l_elem 12 -6*l_elem;
                      6*l_elem 2*l_elem^2 -6*l_elem 4*l_elem^2];

% 全局刚度矩阵
K = zeros(2*n_nodes, 2*n_nodes);

% 组装全局刚度矩阵
elems = load('location of elems.txt');
for i = 1: n_elements
    a = elems(i, 1);
    b = elems(i, 2);
    dof(1) = 2*a-1; % dof是向量
    dof(2) = 2*a;
    dof(3) = 2*b-1;
    dof(4) = 2*b;
    for j = 1: 4
        for m = 1: 4
            K(dof(j), dof(m)) = K(dof(j), dof(m))+ k(j, m);
        end
    end
end

%% 设置边界条件（对角线置1法）
% 边界条件：第1节点和最后1个节点的竖向位移和转角为0
fixed_dofs = [1, 2*n_nodes-1];
for i = 1:length(fixed_dofs)
    dof = fixed_dofs(i);
    K(dof, :) = 0;
    K(:, dof) = 0;
    K(dof, dof) = 1;
end

%% 施加外力
F = zeros(2*n_nodes, 1); % 竖直方向的力和弯矩
n = 9; % 施加力的节点号
F(2*n-1) = -1000; % 竖向力：2n-1；弯矩：2n
F(fixed_dofs) = 0;

%% 计算
U = K\ F;

% 求解单元力（剪力1、弯矩1、剪力2、弯矩2）
Fe = zeros(4, n_elements);
for i = 1: n_elements
    Fe(:, i) = k* U((i-1)*2+1: (i-1)*2+4);
end

% 求解节点力
f = K* U;

% 变形后位置
nodes = load("location of nodes.txt");
new_nodes = nodes(:, end-1: end);
nodes_end = new_nodes;
for i = 1: n_nodes
    % 变形后竖直位置
    nodes_end(i, 2) = nodes_end(i, 2) + U(2*i-1, 1)*1e7;
end

%% 绘制位移
% 绘制初始形状
figure;
plot(nodes(:, 2), nodes(:, 3), 'b-o', 'LineWidth', 2, 'MarkerSize', 5);
hold on;
grid on;
axis equal;
xlabel('X (m)');
ylabel('Y (m)');
title('梁的初始形状和变形后形状');
% 绘制变形后的形状
plot(nodes_end(:, 1), nodes_end(:, 2), 'r-o', 'LineWidth', 2, 'MarkerSize', 5);
% 添加节点编号
for i = 1: n_nodes
    text(nodes(i, 2), nodes(i, 3), sprintf(' %d', i), 'Color', 'b', 'VerticalAlignment', 'bottom');
    text(nodes_end(i, 1), nodes_end(i, 2), sprintf(' %d', i), 'Color', 'r', 'VerticalAlignment', 'top');
end
hold off;

% 绘制弯矩图
figure;
plot(nodes(:, 2), nodes(:, 3), 'b', 'LineWidth', 2, 'MarkerSize', 5);
hold on;
grid on;
axis equal;
M = Fe(2,:);
M(n_elements+1) = -Fe(4,end);
plot(nodes(:, 2), M/1000, 'r', 'LineWidth', 2, 'MarkerSize', 5);
title('弯矩图')
xlabel('X (m)');
ylabel('N.m');
hold off;

% 绘制剪力图
figure;
plot(nodes(:, 2), nodes(:, 3), 'b', 'LineWidth', 2, 'MarkerSize', 5);
hold on;
grid on;
axis equal;
P = Fe(1, :);
P(n_elements+1) = -Fe(3,end);
stairs(nodes(:, 2), P/1000, 'r', 'LineWidth', 2, 'MarkerSize', 5);
title('剪力图')
xlabel('X (m)');
ylabel('N');

hold off;
