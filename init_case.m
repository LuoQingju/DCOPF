%% 初始化算例数据

define_constants; % Defines useful constants for indexing data
mpc = ext2int(mpc0);

baseMVA = mpc.baseMVA;
bus = mpc.bus;
gen = mpc.gen;
branch = mpc.branch;
gencost = mpc.gencost;

nb = size(bus, 1); % 母线数
ng = size(gen, 1); % 机组数
nl = size(branch, 1); % 线路数

slack = find(bus(:, BUS_TYPE) == REF); % 平衡节点

Cg = sparse(gen(:, GEN_BUS), (1:ng)', 1, nb, ng); % connection matrix 母线-发电机关联矩阵
Cl = sparse([branch(:, F_BUS); branch(:, T_BUS)], [(1:nl)'; (1:nl)'], ...
    [ones(nl, 1); -ones(nl, 1)], nb, nl); % connection matrix 母线-线路关联矩阵

% 标幺值 p.u.
Pmin = gen(:, PMIN) ./ baseMVA; % 发电机功率下限
Pmax = gen(:, PMAX) ./ baseMVA; % 发电机功率上限
flow_max = branch(:, RATE_A) / baseMVA; % 线路功率上限
flow_max(flow_max == 0) = Inf; % 修正线路功率上限
Pd = bus(:, PD) ./ baseMVA; % 有功负荷
Gs = bus(:, GS) ./ baseMVA; % Gs, shunt conductance (MW demanded at V = 1.0 p.u.)
Va_ref = bus(slack, VA) ./ 180 .* pi; % 平衡节点电压相角

% 考虑变压器变比
BR_x = branch(:, BR_X); % 支路电抗
tap = branch(:, TAP); % 支路变比
tap(tap == 0) = 1; % 线路等值为变比是1的变压器
BR_x = BR_x .* tap; % 修正支路电抗

% 考虑变压器移相
Pfinj = -(branch(:, SHIFT) * pi / 180) ./ BR_x; %% injected at the from bus ...
% Ptinj = -Pfinj;   %% ... and extracted at the to bus
% Pbusinj = Cl * Pfinj;   %% Pbusinj = Cl * Pfinj; Pbusinj = Cf * Pfinj + Ct * Ptinj;

% 发电成本系数
Qpg = gencost(:, COST) * baseMVA^2; % 二次项发电成本系数
cpg = gencost(:, COST+1) * baseMVA; % 一次项发电成本系数
kpg = gencost(:, COST+2); % 常数项发电成本系数