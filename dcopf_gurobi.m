%% 直流最优潮流 DC Optimal Power Flow

%% GUROBI

% 作者：罗清局
% 邮箱：luoqingju@qq.com

clc
clear

%% Small Transmission System Test Cases
% mpc0 = case6ww;
% mpc0 = case9;
mpc0 = case14;
% mpc0 = case24_ieee_rts;
% mpc0 = case30;
% mpc0 = case_ieee30;
% mpc0 = case39;
% mpc0 = case57;
% mpc0 = case118;
% mpc0 = case145;
% mpc0 = case300;

%% PEGASE European System Test Cases
% mpc0 = case89pegase;
% mpc0 = case1354pegase;
% mpc0 = case2869pegase;
% mpc0 = case9241pegase;
% mpc0 = case13659pegase;

%% RTE French System Test Cases
% mpc0 = case1888rte;
% mpc0 = case1951rte;
% mpc0 = case2848rte;
% mpc0 = case2868rte;
% mpc0 = case6468rte;
% mpc0 = case6470rte;
% mpc0 = case6495rte;
% mpc0 = case6515rte;

%%  Polish System Test Cases
% mpc0 = case2383wp;
% mpc0 = case2736sp;
% mpc0 = case2737sop;
% mpc0 = case2746wop;
% mpc0 = case2746wp;
% mpc0 = case3012wp;
% mpc0 = case3120sp;
% mpc0 = case3375wp;

%%  ACTIV Synthetic Grid Test Cases
% mpc0 = case_ACTIVSg200;
% mpc0 = case_ACTIVSg500;
% mpc0 = case_ACTIVSg2000;
% mpc0 = case_ACTIVSg10k;
% mpc0 = case_ACTIVSg25k;
% mpc0 = case_ACTIVSg70k;

%% 初始化算例数据
init_case;

%% 决策变量索引
nx = 0; % 变量数
idx_x_Pg = nx + (1:ng)'; % 机组出力索引
nx = nx + ng;
idx_x_Pl = nx + (1:nl)'; % 线路有功索引
nx = nx + nl;
idx_x_Va = nx + (1:nb)'; % 电压相角索引
nx = nx + nb;

%% 约束条件索引

ne = 0; % 等式约束数
idx_e_bus = ne + (1:nb)'; % 节点功率平衡约束索引
ne = ne + nb;
idx_e_Pl = ne + (1:nl)'; % 线路功率约束索引
ne = ne + nl;
idx_e_Va_ref = ne + 1; % 平衡节点电压相角约束索引
ne = ne + 1;

ni = 0; % 不等式约束数
% 变量的上下限约束可以直接输入到gurobi中，一般不需要写为不等式约束

%% 约束条件
Ae = sparse(ne, nx); % 等式约束 Ae * x = be
be = zeros(ne, 1); % 等式约束 Ae * x = be
Ai = sparse(ni, nx); % 不等式约束 Ai * x <= bi
bi = zeros(ni, 1); % 不等式约束 Ai * x <= bi

% cons = [cons, Cg * Pg - Cl * Pl == Pd + Gs]; % 节点功率平衡
Ae(idx_e_bus, idx_x_Pg) = Cg;
Ae(idx_e_bus, idx_x_Pl) = -Cl;
be(idx_e_bus) = Pd + Gs;

% cons = [cons, Pl == Cl' * Va ./ BR_x + Pfinj]; % 线路功率约束
Ae(idx_e_Pl, idx_x_Pl) = speye(nl);
Ae(idx_e_Pl, idx_x_Va) = -sparse(1:nl, 1:nl, 1./BR_x, nl, nl) * Cl';
be(idx_e_Pl) = Pfinj;

% cons = [cons, Va(slack) == Va_ref]; % 平衡节点电压相角约束
Ae(idx_e_Va_ref, idx_x_Va(slack)) = 1;
be(idx_e_Va_ref) = Va_ref;

%% 变量上下限
lb = -Inf(nx, 1); % 初始化变量上下限
ub = Inf(nx, 1);

% cons = [cons, Pmin <= Pg <= Pmax]; %#ok<*CHAIN> % 发电机功率上下限
lb(idx_x_Pg) = Pmin; % 发电机功率上下限
ub(idx_x_Pg) = Pmax;

% cons = [cons, -flow_max <= Pl <= flow_max]; % 线路功率上下限
lb(idx_x_Pl) = -flow_max; % 线路功率上下限
ub(idx_x_Pl) = flow_max;

%% 目标函数
% obj = (Pg .* Qpg)' * Pg + cpg' * Pg + sum(kpg); % 发电成本
Q = sparse(nx, nx);
Q(idx_x_Pg, idx_x_Pg) = sparse(1:ng, 1:ng, Qpg, ng, ng);
c = zeros(nx, 1);
c(idx_x_Pg) = cpg;

%% 求解
model.Q = Q;
model.obj = c;
model.A = [Ai; Ae];
model.rhs = [bi; be];
model.sense = [char('<'*ones(ni, 1)); char('='*ones(ne, 1))];
model.vtype = char('c'*ones(nx, 1));
model.modelsense = 'min';
model.lb = lb;
model.ub = ub;
model.objcon = sum(kpg);

result = gurobi(model);
if (strcmp(result.status, 'OPTIMAL') == 1)
    x = result.x;
else
    result.status
    return
end

obj = result.objval;
Pg = x(idx_x_Pg);
Pl = x(idx_x_Pl);
Va = x(idx_x_Va);

%% 与MATPOWER对比
res = rundcopf(mpc0);

disp('系统总成本绝对误差：')
disp(norm(res.f-obj))
disp('系统总成本相对误差：')
disp(norm((res.f - obj)/res.f))

% 可能有多个最优解，所以最优解可能不一样
% disp('发电机出力绝对误差：')
% disp(norm( res.gen(res.gen(:, GEN_STATUS)>0, PG)./mpc0.baseMVA - Pg ))
% disp('线路有功功率绝对误差：')
% disp(norm( res.branch(res.branch(:, BR_STATUS)==1, PF)./mpc0.baseMVA - Pl ))
% disp('电压相角绝对误差：')
% disp(norm( res.bus(:, VA)./180.*pi - Va ))
