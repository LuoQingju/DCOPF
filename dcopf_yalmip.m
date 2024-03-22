%% 直流最优潮流 DC Optimal Power Flow

%% YALMIP

% 作者：罗清局
% 邮箱：luoqingju@qq.com
% 华南理工大学电力学院
% 综合智慧能源系统优化运行与控制团队 ISESOOC 爱思科

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

%% 决策变量
Pg = sdpvar(ng, 1, 'full'); % 机组出力
Pl = sdpvar(nl, 1, 'full'); % 线路有功
Va = sdpvar(nb, 1, 'full'); % 电压相角

%% 约束条件
cons = []; % 初始化约束
cons = [cons, Cg * Pg - Cl * Pl == Pd + Gs]; % 节点功率平衡
cons = [cons, Pl == Cl' * Va ./ BR_x + Pfinj]; % 线路功率约束
cons = [cons, Va(slack) == Va_ref]; % 平衡节点电压相角约束
cons = [cons, Pmin <= Pg <= Pmax]; %#ok<*CHAIN> % 发电机功率上下限
cons = [cons, -flow_max <= Pl <= flow_max]; % 线路功率上下限

%% 目标函数
obj = (Pg .* Qpg)' * Pg + cpg' * Pg + sum(kpg); % 发电成本

%% 求解
ops = sdpsettings('verbose', 2, 'solver', 'gurobi');
sol = optimize(cons, obj, ops);

% 分析错误标志
if sol.problem ~= 0
    sol.info
    yalmiperror(sol.problem)
    return
end
obj = value(obj);
Pg = value(Pg);
Pl = value(Pl);
Va = value(Va);

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
