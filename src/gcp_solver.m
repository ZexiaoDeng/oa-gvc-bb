function [ x, fval, output ] = gcp_solver( GCP )
%   GCP_SOLVER  基于锥形细分的分支定界求解器算法
%   GCP_SOLVER  用于求解一般凹最小化问题 ( general concave minimization problem, GCP ), 形如
%
%       min     f( x )
%       s.t.    x in D
%
%       f( x ) is a concave function
%       D = {x| g(x) <= 0 } is a compact convex set
%       g( (x) is a convex function.    
%
%   GCP_SOLVER  实现了如下算法:
%       (1) 锥形细分 ( conical subdivision )
%       (2) gamma-有效割平面 ( gamma-valid cuts )
%       (3) 外逼近算法( outer approximation )
%       (4) 分支定界框架 ( branch-and-bound )
%   输入:
%       GCP : 一般凹最小化问题
%
%   输出:
%       x       :   最优解
%       fval    :   最优目标函数值
%       output  :   迭代 verbose 信息记录
%
%    see also 
%       全局优化引论, R. Horst, P.M. Pardalos, N.V. Thoai 著, 清华大学出版社, P148
%       Convex Analysis and Global Optimization, 2nd Edition, Hoang Tuy,
%       Spring, 2016, P184~P188
%
% Copyright (C) 2021-2025 Zexiao Deng
%

path = './bt-1.3' ;
addpath( path ) ;

epsilon = 1e-9       ;   % 绝对误差限

D       = GCP.D      ;   % GCP 问题的凸集, mosek 中的 ACC form
oracle  = GCP.oracle ;   % 目标函数
Pk      = GCP.Pk     ;   % 初始覆盖凸集 D 的多面体 Pk

cst     = D.cst   ;      % 凸集 D 的非线性凸约束( 凸集约束 )
subgrad = D.subgrad ;    % 凸集 D 非线性凸约束的次梯度

% ========================
% 初始化
% ========================
% (1) 取在列紧凸集 D 中一个内接的单纯形 S = [ s1, ..., s(n+1) ]
% (2) 取在单纯形 S 中一个内点 x0
S      = GCP.S  ;       % 列紧凸集 D 中一个内接的单纯形
x0     = GCP.x0 ;       % 单纯形 S 中一个内点 x0

% (3) 计算从 x0 发出到单纯形顶点 si 的极射线
%     与列紧凸集 D 的边界交点为 sBari
sBar  = zeros( size( x0, 1 ), size( S, 2 ) ) ;
Theta = zeros( size( S, 2 ), 1 ) ;
fsBar = zeros( size( S, 2 ), 1 ) ;

for idx = 1: size( S, 2 )
    ui      = S( :, idx ) - x0 ;
    options = optimoptions( 'fmincon'  , ...
                            'Display'  ,'off', ...
                            'Algorithm','sqp-legacy' ) ;
    thetai  = fmincon( @(t) -t    , 0                   , ...
                       D.Aineq*ui , D.bineq - D.Aineq*x0, ...
                       []         , []                  , ...
                       []         , []                  , ...
                       @fcst      , options             ) ;   % 迭代初值为 0
   
    sBari   = x0 + thetai*ui ;

    sBar(    :, idx ) = sBari  ;
    Theta( idx, 1   ) = thetai ;
    fsBar( idx, 1   ) = feval( oracle, sBari ) ;
end

plot( sBar( 1, : ), sBar( 2, : ), 'ro' )

% (4) 令 xBar1 in argmin{ f(sBar1), ..., f(sBarn+1) }
%        gamma1 = f( xBar1 )
%     为当前最好可行解: ( xBar )
%     和当前最好解目标值 : ( uppper_bound )
%       ( 局部上边界 )
[ fxBar, idxopt ] = min( fsBar ) ;
xBar        = sBar( :, idxopt ) ;       % 当前最好解
upper_bound = fxBar             ;       % gamma, 当前最好解目标值

% ================================
% (5) 计算以点 x0为起点,
%     过 sBarj, j != i 的 n 条边构成的锥 Mi
% (6) 候选活跃节点 List 配置
candidate_list = cell( size( sBar, 2 ), 1 ) ;
for idx = 1: size( sBar, 2 )
    
    M.V = x0 ;                                  % 锥 M 的顶点
    if idx <= size( sBar, 2 ) - 1               % 锥 M 的 n 个方向
        M.D = sBar( :, idx: idx + 1 ) - x0  ;
    else
        M.D = [ sBar( :, idx ), sBar( :, 1 ) ] - x0 ;
    end
    
    node.model = GCP ;
    node.M     = M   ;
    node.Pk    = Pk  ;          % 当前多面体 Pk
    
    candidate_list( idx )   = { node } ;    % 候选活跃节点集合
end

% 记录使用
LB = [] ;
UB = [] ;
k  = 0 ;

while ( ~isempty( candidate_list ) && k <= 50 )
    
    % 选择节点, 深度优先搜索, 后进先出
    [ node, candidate_list ] = node_choice_DFS( candidate_list ) ;
    
    % =========================
    % 定界操作( bounding )
    % =========================
    M = node.M ;
    
    % 通过 gamma 扩张计算 thetai
    Theta = zeros( size( M.D, 2 ), 1              ) ;
    Z     = zeros( size( M.V, 1 ), size( M.D, 2 ) ) ;
    for idx = 1: size( M.D, 2 )
        [ zi, thetai ] = gamma_extension( upper_bound  , ...
                                          M.V          , ...
                                          M.D( :, idx ), ...
                                          oracle       ) ;
        Z( :      , idx ) = zi     ;	% 锥 M 与集合 C = {x| f(x) > gamma } 边界的交点集合
        Theta( idx, 1   ) = thetai ;    % 扩张倍数 theta
    end
    
    % 求解线性规划问题 LP( Pk, M )
    % M : x = x0 + U*t >= 0 ;
    % Pk: A*x <= b
    n      = size( x0, 1 ) ;
    ATilde = [ node.Pk.Aineq ; -eye( n ) ] ;                    % LP( Pk, M ) 的系统矩阵
    bTilde = [ node.Pk.bineq - node.Pk.Aineq*x0 ; x0 ; ] ;      % LP( Pk, M ) 的右向量
    U      = M.D ;           % 锥 M 的极方向矩阵 U = [ u1, ..., un ]
    
    [ tBar_M , ...           % tBar_M = [ t1 ; ... ; tn ] ;
      mu_M   , ...           % mu( M )
      node.output ] = lpDM_solver( Theta , ...
                                   ATilde, ...
                                   bTilde, ...
                                   U    ) ;
    
    if ( node.output.exitflag <= 0 )
        % 剪枝操作( pruning )
        % 节点的问题是一个不可行问题
        lower_bound = inf ;
        fprintf( 'prune by infeasiblity!\n' ) ;
        continue ;
    end
    
    % 计算 f( Pk cap M ) 的下边界 beta_M = lower_bound
    if mu_M <= 1        
        % 剪枝操作( pruning )
        fprintf( 'prune by mu_M <= 1!\n' ) ;
        continue ;
    else                % mu( M ) > 1
        fY = zeros( size( Theta, 1 ), 1 ) ;
        for idx = 1: size( Theta, 1 )
            thetai = Theta( idx ) ;
            ui     = U( :, idx ) ;
            yi     = x0 + mu_M*thetai*ui ;
            fyi    = feval( oracle, yi ) ;
            fY( idx, 1 ) = fyi ;
        end
        lower_bound = min( fY ) ;
    end
    
    if ( lower_bound >= upper_bound - epsilon )
        % 剪枝操作( pruning )
        % 节点的下边界大于当前最好解( 全局上边界 ), 剪枝
        fprintf( 'prune by bound!\n' ) ;
        continue ;
    end
    
    % ====================
    % 停止条件 2
    % ====================
    omega_M  = x0 + U*tBar_M  ;
    fomega_M = feval( oracle, omega_M ) ;
    
    plot( omega_M(1), omega_M(2), 'r*', 'LineWidth', 2 )
    
    if ( ( max( D.Aineq*omega_M - D.bineq ) <= 0 ) && ...
         ( feval( cst, omega_M ) <= 0            )      )
        % 判断点 omega_M 是否在凸集 D 内
        if fomega_M < upper_bound
            % 更新上界( 当前最好解 )操作
            xBar        = omega_M ;
            upper_bound = fomega_M ;
        end
%         break ;
        continue ;
    end
    
    % ====================================
    % 分支操作( branching )
    % ====================================
    % (1) 单纯 omega-细分
    % 以 x0 为起点, 通过 omega( Mk )点, 进行 omega-细分
    sub_M = conical_subdivision( x0           , ...
                                 omega_M - x0, ...
                                 tBar_M       , ...
                                 U            ) ;
    k     = k + 1 ;
    
    % ============================================================
    % 外逼近操作( branching )
    % 割掉 omegaBar_M
    % Pk+1 = Pk cap { x| <pk, x - omegaBar_M> + g( omegaBar_M ) <= 0 }
    % ============================================================
    ui      = omega_M - x0 ;
    options = optimoptions( 'fmincon'  , ...
                            'Display'  ,'off', ...
                            'Algorithm','sqp-legacy' ) ;
    thetai  = fmincon( @(t) -t   , 0                   , ...
                       D.Aineq*ui , D.bineq - D.Aineq*x0, ...
                       []        , []                  , ...
                       []        , []                  , ...
                       @fcst     , options         ) ;   % 迭代初值为 0
   
    omegaBar_M = x0 + thetai*ui ;
    
    plot( omegaBar_M(1), omegaBar_M(2), 'r*', 'LineWidth', 2 )
    
    if thetai >= 1
        pk          = feval( subgrad, omegaBar_M ) ;
        gomegaBar_M = feval( cst    , omegaBar_M ) ;
        Pk.Aineq    = [ node.Pk.Aineq ; ...
                        pk'           ; ] ; 
        Pk.bineq    = [ node.Pk.bineq                ; ...
                        pk'*omegaBar_M - gomegaBar_M ; ] ;
    else
        flag = ( D.Aineq*omegaBar_M - D.bineq >= 0 ) ;
        Pk.Aineq    = [ node.Pk.Aineq ; ...
                        D.Aineq( flag, : ) ; ] ; 
        Pk.bineq    = [ node.Pk.bineq                ; ...
                        D.bineq( flag ) ; ] ;
    end
    
    node.Pk = Pk ;          % 割掉 omegaBar_M
    
    for idx = 1: length( sub_M )
        
        len = length( candidate_list ) ;
        
        child_node   = node ;
        child_node.M = sub_M{ idx } ;
        
        P   = eval( polyh( sub_M{idx}, 'v' ) ) ;
        opt.color = [ 0, 0, idx+1 ]*0.1 ;
        plot( P, opt ) ;
        axis equal ;
        grid on
        hold on
        
        candidate_list( len + 1 ) = { child_node } ;
    end
    
    LB = [ LB ; lower_bound ] ;
    UB = [ UB ; upper_bound ] ;
    
    
end

LB = [ LB ; lower_bound ] ;
UB = [ UB ; upper_bound ] ;

if abs( lower_bound - upper_bound ) <= 1e-7
    fprintf( 'Optimal solution found!\n' ) ;
    output.message = 'Optimal solution found!' ;
end

x        = xBar ;
fval     = upper_bound ;
output   = node.output ;
output.k = k ;

figure
plot( LB, '-^' ), hold on ;
plot( UB, '-s' ) ;

% =======================================================
%   max     t
%   s.t.    cst( x + t*d ) >= 0
%  单值函数求最优目标函数
% (1) 单峰函数
% (2) 一维精确线搜索( 0.618法, 黄金分割法, Fibonacci法杖)
% (3) 函数逼近法( 牛顿迭代法, 割线法, 抛物线法 )
% (4) 托架技术( bracketing-bisection technique )
% ========================================================
function [ c, ceq ] = fcst( t )
    % 复合仿射映射 g(x) = f( A*x + b )是保凸的
    c   = feval( cst, x0 + t*ui ) ;    % f( x + t*d ) >= 0
    ceq = [] ;
    
end

end

% =======================
% 节点选择策略
% =======================
function [ node, candidate_list ] = node_choice_DFS( candidate_list )
    % 选择节点策略, 深度优先搜索, 后进先出
    len                   = length( candidate_list ) ;
    node                  = candidate_list{ len } ;
    candidate_list( len ) = []                    ;
    return ;
end



