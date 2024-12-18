% ============================================
% 全局优化引论
% R. Horst, P.M. Pardalos, N.V. Thoai 著
% 黄红选 译
% 梁治安 校
% P201
% =============================================
%   min f( x1, x2 ) = -( x1 - 20 )^2 - ( x2 - 10 )^2
%   s.t. x2 - x1/2 <= 10
%        ( x2 - 10 )^2 + ( x1 )^2 <= 500
%        x1 >= 0, x2 >= 0
%
function test_gcp_solver

clc ;
clear ;
close all ;

path = './bt-1.3' ;
addpath( path ) ;

[ rcode, res ] = mosekopt( 'symbcon echo(0)' ) ;

% ===========================
% 约束一: 线性约束
% ===========================
Aineq = [ -1/2,  1  ; ...
          -1  ,  0  ; ...    % lb
           0  , -1  ; ...
           1  ,  0  ; ...    % ub
           0  ,  1  ; ] ;
bineq = [ 10 ; ...
           0 ; ...
           0 ; ...
           30 ; ...
           30 ; ] ;
       
rep.B = Aineq ;
rep.b = bineq ;

P   = eval( polyh( rep, 'h' ) ) ; % 求出多胞体 P 的 H-rep, V-rep, P-rep   
CH  = vrep( P ) ;                 % 获取多胞体 P 的 V-rep
Adj = adj( P ) ;                  % 获取顶点对应的链表( 邻接表表示形式 )

plot( P ) ;
axis equal ;
grid on ;
hold on ;

% ==========================
% 约束二: 凸约束
% ==========================
r = sqrt( 500 ) ;       % 半径
a = 0 ;                 % 横坐标
b = 10 ;                % 纵坐标

para = [ a - r, b - r,  2*r, 2*r ] ;
rectangle( 'Position' , para   , ...
           'Curvature', [ 1 1 ] ) ;

x0 = [ 8  ; 4  ; ] ;
V1 = [ 5  ; 5  ; ] ;
V2 = [ 10 ; 2  ; ] ;
V3 = [ 15 ; 3  ; ] ;
V4 = [ 18 ; 15 ; ] ;
V5 = [ 10 ; 10 ; ] ;
S  = [ V1, V2, V3, V4, V5 ] ;

% =============================================
% 输入格式还需要调整一下?????????????????????
% mosek 的 ACC 输入形式
% 引入变量 t
% t = sqrt( 500 )
% =============================================
% 线性锥
% ===================
% x2 - x1/2 <= 10
% x1 >= 0, x2 >= 0
D.Aineq = [ -1/2, 1, 0 ; ] ;
D.bineq = 10 ;
D.lb    = zeros( 3, 1 ) ;
% D.ub    = 30*ones( 3, 1 ) ;
% D.ub    = 30*ones( 3, 1 ) ;

% ====================================
% 二阶锥( quadratic cone, QC )
% ====================================
% ( x2 - 10 )^2 + ( x1 )^2 <= 500
% { x1, x2 - 10, sqrt( 500 ) } <In> lorentz(2)  % cvx form
QC.f = sparse( [ 0, 0, 0 ; ...
                    1, 0, 0 ; ...
                    0, 1, 0 ; ] ) ;
QC.g = [  sqrt( 500 )  ; ...
             0            ; ...
            -10           ; ] ;
% QC.c = [ res.symbcon.MSK_DOMAIN_QUADRATIC_CONE 3 ] ;    % 三维二阶锥 Q^3

QC.c = [ 4 3 ] ;    % 三维二阶锥 Q^3

% ==========================================
% 旋转二阶锥( rotated quadratic cone, RQC )
% ==========================================
RQC.f = [] ;
RQC.g = [] ;
RQC.c = [] ;

% ==================================
% 指数锥( exponential cone, EXPC )
% ==================================
EXPC.f = [] ;
EXPC.g = [] ;
EXPC.c = [] ;

D.f    = [ QC.f ; RQC.f ; EXPC.f ; ] ;
D.g    = [ QC.g ; RQC.g ; EXPC.g ; ] ;
D.accs = [ QC.c ; RQC.c ; EXPC.c ; ] ;

GCP.D   = D  ;          % 凸集 D
GCP.S   = S  ;          % 凸集 D 的内接单纯形 S
GCP.x0  = x0 ;          % 内接单纯形 S 的内点 x0 in int( S )


rep.V = S ;
P   = eval( polyh( rep, 'v' ) ) ; % 求出多胞体 P 的 H-rep, V-rep, P-rep
opt.color = [ 0.1, 1, 0 ] ;
plot( P, opt ) ;
plot( x0(1), x0(2), 'rs', 'LineWidth', 2 ) ;

% ======================
% 求内接交点 sBari
% ======================
% 计算以 x0 为起点, 到通过 si, i = 1, 2, ..., 4 的射线与 D 的交点
t0      = 0 ;       % 初值 t0
x       = x0 ;
d       = V1 - x0 ;
options = optimoptions( 'fmincon'  , ...
                        'Display'  ,'off', ...
                        'Algorithm','sqp-legacy' ) ;
lambda  = fmincon( @(t) -t    , t0             , ...
                   Aineq*d , bineq - Aineq*x, ...
                   []      , []             , ...
                   []      , []             , ...
                   @nonlcon, options         ) ;
sBar1   = x + lambda*d

d       = V2 - x0 ;
lambda  = fmincon( @fun    , t0             , ...
                   Aineq*d , bineq - Aineq*x, ...
                   []      , []             , ...
                   []      , []             , ...
                   @nonlcon, options         ) ;
sBar2   = x + lambda*d

d       = V3 - x0 ;
lambda  = fmincon( @fun    , t0             , ...
                   Aineq*d , bineq - Aineq*x, ...
                   []      , []             , ...
                   []      , []             , ...
                   @nonlcon, options         ) ;
sBar3   = x + lambda*d

d       = V4 - x0 ;
lambda  = fmincon( @fun    , t0             , ...
                   Aineq*d , bineq - Aineq*x, ...
                   []      , []             , ...
                   []      , []             , ...
                   @nonlcon, options         ) ;
sBar4   = x + lambda*d

d       = V5 - x0 ;
lambda  = fmincon( @fun    , t0             , ...
                   Aineq*d , bineq - Aineq*x, ...
                   []      , []             , ...
                   []      , []             , ...
                   @nonlcon, options         ) ;
sBar5   = x + lambda*d

SBar  = [ sBar1, sBar2, sBar3, sBar4, sBar5 ] ;

rep.V = SBar ;
P     = eval( polyh( rep, 'v' ) ) ; % 求出多胞体 P 的 H-rep, V-rep, P-rep
opt.color = [ 0.1, 0.8, 1 ] ;
plot( P, opt ) ;

plot( [ x0(1), sBar1(1) ], [ x0(2), sBar1(2) ] , '-*' ) ;
plot( [ x0(1), sBar2(1) ], [ x0(2), sBar2(2) ] , '-*' ) ;
plot( [ x0(1), sBar3(1) ], [ x0(2), sBar3(2) ] , '-*' ) ;
plot( [ x0(1), sBar4(1) ], [ x0(2), sBar4(2) ] , '-*' ) ;
plot( [ x0(1), sBar5(1) ], [ x0(2), sBar5(2) ] , '-*' ) ;

% (4) 令 xBar1 in argmin{ f(sBar1), ..., f(sBarn+1) }
%        gamma1 = f( xBar1 )
%     为当前最好可行解
fsBar = zeros( 5, 1 ) ;
for idx = 1: 5
    fsBar( idx, 1 ) = feval( @oracle, SBar( :, idx ) ) ;
end
[ fxBar1, idxopt ] = min( fsBar ) ;
xBar1  = SBar( :, idxopt )
gamma1 = fxBar1

% =========================
% (5) 计算以点 x0为起点
%     过 sBarj, 
%     极方向为 sBarj - x0, 
%     且 j != i 
%     的 n 条边构成的锥 Mi
% ==========================
% lambda  = [ 1 ; 1 ; 1 ; 1 ; 1 ; ] ;
% U       = [ SBar - x0 ]
% n       = size( U, 2 ) ;
% flag    = lambda > 1e-8 ;
% M       = cell( sum( flag ), 1 ) ;          % 锥集 M
% counter = 0 ;                               % 累加器
% Mi.V    = x0 ;                              % 锥 Mi 顶点
% for idx = 1: n
%     if flag( idx )
%         Mi.D = U ;                          % 锥 Mi 的 n 条边
%         Mi.D( :, idx ) = [] ;               % w = [] ;            
%         M( counter + 1, 1 ) = { Mi } ;
%         counter = counter + 1 ;
%     else
%         continue ;
%     end
% end

% ============================
% 设置初始包覆 D 的多面体 P1
% ============================
Aineq = [ -1  ,  0  ; ...    % lb
           0  , -1  ; ...
           1  ,  0  ; ...    % ub
           0  ,  1  ; ] ;
bineq = [  5  ; ...
           5  ; ...
           25 ; ...
           25 ; ] ;

rep.B = Aineq ;
rep.b = bineq ;
P   = eval( polyh( rep, 'h' ) ) ; % 求出多胞体 P 的 H-rep, V-rep, P-rep 
opt.color = [ 0.1, 0.8, 0.5 ] ;
plot( P, opt ) ;

% ============================
% 对每个锥 M
% 基 [ u1, u2, ..., un ]
% 计算点 x0 + thetai*ui
% ============================
U  = [ SBar - x0 ]
M1.V = x0
M1.D = U( :, 1: 2 )

[ y1, theta1 ] = gamma_extension( gamma1, M1.V, M1.D( :, 1 ), @oracle )
[ y2, theta2 ] = gamma_extension( gamma1, M1.V, M1.D( :, 2 ), @oracle )

Theta  = [ theta1 ; theta2 ; ] ;
n      = 2 ;
ATilde = [ Aineq ; -eye( n ) ] ;            % LP( P1, M ) 的系统矩阵
bTilde = [ bineq - Aineq*x0 ; x0 ; ] ;      % LP( P1, M ) 的右向量
U      = M1.D ;          % 锥 M 的极方向矩阵 U = [ u1, ..., un ]
[ tBar_M , ...           % tBar_M = [ t1 ; ... ; tn ] ;
  mu_M   , ...           % mu( M )
  output ] = lpDM_solver( Theta , ...
                               ATilde, ...
                               bTilde, ...
                               U      )

omega_M  = x0 + U*tBar_M
fomega_M = oracle( omega_M )
plot( omega_M(1), omega_M(2), 'r*', 'LineWidth', 2 )

% ============================
U  = [ SBar - x0 ]
M1.V = x0
M1.D = U( :, 2: 3 )

[ y1, theta1 ] = gamma_extension( gamma1, M1.V, M1.D( :, 1 ), @oracle )
[ y2, theta2 ] = gamma_extension( gamma1, M1.V, M1.D( :, 2 ), @oracle )

Theta  = [ theta1 ; theta2 ; ] ;
n      = 2 ;
ATilde = [ Aineq ; -eye( n ) ] ;            % LP( P1, M ) 的系统矩阵
bTilde = [ bineq - Aineq*x0 ; x0 ; ] ;      % LP( P1, M ) 的右向量
U      = M1.D ;          % 锥 M 的极方向矩阵 U = [ u1, ..., un ]
[ tBar_M , ...           % tBar_M = [ t1 ; ... ; tn ] ;
  mu_M   , ...           % mu( M )
  output ] = lpDM_solver( Theta , ...
                               ATilde, ...
                               bTilde, ...
                               U      )

omega_M = x0 + U*tBar_M
fomega_M = oracle( omega_M )
plot( omega_M(1), omega_M(2), 'r*', 'LineWidth', 2 )

% ============================
U  = [ SBar - x0 ]
M1.V = x0
M1.D = U( :, 3: 4 )

[ y1, theta1 ] = gamma_extension( gamma1, M1.V, M1.D( :, 1 ), @oracle )
[ y2, theta2 ] = gamma_extension( gamma1, M1.V, M1.D( :, 2 ), @oracle )

Theta  = [ theta1 ; theta2 ; ] ;
n      = 2 ;
ATilde = [ Aineq ; -eye( n ) ] ;            % LP( P1, M ) 的系统矩阵
bTilde = [ bineq - Aineq*x0 ; x0 ; ] ;      % LP( P1, M ) 的右向量
U      = M1.D ;          % 锥 M 的极方向矩阵 U = [ u1, ..., un ]
[ tBar_M , ...           % tBar_M = [ t1 ; ... ; tn ] ;
  mu_M   , ...           % mu( M )
  output ] = lpDM_solver( Theta , ...
                               ATilde, ...
                               bTilde, ...
                               U      )

omega_M = x0 + U*tBar_M
fomega_M = oracle( omega_M )
plot( omega_M(1), omega_M(2), 'r*', 'LineWidth', 2 )

% ============================
U  = [ SBar - x0 ]
M1.V = x0
M1.D = U( :, 4: 5 )

[ y1, theta1 ] = gamma_extension( gamma1, M1.V, M1.D( :, 1 ), @oracle )
[ y2, theta2 ] = gamma_extension( gamma1, M1.V, M1.D( :, 2 ), @oracle )

Theta  = [ theta1 ; theta2 ; ] ;
n      = 2 ;
ATilde = [ Aineq ; -eye( n ) ] ;            % LP( P1, M ) 的系统矩阵
bTilde = [ bineq - Aineq*x0 ; x0 ; ] ;      % LP( P1, M ) 的右向量
U      = M1.D ;          % 锥 M 的极方向矩阵 U = [ u1, ..., un ]
[ tBar_M , ...           % tBar_M = [ t1 ; ... ; tn ] ;
  mu_M   , ...           % mu( M )
  output ] = lpDM_solver( Theta , ...
                               ATilde, ...
                               bTilde, ...
                               U      )

omega_M = x0 + U*tBar_M
fomega_M = oracle( omega_M )
plot( omega_M(1), omega_M(2), 'r*', 'LineWidth', 2 )

% ============================
U  = [ SBar - x0 ]
M1.V = x0
M1.D = [ U( :, 5 ), U( :, 1 ) ]

[ y1, theta1 ] = gamma_extension( gamma1, M1.V, M1.D( :, 1 ), @oracle )
[ y2, theta2 ] = gamma_extension( gamma1, M1.V, M1.D( :, 2 ), @oracle )

Theta  = [ theta1 ; theta2 ; ] ;
n      = 2 ;
ATilde = [ Aineq ; -eye( n ) ] ;            % LP( P1, M ) 的系统矩阵
bTilde = [ bineq - Aineq*x0 ; x0 ; ] ;      % LP( P1, M ) 的右向量
U      = M1.D ;          % 锥 M 的极方向矩阵 U = [ u1, ..., un ]
[ tBar_M , ...           % tBar_M = [ t1 ; ... ; tn ] ;
  mu_M   , ...           % mu( M )
  output ] = lpDM_solver( Theta , ...
                               ATilde, ...
                               bTilde, ...
                               U      )

omega_M = x0 + U*tBar_M
fomega_M = oracle( omega_M )
plot( omega_M(1), omega_M(2), 'r*', 'LineWidth', 2 )










% =======================================================
%   max     t
%   s.t.    fcst( x + t*d ) >= 0
%  单值函数求最优目标函数
% (1) 单峰函数
% (2) 一维精确线搜索( 0.618法, 黄金分割法, Fibonacci法杖)
% (3) 函数逼近法( 牛顿迭代法, 割线法, 抛物线法 )
% (4) 托架技术( bracketing-bisection technique )
% ========================================================
function f = fun( t )
    f = -t ;        % min -t
end
function [ c, ceq ] = nonlcon( t )
    % 复合仿射映射 g(x) = f( A*x + b )是保凸的
    c   = feval( @fcst, x + t*d ) ;    % f( x + t*d ) >= 0
    ceq = [] ;
    % 约束函数
    function f = fcst( x )
        f = x(1)^2 + ( x(2) - 10 )^2 - 500 ;
    end
end

% ==============================
% 目标函数
% ==============================
function f = oracle( x )
    % 目标函数
    % f( x1, x2 ) = -( x1 - 20 )^2 - ( x2 - 10 )^2
    f = -( x(1) - 20 )^2 - ( x(2) - 10 )^2 ;
end












end



