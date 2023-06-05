#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;    
using namespace Eigen;

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}
TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

//define factorial function, input i, output i!
int TrajectoryGeneratorWaypoint::Factorial(int x)
{
    int fac = 1;
    for(int i = x; i > 0; i--)
        fac = fac * i;
    return fac;
}

MatrixXd TrajectoryGeneratorWaypoint::getQ(int &num_seg, const int &d_order, int &p_num1d, const VectorXd &time){
    MatrixXd Q = MatrixXd::Zero(p_num1d * num_seg, p_num1d * num_seg);

    for (int i = 0; i < num_seg; i++){
        MatrixXd Q_i = MatrixXd::Zero(p_num1d, p_num1d); // 该段轨迹中每一个的参数
        double T = time(i);
        for (int j = d_order; j < p_num1d; j++){
            for (int k = d_order; k < p_num1d; k++){
                Q_i(j, k) = Factorial(j) / Factorial(j-d_order) * Factorial(k) / Factorial(k-d_order) * pow(T, j+k-7) / (j + k - 7);
            }
        }
        Q.block(i*p_num1d, i*p_num1d, p_num1d, p_num1d);
    }
    return Q;
}

// 将多项式参数映射到实际物理含义的变量中
MatrixXd TrajectoryGeneratorWaypoint::getM(int &num_seg, const int &num_order, int &p_num1d, const VectorXd &time, MatrixXd &coeff){
    MatrixXd M = MatrixXd::Zero(p_num1d*num_seg, p_num1d*num_seg);
    for (int i = 0; i < num_seg; i++){
        double t = time(i);
        MatrixXd M_ = MatrixXd(p_num1d, p_num1d);
        // 初始时刻t=0，只有包含没有变量t的项
        for (int j = 0; j < p_num1d; j++){
            M_(j, j) = coeff(j, j);
        }
        // 对后续每一个点的计算
        for (int j = 0; j < num_order; j++){
            for (int k = 0; k < p_num1d; k++){
                if (j == k) M_(j+num_order, k) = coeff(j, k);
                else M_(j+num_order, k) = coeff(j, k) * pow(t, k-1);
            }
        }
        M.block(i*p_num1d, i*p_num1d, p_num1d, p_num1d) = M_;

    }
    return M;
}

MatrixXd TrajectoryGeneratorWaypoint::getCt(int &num_seg, const int &num_order){
    int ct_rows = num_order*2*num_seg;
    int ct_cols = num_order*2*num_seg - (num_seg-1)*num_order;
    MatrixXd Ct = MatrixXd::Zero(ct_rows, ct_cols);

    // build the hash table
    vector<int> table;
    for (int i = 0; i < num_seg; i++){
        for (int t = 0; t < 2; t++){
            for (int j = 0; j < num_order; j++){
                table.push_back(100*i + t*10 + j);
            }
        }
    }

    // 初始化位置信息，分为开头、末尾、中间部分
    int t = 0;      
    int seg = 0;
    int col = 0, row;
    int val;
    int order;

    // 开始节点
    for (int i = 0; i < num_order; i++){
        val = 100 * seg + t*10 + i;
        auto iter = find(table.begin(), table.end(), val);
        row = distance(table.begin(), iter);
        Ct(row, col) = 1;
        col++;
    }

    // 中间节点
    t = 1;
    order = 0; // 中间位置只有0阶已知
    for (int i = 0; i < num_seg-1; i++){
        // 作为一段的终点
        val = 100*i + 10*t + order;
        auto iter = find(table.begin(), table.end(), val);
        row = distance(table.begin(), iter);
        Ct(row, col) = 1;
        // 作为一段的起点
        val = 100*i + 10*(t-1) + order;
        iter = find(table.begin(), table.end(), val);
        row = distance(table.begin(), iter);
        Ct(row, col) = 1;

        col++;
    }

    // 末尾节点
    seg = num_seg - 1; 
    t = 1;
    for (int i = 0; i < num_order; i++){
        val = 100*seg + 10*1 + i;
        auto iter = find(table.begin(), table.end(), val);
        row = distance(table.begin(), iter);
        Ct(row, col) = 1;
        col++;
    }

    // 连续性约束
    t = 1;
    for (int i = 0; i < num_seg-1; i++){
        // 从一阶速度开始
        for (int j = 1; j < num_order; j++){
            val = 100*i + 10*t + j;
            auto iter = find(table.begin(), table.end(), val);
            row = distance(table.begin(), iter);
            Ct(row, col) = 1;    

            val = 100*(i+1) + 10*(t-1) + j;
            iter = find(table.begin(), table.end(), val);
            row = distance(table.begin(), iter);
            Ct(row, col) = 1;
            col++;
        }
    }
    return Ct;
}

/*
    STEP 2: Learn the "Closed-form solution to minimum snap" in L5, then finish this PolyQPGeneration function
    output      MatrixXd PolyCoeff(num_seg, 3 * p_num1d);   // position(x,y,z), so we need (3 * p_num1d) coefficients
*/

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const int d_order,                    // the order of derivative 控制的阶数(p, v, a, j)四阶
            const Eigen::MatrixXd &Path,          // waypoints coordinates (3d)
            const Eigen::MatrixXd &Vel,           // boundary velocity
            const Eigen::MatrixXd &Acc,           // boundary acceleration
            const Eigen::VectorXd &Time)          // time allocation in each segment
{
    // enforce initial and final velocity and accleration, for higher order derivatives, just assume them be 0;
    int p_order   = 2 * d_order - 1;              // the order of polynomial 多项式的次数
    int p_num1d   = p_order + 1;                  // the number of variables in each segment, 1 is the 0 order

    int num_seg = Time.size();                          // the number of segments, each segments
    MatrixXd PolyCoeff = MatrixXd::Zero(num_seg, 3 * p_num1d);           // position(x,y,z), so we need (3 * p_num1d) coefficients
    VectorXd Px(p_num1d * num_seg), Py(p_num1d * num_seg), Pz(p_num1d * num_seg);

    /*   Produce Mapping Matrix A to the entire trajectory, A is a mapping matrix that maps polynomial coefficients to derivatives.   */
    // get Q
    MatrixXd Q = getQ(num_seg, d_order, p_num1d, Time);

    // get M
    MatrixXd coeff = MatrixXd::Zero(p_num1d * num_seg, p_num1d * num_seg);
    coeff << 1,  1,  1,  1,  1,  1,  1,  1,
             0,  1,  2,  3,  4,  5,  6,  7,
             0,  0,  2,  6,  12, 20, 30, 42,
             0,  0,  0,  6,  24, 60, 120,210;
    MatrixXd M = getM(num_seg, d_order, p_num1d, Time, coeff);
    
    // get Ct
    MatrixXd Ct = getCt(num_seg, d_order);
    
    MatrixXd C = Ct.transpose();
    MatrixXd M_inv = M.inverse();
    MatrixXd M_inv_t = M_inv.transpose();
    MatrixXd R = C * M_inv_t * Q * M_inv * Ct;

    // 计算三个轴的系数
    int num_d_F = 2 * d_order + num_seg - 1;
    int num_d_P = (num_seg - 1) * (d_order - 1);

    MatrixXd R_pp = R.bottomRightCorner(num_d_P, num_d_P);
    MatrixXd R_fp = R.topRightCorner(num_d_F, num_d_P);

    for (int dim = 0; dim < 3; dim++){
        VectorXd waypoints = Path.col(dim);
        VectorXd d_F = VectorXd::Zero(num_d_F);

        d_F(0) = waypoints(0);
        for (int i = 0; i < num_seg - 1; i++){
            d_F(d_order+i) = waypoints(i+1);
        }
        d_F(d_order + num_seg - 1) = waypoints(num_seg); // 尝试是否可以d_F(d_F.rows()-1) = waypoints(waypoints.rows()-1)

        VectorXd d_P = -R_pp.inverse() * R_fp.transpose() * d_F;
        VectorXd d(d_F.rows() + d_P.rows()); 
        d << d_F, d_P;

        VectorXd poly_coef_1d = M.inverse() * Ct * d;
        MatrixXd poly_coef_1d_t = poly_coef_1d.transpose();    

        for(int k = 0; k < num_seg; k++){
            PolyCoeff.block(k, dim*p_num1d, 1, p_num1d) = poly_coef_1d_t.block(0,k*p_num1d, 1, p_num1d);
        }
    }

    return PolyCoeff;
}
