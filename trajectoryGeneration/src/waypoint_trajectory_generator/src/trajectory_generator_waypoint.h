#ifndef _TRAJECTORY_GENERATOR_WAYPOINT_H_
#define _TRAJECTORY_GENERATOR_WAYPOINT_H_

#include <Eigen/Eigen>
#include <vector>

#include "OsqpEigen/OsqpEigen.h"

class TrajectoryGeneratorWaypoint {
    private:
		double _qp_cost;
		Eigen::MatrixXd _Q;
		Eigen::VectorXd _Px, _Py, _Pz;
    public:
        TrajectoryGeneratorWaypoint();

        ~TrajectoryGeneratorWaypoint();

        // closed-form: 解析解
        Eigen::MatrixXd PolyQPGeneration(
            const int order,
            const Eigen::MatrixXd &Path,
            const Eigen::MatrixXd &Vel,
            const Eigen::MatrixXd &Acc,
            const Eigen::VectorXd &Time);
        // 数值解
        Eigen::MatrixXd PolyQPGeneration(
            const int order,
            const Eigen::MatrixXd &Path,
            const Eigen::MatrixXd &Vel,
            const Eigen::MatrixXd &Acc,
            const Eigen::VectorXd &Time,
            const OsqpEigen::Solver &solve);
        
        Eigen::MatrixXd getQ(
            int &num_seg, 
            const int &num_order, 
            int &p_num1d, 
            const Eigen::VectorXd &time);
        Eigen::MatrixXd getM(
            int &num_seg, 
            const int &num_order,             
            int &p_num1d, 
            const Eigen::VectorXd &time,
            Eigen::MatrixXd &coeff);
        Eigen::MatrixXd getCt(int &num_seg, const int &num_order);

        int Factorial(int x);
};
        

#endif
