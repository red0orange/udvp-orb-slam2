#ifndef UDVP_H
#define UDVP_H

#include <Eigen/Dense>

#include "Frame.h"
#include "MapPoint.h"

namespace ORB_SLAM2
{

inline double sigmoid(double x)
{
    return 1.0 / (1.0 + exp(-x));
}

double computeUncertainty(Eigen::Matrix4d camera_T, Eigen::Matrix4d estimate_camera_T, Eigen::Vector3d map_point_center, double fx, double fy);  // 对应h函数

double computeObserveError(Eigen::Matrix4d camera_T, Eigen::Vector3d map_point_center, double max_depth, double fov);  // 对应sigma函数

class UDVP
{
public:
    UDVP(const string& strSettingPath);

    cv::Mat optimization(Frame* cur_frmae, std::vector<MapPoint*> local_map_points);

    std::vector<Eigen::Matrix4d> all_optim_camera_T;

private:
    double boundary_function(Eigen::Matrix4d estimate_camera_T, std::vector<MapPoint*> local_map_points, double N_MAX, double max_depth, double fov);  // KKT的约束函数

    float fx;
    float fy;
    float cx;
    float cy;
    cv::Mat K;

    float fov;
    float max_depth;
};
}

#endif