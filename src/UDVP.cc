#include "UDVP.h"
#include "math.h"

#include "Thirdparty/g2o/g2o/core/block_solver.h"
#include "Thirdparty/g2o/g2o/core/optimization_algorithm_levenberg.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_eigen.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"
#include "Thirdparty/g2o/g2o/core/robust_kernel_impl.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_dense.h"
#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"
#include "Thirdparty/g2o/g2o/core/base_multi_edge.h"

#include "Converter.h"


namespace ORB_SLAM2
{

// class UDVPVertex: public g2o::BaseVertex<3, Eigen::Vector3d>
// {
//     public:
//     EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // 字节对齐

//     virtual void setToOriginImpl() // 重置，设定被优化变量的原始值
//     {
//         _estimate << 0,0,0;
//     }

//     virtual void oplusImpl( const double* update ) // 更新
//     {
//         _estimate += Eigen::Vector3d(update);   //update强制类型转换为Vector3d
//     }
//     // 存盘和读盘：留空
//     virtual bool read( istream& in ) {}
//     virtual bool write( ostream& out ) const {}
// };


double computeUncertainty(Eigen::Matrix4d camera_T, Eigen::Matrix4d estimate_camera_T, Eigen::Vector3d map_point_center, double fx, double fy) 
{
    // world axis to camera axis transform
    Eigen::Matrix<double,3,3> rotate_mat = camera_T(Eigen::seq(0, 2), Eigen::seq(0, 2));
    Eigen::Vector3d translation_vector = camera_T(Eigen::seq(0, 2), Eigen::last);
    Eigen::Matrix<double,3,3> estimate_rotate_mat = estimate_camera_T(Eigen::seq(0, 2), Eigen::seq(0, 2));
    Eigen::Vector3d estimate_translation_vector = estimate_camera_T(Eigen::seq(0, 2), Eigen::last);
    // cal world axis camera center
    Eigen::Vector3d camera_center = -rotate_mat.transpose() * translation_vector;
    Eigen::Vector3d estimate_camera_center = -estimate_rotate_mat.transpose() * estimate_translation_vector;

    // **************** 计算delta_depth *********************
    // cosine theorem cal angle
    double a = (map_point_center - camera_center).norm();
    double b = (map_point_center - estimate_camera_center).norm();
    double c = (camera_center - estimate_camera_center).norm();

    double cos_alpha = (b*b + c*c - a*a) / (2*b*c);
    double cos_beta = (a*a + c*c - b*b) / (2*a*c);

    // 注意，函数结果为弧度制
    double delta_beta = atan2(2, (fx + fy));
    double alpha = acos(cos_alpha);
    double beta = acos(cos_beta);

    // 计算最终的公式
    // sin、cos函数也是弧度制，不用转换
    double delta_depth = c * ( sin(beta + delta_beta) / sin(CV_PI - alpha - beta - delta_beta) ) - b;

    // **************** 计算delta_rotate *********************
    // 3D in absolute coordinates
    Eigen::Vector3d Pw = map_point_center;
    Eigen::Vector4d homo_Pw(Pw(0), Pw(1), Pw(2), 1);

    Eigen::Vector4d homo_Pc = estimate_camera_T * homo_Pw;
    Eigen::Vector3d Pc(homo_Pc(0), homo_Pc(1), homo_Pc(2));

    const double &PcX = Pc(0);
    const double &PcY= Pc(1);
    const double &PcZ = Pc(2);

    Eigen::Matrix<double,2,3> tmp_mat_1;
    tmp_mat_1(0,0) = fx / PcZ;
    tmp_mat_1(0,1) = 0;
    tmp_mat_1(0,2) = -( fx * PcX / PcZ * PcZ );
    tmp_mat_1(1,0) = 0;
    tmp_mat_1(1,1) = fy / PcZ;
    tmp_mat_1(1,2) = -( fy * PcY / PcZ * PcZ );

    // 计算第二部分的反对称矩阵
    Eigen::Vector3d tmp_p = rotate_mat * Pw;
    Eigen::Matrix<double, 3, 3> tmp_mat_2;
    tmp_mat_2(0,0) = 0;
    tmp_mat_2(0,1) = -tmp_p(2);
    tmp_mat_2(0,2) = tmp_p(1);
    tmp_mat_2(1,0) = tmp_p(2);
    tmp_mat_2(1,1) = 0;
    tmp_mat_2(1,2) = -tmp_p(0);
    tmp_mat_2(2,0) = -tmp_p(1);
    tmp_mat_2(2,1) = tmp_p(0);
    tmp_mat_2(2,2) = 0;

    // 生成delta_m随机方向向量
    const int N = 999;
    double u = rand() % (N + 1) / (double)(N + 1);
    double v = rand() % (N + 1) / (double)(N + 1);
    double theta = 2 * CV_PI * u;
    double phi = acos(2 * v - 1);

    Eigen::Vector3d unit_rand_vector;
    unit_rand_vector(0) = sin(theta) * sin(phi);
    unit_rand_vector(1) = cos(theta) * sin(phi);
    unit_rand_vector(2) = cos(phi);

    Eigen::Vector2d delta_rotate = -tmp_mat_1 * tmp_mat_2 * unit_rand_vector;

    // **************** 计算uncertainly *********************
    double delta_rotate_euclidean = (delta_rotate.transpose() * delta_rotate)(0);
    double uncertainty = delta_rotate_euclidean * delta_depth * delta_depth;
    return uncertainty;
}

double computeObserveError(Eigen::Matrix4d camera_T, Eigen::Vector3d map_point_center, double max_depth, double fov)
{
    Eigen::Matrix<double,3,3> rotate_mat = camera_T(Eigen::seq(0, 2), Eigen::seq(0, 2));
    Eigen::Vector3d translation_vector = camera_T(Eigen::seq(0, 2), Eigen::last);
    Eigen::Vector3d camera_center = -rotate_mat.transpose() * translation_vector;
    // **************** 计算观测模型 *********************
    // 计算观测部分的误差
    Eigen::Vector3d unit_vector(0, 0, 1);
    Eigen::Vector3d mp_center_c(camera_center(0), camera_center(1), camera_center(2));

    double lambda = acos( (unit_vector.transpose() * mp_center_c)(0) / (unit_vector.norm() * mp_center_c.norm()) );
    double observe_error = sigmoid(max_depth - mp_center_c.norm()) * sigmoid((fov / 2) - abs(lambda));
    return observe_error;
}

// （误差）边的模型    模板参数：观测值维度，类型，连接顶点类型
class UDVPEdge: public g2o::BaseUnaryEdge<1, double, g2o::VertexSE3Expmap>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    UDVPEdge( double fx, double fy, double cx, double cy, double fov, double max_depth,
            Frame* cur_frame,
            MapPoint* map_point,
            double phi, double M, double N, double N_MAX 
            ): BaseUnaryEdge(), 
            fx(fx), fy(fy), cx(cx), cy(cy), fov(fov), max_depth(max_depth),
            cur_frame(cur_frame),
            map_point(map_point),
            phi(phi), M(M), N(N), N_MAX(N_MAX)
    {
    }

    // 计算曲线模型误差
    void computeError()
    {
        // 定义节点的类型，取出n元边的n元
        const g2o::VertexSE3Expmap* pose_node = static_cast<const g2o::VertexSE3Expmap*>(_vertices[0]);

        // 取出估计，即取出变量的含义
        const g2o::SE3Quat pose_se3 = pose_node->estimate();
        // 将估计的pose转换为homo matrix
        Eigen::Matrix<double,4,4> estimate_homo_pose = pose_se3.to_homogeneous_matrix();

        // 取出当前帧的homo matrix
        Eigen::Matrix<double, 4, 4> homo_pose;
        homo_pose(Eigen::seq(0, 2), Eigen::seq(0, 2)) = Converter::toMatrix3d(cur_frame->GetRotation());
        homo_pose(Eigen::seq(0, 2), Eigen::last) = Converter::toVector3d(cur_frame->GetTranslation());
        homo_pose(3, 3) = 1;

        // 取出map_point的world axis center
        Eigen::Vector3d mp_world_center = Converter::toVector3d(map_point->GetWorldPos());

        // 计算uncertainty和observe error
        double uncertainty = computeUncertainty(homo_pose, estimate_homo_pose, mp_world_center, fx, fy);
        double observe_error = computeObserveError(estimate_homo_pose, mp_world_center, max_depth, fov);

        _error(0) = observe_error * (uncertainty + phi - M) - (phi * N_MAX / N) + M; 
    }
    virtual bool read( istream& in ) {}
    virtual bool write( ostream& out ) const {}
public:
    double fx;
    double fy;
    double cx;
    double cy;
    double fov;
    double max_depth;

    double M;  // exploration ratio
    double N;
    double N_MAX; // over-exploration limit num
    double phi;

    Frame* cur_frame;
    MapPoint *map_point;
};

UDVP::UDVP(const string& strSettingsFile)
{
    cv::FileStorage fSettings(strSettingsFile, cv::FileStorage::READ);
    float fx = fSettings["Camera.fx"];
    float fy = fSettings["Camera.fy"];
    float cx = fSettings["Camera.cx"];
    float cy = fSettings["Camera.cy"];
    float fov = fSettings["Camera.fov"];
    int max_depth = fSettings["Camera.max_depth"];

    this->fx = fx;
    this->fy = fy;
    this->cx = cx;
    this->cy = cy;
    this->K = cv::Mat::eye(3, 3, CV_32F);
    K.at<float>(0,0) = fx;
    K.at<float>(1,1) = fy;
    K.at<float>(0,2) = cx;
    K.at<float>(1,2) = cy;
    this->fov = fov;
    this->max_depth = max_depth;
}

double UDVP::boundary_function(Eigen::Matrix4d estimate_camera_T, std::vector<MapPoint*> local_map_points, double N_MAX, double max_depth, double fov)
{
    double tmp_value = 0;
    for (int i = 0; i < local_map_points.size(); i++)
    {
        MapPoint *cur_map_point = local_map_points[i];
        Eigen::Vector3d cur_map_point_center = Converter::toVector3d(cur_map_point->GetWorldPos());
        tmp_value += computeObserveError(estimate_camera_T, cur_map_point_center, max_depth, fov);
    }
    return N_MAX - tmp_value;
}

cv::Mat UDVP::optimization(Frame* cur_frame, std::vector<MapPoint*> local_map_points){
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;
    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolver_6_3::PoseMatrixType>();
    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    // 配置参数
    // 超参数
    double gamma = 0.7;
    double N_t = 50;
    double M = 50;
    double tolerance = 1e-4;  // KKT第二个约束最优化问题的容忍值

    // 动态超参数
    double N = 0;
    Eigen::Vector3d n_w(0, 0, 0);

    Eigen::Vector3d unit_vector(0, 0, 1);
    Eigen::Vector3d cur_camera_center = Converter::toVector3d(cur_frame->GetCameraCenter());
    Eigen::Matrix3d cur_camera_rotation = Converter::toMatrix3d(cur_frame->GetRotation());
    Eigen::Vector3d cur_camera_translation = Converter::toVector3d(cur_frame->GetTranslation());
    Eigen::Matrix4d cur_camera_homo_pose;  // 保留一份homo pose方便计算
    cur_camera_homo_pose(Eigen::seq(0, 2), Eigen::seq(0, 2)) = cur_camera_rotation;
    cur_camera_homo_pose(Eigen::seq(0, 2), Eigen::last) = cur_camera_translation;
    cur_camera_homo_pose(3, 3) = 1;

    for (int i = 0; i < local_map_points.size(); i++)
    {
        // 计算mp的摄像头坐标
        Eigen::Vector3d Pw = Converter::toVector3d(local_map_points[i]->GetWorldPos());
        Eigen::Vector4d homo_Pw(Pw(0), Pw(1), Pw(2), 1);
        Eigen::Vector4d homo_Pc = cur_camera_homo_pose * homo_Pw;
        Eigen::Vector3d Pc(homo_Pc(0), homo_Pc(1), homo_Pc(2));
        // ********************** 计算当前帧姿态N的值 **************************
        // 计算观测error
        double lambda = acos( (unit_vector.transpose() * Pc)(0) / (unit_vector.norm() * Pc.norm()) );
        double observe_error = sigmoid(this->max_depth - Pc.norm()) * sigmoid((fov / 2) - abs(lambda));
        N += observe_error;

        // ********************** 计算优化使用的初始化位姿态 **************************
        n_w += ((Pw - cur_camera_center) / (Pw - cur_camera_center).norm());
    }
    // ********************** 根据N、N_t、gamma得到N_max值 **************************
    double N_MAX = max(gamma * N, N_t);
    // ********************** 根据n_w计算得到初始旋转矩阵 **************************
    Eigen::Matrix4d init_pose; 
    init_pose(Eigen::seq(0, 2), Eigen::last) = cur_camera_translation;  // 平移矩阵固定
    Eigen::AngleAxisd init_rotate_vector(
        acos((n_w.transpose() * -unit_vector)(0) / (n_w.norm() * -unit_vector.norm())),
        ((-n_w.transpose().cross(unit_vector)) / (-n_w.transpose().cross(unit_vector).norm()))
        ); // TODO 不明白它这里求反方向+旋转向量的方法，不确定对不对
    Eigen::Matrix3d init_rotate_matrix = init_rotate_vector.toRotationMatrix();
    init_pose(Eigen::seq(0, 2), Eigen::seq(0, 2)) = init_rotate_matrix;

    g2o::VertexSE3Expmap * pose = new g2o::VertexSE3Expmap();
    pose->setEstimate(Converter::toSE3Quat(init_pose));
    pose->setId(0);
    pose->setFixed(false);
    optimizer.addVertex(pose);

    // 构造优化图
    vector<UDVPEdge*> vpEdgesMono;
    vector<size_t> vnIndexEdgeMono;
    double tmp_phi = 0;
    for (int i = 0; i < local_map_points.size(); i++)
    {
        UDVPEdge* e = new UDVPEdge(fx, fy, cx, cy, fov, max_depth, cur_frame, local_map_points[i], tmp_phi, M, N, N_MAX);

        e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));

        Eigen::Matrix<double, 6, 6> pure_rotate_information_matrix;
        pure_rotate_information_matrix.setIdentity();
        pure_rotate_information_matrix(3, 3) = 1000;
        pure_rotate_information_matrix(4, 4) = 1000;
        pure_rotate_information_matrix(5, 5) = 1000;
        // e->setInformation(pure_rotate_information_matrix);

        // e->setMeasurement(obs);
        // const float invSigma2 = pFrame->mvInvLevelSigma2[kpUn.octave];
        // e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

        // g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
        // e->setRobustKernel(rk);
        // rk->setDelta(deltaMono);

        optimizer.addEdge(e);

        vpEdgesMono.push_back(e);
        vnIndexEdgeMono.push_back(i);
    }

    // ********************* 开始优化，添加图的边，KKT 条件使用 ********************
    Eigen::Matrix4d optim_pose = init_pose;
    // 无约束情况
    if (N <= N_t) {
        double phi = 0;
        for (auto edge : vpEdgesMono)
        {
            edge->phi = phi;
        }
        optimizer.initializeOptimization(0);
        optimizer.optimize(10); // 默认迭代次数为10

        g2o::VertexSE3Expmap* optim_pose_SE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
        g2o::SE3Quat optim_pose_SE3quatv = optim_pose_SE3->estimate();
        optim_pose = optim_pose_SE3quatv.to_homogeneous_matrix();
    }
    // 有约束情况
    else
    {
        int k = 0;
        double epsilon = N / N_MAX;

        double tmp_value = 0;
        for (int i = 0; i < local_map_points.size(); i++) {
            MapPoint *cur_map_point = local_map_points[i];
            Eigen::Vector3d cur_map_point_center = Converter::toVector3d(cur_map_point->GetWorldPos());
            tmp_value += (computeObserveError(init_pose, cur_map_point_center, this->max_depth, this->fov) / N_MAX) + (computeUncertainty(cur_camera_homo_pose, init_pose, cur_map_point_center, this->fx, this->fy) / N);
        }
        double phi = M - tmp_value;

        while (boundary_function(optim_pose, local_map_points, N_MAX, max_depth, fov) > tolerance)
        {
            for (auto edge : vpEdgesMono)
            {
                edge->phi = phi;
            }

            optimizer.initializeOptimization(0);
            optimizer.optimize(10); // 默认迭代次数为10

            k += 1;
            phi += k * epsilon;

            g2o::VertexSE3Expmap* optim_pose_SE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(0));
            g2o::SE3Quat optim_pose_SE3quatv = optim_pose_SE3->estimate();
            optim_pose = optim_pose_SE3quatv.to_homogeneous_matrix();
        }
    }

    return Converter::toCvMat(optim_pose);
    // return cv::Mat::zeros(cv::Size(3, 3), CV_32F);
}
}
