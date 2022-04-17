#include "UDVPDrawer.h"

#include "Converter.h"

namespace ORB_SLAM2
{

UDVPDrawer::UDVPDrawer(UDVP* mudvp, const string &strSettingPath): mudvp(mudvp)
{
    cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);

    mUDVPSize = fSettings["Viewer.UDVPSize"];
    mUDVPLineWidth = fSettings["Viewer.UDVPLineWidth"];
}

void UDVPDrawer::DrawUDVPOptimPose()
{
    const float &w = mUDVPSize;
    const float h = w*0.75;
    const float z = w*0.6;

    for(size_t i=0; i<mudvp->all_optim_camera_T.size(); i++)
    {
        Eigen::Matrix4d optim_pose = mudvp->all_optim_camera_T[i];
        Eigen::Matrix3d eigen_Rcw = optim_pose(Eigen::seq(0, 2), Eigen::seq(0, 2));
        cv::Mat Rcw = Converter::toCvMat(eigen_Rcw);
        Eigen::Vector3d eigen_tcw = optim_pose(Eigen::seq(0, 2), Eigen::last);
        cv::Mat tcw = Converter::toCvMat(eigen_tcw);
        cv::Mat Rwc = Rcw.t();
        cv::Mat twc = -Rwc*tcw;

        cv::Mat Twc = cv::Mat::eye(4,4,CV_32F);
        Rwc.copyTo(Twc.rowRange(0,3).colRange(0,3));
        twc.copyTo(Twc.rowRange(0,3).col(3));
        Twc = Twc.t();

        glPushMatrix();

        glMultMatrixf(Twc.ptr<GLfloat>(0));

        glLineWidth(mUDVPLineWidth);
        glColor3f(0.0f,0.0f,1.0f);
        glBegin(GL_LINES);
        glVertex3f(0,0,0);
        glVertex3f(w,h,z);
        glVertex3f(0,0,0);
        glVertex3f(w,-h,z);
        glVertex3f(0,0,0);
        glVertex3f(-w,-h,z);
        glVertex3f(0,0,0);
        glVertex3f(-w,h,z);

        glVertex3f(w,h,z);
        glVertex3f(w,-h,z);

        glVertex3f(-w,h,z);
        glVertex3f(-w,-h,z);

        glVertex3f(-w,h,z);
        glVertex3f(w,h,z);

        glVertex3f(-w,-h,z);
        glVertex3f(w,-h,z);
        glEnd();

        glPopMatrix();
    }
}

void UDVPDrawer::DrawUDVPInitPose()
{
    const float &w = mUDVPSize;
    const float h = w*0.75;
    const float z = w*0.6;

    for(size_t i=0; i<mudvp->all_init_camera_T.size(); i++)
    {
        Eigen::Matrix4d init_pose = mudvp->all_init_camera_T[i];
        Eigen::Matrix3d eigen_Rcw = init_pose(Eigen::seq(0, 2), Eigen::seq(0, 2));
        cv::Mat Rcw = Converter::toCvMat(eigen_Rcw);
        Eigen::Vector3d eigen_tcw = init_pose(Eigen::seq(0, 2), Eigen::last);
        cv::Mat tcw = Converter::toCvMat(eigen_tcw);
        cv::Mat Rwc = Rcw.t();
        cv::Mat twc = -Rwc*tcw;

        cv::Mat Twc = cv::Mat::eye(4,4,CV_32F);
        Rwc.copyTo(Twc.rowRange(0,3).colRange(0,3));
        twc.copyTo(Twc.rowRange(0,3).col(3));
        Twc = Twc.t();

        glPushMatrix();

        glMultMatrixf(Twc.ptr<GLfloat>(0));

        glLineWidth(mUDVPLineWidth);
        glColor3f(0.0f,1.0f,1.0f);
        glBegin(GL_LINES);
        glVertex3f(0,0,0);
        glVertex3f(w,h,z);
        glVertex3f(0,0,0);
        glVertex3f(w,-h,z);
        glVertex3f(0,0,0);
        glVertex3f(-w,-h,z);
        glVertex3f(0,0,0);
        glVertex3f(-w,h,z);

        glVertex3f(w,h,z);
        glVertex3f(w,-h,z);

        glVertex3f(-w,h,z);
        glVertex3f(-w,-h,z);

        glVertex3f(-w,h,z);
        glVertex3f(w,h,z);

        glVertex3f(-w,-h,z);
        glVertex3f(w,-h,z);
        glEnd();

        glPopMatrix();
    }
}
}
