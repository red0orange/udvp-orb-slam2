#ifndef UDVPDRAWER_H
#define UDVPDRAWER_H

#include "UDVP.h"

#include<pangolin/pangolin.h>

#include<mutex>

namespace ORB_SLAM2
{

class UDVPDrawer
{
public:
    UDVPDrawer(UDVP* mudvp, const string &strSettingPath);

    void DrawUDVPPose();

private:

    std::mutex mMutexCamera;

    UDVP* mudvp;

    float mUDVPSize;
    float mUDVPLineWidth;
};

} //namespace ORB_SLAM

#endif