# Robot_Planning
This is an introductory project about motion planning for mobile robots. The main content includes graph search-based algorithms (A*, jps), sampling-based algorithms (RRT), path planning and trajectory generation considering object dynamic constraints.
### Prerequisities
- Ubuntu 18.04 with ros-melodic
- eigen3
- pcl-1.8
- ompl-1.6.0
- osqp and osqp-eigen
- cmake version > 3.20
### How to use?
    git clone https://github.com/LwXxz/Robot_Planning.git
    cd [the project]
    catkin_make
    source ./devel/setup.bash
    roslaunch [node name] [launch file]
> reference: https://github.com/teamo1996/Motion-plan
