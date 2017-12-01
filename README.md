# Parallel-Programming-Project
Semester project for Applied parallel programming class EE 5351 , fall 2017. 



Introduction
Localization is an important aspect in the Autonomous Vehicles technologies, especially in an Urban
environment. For cars to follow traffic rules and regulations, it should have precise information about its
current lanes and the path to be taken ahead. In order to obtain this information, it needs to localize itself
in a given HD Map. This localization is carried out using Iterative Closest Point algorithm (ICP). It is also
crucial to implement ICP algorithm fast enough to be suitable for real time localization. We intent to
implement a parallel ICP algorithm for localizing a development vehicle in a given HD Map.


Description
Iterative Closest Point (ICP) is an algorithm used to minimize the difference between 2 point clouds. It
gives the rotation and translation required to be applied to one point cloud(source) in order to match it
with the other(reference). The matching parameter used is nothing but a mean square error between the 2
point clouds. So the general process involved in ICP is you apply the transformation R and t to the source
point cloud Ps and calculate the mean square error between the transformed point cloud and the reference
point cloud. The objective is to minimize the mean square error between the 2 point clouds for a certain
Rotation R and translation t.


Objective
We intend to implement this application at a local automotive research company VSI-Labs​ ​, St. Louis
Park facility where we both project members currently are employed as interns. The main objective is to
use HD maps data from the map developer company HERE and carry out Localization of our
development vehicle. In order to obtain precise localization we will be implementing ICP to find the best
fit for our LIDAR data after comparing it with the HD map point cloud data. As the point cloud data is
very large , carrying out ICP on that is very computationally expensive making it difficult to obtain results
in real time. Since for an Autonomous Vehicle, real-time processing is very crucial we intend to use the
GPU to ease this task. For development purposes initially the programming would take place using
NVIDIA GTX 1060 GPU and then ported onto NVIDIA’s Automotive computer DRIVE PX-2 equipped
with 2 Pascal GPU’s if time permits.


Background
We both are currently graduate students in the EE department with background in Robotics, Computer
Vision, Algorithms and Optimization theory. The key skills required to understand and implement ICP is
a basic knowledge of Linear Algebra, Kinematics, understanding of convergence of functions (Calculus)
and Estimation theory.



Resources

● A. Ratter, C. Sammut and M. McGill, "GPU accelerated graph SLAM and occupancy voxel
based ICP for encoder-free mobile robots," 2013 IEEE/RSJ International Conference on
Intelligent Robots and Systems, Tokyo, 2013, pp. 540-547.
● C. Langis, M. Greenspan and G. Godin, "The parallel iterative closest point algorithm,"
Proceedings Third International Conference on 3-D Digital Imaging and Modeling, Quebec City,
Que., 2001, pp. 195-202.
● A. Milstein, M. McGill, T. Wiley, R. Salleh, and C. Sammut, "Occupancy voxel metric based
iterative closest point for position tracking in 3d environments," in ICRA. IEEE, 2011, pp.
4048-4053.
● https://www.coursera.org/learn/robotics-learning/lecture/1jDix/iterative-closest-point
