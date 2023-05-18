#ifndef _NODE_H_
#define _NODE_H_

#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include "backward.hpp"

#define inf 1>>20
#define useTieBreaker 1
struct GridNode;
typedef GridNode* GridNodePtr;
enum distanceTye {Euclidean, Manhattan, Diagnol};

struct GridNode
{     
    int id;        // 1--> open set, -1 --> closed set 
    Eigen::Vector3i dir;   // direction of expanding
    Eigen::Vector3d coord; // index is defferent from coord in data type
    Eigen::Vector3i index; // index in map, index is Discretization coord
	
    double gScore, fScore;
    GridNodePtr cameFrom; // point to front node
    std::multimap<double, GridNodePtr>::iterator nodeMapIt; 

    GridNode(Eigen::Vector3i _index, Eigen::Vector3d _coord){  
      id = 0;
      index = _index;
      coord = _coord;
      dir   = Eigen::Vector3i::Zero();

      gScore = inf;
      fScore = inf;
      cameFrom = NULL;
    }

    GridNode(){};
    ~GridNode(){};
};


#endif
