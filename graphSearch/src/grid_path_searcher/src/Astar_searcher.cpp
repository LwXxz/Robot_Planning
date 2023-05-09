#include "Astar_searcher.h"

using namespace std;
using namespace Eigen;

void AstarPathFinder::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    
    GLX_SIZE = max_x_id;
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
    
    GridNodeMap = new GridNodePtr ** [GLX_SIZE];
    for(int i = 0; i < GLX_SIZE; i++){
        GridNodeMap[i] = new GridNodePtr * [GLY_SIZE];
        for(int j = 0; j < GLY_SIZE; j++){
            GridNodeMap[i][j] = new GridNodePtr [GLZ_SIZE];
            for( int k = 0; k < GLZ_SIZE;k++){
                Vector3i tmpIdx(i,j,k);
                Vector3d pos = gridIndex2coord(tmpIdx);
                GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);
            }
        }
    }
}

void AstarPathFinder::resetGrid(GridNodePtr ptr)
{
    ptr->id = 0;
    ptr->cameFrom = NULL;
    ptr->gScore = inf;
    ptr->fScore = inf;
}

void AstarPathFinder::resetUsedGrids()
{   
    for(int i=0; i < GLX_SIZE ; i++)
        for(int j=0; j < GLY_SIZE ; j++)
            for(int k=0; k < GLZ_SIZE ; k++)
                resetGrid(GridNodeMap[i][j][k]);
}

void AstarPathFinder::setObs(const double coord_x, const double coord_y, const double coord_z)
{   
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      

    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

vector<Vector3d> AstarPathFinder::getVisitedNodes()
{   
    vector<Vector3d> visited_nodes;
    for(int i = 0; i < GLX_SIZE; i++)
        for(int j = 0; j < GLY_SIZE; j++)
            for(int k = 0; k < GLZ_SIZE; k++){   
                //if(GridNodeMap[i][j][k]->id != 0) // visualize all nodes in open and close list
                if(GridNodeMap[i][j][k]->id == -1)  // visualize nodes in close list only
                    visited_nodes.push_back(GridNodeMap[i][j][k]->coord);
            }

    ROS_WARN("visited_nodes size : %d", visited_nodes.size());
    return visited_nodes;
}

// grid to world coordinate system
Vector3d AstarPathFinder::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

// world coordinate system to grid
Vector3i AstarPathFinder::coord2gridIndex(const Vector3d & pt) 
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d AstarPathFinder::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

inline bool AstarPathFinder::isOccupied(const Eigen::Vector3i & index) const
{
    return isOccupied(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isFree(const Eigen::Vector3i & index) const
{
    return isFree(index(0), index(1), index(2));
}

inline bool AstarPathFinder::isOccupied(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return  (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
            (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] == 1));
}

inline bool AstarPathFinder::isFree(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

inline void AstarPathFinder::AstarGetSucc(GridNodePtr currentPtr, vector<GridNodePtr> & neighborPtrSets, vector<double> & edgeCostSets)
{   
    neighborPtrSets.clear();
    edgeCostSets.clear();

    // get the position
    int current_x = currentPtr->index[0];
    int current_y = currentPtr->index[1];
    int current_z = currentPtr->index[2];
    Vector3d currentCoord = currentPtr->coord;
    // from position to the neibor
    for (int x = -1; x <= 1; x++){
        for (int y = -1; y <= 1; y++){
            for (int z = -1; z <= 1; z++){
                if (x == 0 && y == 0 && z == 0){
                    continue;
                }
                // the expand index
                int expand_x = current_x + x;
                int expand_y = current_y + y;
                int expand_z = current_z + z;
                // use index, should compare with the GLX_SIZE...
                if (expand_x < 0 || (expand_x > (GLX_SIZE - 1)) || expand_y < 0 || (expand_y > (GLY_SIZE - 1)) || expand_z < 0 || (expand_z > (GLZ_SIZE - 1))) {
                    continue;
                }
                // whether occupied
                if (isOccupied(expand_x, expand_y, expand_z)){
                    continue;
                }
                // whether in close set
                if (GridNodeMap[expand_x][expand_y][expand_z]->id == -1){
                    continue;
                }
                GridNodePtr tempNode = GridNodeMap[expand_x][expand_y][expand_z];

                neighborPtrSets.push_back(tempNode);
                // edgeCost uesd the coord, Vector3d can be indexed by [] or (), MatrixXd just (). 
                double edgeCost = std::sqrt(pow(tempNode->coord(0) - currentCoord(0), 2) + pow(tempNode->coord(1) - currentCoord(1), 2) + pow(tempNode->coord(2) - currentCoord(2), 2));
                edgeCostSets.push_back(edgeCost);
                
            }
        }
    }
    
}

double AstarPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2)
{
    /* 
    base on world Coordinate System!!!!
    choose possible heuristic function you want
    Manhattan, Euclidean, Diagonal, or 0 (Dijkstra)
    */
    distanceTye D = distanceTye::Manhattan;
    double HeuristicValue;
    switch (D)
    {
    case 0:
    {
        // Euclidean
        HeuristicValue = (node1->coord - node2->coord).norm();
        break;
    }
    case 1:
    {
        // Manhattan
        HeuristicValue = (node1->coord.array() - node2->coord.array()).abs().sum();
        break;
    }
    case 2:
    {
        double dx = std::abs(node1->coord(0) - node2->coord(0) );
        double dy = std::abs(node1->coord(1) - node2->coord(1) );
        double dz = std::abs(node1->coord(2) - node2->coord(2) );
        double min_xyz = std::min({dx, dy, dz});
        HeuristicValue = dx + dy + dz + (std::sqrt(3.0) - 3) * min_xyz;
    }
    default:
        break;
    }

    if (useTieBreaker){
        double p = 1 / 100;
        HeuristicValue = HeuristicValue * (1 + p);
    }
    
    return HeuristicValue;

}

void AstarPathFinder::AstarGraphSearch(Vector3d start_pt, Vector3d end_pt)
{   
    ros::Time time_1 = ros::Time::now();  // record the search time

    //index of start_point and end_point
    Vector3i start_idx = coord2gridIndex(start_pt);
    Vector3i end_idx   = coord2gridIndex(end_pt);
    goalIdx = end_idx;

    //position of start_point and end_point
    start_pt = gridIndex2coord(start_idx);
    end_pt   = gridIndex2coord(end_idx);

    //Initialize the pointers of struct GridNode which represent start node and goal node
    GridNodePtr startPtr = new GridNode(start_idx, start_pt);
    GridNodePtr endPtr   = new GridNode(end_idx,   end_pt);

    //openSet is the open_list implemented through multimap in STL library
    openSet.clear();
    // currentPtr represents the node with lowest f(n) in the open_list
    GridNodePtr currentPtr  = NULL;
    GridNodePtr neighborPtr = NULL;

    //put start node in open set
    startPtr -> gScore = 0;
    startPtr -> fScore = getHeu(startPtr, endPtr);  
    std::cout << startPtr->fScore << std::endl; 
    startPtr -> id = 1; 
    startPtr -> coord = start_pt;
    openSet.insert( make_pair(startPtr -> fScore, startPtr) );

    GridNodeMap[start_idx[0]][start_idx[1]][start_idx[2]] -> id = 1;

    vector<GridNodePtr> neighborPtrSets; // expand queue
    vector<double> edgeCostSets;         // f cost = g + h

    // this is the main loop
    while ( !openSet.empty() ){
        // Multimap can Automatic sorting, one-to-many data structure, openSet is composed of GridNodeMap[x][y][z]
        currentPtr = openSet.begin() -> second;  
        // GridNodeMap[currentPtr->index[0]][currentPtr->index[1]][currentPtr->index[2]]->id = -1; // set to the close set
        openSet.erase(openSet.begin()); // erase from the open set
        
        if(currentPtr->id == -1)
            continue;
        currentPtr->id = -1;

        if (currentPtr == GridNodeMap[currentPtr->index[0]][currentPtr->index[1]][currentPtr->index[2]]){
            ROS_INFO("Map == openSet");  // equal
        } else{
            ROS_INFO("Not equal!!!!");
        }
        
        // if the current node is the goal 
        if( currentPtr->index == goalIdx ){
            ros::Time time_2 = ros::Time::now();
            terminatePtr = currentPtr;
            ROS_WARN("[A*]{sucess}  Time in A*  is %f ms, path cost if %f m", (time_2 - time_1).toSec() * 1000.0, currentPtr->gScore * resolution );            
            return;
        }
        
        //get the succetion
        AstarGetSucc(currentPtr, neighborPtrSets, edgeCostSets);  //STEP 4: finish AstarPathFinder::AstarGetSucc yourself     
        
        ROS_INFO("size is %d.", ((int)neighborPtrSets.size()));
        for(int i = 0; i < (int)neighborPtrSets.size(); i++){
            /*
            *
            Judge if the neigbors have been expanded
            please write your code below
            IMPORTANT NOTE!!!
            neighborPtrSets[i]->id = -1 : unexpanded
            neighborPtrSets[i]->id = 1 : expanded, equal to this node is in close set
            *        
            */
            neighborPtr = neighborPtrSets[i];
            if(neighborPtr -> id == 1){ //discover a new node, which is not in the closed set and open set

                double currentGCost = currentPtr->gScore + edgeCostSets[i];
                if (currentGCost < neighborPtr->gScore){
                    neighborPtr->gScore = currentGCost;
                    neighborPtr->fScore = neighborPtr -> gScore + getHeu(neighborPtr, endPtr);
                    neighborPtr->cameFrom = currentPtr;
                }
                continue;
            }
            else if(neighborPtr->id == 0){ //this node is in open set and need to judge if it needs to update, the "0" should be deleted when you are coding
                neighborPtr->gScore = edgeCostSets[i] + currentPtr->gScore;
                neighborPtr->fScore = neighborPtr->gScore + getHeu(neighborPtr, endPtr);
                neighborPtr->cameFrom = currentPtr;
                openSet.insert(make_pair(neighborPtr->fScore, neighborPtr));
                neighborPtr->id = 1;
                continue;
            }
            else{//this node is in closed set
                continue;
            }
        }      
    }

    //if search fails
    ros::Time time_2 = ros::Time::now();
    if((time_2 - time_1).toSec() > 0.1)
        ROS_WARN("Time consume in Astar path finding is %f", (time_2 - time_1).toSec() );
    if (currentPtr->coord == GridNodeMap[currentPtr->index[0]][currentPtr->index[1]][currentPtr->index[2]]->coord){
        ROS_INFO("Map == openSet");  // equal
    } else{
        ROS_INFO("Not equal!!!!");
    }
}


vector<Vector3d> AstarPathFinder::getPath() 
{   
    vector<Vector3d> path;
    vector<GridNodePtr> gridPath;

    GridNodePtr tempNode = terminatePtr; // use the front the find the path, backtrack
    while(tempNode -> cameFrom != nullptr){
        gridPath.push_back(tempNode);
        tempNode = tempNode->cameFrom;
    }

    for (auto ptr: gridPath)
        path.push_back(ptr->coord);
        
    reverse(path.begin(),path.end());

    return path;
}