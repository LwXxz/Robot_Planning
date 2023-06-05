#include <hw_tool.h>

using namespace std;
using namespace Eigen;

void Homeworktool::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
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
}

void Homeworktool::setObs(const double coord_x, const double coord_y, const double coord_z)
{   
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      
    
    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

bool Homeworktool::isObsFree(const double coord_x, const double coord_y, const double coord_z)
{
    Vector3d pt;
    Vector3i idx;
    
    pt(0) = coord_x;
    pt(1) = coord_y;
    pt(2) = coord_z;
    idx = coord2gridIndex(pt);

    int idx_x = idx(0);
    int idx_y = idx(1);
    int idx_z = idx(2);

    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

Vector3d Homeworktool::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector3i Homeworktool::coord2gridIndex(const Vector3d & pt) 
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d Homeworktool::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

double Homeworktool::OptimalBVP(Eigen::Vector3d _start_position,Eigen::Vector3d _start_velocity,Eigen::Vector3d _target_position)
{
    double optimal_cost = DBL_MAX;
    
    double startx = _start_position[0];
    double starty = _start_position[1];
    double startz = _start_position[2];

    double startxv = _start_velocity[0];
    double startyv = _start_velocity[1];
    double startzv = _start_velocity[2];

    double targetx = _target_position[0];
    double targety = _target_position[1];
    double targetz = _target_position[2];

    double targetxv = 0;
    double targetyv = 0;
    double targetzv = 0;

    Vector4d Fcoeff; // from low frequency
    /*
    F = T + 
        (4.0*startxv**2 + 4.0*startyv**2 + 4.0*startzv**2)/T + 
        (12.0*startx*startxv - 12.0*startxv*targetx + 12.0*starty*startyv - 12.0*startyv*targety + 12.0*startz*startzv - 12.0*startzv*targetz)/T**2 + 
        (12.0*startx**2 - 24.0*startx*targetx + 12.0*starty**2 - 24.0*starty*targety + 12.0*startz**2 - 24.0*startz*targetz + 12.0*targetx**2 + 12.0*targety**2 + 12.0*targetz**2)/T**3
    */
    Fcoeff[0] = 12 * (pow((startx - targetx), 2) + pow((starty - targety), 2) + pow((startz - targetz), 2)); 
    Fcoeff[1] = 12.0 * (startx*startxv - startxv*targetx + starty*startyv - startyv*targety + startz*startzv - startzv*targetz);
    Fcoeff[2] = 4.0 * (pow(startxv, 2) + pow(startyv, 2) + pow(startzv, 2));
    Fcoeff[3] = 1.0;

    Vector4d coeff;
    /*
    F` = 1 + 
        (-4.0*startxv**2 - 4.0*startyv**2 - 4.0*startzv**2)/T**2 + 
        (-24.0*startx*startxv + 24.0*startxv*targetx - 24.0*starty*startyv + 24.0*startyv*targety - 24.0*startz*startzv + 24.0*startzv*targetz)/T**3 + 
        (-36.0*startx**2 + 72.0*startx*targetx - 36.0*starty**2 + 72.0*starty*targety - 36.0*startz**2 + 72.0*startz*targetz - 36.0*targetx**2 - 36.0*targety**2 - 36.0*targetz**2)/T**4
    */
    coeff[0] = -36 * (pow((startx - targetx), 2) + pow((starty - targety), 2) + pow((startz - targetz), 2));
    coeff[1] = 24 * (startx*startxv - startxv*targetx + starty*startyv - startyv*targety + startz*startzv - startzv*targetz);
    coeff[2] = -4.0 * (pow(startxv, 2) + pow(startyv, 2) + pow(startzv, 2));
    coeff[3] = 1.0;

    // Use an adjoint matrix
    Matrix<double, 4, 4> C;
    C << 0, 0, 0, -coeff[0],
        1, 0, 0, -coeff[1],
        0, 1, 0, -coeff[2],
        0, 0, 1, -coeff[3];

    Matrix<complex<double>, Dynamic, Dynamic> eigen_values ;  // The size of the matrix is determined ai runtime.
    eigen_values = C.eigenvalues();

    for (int i = 0; i < eigen_values.size(); i++){
        if (abs(eigen_values(i).imag()) != 0.0 || eigen_values(i).real() < 0.0) continue;
        double T = eigen_values(i).real();

        double real_cost = Fcoeff[0] / pow(T, 3) + Fcoeff[1] / pow(T, 2) + Fcoeff[2] / T + T;

        optimal_cost = optimal_cost < real_cost ? optimal_cost: real_cost;
    }

    // Use the eigen API
    // PolynomialSolver<double, Dynamic> solver;
    // solver.compute(coeff);
    // const PolynomialSolver<double, Dynamic>::RootsType& r = solver.roots();

    // for (int i = 0; i < r.rows(); i++){
    //     if (r[i].imag() != 0.0 || r[i].real() < 0) continue;
    //     double T = r[i].real();
    //     double real_cost = Fcoeff[0] / pow(T, 3) + Fcoeff[1] / pow(T, 2) + Fcoeff[2] / T + T;
    //     optimal_cost = optimal_cost < real_cost ? optimal_cost: real_cost;
    // }

    return optimal_cost;
}
