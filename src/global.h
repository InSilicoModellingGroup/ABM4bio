// =============================================================================
//
//   Copyright (C) 2020-2024 Vasileios Vavourakis (vasvav@gmail.com)
//   All Rights Reserved.
//
//   Licensed under the GNU General Public License v3.0 (the "License").
//   See the LICENSE file provided in this project details the License.
//   You cannot use this file except in compliance with the License.
//
// =============================================================================

// =============================================================================
#ifndef _GLOBAL_H_
#define _GLOBAL_H_
// =============================================================================
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <algorithm>
#include <functional>
#include <random>
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cstddef>
#include <typeinfo>
#include <sys/stat.h>
// =============================================================================
#include "biodynamo.h"
#include "core/util/log.h"
#include "core/resource_manager.h"
#include "core/substance_initializers.h"
#include "core/agent/agent.h"
#include "core/agent/cell_division_event.h"
#include "neuroscience/param.h"
#include "neuroscience/neuroscience.h"
#include "neuroscience/neurite_element.h"
// =============================================================================
namespace bdm {
// =============================================================================
typedef bdm::MathArray<bdm::MathArray<double, 3>, 3>   Double3x3;
typedef bdm::neuroscience::NeuriteElement       TubularElement;
typedef bdm::AgentPointer<bdm::TubularElement>  TubularElementPointer;
// =============================================================================
} // ...end of namespace
// =============================================================================
#ifndef ABORT_
#define ABORT_(_msg_) \
  bdm::Log::Fatal( " *** @file "+std::string(__FILE__)+" @line "+std::to_string(__LINE__), \
                   std::string(_msg_));
#endif
#ifndef ASSERT_
#define ASSERT_(_cond_,_msg_) \
  if (!(_cond_)) \
    bdm::Log::Fatal( " *** @file "+std::string(__FILE__)+" @line "+std::to_string(__LINE__), \
                     std::string(_msg_));
#endif
// =============================================================================
#include "./model_parameters.h"
#include "./ida.h"
// =============================================================================
inline
void status_bar(std::ostream& os, const int x, const int n, int barWidth =50) {
  const double progress = x / static_cast<double>(n);
  const int pos = barWidth * progress;
  os << "[";
  for (int i=0; i<barWidth; i++)
    {
      if      (i <  pos)
        os << "=";
      else if (i == pos)
        os << ">";
      else
        os << " ";
    }
  os << "]";
  os << std::fixed << std::setprecision(3) << (progress*100.0) << "%\r" << std::flush;
  if ( n == x )
    {
      for (int i=0; i<barWidth+50; i++)
        os << ' ';
      os << '\r' << std::flush;
    }
}
inline
void time_status_bar(std::ostream& os, const int x, const int n, double T, int barWidth =50) {
  const double progress = x / static_cast<double>(n);
  const int pos = barWidth * progress;
  os << "[";
  for (int i=0; i<barWidth; i++)
    {
      if      (i <  pos)
        os << "=";
      else if (i == pos)
        os << ">";
      else
        os << " ";
    }
  os << "]";
  os << " " << std::fixed << std::setprecision(3) << (progress*100.0)
     << "% (" << std::scientific << T << ")" << '\r';
  if ( n == x )
    {
      for (int i=0; i<barWidth+50; i++)
        os << ' ';
      os  << '\r';
    }
  os << std::flush;
}
// =============================================================================
inline double pow2(const double& v) { return v*v; }
inline double pow3(const double& v) { return v*pow2(v); }
inline double pow4(const double& v) { return pow2(pow2(v)); }
inline double pow5(const double& v) { return v*pow4(v); }
inline double pow6(const double& v) { return pow2(pow3(v)); }
inline double pow7(const double& v) { return v*pow6(v); }
inline double pow8(const double& v) { return pow2(pow4(v)); }
inline double pow9(const double& v) { return v*pow8(v); }
// =============================================================================
inline int sign(const double r) { return ( r>=0.0 ? +1 : -1 ); }
inline double Macaulay(double (*f)(double), const double x) { return (f(x)<0.0 ? 0.0 : f(x)); }
// =============================================================================
inline double degrees_to_radians(const double& d) { return (bdm::Math::kPi*d/180.0); }
inline double radians_to_degrees(const double& r) { return (180.0*r/bdm::Math::kPi); }
// =============================================================================
inline
double apply_lbound(const double& L, const double& X) { return (X < L ? L : X); }
inline
double apply_ubound(const double& X, const double& U) { return (X > U ? U : X); }
inline
double apply_bounds(const double& L, const double& X, const double& U) { return (X < L ? L : (X > U ? U : X)); }
// =============================================================================
inline
int uniform_distro(const int from, const int to)
{
  std::random_device rd;
  std::default_random_engine rgen(rd());
  std::uniform_int_distribution<int> dist(from, to);
  return dist(rgen);
}
inline
double uniform_distro(const double from, const double to)
{
  std::random_device rd;
  std::default_random_engine rgen(rd());
  std::uniform_real_distribution<double> dist(from, to);
  return dist(rgen);
}
inline
double normal_distro(const double mean, const double stdev)
{
  std::random_device rd;
  std::mt19937 rgen(rd());
  std::normal_distribution<double> dist(mean, stdev);
  return dist(rgen);
}
inline
std::string random_string(int len =8)
{
  std::string str("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
  std::random_device rd;
  std::mt19937 generator(rd());
  std::shuffle(str.begin(), str.end(), generator);
  return str.substr(0, len+1);
}
// =============================================================================
inline
bdm::Double3x3 zero()
{
  bdm::Double3x3 r;
  r[0][0] = r[0][1] = r[0][2] =
  r[1][0] = r[1][1] = r[1][2] =
  r[2][0] = r[2][1] = r[2][2] = 0.0;
  return r;
}
inline
bdm::Double3x3 ones()
{
  bdm::Double3x3 r;
  r[0][0] = r[0][1] = r[0][2] =
  r[1][0] = r[1][1] = r[1][2] =
  r[2][0] = r[2][1] = r[2][2] = 1.0;
  return r;
}
inline
bdm::Double3x3 eye()
{
  bdm::Double3x3 r = zero();
  r[0][0] = r[1][1] = r[2][2] = 1.0;
  return r;
}
inline
bdm::Double3x3 diag(double d0, double d1, double d2)
{
  bdm::Double3x3 r = zero();
  r[0][0] = d0; r[1][1] = d1; r[2][2] = d2;
  return r;
}
inline
bdm::Double3x3 tensor(const bdm::Double3& u, const bdm::Double3& v)
{
  bdm::Double3x3 r;
  r[0][0] = u[0]*v[0]; r[0][1] = u[0]*v[1]; r[0][2] = u[0]*v[2];
  r[1][0] = u[1]*v[0]; r[1][1] = u[1]*v[1]; r[1][2] = u[1]*v[2];
  r[2][0] = u[2]*v[0]; r[2][1] = u[2]*v[1]; r[2][2] = u[2]*v[2];
  return r;
}
inline
bdm::Double3x3 tensor(const bdm::Double3& u, const bdm::Double3& v, double l)
{
  bdm::Double3x3 r;
  r[0][0] = u[0]*v[0]*l; r[0][1] = u[0]*v[1]*l; r[0][2] = u[0]*v[2]*l;
  r[1][0] = u[1]*v[0]*l; r[1][1] = u[1]*v[1]*l; r[1][2] = u[1]*v[2]*l;
  r[2][0] = u[2]*v[0]*l; r[2][1] = u[2]*v[1]*l; r[2][2] = u[2]*v[2]*l;
  return r;
}
inline
bdm::Double3x3 add(const bdm::Double3x3& a, const bdm::Double3x3& b)
{
  bdm::Double3x3 r;
  r[0][0] = a[0][0]+b[0][0]; r[0][1] = a[0][1]+b[0][1]; r[0][2] = a[0][2]+b[0][2];
  r[1][0] = a[1][0]+b[1][0]; r[1][1] = a[1][1]+b[1][1]; r[1][2] = a[1][2]+b[1][2];
  r[2][0] = a[2][0]+b[2][0]; r[2][1] = a[2][1]+b[2][1]; r[2][2] = a[2][2]+b[2][2];
  return r;
}
// =============================================================================
inline
double L2norm(const bdm::Double3& v)
{
  return sqrt(pow2(v[0])+pow2(v[1])+pow2(v[2]));
}
inline
bool normalize(const bdm::Double3& v, bdm::Double3& r)
{
  double v_mag = L2norm(v);
  if (v_mag>0.0)
    {
      r = (v*(1.0/v_mag));
      return true;
    }
  r = bdm::Double3();
  return false;
}
inline
bdm::Double3 normalize(const bdm::Double3& v)
{
  bdm::Double3 r;
  normalize(v, r);
  return r;
}
inline
double relative_angle(const bdm::Double3& u, const bdm::Double3& v)
{
  return acos((u*v)/L2norm(u)/L2norm(v));
}
inline
bdm::Double3 cross(const bdm::Double3& u, const bdm::Double3& v)
{
  return { (u[1]*v[2]-u[2]*v[1]) ,
          -(u[0]*v[2]-u[2]*v[0]) ,
           (u[0]*v[1]-u[1]*v[0]) };
}
inline
bdm::Double3 project_to_plane(const bdm::Double3& normal, const bdm::Double3& v)
{
  const double e2 = pow2(L2norm(normal));
  bdm::Double3 v_X_e = cross(v, normal),
               e_X_v_X_e = cross(normal, v_X_e);
  return e_X_v_X_e / e2;
}
inline
bool line_intersects_plane(const bdm::Double3& l0, const bdm::Double3& l1,
                           const bdm::Double3& N, const bdm::Double3& c,
                           bdm::Double3& i, double tol = 1.0e-5)
{
  bdm::Double3 u = l1 - l0;
  u.Normalize();
  //
  const bdm::Double3 n = N * (1.0/L2norm(N));
  //
  if (fabs(n*u)<=tol) return false;
  // see: https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
  const bdm::Double3 v = c - l0;
  //
  const double v_DOT_n = (v*n) / (u*n);
  //
  i = l0 + u * v_DOT_n;
  //
  return true;
}
inline
bdm::Double3 project_to_line(const bdm::Double3& l0, const bdm::Double3& l1,
                             const bdm::Double3& p)
{
  bdm::Double3 u = l1 - l0;
  u.Normalize();
  //
  const double u_DOT_p_l0 = u * (p - l0);
  //
  return (l0 + u * u_DOT_p_l0);
}
inline
bdm::Double3 project_to_plane(const bdm::Double3& N, const bdm::Double3& c,
                              const bdm::Double3& p, double tol = 1.0e-5)
{
  const bdm::Double3 n = N * (1.0/L2norm(N));
  //
  if (fabs(n*p)<=tol) return p;
  // see: https://exceptionshub.com/how-to-project-a-point-onto-a-plane-in-3d.html
  const bdm::Double3 v = p - c;
  const double d = n * v;
  //
  return (p - n * d);
}
inline
double triangle_area(const bdm::Double3& A, const bdm::Double3& B, const bdm::Double3& C)
{
  const bdm::Double3 AB = B - A,
                     AC = C - A;
  // see: https://en.wikipedia.org/wiki/Cross_product#Matrix_notation
  return (0.5*L2norm(cross(AB, AC))); // ...(half) cross product formula
}
inline
bool is_inside_segment (const bdm::Double3& s0, const bdm::Double3& s1,
                        const bdm::Double3& p)
{
  bdm::Double3 u = s1 - s0;
  const double L = L2norm(u);
  u.Normalize();
  //
  const double u_DOT_p_s0 = u * (p - s0);
  //
  if ( u_DOT_p_s0<0.0 || u_DOT_p_s0>L ) return false;
  return true;
}
inline
bool is_inside_triangle(const bdm::Double3& A, const bdm::Double3& B, const bdm::Double3& C,
                        const bdm::Double3& p, double tol = 1.0e-5)
{
  const double area1 = triangle_area(p, B, C),
               area2 = triangle_area(A, p, C),
               area3 = triangle_area(A, B, p);
  // when all 3 triangles are forming the whole area of triangle ABC
  return (fabs(triangle_area(A, B, C)-(area1+area2+area3))<tol);
}
inline
bool is_inside_cube(const bdm::Double3& A, const bdm::Double3& B, const bdm::Double3& C, const bdm::Double3& D,
                    const bdm::Double3& E, const bdm::Double3& F, const bdm::Double3& G, const bdm::Double3& H,
                    const bdm::Double3& p)
{
  const bdm::Double3 M = (A+G) * 0.5, // box center
                     V = p - M; // direction vector from point to center
  bdm::Double3 Xl = B - A, // local axes of the box
               Yl = D - A,
               Zl = E - A;
  const double L_X = L2norm(Xl),
               L_Y = L2norm(Yl),
               L_Z = L2norm(Zl);
  Xl /= L_X;
  Yl /= L_Y;
  Zl /= L_Z;
  const double pX = 2.0 * fabs(V*Xl),
               pY = 2.0 * fabs(V*Yl),
               pZ = 2.0 * fabs(V*Zl);
  return (pX<=L_X && pY<=L_Y && pZ<=L_Z);
}
inline
bool is_inside_sphere(const bdm::Double3& C, const double& R,
                      const bdm::Double3& p)
{
  const bdm::Double3 V = p - C; // direction vector from point to center
  const double r = L2norm(V);
  return (r<=R);
}
inline
bdm::Double3 rotate(const bdm::Double3& v, double theta_x, double theta_y, double theta_z)
{
  const double v_x = v[0],
               v_y = v[1],
               v_z = v[2];
  const double S_x = sin(theta_x), C_x = cos(theta_x),
               S_y = sin(theta_y), C_y = cos(theta_y),
               S_z = sin(theta_z), C_z = cos(theta_z);
  bdm::Double3 r;
  r[0] = v_z*(S_x*S_z + C_x*C_z*S_y) - v_y*(C_x*S_z - C_z*S_x*S_y) + C_y*C_z*v_x;
  r[1] = v_y*(C_x*C_z + S_x*S_y*S_z) - v_z*(C_z*S_x - C_x*S_y*S_z) + C_y*S_z*v_x;
  r[2] = C_x*C_y*v_z - S_y*v_x + C_y*S_x*v_y;
  return r;
}
inline
bdm::Double3x3 rotate(const bdm::Double3x3& m, double theta_x, double theta_y, double theta_z)
{
  const double m_xx = m[0][0], m_xy = m[0][1], m_xz = m[0][2],
               m_yx = m[1][0], m_yy = m[1][1], m_yz = m[1][2],
               m_zx = m[2][0], m_zy = m[2][1], m_zz = m[2][2];
  const double S_x = sin(theta_x), C_x = cos(theta_x),
               S_y = sin(theta_y), C_y = cos(theta_y),
               S_z = sin(theta_z), C_z = cos(theta_z);
  bdm::Double3x3 r;
  r[0][0] = (S_x*S_z + C_x*C_z*S_y)*(m_zz*(S_x*S_z + C_x*C_z*S_y) - m_yz*(C_x*S_z - C_z*S_x*S_y) + C_y*C_z*m_xz) - (C_x*S_z - C_z*S_x*S_y)*(m_zy*(S_x*S_z + C_x*C_z*S_y) - m_yy*(C_x*S_z - C_z*S_x*S_y) + C_y*C_z*m_xy) + C_y*C_z*(m_zx*(S_x*S_z + C_x*C_z*S_y) - m_yx*(C_x*S_z - C_z*S_x*S_y) + C_y*C_z*m_xx);
  r[0][1] = (C_x*C_z + S_x*S_y*S_z)*(m_zy*(S_x*S_z + C_x*C_z*S_y) - m_yy*(C_x*S_z - C_z*S_x*S_y) + C_y*C_z*m_xy) - (C_z*S_x - C_x*S_y*S_z)*(m_zz*(S_x*S_z + C_x*C_z*S_y) - m_yz*(C_x*S_z - C_z*S_x*S_y) + C_y*C_z*m_xz) + C_y*S_z*(m_zx*(S_x*S_z + C_x*C_z*S_y) - m_yx*(C_x*S_z - C_z*S_x*S_y) + C_y*C_z*m_xx);
  r[0][2] = C_x*C_y*(m_zz*(S_x*S_z + C_x*C_z*S_y) - m_yz*(C_x*S_z - C_z*S_x*S_y) + C_y*C_z*m_xz) - S_y*(m_zx*(S_x*S_z + C_x*C_z*S_y) - m_yx*(C_x*S_z - C_z*S_x*S_y) + C_y*C_z*m_xx) + C_y*S_x*(m_zy*(S_x*S_z + C_x*C_z*S_y) - m_yy*(C_x*S_z - C_z*S_x*S_y) + C_y*C_z*m_xy);
  r[1][0] = (S_x*S_z + C_x*C_z*S_y)*(m_yz*(C_x*C_z + S_x*S_y*S_z) - m_zz*(C_z*S_x - C_x*S_y*S_z) + C_y*S_z*m_xz) - (C_x*S_z - C_z*S_x*S_y)*(m_yy*(C_x*C_z + S_x*S_y*S_z) - m_zy*(C_z*S_x - C_x*S_y*S_z) + C_y*S_z*m_xy) + C_y*C_z*(m_yx*(C_x*C_z + S_x*S_y*S_z) - m_zx*(C_z*S_x - C_x*S_y*S_z) + C_y*S_z*m_xx);
  r[1][1] = (C_x*C_z + S_x*S_y*S_z)*(m_yy*(C_x*C_z + S_x*S_y*S_z) - m_zy*(C_z*S_x - C_x*S_y*S_z) + C_y*S_z*m_xy) - (C_z*S_x - C_x*S_y*S_z)*(m_yz*(C_x*C_z + S_x*S_y*S_z) - m_zz*(C_z*S_x - C_x*S_y*S_z) + C_y*S_z*m_xz) + C_y*S_z*(m_yx*(C_x*C_z + S_x*S_y*S_z) - m_zx*(C_z*S_x - C_x*S_y*S_z) + C_y*S_z*m_xx);
  r[1][2] = C_x*C_y*(m_yz*(C_x*C_z + S_x*S_y*S_z) - m_zz*(C_z*S_x - C_x*S_y*S_z) + C_y*S_z*m_xz) - S_y*(m_yx*(C_x*C_z + S_x*S_y*S_z) - m_zx*(C_z*S_x - C_x*S_y*S_z) + C_y*S_z*m_xx) + C_y*S_x*(m_yy*(C_x*C_z + S_x*S_y*S_z) - m_zy*(C_z*S_x - C_x*S_y*S_z) + C_y*S_z*m_xy);
  r[2][0] = (S_x*S_z + C_x*C_z*S_y)*(C_x*C_y*m_zz - S_y*m_xz + C_y*S_x*m_yz) - (C_x*S_z - C_z*S_x*S_y)*(C_x*C_y*m_zy - S_y*m_xy + C_y*S_x*m_yy) + C_y*C_z*(C_x*C_y*m_zx - S_y*m_xx + C_y*S_x*m_yx);
  r[2][1] = (C_x*C_z + S_x*S_y*S_z)*(C_x*C_y*m_zy - S_y*m_xy + C_y*S_x*m_yy) - (C_z*S_x - C_x*S_y*S_z)*(C_x*C_y*m_zz - S_y*m_xz + C_y*S_x*m_yz) + C_y*S_z*(C_x*C_y*m_zx - S_y*m_xx + C_y*S_x*m_yx);
  r[2][2] = C_x*C_y*(C_x*C_y*m_zz - S_y*m_xz + C_y*S_x*m_yz) - S_y*(C_x*C_y*m_zx - S_y*m_xx + C_y*S_x*m_yx) + C_y*S_x*(C_x*C_y*m_zy - S_y*m_xy + C_y*S_x*m_yy);
  return r;
}
inline
bdm::Double3 rotate_about_axis(bdm::Double3 axis, double angle, const bdm::Double3& X)
{
  axis.Normalize();
  // Euler-Rodrigues formula:
  //    https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula
  const double sin_phi = sin(0.5*angle), cos_phi = cos(0.5*angle),
               a = cos_phi,
               b = axis[0]*sin_phi,
               c = axis[1]*sin_phi,
               d = axis[2]*sin_phi;
  bdm::Double3x3 r;
  r[0][0] = a*a+b*b-c*c-d*d;
  r[0][1] = 2.0*(b*c-a*d);
  r[0][2] = 2.0*(b*d+a*c);
  r[1][0] = 2.0*(b*c+a*d);
  r[1][1] = a*a-b*b+c*c-d*d;
  r[1][2] = 2.0*(c*d-a*b);
  r[2][0] = 2.0*(b*d-a*c);
  r[2][1] = 2.0*(c*d+a*b);
  r[2][2] = a*a-b*b-c*c+d*d;
  return { r[0][0]*X[0]+r[0][1]*X[1]+r[0][2]*X[2] ,
           r[1][0]*X[0]+r[1][1]*X[1]+r[1][2]*X[2] ,
           r[2][0]*X[0]+r[2][1]*X[1]+r[2][2]*X[2] };
}
inline
void check_agent_position_in_domain(const double D_min, const double D_max,
                                    const bdm::Double3& X,
                                    std::map<double, bdm::Double3>& map_proj)
{
  std::vector<bdm::Double3> N(6);
  N[0] = { +1.0,  0.0,  0.0 };
  N[1] = { -1.0,  0.0,  0.0 };
  N[2] = { 0.0,  +1.0,  0.0 };
  N[3] = { 0.0,  -1.0,  0.0 };
  N[4] = { 0.0,  0.0,  +1.0 };
  N[5] = { 0.0,  0.0,  -1.0 };
  std::vector<bdm::Double3> c(6);
  c[0] = { D_min, 0.0,   0.0   };
  c[1] = { D_max, 0.0,   0.0   };
  c[2] = { 0.0,   D_min, 0.0   };
  c[3] = { 0.0,   D_max, 0.0   };
  c[4] = { 0.0,   0.0,   D_min };
  c[5] = { 0.0,   0.0,   D_max };
  map_proj.clear();
  for (unsigned int s=0; s<6; s++)
    {
      bdm::Double3 p, d;
      {
        const bdm::Double3 v = X - c[s];
        //
        d = N[s] * (N[s]*v);
        p = X - d;
      }
      if ((d*N[s])<=0.0) map_proj.insert(std::make_pair(0.0,       p));
      else               map_proj.insert(std::make_pair(L2norm(d), p));
    }
}
inline
bool check_agent_position_in_domain(const double D_min, const double D_max,
                                    const bdm::Double3& X, const double& tol)
{
  std::map<double, bdm::Double3> map_proj;
  check_agent_position_in_domain(D_min, D_max, X, map_proj);
  //
  for (auto proj : map_proj)
    {
      const double distance = proj.first;
      // check if agent position is very close to simulation bounds
      if (distance<=tol) return false;
    }
  return true;
}
// =============================================================================
inline
std::vector<std::string> extract_words_vector(const std::string& Text)
{
  std::vector<std::string> Words;
  std::stringstream ss(Text);
  std::string Buf;
  while (ss >> Buf)
    {
      Words.push_back(Buf);
    }
  return Words;
}
// =============================================================================
inline
double flat_function(const double* param, const double t)
{
  const double amp = param[1]; // amplitude of time-function
  return amp;
}
inline
double linear_function(const double* param, const double t)
{
  const double amp = param[1]; // amplitude of time-function
  return amp * t;
}
inline
double step_function(const double* param, const double t)
{
  const double amp = param[1]; // amplitude of time-function
  const double ti  = param[2]; // initialization time
  const double td  = param[3]; // duration time for flat response
  if      (t<=ti)    return 0.0;
  else if (t<=ti+td) return amp;
  else               return 0.0;
}
inline
double ramp_function(const double* param, const double t)
{
  const double amp = param[1]; // amplitude of time-function
  const double ti  = param[2]; // initialization time
  const double ta  = param[3]; // arrival (max amp) time
  const double td  = param[4]; // duration time for flat response
  if      (t<=ti)       return 0.0;
  else if (t<=ti+ta)    return amp * (t-ti)/ta;
  else if (t<=ti+td)    return amp;
  else if (t<=ti+td+ta) return amp * (td+ta-(t-ti))/ta;
  else                  return 0.0;
}
inline
double gaussian_function(const double* param, const double t)
{
  const double amp = param[1]; // amplitude of time-function
  const double tp  = param[2]; // peak time
  const double rms = param[3]; // gaussian RMS width
  return amp * exp(-0.5*pow2((t-tp)/rms));
}
inline
double logistic_function(const double* param, const double t)
{
  const double amp = param[1]; // amplitude of time-function
  const double tm0 = param[2]; // mid-point time (increase)
  const double tm1 = param[3]; // mid-point time (decrease)
  const double cs  = param[4]; // curve steepness
  if (t<0.5*(tm0+tm1)) return amp / (1.0+exp(-cs*(t-tm0)));
  else                 return amp / (1.0+exp(+cs*(t-tm1)));
}
inline
double linear_periodic_function(const double* param, const double t)
{
  const double period = param[2]; // period
  const int N = std::floor(t/period); // total number of periods run
  return linear_function(param, t-N*period);
}
inline
double step_periodic_function(const double* param, const double t)
{
  const double period = param[4]; // period
  const int N = std::floor(t/period); // total number of periods run
  return step_function(param, t-N*period);
}
inline
double ramp_periodic_function(const double* param, const double t)
{
  const double period = param[5]; // period
  const int N = std::floor(t/period); // total number of periods run
  return ramp_function(param, t-N*period);
}
inline
double gaussian_periodic_function(const double* param, const double t)
{
  const double period = param[4]; // period
  const int N = std::floor(t/period); // total number of periods run
  return gaussian_function(param, t-N*period);
}
inline
double logistic_periodic_function(const double* param, const double t)
{
  const double period = param[5]; // period
  const int N = std::floor(t/period); // total number of periods run
  return logistic_function(param, t-N*period);
}
// =============================================================================
#endif // _GLOBAL_H_
// =============================================================================
