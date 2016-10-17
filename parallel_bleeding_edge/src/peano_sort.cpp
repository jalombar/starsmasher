#include <vector>
#include <cmath>
#include <algorithm>

#include "peano.h"

struct cmp_peanokey_index{
	bool operator() (const peano_hilbert::peano_struct &lhs, const peano_hilbert::peano_struct &rhs){
		return lhs.key < rhs.key;
	}
};

  template< class T >
void reorder(const std::vector<int> &order, std::vector<T> &v)
{
  const std::vector<T> v_orig(v);
  const int order_size = order.size();
  for (int i = 0; i < order_size; i++)
    v[i] = v_orig[order[i]];
}

struct Particle
{
  typedef std::vector<Particle> Vector;
  double x, y, z;
  double vx, vy, vz;
  double hp, u;
};
void peano_sort(
    int local_n,
    double *p_x,  double *p_y,  double *p_z,
    double *p_vx, double *p_vy, double *p_vz,
    double *p_hp, double *p_u)
{
  if (local_n == 0) return;

  Particle::Vector ptcl_local(local_n);
    
  double xmin = +HUGE, ymin = +HUGE, zmin = +HUGE;
  double xmax = -HUGE, ymax = -HUGE, zmax = -HUGE;
  for (int i = 0; i < local_n; i++)
  {
    ptcl_local[i].x = p_x[i];
    ptcl_local[i].y = p_y[i];
    ptcl_local[i].z = p_z[i];
    
    ptcl_local[i].vx = p_vx[i];
    ptcl_local[i].vy = p_vy[i];
    ptcl_local[i].vz = p_vz[i];

    ptcl_local[i].hp = p_hp[i];
    ptcl_local[i].u  = p_u [i];

    xmin = std::min(xmin, p_x[i]);
    ymin = std::min(ymin, p_y[i]);
    zmin = std::min(zmin, p_z[i]);
    
    xmax = std::max(xmax, p_x[i]);
    ymax = std::max(ymax, p_y[i]);
    zmax = std::max(zmax, p_z[i]);
  }

  const double xsize = xmax - xmin;
  const double ysize = ymax - ymin;
  const double zsize = zmax - zmin;

  static std::vector<int> order;
  order.resize(local_n);
  order.clear();

  {
    static std::vector<peano_hilbert::peano_struct> keys;
    keys.resize(local_n);
    keys.clear();


    const float size = std::max(xsize, std::max(ysize, zsize));

    const float domain_fac = 1.0f / size * (((peano_hilbert::peanokey)1) << (BITS_PER_DIMENSION));

    for (int i = 0; i < (int)local_n; i++) 
    {
      const int x = (int)((ptcl_local[i].x - xmin) * domain_fac);
      const int y = (int)((ptcl_local[i].y - ymin) * domain_fac);
      const int z = (int)((ptcl_local[i].z - zmin) * domain_fac);
      keys[i].idx = i;
      keys[i].key = peano_hilbert::peano_hilbert_key(x, y, z, BITS_PER_DIMENSION);
    }

    std::sort(keys.begin(), keys.end(), cmp_peanokey_index());

    for (int i = 0; i < (int)local_n; i++)
      order[i] = keys[i].idx;
  }

  reorder(order, ptcl_local);
  
  for (int i = 0; i < local_n; i++)
  {
    p_x[i] = ptcl_local[i].x;
    p_y[i] = ptcl_local[i].y;
    p_z[i] = ptcl_local[i].z;
    
    p_vx[i] = ptcl_local[i].vx;
    p_vy[i] = ptcl_local[i].vy;
    p_vz[i] = ptcl_local[i].vz;
    
    p_hp[i] = ptcl_local[i].hp;
    p_u[i] = ptcl_local [i].u;
  }
}

extern "C"
{
  void peano_sort_(int *n, 
      double *x, double *y, double *z, 
      double *vx, double *vy, double *vz, 
      double *hp, double *u)
  {
    peano_sort(*n, x,y,z,  vx,vy,vz, hp,u);
  }
}


