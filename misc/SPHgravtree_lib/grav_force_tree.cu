#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "libSequoia/include/sequoiaInterface.h"
#include "vector3.h"

struct Particle
{
  typedef std::vector<Particle> Vector;
  vec3   pos;
  double mass, h, hflag;
  Particle() {}
  Particle(const vec3 &_pos, double _mass, const double _h, const double _hf) :
    pos(_pos), mass(_mass), h(_h), hflag(_hf) {}
};

struct Force
{
  typedef std::vector<Force> Vector;
  vec3 acc;
  double pot;
  Force() {}
  Force(const vec3 &_acc, const double _pot = 0) : acc(_acc), pot(_pot) {}
};

struct SPHgrav_tree
{
  Particle::Vector ptcl;
  Force   ::Vector force;
  
  cudaEvent_t start, stop;

  my_dev::context devContext;
    
  /* Input */
  my_dev::dev_mem<real4> j_bodies_pos;
  my_dev::dev_mem<real4> j_bodies_h  ;
  my_dev::dev_mem<int>   j_bodies_IDs;

  my_dev::dev_mem<real4> i_bodies_pos;
  my_dev::dev_mem<real4> i_bodies_h  ;
  my_dev::dev_mem<int>   i_bodies_IDs;

  /* Output */
  my_dev::dev_mem<real4> bodies_acc;
  my_dev::dev_mem<real>  bodies_ds2;
  my_dev::dev_mem<int>   bodies_ngb;



  int nibeg, ni;
  int ndevice;
  std::vector<bool> can_use_device;
  char *argv[1];
  bool init_context;
  SPHgrav_tree(const int rank = 0)
  {
    init_context = false;
    assert(cudaGetDeviceCount(&ndevice) == 0);
    if (rank == 0)
      fprintf(stderr, " SPHgrav found %d CUDA devices \n", ndevice);
    assert(ndevice > 0);
    can_use_device = std::vector<bool>(ndevice, false);
    int no_supported = 0;
    for (int dev = 0; dev < ndevice; dev++)
    {
      cudaDeviceProp p;
      assert(cudaGetDeviceProperties(&p, dev) == cudaSuccess);
      const bool supported = p.major > 1 || (p.major == 1 && p.minor >= 3);
      if (rank == 0)
        fprintf(stderr,"  Device= %d: %s computeMode= %d computeCapability= %d_%d  supported= %s\n", 
            dev, p.name, p.computeMode, p.major, p.minor, supported ? "YES" : "NO");
      no_supported += supported ? 1 : 0;
      can_use_device[dev] = true;
    }
    assert(no_supported > 0);

    argv[0] = new char[256];
    sprintf(argv[0], "./");
#if 1
    sprintf(argv[0], "./src/SPHgravtree_lib/libSequoia/");
#endif
  }
  ~SPHgrav_tree()
  {
    delete[] argv[0];
    sequoia_cleanup();
  }

  void initDevice(const int device, const double theta = 0.3)
  {
    assert(!init_context);
    init_context = true;
    assert(device < ndevice);
    assert(can_use_device[device]);
    fprintf(stderr, " >>> Setting device %d , theta= %g <<< \n", device, theta);
    devContext =  sequoia_init(argv, device, theta, 0.0);
   
    /* Input */ 
    j_bodies_pos.setContext(devContext);    //Bodies positions
    j_bodies_h  .setContext(devContext);    //Bodies range
    j_bodies_IDs.setContext(devContext);    //Bodies idx
    
    i_bodies_pos.setContext(devContext);    //Bodies positions
    i_bodies_h  .setContext(devContext);    //Bodies range
    i_bodies_IDs.setContext(devContext);    //Bodies idx
    
    /* Output */
    bodies_acc.setContext(devContext);    //Bodies Accelerations
    bodies_ds2.setContext(devContext);    //Bodies distance to nearest neighbour squared
    bodies_ngb.setContext(devContext);    //Bodies nearest neighbour
  }

  void first_half(const int nibeg, const int ni)
  {
    assert(init_context);
    cudaEventCreate( &start );
    cudaEventCreate( &stop  );
    cudaEventRecord( start, 0 );

    this->nibeg = nibeg;
    this->ni    = ni;
    force.resize(ni);

    int nj = ptcl.size();

    /* Input */
    j_bodies_pos.cmalloc(nj);    //Bodies positions
    j_bodies_h  .cmalloc(nj);    //Bodies range
    j_bodies_IDs.cmalloc(nj);    //Bodies idx

    //Jeroen changed this to size ni, since its possible
    //that ni is not as large as nj in parallel code
    i_bodies_pos.cmalloc(ni);    //Bodies positions
    i_bodies_h  .cmalloc(ni);    //Bodies range
    i_bodies_IDs.cmalloc(ni);    //Bodies idx

    /* Output */
    bodies_acc.cmalloc(ni);    //Bodies Accelerations
    bodies_ds2.cmalloc(ni);    //Bodies distance to nearest neighbour squared
    bodies_ngb.cmalloc(ni);    //Bodies nearest neighbour

    for (int i = 0; i < nj; i++)
    {
      j_bodies_pos[i].x = ptcl[i].pos.x;
      j_bodies_pos[i].y = ptcl[i].pos.y;
      j_bodies_pos[i].z = ptcl[i].pos.z;
      j_bodies_pos[i].w = ptcl[i].mass;

      j_bodies_h[i].x =     ptcl[i].h;
      j_bodies_h[i].y = 1.0/ptcl[i].h;
      j_bodies_h[i].z =     ptcl[i].hflag;

      j_bodies_IDs[i] = i;
    }

    j_bodies_pos.h2d();
    j_bodies_h  .h2d();
    j_bodies_IDs.h2d();

    for (int i = 0; i < ni; i++)
    {
      const int j = i + nibeg;

      i_bodies_pos[i].x = ptcl[j].pos.x; //Jeroen, first i was j, copy paste error
      i_bodies_pos[i].y = ptcl[j].pos.y;
      i_bodies_pos[i].z = ptcl[j].pos.z;
      i_bodies_pos[i].w = ptcl[j].mass;

      i_bodies_h[i].x =     ptcl[j].h;
      i_bodies_h[i].y = 1.0/ptcl[j].h;
      i_bodies_h[i].z =     ptcl[j].hflag;

      i_bodies_IDs[i] = j;
    }

    i_bodies_pos.h2d();
    i_bodies_h  .h2d();
    i_bodies_IDs.h2d();

    sequoia_setParticlesAndGetGravity(
        j_bodies_pos,      //Positions J-particles
        j_bodies_h,       
        j_bodies_IDs,      //Particle IDs J-particles
        nj,         //Number of J-particles
        i_bodies_pos,      //Positions I-particles (Can be the same or different than J-particles)
        i_bodies_h,       
        i_bodies_IDs,      //Particle IDs J-particles (Can be the same or different than J-particles)
        ni,                  //Number of I-particles (Can be the same or different than J-particles)
        true,              //Do we need to sort J-particles?
        true,              //Do we need to sort I-particles? Can be false if i-bodies are the same as j-bodies           
        bodies_acc,        //OUT  Accelerations for I-particles
        bodies_ds2,        //OUT  min distance squared for I-particles
        bodies_ngb);       //OUT  J-ID of the nearest neighbour for I-particles

  }

  float last_half()
  {
    bodies_acc.d2h();
    i_bodies_IDs.d2h();

    for (int i = 0; i < ni; i++)
    {
      const int idx = i_bodies_IDs[i];
      assert(idx >= nibeg);
      assert(idx < ni+nibeg);
      //force[idx] = Force(vec3(bodies_acc[i].x, bodies_acc[i].y, bodies_acc[i].z), bodies_acc[i].w);
      //Jeroen, idx-nibeg otherwise there is no valid memory location if nibeg > 0 and all kind of bad
      //things start to happen
      force[idx-nibeg] = Force(vec3(bodies_acc[i].x, bodies_acc[i].y, bodies_acc[i].z), bodies_acc[i].w);
    }

    cudaEventRecord( stop, 0 );
    cudaThreadSynchronize();
    float elapsed_time_ms;
    cudaEventElapsedTime( &elapsed_time_ms, start, stop );
    
    /* Input */
    j_bodies_pos.free_mem();
    j_bodies_h  .free_mem();
    j_bodies_IDs.free_mem();

    i_bodies_pos.free_mem();
    i_bodies_h  .free_mem();
    i_bodies_IDs.free_mem();

    /* Output */
    bodies_acc.free_mem();
    bodies_ds2.free_mem();
    bodies_ngb.free_mem();


    return elapsed_time_ms/1000.0;
  }
};

namespace SPHgrav
{
  SPHgrav_tree grav;

  void firsthalf_grav_force(
      const int nj, const int nibeg, const int ni,
      double *px, double *py, double *pz,
      double *mass, double *h2)
  {
    grav.ptcl.resize(nj);

    assert(nibeg >= 0);

    for (int i = 0; i < nj; i++)
    {
      assert(h2  [i] != 0.0);
      assert(mass[i] >  0.0);
      grav.ptcl[i] = Particle(vec3(px[i], py[i], pz[i]), mass[i], 
          std::sqrt(h2[i] > 0.0 ? h2[i] : -h2[i]),
          h2[i] > 0.0 ? 1.0 : 0.0);
    }

    grav.first_half(nibeg, ni);
  }

  void lasthalf_grav_force(double *ax, double *ay, double *az, double *pot)
  {
    const int ni  = grav.ni;
    const float dt = grav.last_half();

    for (int i = 0; i < ni; i++)
    {
      const int j = i + grav.nibeg;
      ax [j] = grav.force[i].acc.x;
      ay [j] = grav.force[i].acc.y;
      az [j] = grav.force[i].acc.z;
      pot[j] = grav.force[i].pot;
    }

#if 1
    fprintf(stderr, " >>> SPHgrav_lib took %g sec \n", dt);
#endif
  }

  void initDevice(const int device, const double theta)
  {
    grav.initDevice(device, theta);   /* opening angle, 0.3 or 0.4 should be safe */
  }
}

extern "C"
{
  void firsthalf_grav_forces_(
      int *n, int *n_lower, int *my_length, 
      double *px, double *py, double *pz,
      double *mass, double *range2, int *q_)
  {
    SPHgrav::firsthalf_grav_force(*n, *n_lower - 1, *my_length, px, py, pz, mass, range2);
  }

  void lasthalf_grav_forces_(
      int *n, 
      double *ax, double *ay, double *az,
      double *pot,int *myrank) 
  {
    SPHgrav::lasthalf_grav_force(ax, ay, az, pot);
  }

  void gpu_init_dev_(int *myrank, double *theta)
  {
    SPHgrav::initDevice(*myrank, *theta);
  }
}

