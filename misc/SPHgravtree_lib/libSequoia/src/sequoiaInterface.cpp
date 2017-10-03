#include "../include/sequoiaInterface.h"
#include "../include/octree.h"


//The memory counters
long long my_dev::base_mem::currentMemUsage;
long long my_dev::base_mem::maxMemUsage;

//The Bonsai Octree class

octree *sequoia;

bool initFlag = false;
ofstream logFile;
  
  
//extern "C" {


my_dev::context & sequoia_init(char** argv, int device, const float _theta, const float eps)
{  
  sequoia = new octree(argv, device, _theta, eps);
   
  //Used for profiler
  char *gpu_prof_log;
  gpu_prof_log=getenv("CUDA_PROFILE_LOG");
  if(gpu_prof_log){
    char tmp[50];
    sprintf(tmp,"process%d_%s",0,gpu_prof_log);
    setenv("CUDA_PROFILE_LOG",tmp,1);
  }

  string logFileName    = "gpuLog.log";
  logFile.open(logFileName.c_str());
    
  sequoia->set_context(logFile, false); //Do logging to file and enable timing (false = enabled)
  
  //Load the device kernels
  sequoia->load_kernels();
 
  initFlag = true;
  
  return sequoia->getDevContext();
}


int sequoia_cleanup()
{
  assert(initFlag);
  
  
  delete sequoia;
  sequoia = NULL;
  
  
  
  initFlag = false;
  return 0;
}

int sequoia_sortBodies(my_dev::dev_mem<real4>  &bodies_pos, my_dev::dev_mem<uint> &permutation, int n_bodies)
{
  sequoia->sort_bodies(sequoia->localTree, bodies_pos, permutation, n_bodies);  
  return 0; 
}

int sequoia_reorderReal4(my_dev::dev_mem<real4>  &data, my_dev::dev_mem<uint> &permutation, int n_items)
{
  sequoia->reorder_dataR4(data, permutation, n_items);  
  return 0; 
}

int sequoia_reorderReal2(my_dev::dev_mem<real2>  &data, my_dev::dev_mem<uint> &permutation, int n_items)
{
  sequoia->reorder_dataR2(data, permutation, n_items);  
  return 0; 
}

int sequoia_reorderInt1(my_dev::dev_mem<int>  &data, my_dev::dev_mem<uint> &permutation, int n_items)
{
  sequoia->reorder_dataI1(data, permutation, n_items);  
  return 0; 
}

int sequoia_buildTreeStructure(my_dev::dev_mem<real4>  &bodies_pos, int n_bodies)
{
  sequoia->build(sequoia->localTree, bodies_pos, n_bodies);
  sequoia->compute_properties(sequoia->localTree, bodies_pos, n_bodies);
  return 0;
}

int sequoia_computeTreeProperties(my_dev::dev_mem<real4>  &bodies_pos, int n_bodies)
{
    sequoia->compute_properties(sequoia->localTree, bodies_pos, n_bodies);
    return 0;
}

int sequoia_createGroups(my_dev::dev_mem<real4> &bodies_pos, int n_bodies)
{
  sequoia->createGroups(sequoia->localTree, bodies_pos, n_bodies);
  return 0;
}

int sequoia_computeGravity(my_dev::dev_mem<real4> &j_bodies_pos,                           
                           my_dev::dev_mem<real4> &j_bodies_h,
                           my_dev::dev_mem<int>   &j_bodies_IDs, 
                           my_dev::dev_mem<real4> &i_bodies_pos, 
                           my_dev::dev_mem<real4> &i_bodies_h,
                           my_dev::dev_mem<int>   &i_bodies_IDs, 
                           my_dev::dev_mem<real4> &i_bodies_acc, 
                           my_dev::dev_mem<real>  &i_bodies_ds2, 
                           my_dev::dev_mem<int>   &i_bodies_ngb, 
                           int                     n_i_bodies)
{
  
  sequoia->approximate_gravity(sequoia->localTree, 
                               j_bodies_pos, j_bodies_h, j_bodies_IDs,
                               i_bodies_pos, i_bodies_h, i_bodies_IDs, n_i_bodies, 
                               i_bodies_acc, i_bodies_ds2, i_bodies_ngb);
  
  return 0;
}

int sequoia_setParticlesAndGetGravity( my_dev::dev_mem<real4> &j_bodies_pos,      //Positions J-particles
                                       my_dev::dev_mem<real4> &j_bodies_h,
                                       my_dev::dev_mem<int>   &j_bodies_IDs,      //Particle IDs J-particles
                                       int                    n_j_bodies,         //Number of J-particles
                                       my_dev::dev_mem<real4> &i_bodies_pos,      //Positions I-particles 
                                       my_dev::dev_mem<real4> &i_bodies_h,
                                       my_dev::dev_mem<int>   &i_bodies_IDs,      //Particle IDs J-particles
                                       int                    n_i_bodies,         //Number of I-particles
                                       bool                   sortJBodies,        //Do we need to sort J-particles?
                                       bool                   sortIBodies,        //Do we need to sort I-particles?           
                                       my_dev::dev_mem<real4> &i_bodies_acc,      //OUT  Accelerations for I-particles
                                       my_dev::dev_mem<real>  &i_bodies_ds2,      //OUT  min distance squared for I-particles
                                       my_dev::dev_mem<int>   &i_bodies_ngb)      //OUT  J-ID of the nearest neighbour for I-particles
{
  assert(initFlag);
  
  //If required sort J particles
  if(sortJBodies)
  {
    my_dev::dev_mem<uint> permutation(sequoia->getDevContext(), n_j_bodies);
    
    sequoia_sortBodies  (j_bodies_pos, permutation, n_j_bodies);
    sequoia_reorderReal4(j_bodies_pos, permutation, n_j_bodies);
    sequoia_reorderReal4(j_bodies_h,   permutation, n_j_bodies);
    sequoia_reorderInt1 (j_bodies_IDs, permutation, n_j_bodies);
  }
  //If required sort I particles
  if(sortIBodies)
  {
    my_dev::dev_mem<uint> permutation(sequoia->getDevContext(), n_i_bodies);
    
    sequoia_sortBodies  (i_bodies_pos, permutation, n_i_bodies);
    sequoia_reorderReal4(i_bodies_pos, permutation, n_i_bodies);
    sequoia_reorderReal4(i_bodies_h,   permutation, n_i_bodies);
    sequoia_reorderInt1 (i_bodies_IDs, permutation, n_i_bodies);
  }
  
  //Build the tree-structure
  sequoia_buildTreeStructure(j_bodies_pos, n_j_bodies);
  //Create the groups that walk the tree
  sequoia_createGroups(i_bodies_pos, n_i_bodies);
   
  //Finally compute the gravity and get the nearest neighbour + distance 
  sequoia_computeGravity(j_bodies_pos, j_bodies_h, j_bodies_IDs, 
                         i_bodies_pos, i_bodies_h, i_bodies_IDs,
                         i_bodies_acc, i_bodies_ds2, i_bodies_ngb, 
                         n_i_bodies);  
  
  
  return 0;
}


int sequoia_get_ngb(my_dev::dev_mem<real4> &j_bodies_pos,                           
                    my_dev::dev_mem<int>   &j_bodies_IDs, 
                    my_dev::dev_mem<real4> &i_bodies_pos, 
                    my_dev::dev_mem<int>   &i_bodies_IDs, 
                    my_dev::dev_mem<int>   &i_bodies_ngb)
{
  assert(initFlag);
    
  sequoia->get_ngb(sequoia->localTree, 
                   j_bodies_pos, j_bodies_IDs,
                   i_bodies_pos, i_bodies_IDs,
                   i_bodies_ngb);
  
  return 0;
}




  
//} //end extern

