#ifndef _CUVEC3_H_
#define _CUVEC3_H_

template<class REAL>
struct cuVec3
{
  REAL x, y, z;
  __host__ __device__ cuVec3() {}
  __host__ __device__ cuVec3(const REAL v) : x(v), y(v), z(v) {}
  __host__ __device__ cuVec3(const REAL _x, const REAL _y, const REAL _z) : x(_x), y(_y), z(_z) {}
  
  __host__ __device__ cuVec3 operator=(const cuVec3<float> v) {x = v.x; y = v.y; z = v.z; return *this;};
  __host__ __device__ cuVec3 operator=(const cuVec3<double > v) {x = v.x; y = v.y; z = v.z; return *this;};
  

  __host__ __device__ REAL   operator*(const cuVec3<REAL> v) const {return       (x*v.x + y*v.y + z*v.z);}
  __host__ __device__ cuVec3 operator*(const        REAL  v) const {return cuVec3(x*v,   y*v,   z*v);}
//  __host__ __device__ cuVec3 operator+(const cuVec3<REAL> v) const {return cuVec3(x+v.x, y+v.y, z+v.z);}
  __host__ __device__ cuVec3 operator-(const cuVec3<REAL> v) const {return cuVec3(x-v.x, y-v.y, z-v.z);}
  __host__ __device__ cuVec3 operator%(const cuVec3<REAL> v) const {return cuVec3(x*v.y - y*v.x, y*v.z-z*v.y, z*v.x - x*v.z);}
  __host__ __device__ cuVec3 operator-() const {return cuVec3(-x, -y, -z);}
  
  __host__ __device__ cuVec3 operator+(const cuVec3<float> v) const {return cuVec3(x+v.x, y+v.y, z+v.z);}
  __host__ __device__ cuVec3 operator+(const cuVec3<double > v) const {return cuVec3(x+v.x, y+v.y, z+v.z);}


	__host__ __device__ cuVec3 operator += (const cuVec3<REAL> v)
  {
		*this = *this + v;
		return *this;
	}

	__host__ __device__ cuVec3 operator -= (const cuVec3<REAL> v)
  {
		*this = *this - v;
		return *this;
	}
	__host__ __device__ cuVec3 operator *= (const REAL s)
  {
		*this = *this * s;
		return *this;
	}
  __host__ __device__ friend cuVec3 operator * (const REAL s ,const cuVec3<REAL> v)
  {
    return v*s;
  }


  __host__ __device__ REAL norm2() const {return (*this)*(*this);};
};

#endif /* _CUVEC3_H_ */
