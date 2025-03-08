#ifndef _PVECTOR_
#define _PVECTOR_

#include<initializer_list> // for vector initialization
#include<iostream> // input/output
#include<cmath>
#include<string> // strings
#include "./randnumgen.hpp"

inline randnumgen rng;

//pvector class
template <typename ntype, int NT>
class pvector
{
  ntype v[NT]; // private member
public:

  pvector() // void constructor
    {
      for(int i = 0; i < NT; i++)
        {
          v[i]=0;
        }
    }

  // overloading of constructor for handling curly brace initialization
  pvector(std::initializer_list<ntype> list)
    {
      int c=0;
      for (ntype el: list)
        {
          if (c < NT)
            {
              v[c] = el;
            }
          c++;
        }
      for (;c < NT; c++) // if list length is less than NT,
                        // following elements initialized to 0
        {
          v[c]=0.0;
        }
    }


  ~pvector() // destructor
    = default;

    //show method
    void show(const std::string& s="") const
    {
      std::cout << s << "(";
      for (int i=0; i < NT; i++)
        {
          std::cout << v[i];
          if (i < NT-1) 
            std::cout << ",";
        }
      std::cout << ")\n";
    }

  //overloading = operator
  pvector& operator=(const pvector& v2)
    {
      for (int i=0; i < NT; i++)
        {

          (*this).v[i] = v2.v[i];
          // equivalently:
          // (*this).v[i] = v2.v[i];
        }
      return(*this);
    }


  // overloading operator+
  pvector operator+(const pvector& v2) const
    {
      pvector vs;
      for (int i=0; i < NT; i++)
        {
          vs.v[i] = v[i] + v2.v[i];
        }
      return vs;
    } 

  // overloading operator -
  pvector operator-(const pvector& v2) const
    {
      pvector vs;
      for (int i=0; i < NT; i++)
        {
          vs.v[i] = v[i] - v2.v[i];
        }
      return vs;
    } 


  // sum based on overloaded + operator
  pvector sum(const pvector& v2) const
    {
      return (*this)+v2;
    }


  // get
  ntype get(int i) const
    {
      return v[i];
    }

  // set
  ntype set(int i, ntype val) 
    {
      return v[i]=val;
    }

  //overloading operator +=
  pvector& operator+=(const pvector& v2)
    {
      for (int i=0; i < NT; i++)
        {
          v[i] += v2.v[i];
        }
      return (*this);
    } 

  //overloading  operator -=
  pvector& operator-=(const pvector& v2)
    {
      for (int i=0; i < NT; i++)
        {
          v[i] -= v2.v[i];
        }
      return (*this);
    } 

  //overloading operator () (not modifying the object)
  ntype operator()(int idx) const
    {
      return v[idx]; 
    }

  //overloading operator () (modifying the object)
  ntype& operator()(int idx)
    {
      return v[idx]; 
    }

  // scalar product
  ntype operator*(const pvector& vec) const
    {
      ntype sp=0;
      for (int i=0; i < NT; i++)
        sp += v[i]*vec.v[i];
      return sp;
    }

  //norm of the vector
  ntype norm() const
    {
      return sqrt((*this)*(*this));
    }
    
  // vector times scalar 
  pvector operator*(ntype s) const
    {
      pvector vt;
      for (int i=0; i < NT; i++)
        vt.v[i] = v[i]*s;
      return vt;
    }

  // vector divided by scalar 
  pvector operator/(ntype s) const
    {
      pvector vt;
      for (int i=0; i < NT; i++)
        vt.v[i] = v[i]/s;
      return vt;
    }

  // multiply by scalar and assign result to vector
  pvector& operator *=(ntype s)
    {
      for (int i=0; i < NT; i++)
        v[i] *= s;
      return (*this);
    }

  // divide by scalar and assign result to vector
  pvector& operator /=(ntype s)
    {
      for (int i=0; i < NT; i++)
        v[i] /= s;
      return (*this);
    }


  //divide by vector (with all components different from zero) and assign result to vector
  pvector& operator /(pvector v2)
  {
    for (int i=0; i < NT; i++)
    {
      v[i] /= v2.v[i];
    }
    return (*this);
  }

  //overloading of == operator
  bool operator==(pvector vec)
    {
      for (int i=0; i < NT; i++)
        {
          if (v[i] != vec.v[i])
            return false;
        }
      return true;
    }

  //vector product
  pvector operator^(const pvector& vec) const
    {
      if (NT==3)
        {
          pvector vt;
          vt.v[0] = v[1]*vec.v[2]-v[2]*vec.v[1];
          vt.v[1] = v[2]*vec.v[0]-v[0]*vec.v[2];
          vt.v[2] = v[0]*vec.v[1]-v[1]*vec.v[0];
          return vt;
        }
      else
        {
          std::cout << "Cross product not defined\n";
          exit(1);
        }
    }

  // scalar time vector
  friend pvector operator*(ntype s, const pvector& vec)
    {
      pvector vt;
      for (int i=0; i < NT; i++)
        {
          vt.v[i] = s*vec.v[i];
        }
      return vt;
    }

  // rint() method
  pvector rint()
    {
      pvector vt;
      for (auto i=0; i < NT; i++)
        {
          vt.v[i] = rint(v[i]);
        }
      return vt;
    }

  //overloading of <<
  friend std::ostream& operator<<(std::ostream& os, const pvector& vec)
    {
      os << "(";
      for (int i=0; i < NT; i++)
        {
          os << vec.v[i];
          if (i < NT-1)
           os << ","; 
        }
      os << ")";
      return os;
    }

  // component-wise multiplication of two vectors
  pvector mulcw(const pvector& vec) const
    {
      pvector vt;
      for (int i=0; i < NT; i++)
        {
          vt(i) = v[i]*vec(i);
        }
      return vt;
    }
 
  // component-wise division of two vectors
  pvector divcw(const pvector& vec) const
    {
      pvector vt;
      for (int i=0; i < NT; i++)
        {
          vt(i) = v[i]/vec(i);
        }
      return vt;
    }

  //assign to the vector random position in the box
  pvector& random(const ntype& L)
    {
      for (int i=0; i < NT; i++)
        {
          v[i] = (rng.ranf()-0.5)*L; // assign a random value in [-L/2,L/2]
        }
      return (*this);
    }

  //assign to the vector a random orientation on the unit sphere
  void random_orient()
    {
      ntype rS, S, V1, V2;
      if (NT==3)
        {
          do
            {
              V1 = 2.0*rng.ranf()-1.0;
              V2 = 2.0*rng.ranf()-1.0;
              S = V1*V1+V2*V2;
            }
          while (S >= 1.0);
          rS = sqrt(1.0-S);
#if 0
          v[0] = 2.0*rS*V1;
          v[1] = 2.0*rS*V2;
          v[2] = 1.0-2.0*S;
#else
          (*this) = {2.0*rS*V1, 2.0*rS*V2, 1.0-2.0*S};
#endif

        }
      else
        {
          std::cout << "[random_orient] Only 3D vectors are supported\n";
        }
    }


};

// we make rint() generic function (function template)
template<typename ntype, int NT>
pvector<ntype,NT> rint(const pvector<ntype,NT>& vec)
{
  pvector<ntype,NT> vt;

  for (int i=0; i < NT; i++)
    {
      vt(i) = rint(vec(i));
    }
  return vt;
}
#if 0
template<typename ntype, int NT>
std::ostream& operator<<(std::ostream& os, const pvector<ntype,NT>& vec)
    {
      os << "(";
      for (int i=0; i < NT; i++)
        {
          os << vec(i);
          if (i < NT-1)
           os << ","; 
        }
      os << ")";
      return os;
    }
 
#endif
#if 0
// non friend implementation of scalar time vector
// scalar time vector 
template<typename ntype, int NT>
pvector<ntype,NT> operator*(ntype s, const pvector<ntype, NT>& vec)
{
  pvector<ntype,NT> vt;
  for (int i=0; i < NT; i++)
    {
      vt(i) = s*vec(i);
    }
  return vt;
}
#endif

// predefined 3d vector of doubles
using pvec3d=pvector<double,3>;
#endif
