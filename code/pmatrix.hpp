#ifndef _PMATRIX_
#define _PMATRIX_ 

#include "./pvector.hpp"
template <typename ntype, int NT>
class pmatrixq
{
  ntype m[NT][NT];
public:

  //destructor
  ~pmatrixq()
    {

    }
    //trivial constructor
  pmatrixq()
    {
      for (int i=0; i < NT; i++)
        for (int j=0; j < NT; j++)
          m[i][j]=0;
    }

    //constructor with initializer list
  pmatrixq(std::initializer_list<ntype> list)
    {
      int cc=0;
      for (auto el: list)
        {
          m[cc/NT][cc%NT] = el;
          cc++;  
        }
    } 

  // overloading of () operator to set (i,j)-th element
  ntype& operator()(int i, int j)
    {
      return m[i][j];
     }
  
  // overloading of () operator to get (i,j)-th element
  ntype operator()(int i, int j) const
    {
      return m[i][j];
     }

  // << overloading
  friend std::ostream& operator<<(std::ostream& os, const pmatrixq& m)
    {
      // write your code here
      int i, j;
      os << "{";
      for (i=0; i < NT; i++)
        {
          os << "{";
          for (j=0; j < NT; j++)
            {
              os << m(i,j);
              if (j < NT-1)
                os << ",";
            }
          os << "}";
          if (i < NT-1)
            os << ",";
        }
      os << "}";
      return os;
    }

    // overloading operator +
  pmatrixq operator+(const pmatrixq &m2) const
    {
      pmatrixq mt;

      for (int i=0; i < NT; i++)
          for (int j=0; j < NT; j++)
              mt(i,j) = m[i][j]+m2.m[i][j];
      return mt;
    }


    // overloading operator -
  pmatrixq operator-(const pmatrixq &m2) const
  {
    pmatrixq mt;

    for (int i=0; i < NT; i++)
      for (int j=0; j < NT; j++)
        mt(i,j) = m[i][j]-m2.m[i][j];
    return mt;
  }

      //overloading operator +=
   pmatrixq& operator+=(const pmatrixq &m2)
    {
      return (*this = *this + m2);
    }
  
  //overloading operator -=
  pmatrixq& operator-=(const pmatrixq &m2)
    {
      return (*this = *this - m2);
    }


    //overloading operator *= a scalar
  pmatrixq& operator*=(const ntype& s)
    {
      for(int i=0; i < NT; i++)
          for (int j=0; j < NT; j++)
              m[i][j] *= s;
      return (*this);
    } 


  // matrix times scalar
  pmatrixq operator*(const ntype& s) const
    {
      pmatrixq mt;
      for(int i=0; i < NT; i++)
          for (int j=0; j < NT; j++)
              mt(i,j) = m[i][j]*s;

      return mt;
    }   
 
  // scalar times matrix 
  friend pmatrixq operator*(const ntype& s, const pmatrixq &m2)
    {
      pmatrixq mt = m2*s;

      return mt;
    }


  // matrix divided by scalar
  pmatrixq operator/(const ntype& s) const
  {
    pmatrixq mt;
    for(int i=0; i < NT; i++)
      for (int j=0; j < NT; j++)
        mt(i,j) = m[i][j]/s;

    return mt;
  }

  //overloading operator /= a scalar
  pmatrixq& operator/=(const ntype& s)
  {
    for(int i=0; i < NT; i++)
      for (int j=0; j < NT; j++)
        m[i][j] /= s;
    return (*this);
  }

  // matrix times vector
  pvector<ntype,NT> operator*(const pvector<ntype,NT> &v2) const 
    {
      pvector<ntype,NT> vt;
      for (int i=0; i < NT; i++)
          for (int j=0; j < NT; j++)
              vt(i) += m[i][j]*v2(j);

      return vt;
    }
 
  // vector times matrix
  friend pvector<ntype,NT> operator*(const pvector<ntype,NT>& v1, const pmatrixq& m2)
    {
      pvector<ntype,NT> vt;
      for (int i=0; i < NT; i++)
          for (int j=0; j < NT; j++)
              vt(i) += v1(j)*m2.m[j][i];

      return vt;
    }
 
  // matrix times matrix
  pmatrixq operator*(const pmatrixq &m2) const
    {
      pmatrixq mt;

      for (int i=0; i < NT; i++)
          for (int j=0; j < NT; j++)
              for (int l=0; l < NT; l++)
                  mt(i,j) += m[i][l]*m2.m[l][j];

      return mt;
    }


  //overloading operator *= a matrix
  pmatrixq& operator*=(const pmatrixq &m2)
  {
    return (*this = *this * m2);
  }
  
  //transpose
  pmatrixq transpose(void)
    {
      pmatrixq mt;
      for (int i=0; i < NT; i++)
          for (int j=0; j < NT; j++)
              mt(i,j) = m[j][i];

      return mt;
    }

    //overloading operator ==
  bool operator==(const pmatrixq<ntype,NT>& m2) const
    {
      for (int i=0; i < NT; i++)
          for (int j=0; j < NT; j++)
              if (m2.m[i][j]!=m[i][j])
                  return 0;

      return 1;
    }
};

#if 0
// =============================================================
  template<typename ntype, int NT>
  std::ostream& operator<<(std::ostream& os, const pmatrixq<ntype,NT>& m)
    {
      // write your code here
      int i, j;
      os << "{";
      for (i=0; i < NT; i++)
        {
          os << "{";
          for (j=0; j < NT; j++)
            {
              os << m(i,j);
              if (j < NT-1)
                os << ",";
            }
          os << "}";
          if (i < NT-1)
            os << ",";
        }
      os << "}";
      return os;
    }
#endif
#endif
