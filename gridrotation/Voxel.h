#ifndef VOXEL_H
#define VOXEL_H

#include <cassert>
#include <cstdlib>
#include <functional>

#define GLEW_STATIC 1

#include "BasicPoint.h"

class Voxel
{
    int _i , _j , _k;

public:
    Voxel() : _i(0) , _j(0) , _k(0){}
    Voxel( int x , int y , int z  ) : _i(x) , _j(y) , _k(z) {}

    inline int i() const {return _i;}
    inline int j() const {return _j;}
    inline int k() const {return _k;}

    inline void setI(int i) {_i = i;}
    inline void setJ(int j) {_j = j;}
    inline void setK(int k) {_k = k;}

    inline bool contains (int x) const { return (_i == x || _j == x || _k == x); }


    int operator * ( Voxel const & other ) const
    {
        return _i * other._i + _j * other._j + _k * other._k;
    }
    int operator [] (unsigned int c) const
    {
        if( c == 0 )
            return _i;
        if( c == 1 )
            return _j;
        if( c == 2 )
            return _k;

        assert( 0  &&  "Give a index between 0 and 2 as a coordinate for Voxel" );
        return -1;
    }

    int & operator [] (unsigned int c)
    {
        if( c == 0 )
            return _i;
        if( c == 1 )
            return _j;
        if( c == 2 )
            return _k;

        assert( 0  &&  "Give a index between 0 and 2 as a coordinate for Voxel" );
    }

    Voxel operator % ( Voxel const & other ) const
    {
        return Voxel( _j * other._k - _k * other._j , _k * other._i - _i * other._k , _i * other._j - _j * other._i );
    }

    Voxel operator + ( Voxel const & other ) const
    {
        return Voxel( _i + other._i , _j + other._j , _k + other._k );
    }
    Voxel operator - ( Voxel const & other ) const
    {
        return Voxel( _i - other._i , _j - other._j , _k - other._k );
    }
    void operator += ( Voxel const & other )
                     {
        _i += other._i;
        _j += other._j;
        _k += other._k;
    }
    void operator -= ( Voxel const & other )
                     {
        _i -= other._i;
        _j -= other._j;
        _k -= other._k;
    }
    void operator /= (float m)
                     {
        _i /= m;
        _j /= m;
        _k /= m;
    }
    void operator /= (double m)
                     {
        _i /= m;
        _j /= m;
        _k /= m;
    }

    int sqrnorm() const
    {
        return _i * _i + _j * _j + _k * _k;
    }

    float norm() const
    {
        return sqrt( sqrnorm() );
    }



    double normDoublePrecision()
    {
        return sqrt( (double)(_i) * (double)(_i) + (double)(_i) * (double)(_j) + (double)(_k) * (double)(_k) );
    }


    void operator << ( Voxel const & p )
    {
        _i = std::min( _i , p[0] );
        _j = std::min( _j , p[1] );
        _k = std::min( _k , p[2] );
    }
    void operator >> ( Voxel const & p )
    {
        _i = std::max( _i , p[0] );
        _j = std::max( _j , p[1] );
        _k = std::max( _k , p[2] );
    }

    static inline Voxel min ( Voxel const & p , Voxel const & p2 )
    {
        return Voxel( std::min( p2[0] , p[0] ),
                           std::min( p2[1] , p[1] ),
                           std::min( p2[2] , p[2] ) );
    }
    static inline Voxel max ( Voxel const & p , Voxel const & p2 )
    {
        return Voxel( std::max( p2[0] , p[0] ),
                           std::max( p2[1] , p[1] ),
                           std::max( p2[2] , p[2] ) );
    }
};

inline BasicPoint operator * (float m , Voxel const & p)
{
    return BasicPoint( m * p[0] , m * p[1] , m * p[2] );
}
inline BasicPoint operator * (double m , Voxel const & p)
{
    return BasicPoint( m * p[0] , m * p[1] , m * p[2] );
}
inline BasicPoint operator * (Voxel const & p , float m)
{
    return BasicPoint( m * p[0] , m * p[1] , m * p[2] );
}
inline BasicPoint operator * (Voxel const & p , double m)
{
    return BasicPoint( m * p[0] , m * p[1] , m * p[2] );
}

inline std::ostream & operator << (std::ostream & s , Voxel const & p)
{
    s << p[0] << " " << p[1] << " " << p[2] << " ";
//    if(p.boundary) s << " boundary";
    return s;
}

inline bool operator< (Voxel const & a, Voxel const & b) {
    return (a[0] < b[0] && a[1] < b[1] && a[2] < b[2]);
}

inline bool operator== (Voxel const & p1, Voxel const & p2) {
    return (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]);
}

inline bool operator>= (Voxel const & a, Voxel const & b) {
    return (a[0] >= b[0] || a[1] >= b[1] || a[2] >= b[2]);
}

inline bool operator!= (Voxel const & p1, Voxel const & p2) {
    return (p1[0] != p2[0] || p1[1] != p2[1] || p1[2] != p2[2]);
}

inline Voxel operator / (Voxel const & p , float m)
{
    return Voxel( p[0] / m , p[1] / m , p[2] / m );
}
inline Voxel operator / (Voxel const & p , double m)
{
    return Voxel( p[0] / m , p[1] / m , p[2] / m );
}

inline Voxel cross( Voxel const & v1 , Voxel const & v2 )
{
    return v1 % v2;
}

inline int dot( Voxel const & v1 , Voxel const & v2 )
{
    return v1 * v2;
}


namespace VoxelM
{
    inline Voxel random( int range )
    {
        return Voxel( range * (rand()) / ( RAND_MAX ) , range * (rand()) / ( RAND_MAX ) , range * (rand()) / ( RAND_MAX ) );
    }
}

struct compareVoxel {
    inline bool operator()(const Voxel p1, const Voxel p2) const {
        return p1 == std::min(p1, p2);
    }
};

typedef std::map<Voxel, int, compareVoxel> VoxelMapIndex;
#endif // Voxel_H


