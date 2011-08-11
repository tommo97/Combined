/*
This file is part of the Combined Wake Modelling Code Version 1.0

VTM Code Copyright Tom McCombes 2009
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev:: 2                $:  Revision of last commit
$Author:: tom           $:  Author of last commit
$Date:: 2009-10-28 20:1#$:  Date of last commit

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef TEST_MODE
#include "types.hpp"
#endif

#ifndef VECT34_HPP
#define VECT34_HPP
/**************************************************************/
class Vect3 
{
    public:
	REAL x, y, z;
	Vect3(REAL ix, REAL iy, REAL iz) : x(ix), y(iy), z(iz) {}; 
	Vect3() : x(0.0), y(0.0), z(0.0) {};
	Vect3(REAL A) : x(A), y(A), z(A) {};
	Vect3(REAL A[3]) : x(A[0]), y(A[1]), z(A[2]) {};
    Vect3 (const Vect3&A) : x(A.x), y(A.y), z(A.z) {};
	Vect3 Cross(const Vect3 B) {return Vect3(y*B.z - z*B.y, z*B.x - x*B.z, x*B.y - y*B.x);}
	REAL Dot(const Vect3 B) {return x*B.x + y*B.y + z*B.z;}
	REAL Mag() {return sqrt(x*x + y*y + z*z);}
	//	Addition operators
	inline Vect3 operator + (const Vect3 B) {return Vect3(x+B.x, y+B.y, z+B.z);}
	inline friend Vect3 operator + (const Vect3 A, const REAL B ) {return Vect3(A.x+B,A.y+B,A.z+B);}
	inline friend Vect3 operator + (const REAL A, const Vect3 B ) {return Vect3(B.x+A,B.y+A,B.z+A);}
	inline void operator += (const Vect3 B) {x += B.x, y += B.y, z += B.z;}
	inline friend void operator += (Vect3 A, const REAL B ) {A.x+=B, A.y+=B, A.z+=B;}
	inline friend void operator += (const REAL A, Vect3 B ) {B.x+=A, B.y+=A, B.z+=A;}
	//	Subtraction operators
        inline Vect3 operator - () {return Vect3(-x, -y, -z);}
	inline Vect3 operator - (const Vect3 B) {return Vect3(x-B.x,y-B.y,z-B.z);}
	inline friend Vect3 operator - (const Vect3 A, const REAL B ) {return Vect3(A.x-B,A.y-B,A.z-B);}
	inline friend Vect3 operator - (const REAL A, const Vect3 B ) {return Vect3(B.x-A,B.y-A,B.z-A);}
	inline void operator -= (const Vect3 B) {x -= B.x, y -= B.y, z -= B.z;}
	inline friend void operator -= (Vect3 A, const REAL B ) {A.x-=B, A.y-=B, A.z-=B;}
	inline friend void operator -= (const REAL A, Vect3 B ) {B.x-=A, B.y-=A, B.z-=A;}
	//	Multiplication operators
	inline Vect3 operator * (const Vect3 B) {return Vect3(x*B.x,y*B.y,z*B.z);}
	inline friend Vect3 operator * (const Vect3 A, const REAL B ) {return Vect3(A.x*B,A.y*B,A.z*B);}
	inline friend Vect3 operator * (const REAL A, const Vect3 B ) {return Vect3(B.x*A,B.y*A,B.z*A);}
	inline void operator *= (const Vect3 B) {x *= B.x, y *= B.y, z *= B.z;}
	inline friend void operator *= (Vect3 A, const REAL B ) {A.x*=B, A.y*=B, A.z*=B;}
	inline friend void operator *= (const REAL A, Vect3 B ) {B.x*=A, B.y*=A, B.z*=A;}
	//	Division operators
	inline Vect3 operator / (const Vect3 B) {return Vect3(x/B.x,y/B.y,z/B.z);}
	inline friend Vect3 operator / (const Vect3 A, const REAL B ) {return Vect3(A.x/B,A.y/B,A.z/B);}
	inline friend Vect3 operator / (const REAL A, const Vect3 B ) {return Vect3(B.x/A,B.y/A,B.z/A);}
	inline void operator /= (const Vect3 B) {x /= B.x, y /= B.y, z /= B.z;}
	//inline friend void operator /= (Vect3 A, const REAL B ) {A.x/=B, A.y/=B, A.z/=B;}
	//inline friend void operator /= (const REAL A, Vect3 B ) {B.x/=A, B.y/=A, B.z/=A;}

	//	Relational operators
        inline friend REAL max (const Vect3 A) {return max(max(A.x,A.y),A.z);}
        inline friend Vect3 max(const Vect3 A, const Vect3 B) {return Vect3(max(A.x,B.x), max(A.y,B.y), max(A.z,B.z));}
        inline friend Vect3 min(const Vect3 A, const Vect3 B) {return Vect3(min(A.x,B.x), min(A.y,B.y), min(A.z,B.z));}
	inline Vect3 operator >  (const REAL B) 	{return Vect3((REAL) x>B, 	(REAL) y>B, 	(REAL) z>B);}
	inline Vect3 operator <  (const REAL B) 	{return Vect3((REAL) x<B, 	(REAL) y<B, 	(REAL) z<B);}
	inline Vect3 operator >= (const REAL B) 	{return Vect3((REAL) x>=B, 	(REAL) y>=B, 	(REAL) z>=B);}
	inline Vect3 operator <= (const REAL B) 	{return Vect3((REAL) x<=B, 	(REAL) y<=B, 	(REAL) z<=B);}
	inline Vect3 operator >  (const Vect3 B) 	{return Vect3((REAL) x>B.x, 	(REAL) y>B.x, 	(REAL) z>B.x);}
	inline Vect3 operator <  (const Vect3 B) 	{return Vect3((REAL) x<B.x, 	(REAL) y<B.x, 	(REAL) z<B.x);}
	inline Vect3 operator >= (const Vect3 B) 	{return Vect3((REAL) x>=B.x, 	(REAL) y>=B.x, 	(REAL) z>=B.x);}
	inline Vect3 operator <= (const Vect3 B) 	{return Vect3((REAL) x<=B.x, 	(REAL) y<=B.x, 	(REAL) z<=B.x);}
	inline bool operator ==  (const Vect3 B) 	{return ((x==B.x) && 	(y==B.y) && 	(z==B.z));}
	inline bool operator !=  (const Vect3 B) 	{return ((x!=B.x) || 	(y!=B.y) || 	(z!=B.z));}
	/*    void operator = (REAL ix, REAL iy, REAL iz) {x = ix; y = iy; z = iz;}*/
	//	Floating point type operators
	inline friend Vect3 fabs(const Vect3 B) {return Vect3(fabs(B.x), fabs(B.y), fabs(B.z));}
	//	An assignemnt operator
	inline void operator = (const REAL B) {x = y = z = B;}
#ifndef SIGN
#define SIGN(a) ((a) < 0.0 ? -1 : 1)
#endif
	inline friend Vect3 sign(const Vect3 B) {return Vect3(SIGN(B.x), SIGN(B.y), SIGN(B.z));}
	//	Some other operators
	inline friend ostream& operator << (ostream& os, const Vect3 in) {return os << in.x << " " << in.y << " " << in.z;}
	inline friend Vect3 floor(const Vect3 B) {return Vect3(floor(B.x), floor(B.y), floor(B.z));}
        inline friend Vect3 ceil(const Vect3 B) {return Vect3(ceil(B.x), ceil(B.y), ceil(B.z));}
	inline friend Vect3 sqrt(const Vect3 B) {return Vect3(sqrt(B.x), sqrt(B.y), sqrt(B.z));}
	inline friend REAL pow(const Vect3 B, int a, int b, int c) {return pow(B.x,a)*pow(B.y,b)*pow(B.z,c);}
	inline void assign(REAL a, REAL b, REAL c) {x = a, y = b, z = c;} 
	inline void assign(REAL a[]) {x = a[0], y = a[1], z = a[2];} 
	inline REAL& operator [] (const int index) throw (const char *)
	{
	    if (index==0) return x;
	    else if (index==1) return y;
	    else if (index==2) return z;
	    else throw "Vect3 Bounds Violation";
	}
        ~Vect3() {};
};

inline Vect3 VectMultMatrix(Vect3 M[], Vect3 Vin)
{
    return Vect3(M[0].Dot(Vin), M[1].Dot(Vin), M[2].Dot(Vin));
}
inline Vect3 VectMultMatrixTranspose(Vect3 M[], Vect3 Vin)
{
    return Vect3(M[0].x * Vin[0] + M[1].x * Vin[1] + M[2].x * Vin[2],
		 M[0].y * Vin[0] + M[1].y * Vin[1] + M[2].y * Vin[2],
		 M[0].z * Vin[0] + M[1].z * Vin[1] + M[2].z * Vin[2]);
}
inline void VectMultMatrix(Vect3 M[3], Vect3 Vin, Vect3 &Vout) 
{
Vout.x = M[0].Dot(Vin);
Vout.y = M[1].Dot(Vin);
Vout.z = M[2].Dot(Vin);
}

inline void VectMultMatrixPlus(Vect3 M[3], Vect3 Vin, Vect3 &Vout)
{
Vout.x += M[0].Dot(Vin);
Vout.y += M[1].Dot(Vin);
Vout.z += M[2].Dot(Vin);
}

inline void VectMultMatrixTransposePlus(Vect3 M[3], Vect3 Vin, Vect3 &Vout) 
{
Vout.x += M[0].x * Vin.x + M[1].x * Vin.y + M[2].x * Vin.z;
Vout.y += M[0].y * Vin.x + M[1].y * Vin.y + M[2].y * Vin.z;
Vout.z += M[0].z * Vin.x + M[1].z * Vin.y + M[2].z * Vin.z;
}

inline void VectMultMatrixTranspose(Vect3 M[3], Vect3 Vin, Vect3 &Vout)
{
Vout.x = M[0].x * Vin.x + M[1].x * Vin.y + M[2].x * Vin.z;
Vout.y = M[0].y * Vin.x + M[1].y * Vin.y + M[2].y * Vin.z;
Vout.z = M[0].z * Vin.x + M[1].z * Vin.y + M[2].z * Vin.z;
}
/**************************************************************/
class Vect4 
{
    public:
    REAL a, b, c, d;
    Vect4(REAL ix, REAL iy, REAL iz, REAL it) : a(ix), b(iy), c(iz), d(it) {}; 
    Vect4() : a(0.), b(0.), c(0.), d(0.) {};
    //	Addition operators
    inline Vect4 operator + (const Vect4 B) {return Vect4(a+B.a, b+B.b, c+B.c, d+B.d);}
    inline friend Vect4 operator + (const Vect4 A, const REAL B ) {return Vect4(A.a+B,A.b+B,A.c+B,A.d+B);}
    inline friend Vect4 operator + (const REAL A, const Vect4 B ) {return Vect4(B.a+A,B.b+A,B.c+A,B.d+A);}
    inline void operator += (const Vect4 B) {a += B.a, b += B.b, c += B.c, d += B.d;}
    inline friend void operator += (Vect4 A, const REAL B ) {A.a+=B, A.b+=B, A.c+=B, A.d+=B;}
    inline friend void operator += (const REAL A, Vect4 B ) {B.a+=A, B.b+=A, B.c+=A, B.d+=A;}
    //	Subtraction operators	
    inline Vect4 operator - (const Vect4 B) {return Vect4(a-B.a, b-B.b, c-B.c, d-B.d);}
    inline friend Vect4 operator - (const Vect4 A, const REAL B ) {return Vect4(A.a-B,A.b-B,A.c-B,A.d-B);}
    inline friend Vect4 operator - (const REAL A, const Vect4 B ) {return Vect4(B.a-A,B.b-A,B.c-A,B.d-A);}
    inline void operator -= (const Vect4 B) {a -= B.a, b -= B.b, c -= B.c, d -= B.d;}
    inline friend void operator -= (Vect4 A, const REAL B ) {A.a-=B, A.b-=B, A.c-=B, A.d-=B;}
    inline friend void operator -= (const REAL A, Vect4 B ) {B.a-=A, B.b-=A, B.c-=A, B.d-=A;}
    //	Multiplication operators		
    inline Vect4 operator * (const Vect4 B) {return Vect4(a*B.a, b*B.b, c*B.c, d*B.d);}
    inline friend Vect4 operator * (const Vect4 A, const REAL B ) {return Vect4(A.a*B,A.b*B,A.c*B,A.d*B);}
    inline friend Vect4 operator * (const REAL A, const Vect4 B ) {return Vect4(B.a*A,B.b*A,B.c*A,B.d*A);}
    inline void operator *= (const Vect4 B) {a *= B.a, b *= B.b, c *= B.c, d *= B.d;}
    inline friend void operator *= (Vect4 A, const REAL B ) {A.a*=B, A.b*=B, A.c*=B, A.d*=B;}
    inline friend void operator *= (const REAL A, Vect4 B ) {B.a*=A, B.b*=A, B.c*=A, B.d*=A;}
    //	Division operators		
    inline Vect4 operator / (const Vect4 B) {return Vect4(a/B.a, b/B.b, c/B.c, d/B.d);}
    inline friend Vect4 operator / (const Vect4 A, const REAL B ) {return Vect4(A.a/B,A.b/B,A.c/B,A.d/B);}
    inline friend Vect4 operator / (const REAL A, const Vect4 B ) {return Vect4(B.a/A,B.b/A,B.c/A,B.d/A);}
    inline void operator /= (const Vect4 B) {a /= B.a, b /= B.b, c /= B.c, d /= B.d;}
    inline friend void operator /= (Vect4 A, const REAL B ) {A.a/=B, A.b/=B, A.c/=B, A.d/=B;}
    inline friend void operator /= (const REAL A, Vect4 B ) {B.a/=A, B.b/=A, B.c/=A, B.d/=A;}
    inline friend Vect4 atan2(const Vect4 A, const Vect4 B) {return Vect4(atan2(A.a,B.a), atan2(A.b,B.b), atan2(A.c,B.c), atan2(A.d,B.d));}
    inline friend Vect4 log(const Vect4 A) {return Vect4(log(A.a), log(A.b), log(A.c), log(A.d));}
    inline friend Vect4 permute_1(const Vect4 A) {return Vect4(A.b, A.c, A.d, A.a);}
    inline friend REAL sum(const Vect4 A) {return A.a+A.b+A.c+A.d;}
//	Some relational operators	
    bool operator == (const Vect4 B) {return ((a==B.a) && (b==B.b) && (c==B.c) && (d==B.d));}
    bool operator != (const Vect4 B) {return ((a!=B.a) || (b!=B.b) || (c!=B.c) || (d!=B.d));}
    friend ostream& operator << (ostream& os, const Vect4 in) {return os << in.a << " " << in.b << " " << in.c << " "  << in.d << endl;}
    inline void assign(REAL i, REAL j, REAL k, REAL l) {a = i, b = j, c = k, d = l;} 
    inline void assign(REAL A[]) {a = A[0], b = A[1], c = A[2], d = A[3];} 
    inline friend Vect4 sqrt(const Vect4 B) {return Vect4(sqrt(B.a), sqrt(B.b), sqrt(B.c), sqrt(B.d));}
    REAL& operator [] (const int index) throw (const char *)
    {
	if (index==0) return a;
	else if (index==1) return b;
	else if (index==2) return c;
	else if (index==3) return d;
	else throw "Vect4 Bounds Violation";
    }
    ~Vect4() {};
};


#endif
