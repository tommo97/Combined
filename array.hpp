/*
This file is part of the Combined Wake Modelling Code Version 1.0

VTM Code Copyright Tom McCombes 2009
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velociy vorticity form


$Rev::                  $:  Revision of last commit
$Author::               $:  Author of last commit
$Date::                 $:  Date of last commit

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



#ifndef ARRAY_INCL
#define ARRAY_INCL
#ifndef TEST_MODE
#include "includes.hpp"
#include "types.hpp"
#endif
#ifdef USE_ARRAY

template <class T>
class Array {
public:

    Array() : length(0), data_(NULL) {
    };

    Array(int n) : length(n), data_(new T[length])
    {
#ifndef ARRAY_NO_CHECK
        if (!data_) throw MEMFAIL;
#endif
    };

    Array(int n, const T &t) : length(n), data_(new T[length])
    {
#ifndef ARRAY_NO_CHECK
        if (!data_) throw MEMFAIL;
#endif
        for (unsigned int i = 0; i < length; ++i) data_[i] = t;
    };

    Array(const Array <T> &a) : length(a.length), data_(new T[length])
    {
#ifndef ARRAY_NO_CHECK
        if (!data_) throw MEMFAIL;
#endif
//        memcpy (data_,a.data_, length * sizeof(T));
        for (unsigned int i = 0; i < length; ++i) data_[i] = a.data_[i];
    };
    
    ~Array();
    Array & operator=(const Array&);
    Array & operator=(const T&);
    void clear();
    Array& push_back(const T&);
    Array& push_front(const T&);
    Array& pop_front();
    Array& pop_back();

    int size() {
        return length;
    }
#ifndef TEST_MODE
    void print() {
        for (unsigned int i = 0; i < length; ++i) if (WRITE_TO_SCREEN) std::cout << data_[i] << " ";
        if (WRITE_TO_SCREEN) std::cout << std::endl;
    };
#endif
    Array& allocate(int i);
    Array& assign(int, const T&);

    class BadSize {
    };
    //      Access methods to get the (i,j,k,l) element:
    T & operator() (int ind);
    T & front();
    T & back();
    const T & operator() (int ind) const;
#ifndef ARRAY_NO_CHECK

    T & operator[] (unsigned int i) {
        if (i >= length) {
            printf("Array Bounds Violation: Attempted to access element %d, array length %d - ", i, length);

            throw BoundsViolation();
        }
        return data_[i];
    }
#else

    T & operator[] (int i) {
        return data_[i];
    }
#endif
    // Exception enum

    enum exception {
        MEMFAIL
    };
    // These throw a BoundsViolation object if i or j is too big

    class BoundsViolation {
    };

    class WrongSizeAssign {
    };
private:
    unsigned int length;
    T* data_;
};

/************************************************/

/*  Destructor  */
template <class T> Array<T>::~Array() {
    clear();
}
/************************************************/

/*  Exact Copy  */
template <class T> Array<T>& Array<T>::operator =(const Array &a) {

    if (this != &a) // in case somebody tries assign array to itself
    {
        if (data_) delete[] data_;
        length = a.length;

        data_ = new T[length];

        for (unsigned int i = 0; i < length; ++i) data_[i] = a.data_[i];
    }
    return *this;
}
/************************************************/

/*  Set to values  */
template <class T> Array<T>& Array<T>::operator =(const T &a) {
    for (unsigned int i = 0; i < length; i++) data_[i] = a;
    return *this;
}

/************************************************/
template <class T> void Array<T>::clear() // clear array memory
{
    if (data_) delete[] data_;
    data_ = NULL;
    length = 0;
}
/************************************************/
template <class T> T & Array<T>::front(){
        if (data_)
            return data_[0];
        else
            throw BoundsViolation();
}
/************************************************/
template <class T> T & Array<T>::back(){
        if (data_)
            return data_[length-1];
        else
            throw BoundsViolation();
}
/************************************************/
template <class T> Array<T>& Array<T>::push_back(const T &insrt) {
    T* data_new = new T [length + 1];

    if (data_) {
        for (unsigned int i = 0; i < length; ++i) data_new[i] = data_[i];
        delete[] data_;
    }
    data_new[length] = insrt;
    data_ = data_new;
    length += 1;
    return *this;
}
/************************************************/
template <class T> Array<T>& Array<T>::push_front(const T &insrt) {
    T* data_new = new T [length + 1];

    if (data_) {
        for (unsigned int i = 0; i < length; ++i) data_new[i+1] = data_[i];
        delete[] data_;
    }
    data_new[0] = insrt;
    data_ = data_new;
    length += 1;
    return *this;
}
/************************************************/
template <class T> Array<T>& Array<T>::pop_back() {

#ifndef ARRAY_NO_CHECK
    if (!data_) throw BadSize();
#endif
        if (length == 1) {
        clear();
        return *this;
    }
    T* data_new = new T [length - 1];
    if (data_) {
        for (unsigned int i = 0; i < length - 1; ++i) data_new[i] = data_[i];
        delete[] data_;
    }

    data_ = data_new;
    length -= 1;
    return *this;
}

/************************************************/
template <class T> Array<T>& Array<T>::pop_front() {
#ifndef ARRAY_NO_CHECK
    if (!data_) throw BadSize();
#endif
    if (length == 1) {
        clear();
        return *this;
    }
    T* data_new = new T [length - 1];

    if (data_) {
        for (unsigned int i = 1; i < length; ++i) data_new[i - 1] = data_[i];
        delete[] data_;
    }

    data_ = data_new;
    length -= 1;
    return *this;
}

/************************************************/
template <class T> Array<T>& Array<T>::allocate(int l) {
    if (data_) delete[] data_;
#ifndef ARRAY_NO_CHECK
    if (l == 0) throw BadSize();
#endif
    data_ = new T[l];
    length = l;
#ifndef ARRAY_NO_CHECK
    if (!data_) throw MEMFAIL;
#endif
    return *this;
}

/************************************************/
template <class T> Array<T>& Array<T>::assign(int l, const T &asignee) {
    clear();

#ifndef ARRAY_NO_CHECK
    if (l == 0) throw BadSize();
#endif
    data_ = new T[l];
    length = l;
#ifndef ARRAY_NO_CHECK
    if (!data_) throw MEMFAIL;
#endif
    for (int i = 0; i < l; ++i) data_[i] = asignee;

    return *this;
}
#else

#define Array vector

#endif

/**************************************************************/
template <class T>
class JaggedArray {
private:
    Array < Array < Array <T> > > data_;
public:

    JaggedArray() {
    };
    JaggedArray(int m);
    JaggedArray(int k1, int k2, int k3);

    ~JaggedArray() {
        data_.clear();
    };

    JaggedArray <T> & operator =(const JaggedArray &a) {
        data_ = a.data_;
        return *this;
    };

    void operator =(const T X) {
        if (type == 0) {
            for (int i = 0; i < M; ++i)
                for (int j = 0; i + j < M; ++j)
                    for (int k = 0; i + j + k < M; ++k)
                        data_[i][j][k] = X;
        } else {
            for (int n1 = 0; n1 <= K1; ++n1)
                for (int n2 = 0; n2 <= K2; ++n2)
                    for (int n3 = 0; n3 <= K3; ++n3)
                        data_[n1][n2][n3] = X;
        }
    };

    Array < Array < T > > & operator [] (unsigned int i) {
        return data_[i];
    };

    bool type;
    int M, K1, K2, K3;

};

template <class T> JaggedArray<T>::JaggedArray(int m) : M(m) {
    type = 0;
    data_.assign(m, Array <Array <T> > ());
    for (int k1 = 0; k1 < m; ++k1) {
        data_[k1].assign(m - k1, Array <T > ());
        for (int k2 = 0; k1 + k2 < m; ++k2)
            data_[k1][k2].assign(m - k1 - k2, T());
    }
}

template <class T> JaggedArray<T>::JaggedArray(int k1, int k2, int k3) : K1(k1), K2(k2), K3(k3) {
    type = 1;
    data_.assign(k1 + 1, Array <Array <T> > ());
    for (int n1 = 0; n1 <= k1; ++n1) {
        data_[n1].assign(k2 + 1, Array <T > ());
        for (int n2 = 0; n2 <= k2; ++n2)
            data_[n1][n2].assign(k3 + 1, T());
    }
}

/**************************************************************/
template <class C>
class NeighbSet {
public:

    C N, S, E, W, T, B;

    NeighbSet() {
    };

    NeighbSet(const NeighbSet& P) {
        N = P.N;
        S = P.S;
        E = P.E;
        W = P.W;
        T = P.T;
        B = P.B;
    }

    NeighbSet(const C &P) {
        N = S = E = W = T = B = P;
    }

    C & operator[] (unsigned int i) {
        switch (i) {
            case (0) : return N;
            case (1) : return S;
            case (2) : return E;
            case (3) : return W;
            case (4) : return T;
            case (5) : return B;
            default : throw BoundsViolation();
        }
    }

    class BoundsViolation {
    };
};
/**************************************************************/
template <class C>
class PanelNeighbSet {
public:

    C T, B, L, R;

    PanelNeighbSet() {
    };

    PanelNeighbSet(const PanelNeighbSet& P) {
        T = P.T;
        B = P.B;
        L = P.L;
        R = P.R;
    }

    PanelNeighbSet(const C &P) {
        T = B = L = R = P;
    }

    C & operator[] (unsigned int i) {
        switch (i) {
            case (0) : return L;
            case (1) : return T;
            case (2) : return R;
            case (3) : return B;
            default : throw BoundsViolation();
        }
    }

    class BoundsViolation {
    };
};

#endif /* if !defined(ARRAY_INCL) */
