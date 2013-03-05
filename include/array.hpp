/*
This file is part of the Combined Wake Modelling Code Version 1.0

V3D Code Copyright Tom McCombes 2013
This code solves the 3D unsteady incompressible
Navier-Stokes equations in velocity vorticity form

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

#define USE_ARRAY

#ifndef ARRAY_INCL
#define ARRAY_INCL
#ifndef TEST_MODE
#include "includes.hpp"
//#include "types.hpp"
#endif
#ifdef USE_ARRAY

template <class T>
class Array {
public:

    Array() : length(0), data_(NULL) {
    };

    Array(int n) : length(n), data_(new T[length]) {
#ifndef ARRAY_NO_CHECK
        if (!data_) throw MEMFAIL;
#endif
    };


#ifdef ARRAY_NO_CHECK

    Array(int n, const T &t) : length(n), data_(new T[length]) {
#else

    Array(int n, const T &t) {

        if (n < 0) {

            printf("Attempted to create an array of length %d. Aborting. ", n);
            throw BadSize();
        }
        length = n;
        data_ = new T[length];
        if (!data_) throw MEMFAIL;
#endif
        for (int i = 0; i < length; ++i) data_[i] = t;
    };

    Array(const Array <T> &a) : length(a.length), data_(new T[length]) {
#ifndef ARRAY_NO_CHECK
        if (!data_) throw MEMFAIL;
#endif
        //        memcpy (data_,a.data_, length * sizeof(T));
        for (int i = 0; i < length; ++i) data_[i] = a.data_[i];
    };

    ~Array();



    Array & operator=(const Array&);
    Array & operator=(const T&);
    std::vector <T> & ReturnVector(){
        return std::vector <T> (data_,data_ + length);
    };

    void clear();
    
    static int partition(Array <T> &y, int f, int l) {
        int up, down;
        T temp, piv = y[f];
        up = f;
        down = l;
        goto partLS;
        do {
            temp = y[up];
            y[up] = y[down];
            y[down] = temp;
partLS:
            while (y[up] <= piv && up < l) {
                up++;
            }
            while (y[down] > piv && down > f) {
                down--;
            }
        } while (down > up);
        y[f] = y[down];
        y[down] = piv;
        return down;
    };

    static void quicksort1(Array <T> &x, int first, int last) {
        int pivIndex = 0;
        if (first < last) {
            pivIndex = Array::partition(x, first, last);
            Array::quicksort1(x, first, (pivIndex - 1));
            Array::quicksort1(x, (pivIndex + 1), last);
        }
    };

    static void QuickSortA(Array <T> &In) {
        Array::quicksort1(In, 0, In.size() - 1);
    };

    static void QuickSortB(Array <T> &In) {
        Array::quicksort2(In, 0, In.size() - 1);
    };
    
    static void quicksort2(Array <T> &list, int beg, int end) {
        T piv;
        T tmp;

        int l, r, p;

        while (beg < end) // This while loop will avoid the second recursive call
        {
            l = beg;
            p = (beg + end) / 2;
            r = end;

            piv = list[p];

            while (1) {
                while ((l <= r) && (list[l]  <= piv)) l++;
                while ((l <= r) && (list[r]  > piv)) r--;

                if (l > r) break;

                tmp = list[l];
                list[l] = list[r];
                list[r] = tmp;

                if (p == r) p = l;

                l++;
                r--;
            }

            list[p] = list[r];
            list[r] = piv;
            r--;

            // Recursion on the shorter side & loop (with new indexes) on the longer
            if ((r - beg)<(end - l)) {
                quicksort2(list, beg, r);
                beg = l;
            } else {
                quicksort2(list, l, end);
                end = r;
            }
        }
    };
    

    void pop(int);
    Array& push_back(const T&);
    Array& push_back(const Array <T>&);
    Array& push_front(const T&);
    Array& pop_front();
    Array& pop_back();

    inline friend std::ostream & operator <<(std::ostream& os, const Array &in) {
        for (int i = 0; i < in.length - 1; ++i)
            os << in.data_[i] << "\t ";

        return os << in.data_[in.length - 1] << "\n";
    };

    int size() {
        return (int) length;
    }
#ifndef TEST_MODE

    void print() {
        for (int i = 0; i < length; ++i) std::cout << data_[i] << " ";
        std::cout << std::endl;
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

    T & operator[] (int i) const {
        if (i >= length || i < 0) {
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


    inline friend Array <T> sin(const Array<T> A) {
        Array <T> L(A.length);
        for (int i = 0; i < A.length; ++i)
            L.data_[i] = sin(A.data_[i]);

        return L;
    }

    inline friend T max(const Array<T> A) {
        T L = A.data_[0];
        for (int i = 1; i < A.length; ++i)
            L = max(L, A.data_[i]);

        return L;
    }

    inline friend T min(const Array<T> A) {
        T L = A.data_[0];
        for (int i = 1; i < A.length; ++i)
            L = min(L, A.data_[i]);

        return L;
    }

    inline friend Array <T> cos(const Array<T> A) {
        Array <T> L(A.length);
        for (int i = 0; i < A.length; ++i)
            L.data_[i] = cos(A.data_[i]);

        return L;
    }

    inline friend Array <T> tan(const Array<T> A) {
        Array <T> L(A.length);
        for (int i = 0; i < A.length; ++i)
            L.data_[i] = tan(A.data_[i]);

        return L;
    }
    
    inline friend Array <T> fabs(const Array<T> A) {
        Array <T> L(A.length);
        for (int i = 0; i < A.length; ++i)
            L.data_[i] = fabs(A.data_[i]);

        return L;
    }

    inline friend Array <T> sqrt(const Array<T> A) {
        Array <T> L(A.length);
        for (int i = 0; i < A.length; ++i)
            L.data_[i] = sqrt(A.data_[i]);

        return L;
    }

    inline friend Array <T> atan2(const Array<T> A, const Array<T> B) {
        Array <T> L(A.length);
        for (int i = 0; i < A.length; ++i)
            L.data_[i] = atan2(A.data_[i], B.data_[i]);

        return L;
    }

    inline friend Array <T> atan(const Array<T> A, const Array<T> B) {
        Array <T> L(A.length);
        for (int i = 0; i < A.length; ++i)
            L.data_[i] = atan(A.data_[i] / B.data_[i]);

        return L;
    }

    inline friend Array <T> log(const Array<T> A) {
        Array <T> L(A.length);
        for (int i = 0; i < A.length; ++i)
            L.data_[i] = log(A.data_[i]);

        return L;
    }

    //  Some addition operators

    class WrongSizeSum {
    };

    template <class U> Array operator+(const U &a) {
        Array b(length);
        for (int i = 0; i < length; i++)
            b.data_[i] = data_[i] + a;
        return b;
    }

    //     template <class U> friend Array operator+(const Array A, const U B ) {return A + B;}

    template <class U> friend Array operator+(const U &A, const Array &B) {
        Array c(B.length);
        for (int i = 0; i < B.length; i++)
            c.data_[i] = B.data_[i] + A;
        return c;
    }

    Array operator+(const Array &a) {
#ifndef ARRAY_NO_CHECK
        if (length != a.length) {
            std::cout << "Size of *this array: " << length << ", size of input array: " << a.length << std::endl;
            std::cout << a << std::endl;
            throw WrongSizeSum();
        }
#endif
        Array b(length);
        for (int i = 0; i < length; i++)
            b.data_[i] = data_[i] + a.data_[i];
        return b;
    }

    inline void operator +=(const Array &B) {
#ifndef ARRAY_NO_CHECK
        if (length != B.length) throw WrongSizeSum();
#endif
        for (int i = 0; i < length; ++i)
            data_[i] += B.data_[i];
    }

    template <class U> inline friend void operator +=(Array &A, const U &B) {
        for (int i = 0; i < A.length; ++i)
            A.data_[i] += B;
    }

    template <class U> inline friend void operator +=(const U &A, Array &B) {
        for (int i = 0; i < B.length; ++i)
            B.data_[i] += A;
    }


    //  Some subtraction operators

    class WrongSizeSubtract {
    };

    template <class U> Array operator-(const U &a) {
        Array b(length);
        for (int i = 0; i < length; i++)
            b.data_[i] = a - data_[i];
        return b;
    }

    template <class U> friend Array operator-(const U &A, const Array &B) {
        Array c(B.length);
        for (int i = 0; i < B.length; i++)
            c.data_[i] = B.data_[i] - A;
        return c;
    }

    Array operator-(const Array &a) {
#ifndef ARRAY_NO_CHECK
        if (length != a.length) {
            std::cout << "Size of *this array: " << length << ", size of input array: " << a.length << std::endl;
            std::cout << a << std::endl;
            throw WrongSizeSubtract();
        }
#endif
        Array b(length);
        for (int i = 0; i < length; i++)
            b.data_[i] = data_[i] - a.data_[i];
        return b;
    }

    inline void operator -=(const Array &B) {
#ifndef ARRAY_NO_CHECK
        if (length != B.length) throw WrongSizeSubtract();
#endif
        for (int i = 0; i < length; ++i)
            data_[i] -= B.data_[i];
    }

    template <class U> inline friend void operator -=(Array &A, const U &B) {
        for (int i = 0; i < A.length; ++i)
            A.data_[i] -= B;
    }

    template <class U> inline friend void operator -=(const U &A, Array &B) {
        for (int i = 0; i < B.length; ++i)
            B.data_[i] -= A;
    }

    //  Some multiplication operators

    class WrongSizeMultiply {
    };

    template <class U> Array operator*(const U &a) {
        Array b(length);
        for (int i = 0; i < length; i++)
            b.data_[i] = data_[i] * a;
        return b;
    }

    template <class U> friend Array operator*(const U &A, const Array &B) {
        Array c(B.length);
        for (int i = 0; i < B.length; i++)
            c.data_[i] = B.data_[i] * A;
        return c;
    }

    Array operator*(const Array &a) {
#ifndef ARRAY_NO_CHECK
        if (length != a.length) {
            std::cout << "Size of *this array: " << length << ", size of input array: " << a.length << std::endl;
            std::cout << a << std::endl;
            throw WrongSizeMultiply();
        }
#endif
        Array b(length);
        for (int i = 0; i < length; i++)
            b.data_[i] = data_[i] * a.data_[i];
        return b;
    }

    inline void operator *=(const Array &B) {
#ifndef ARRAY_NO_CHECK
        if (length != B.length) throw WrongSizeMultiply();
#endif
        for (int i = 0; i < length; ++i)
            data_[i] *= B.data_[i];
    }

    template <class U> inline friend void operator *=(Array &A, const U &B) {
        for (int i = 0; i < A.length; ++i)
            A.data_[i] *= B;
    }

    template <class U> inline friend void operator *=(const U &A, Array &B) {
        for (int i = 0; i < B.length; ++i)
            B.data_[i] *= A;
    }

    //  Some division operators

    class WrongSizeDivide {
    };

    template <class U> Array operator/(const U &a) {
        Array b(length);
        for (int i = 0; i < length; i++)
            b.data_[i] = data_[i] / a;
        return b;
    }

    template <class U> friend Array operator/(const U &A, const Array &B) {
        Array c(B.length);
        for (int i = 0; i < B.length; i++)
            c.data_[i] = A / B.data_[i];
        return c;
    }

    Array operator/(const Array &a) {
#ifndef ARRAY_NO_CHECK
        if (length != a.length) {
            std::cout << "Size of /this array: " << length << ", size of input array: " << a.length << std::endl;
            std::cout << a << std::endl;
            throw WrongSizeDivide();
        }
#endif
        Array b(length);
        for (int i = 0; i < length; i++)
            b.data_[i] = data_[i] / a.data_[i];
        return b;
    }

    inline void operator /=(const Array &B) {
#ifndef ARRAY_NO_CHECK
        if (length != B.length) throw WrongSizeDivide();
#endif
        for (int i = 0; i < length; ++i)
            data_[i] /= B.data_[i];
    }

    template <class U> inline friend void operator /=(Array &A, const U &B) {
        for (int i = 0; i < A.length; ++i)
            A.data_[i] /= B;
    }


private:
    int length;
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

        for (int i = 0; i < length; ++i) data_[i] = a.data_[i];
    }
    return *this;
}
/************************************************/

/*  Pop out a single element  */
template <class T> void Array<T>::pop(int p) {
    {
#ifndef ARRAY_NO_CHECK
        if (p >= length || length == 1) {
            printf("Array Bounds Violation: Attempted to pop out element %d, or array too short: array length %d - ", p, length);
            throw BoundsViolation();
        }
#endif
        printf("L %d ", length);
        T* data_new = new T [length - 1];
        memcpy(data_new, data_, p * sizeof (T));
        memcpy(data_new + p, data_ + p + 1, (length - p - 1) * sizeof (T));
        delete[] data_;
        data_ = data_new;
        length--;

    }


}


/************************************************/

/*  Set to values  */
template <class T> Array<T>& Array<T>::operator =(const T &a) {
    for (int i = 0; i < length; i++) data_[i] = a;
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
template <class T> T & Array<T>::front() {
    if (data_)
        return data_[0];
    else
        throw BoundsViolation();
}

/************************************************/
template <class T> T & Array<T>::back() {
    if (data_)
        return data_[length - 1];
    else
        throw BoundsViolation();
}

/************************************************/
template <class T> Array<T>& Array<T>::push_back(const T &insrt) {
    T* data_new = new T [length + 1];

    if (data_) {
        for (int i = 0; i < length; ++i) data_new[i] = data_[i];
        delete[] data_;
    }
    data_new[length] = insrt;
    data_ = data_new;
    length += 1;
    return *this;
}
/************************************************/
template <class T> Array<T>& Array<T>::push_back(const Array <T> &insrt) {
    T* data_new = new T [length + insrt.length];

    int cnt = 0;
    if (data_) {
        for (int i = 0; i < length; ++i, ++cnt) data_new[i] = data_[i];
        delete[] data_;
    }
    if (insrt.data_) {
        for (int i = 0; i < insrt.length; ++i, ++cnt) data_new[cnt] = insrt.data_[i];
    }

    data_ = data_new;
    length = length + insrt.length;
    return *this;
}
/************************************************/
template <class T> Array<T>& Array<T>::push_front(const T &insrt) {
    T* data_new = new T [length + 1];

    if (data_) {
        for (int i = 0; i < length; ++i) data_new[i + 1] = data_[i];
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
        for (int i = 0; i < length - 1; ++i) data_new[i] = data_[i];
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
        for (int i = 1; i < length; ++i) data_new[i - 1] = data_[i];
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
    if (l < 0) {
        printf("Attempted to allocate an array of length %d. Aborting.", l);
        throw BadSize();
    }
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
    if (l < 0) {
        printf("Attempted to assign an array of length %d. Aborting.", l);
        throw BadSize();
    }
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

    Array < Array < T > > & operator [] (int i) {
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

    C & operator[] (int i) {
        switch (i) {
            case (0): return N;
            case (1): return S;
            case (2): return E;
            case (3): return W;
            case (4): return T;
            case (5): return B;
            default: throw BoundsViolation();
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

    C & operator[] (int i) {
        switch (i) {
            case (0): return B;
            case (1): return R;
            case (2): return T;
            case (3): return L;
            default: throw BoundsViolation();
        }
    }

    class BoundsViolation {
    };
};

/**************************************************************/
class Vect3;

class MATLABOutputStruct {
public:
    Array < Array < int > > Int1DArrays;
    Array < std::string > Int1DArrayStrings;

    Array < Array < double > > Double1DArrays;
    Array < std::string > Double1DArrayStrings;

    Array < Array < Vect3 > > Vect1DArrays;
    Array < std::string > Vect1DArrayStrings;

    Array < Array < Array < int> > > Int2DArrays;
    Array < std::string > Int2DArrayStrings;

    Array < Array < Array <double> > > Double2DArrays;
    Array < std::string > Double2DArrayStrings;

    Array < Array < Array <Vect3> > > Vect2DArrays;
    Array < std::string > Vect2DArrayStrings;

    
};
/**************************************************************/
#endif /* if !defined(ARRAY_INCL) */
