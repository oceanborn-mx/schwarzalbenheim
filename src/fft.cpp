///////////////////////////////////////////////////////////////////////////////
//    FFT: Fast Fourier Transform                                             
//    Copyright (C) 2015  Daniel A Razo M (da.razo@yahoo.com)             
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////

// TODO: * clean-up the code.
//       * make generic abount the number of FFT points, always in power of 2
//       * remove the hard-coded values.
//       * try to optimize the ** SWAP ** function.
//       * use c++ 14 initializers, constructors, etc.
//       * your contribution is welcome

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <iomanip>

#define PI      M_PI        // pi to machine precision, defined in math.h
#define TWOPI   (2.0 * PI)
#define NBITS   4           // WA: hard-coded value to number of bits
#define _DEBUG_

using namespace std;

const unsigned int NN       = 16;           // Number of points of the FFT
//const unsigned int NBITS    = log2(NN);     // Number of bits to fit the number of FFT points
const unsigned int NSTAGES  = log2(NN);     // Number of stages to complete the FFT

// struct to hold the number of bits needed 
// to fit the number of FFT points
typedef struct BITS 
{
    unsigned int bits : NBITS;  // hard-coded value needed, otherwise, copiling error
} BITS;

// struct to hold the number of bits needed 
// to fit the twiddle factor
typedef struct BITS_TWIDDLE
{
    unsigned int bits : 3;
} BITS_TWIDDLE;

// struct to define the layout of a complex number
typedef struct COMPLEX_NUMBER
{
    union
    {
        long double complexN;   // complex number

        struct
        {
            double re;          // real
            double im;          // imaginary
        };
    };
} COMPLEX_NUMBER;

void fourier(COMPLEX_NUMBER data[], int nn, int isign)
{
    // variables to perform the bit reversal algorithm
    BITS u, _au, _bu, _cu, _xu, _yu, _u;
    BITS v, _av, _bv, _cv, _xv, _yv, _v;

    BITS uaux, vaux, array[16];

    BITS a, b, c, x, y, r;
    BITS aa, bb, cc, xx, yy, s;

    // bits for twiddle factor counter
    BITS_TWIDDLE _w;

    //int array[16];

    // butterfly
    double theta, wr, wi;
    double temprA, tempiA;
    double temprB, tempiB;
    double auxrB, auxiB;
    double temp1, temp2;
    int k;

    // twiddle factor
	wr = 1.0;
	wi = 0.0;

#ifdef _DEBUG_
    cout << "NSTAGES: " << NSTAGES << endl;
#endif

    for (int i = 0; i < NSTAGES; ++i)
    {
        // bit counter for twiddle factor
        k = 0;

#ifdef _DEBUG_ 
        //cout << "### i: " << i << endl;
#endif
        for (int j = 0; j < NN - 1; j += 2)
        {
            // twiddle factor bits 
            _w.bits = k;
            _w.bits <<= i;
            ++k;
        
            // twiddle factor
            theta = (isign * TWOPI * _w.bits) / NN;     // angle
            wr = cos(theta);    // real
            wi = sin(theta);    // imaginary

                                    
            // bit reversal algorithm
            u.bits = j;
            v.bits = j + 1;

            // auxiliar counter for **swap** function
            // initialize
            if (0 == i)
            {
                uaux.bits = u.bits;
                vaux.bits = v.bits;
                array[j].bits = u.bits;
                array[j + 1].bits = v.bits;
            }
            
            // hold the actual value
            uaux.bits = uaux.bits;
            vaux.bits = vaux.bits;
            array[j].bits = array[j].bits;
            array[j + 1].bits = array[j + 1].bits;

            a.bits = array[j].bits & 0x1;
            b.bits = array[j].bits & (0xF >> i);
            c.bits = array[j].bits & ((~0xF) >> i);
            x.bits = a.bits << (NSTAGES - 1 - i);
            y.bits = b.bits >> 1;

            aa.bits = array[j+1].bits & 0x1;
            bb.bits = array[j+1].bits & (0xF >> i);
            cc.bits = array[j+1].bits & ((~0xF) >> i);
            xx.bits = aa.bits << (NSTAGES - 1 - i);
            yy.bits = bb.bits >> 1;

            r.bits = x.bits | y.bits | c.bits;
            s.bits = xx.bits | yy.bits | cc.bits;


            array[j].bits = r.bits;
            array[j+1].bits = s.bits;

            // this first approach was intended to do what is actually
            // perform in the code above (see the README file), but, 
            // it actually doing what the butterfly architecture
            // is intended to be, this initial bug is in fact a feature,
            // is a kind of penicillin :)
            _au.bits = u.bits & 0x1;
            _bu.bits = u.bits & (0xF >> i);
            _cu.bits = u.bits & ((~0xF) >> i);
            _xu.bits = _au.bits << (NSTAGES - 1 - i);
            _yu.bits = _bu.bits >> 1;

            _av.bits = v.bits & 0x1;
            _bv.bits = v.bits & (0xF >> i);
            _cv.bits = v.bits & ((~0xF) >> i);
            _xv.bits = _av.bits << (NSTAGES - 1 - i);
            _yv.bits = _bv.bits >> 1;
   
            _u.bits = _xu.bits | _yu.bits | _cu.bits;
            _v.bits = _xv.bits | _yv.bits | _cv.bits;


            // butterfly
            temprA = data[_u.bits].re + data[_v.bits].re;   // A = a + b
            tempiA = data[_u.bits].im + data[_v.bits].im;

            auxrB = data[_u.bits].re - data[_v.bits].re;    // (a - b)
            auxiB = data[_u.bits].im - data[_v.bits].im;

            temprB = auxrB * wr - auxiB * wi;               // B = (a - b) * W
            tempiB = auxiB * wr + auxrB * wi;

            data[_u.bits].re = temprA;
            data[_u.bits].im = tempiA;

            data[_v.bits].re = temprB;
            data[_v.bits].im = tempiB;
            
#ifdef _DEBUG_ 
            //cout << "Stage: " << i << endl;
            //cout << "_u.bits = j: " << j << endl;
            //cout << "_w.bits: " << _w.bits << endl;
            
            //cout << "_u.bits (k)    : " << _u.bits << endl;
            //cout << "_v.bits (k + 1): " << _v.bits << endl;

            //cout << "uaux.bits (k)    : " << uaux.bits << endl;
            //cout << "vaux.bits (k + 1): " << vaux.bits << endl;

            //cout << "array[" << j <<     "].bits (k)    : " << array[j].bits << endl;
            //cout << "array[" << j + 1 << "].bits (k + 1): " << array[j+1].bits << endl;

            //cout << "dataRe[" << _u.bits << "] = " << data[_u.bits].re << endl;
            //cout << "dataIm[" << _u.bits << "] = " << data[_u.bits].im << endl;

            //cout << "dataRe[" << _v.bits << "] = " << data[_v.bits].re << endl;
            //cout << "dataIm[" << _v.bits << "] = " << data[_v.bits].im << endl;
            //
            // alternate
            //cout << "dataRe[" << r.bits << "] = " << data[r.bits].re << endl;
            //cout << "dataIm[" << r.bits << "] = " << data[r.bits].im << endl;

            //cout << "dataRe[" << s.bits << "] = " << data[s.bits].re << endl;
            //cout << "dataIm[" << s.bits << "] = " << data[s.bits].im << endl;
#endif
        }   // end for
    }   // end for

    // swap
    for (int jj = 1; jj < (NN / 4); ++jj)
    {
        temp1 = data[array[jj].bits].re;
        temp2 = data[array[jj].bits].im;

        data[array[jj].bits].re = data[jj].re;
        data[array[jj].bits].im = data[jj].im;

        data[jj].re = temp1;
        data[jj].im = temp2;
    }

    for (int tt = NN / 4 + 1; tt < (NN / 2); ++tt)
    {
        temp1 = data[array[tt].bits].re;
        temp2 = data[array[tt].bits].im;

        data[array[tt].bits].re = data[tt].re;
        data[array[tt].bits].im = data[tt].im;

        data[tt].re = temp1;
        data[tt].im = temp2;
    }

    for (int kk = NN / 2 + NN / 4 - 1; kk < (NN * 3 / 4); ++kk)
    {
        temp1 = data[array[kk].bits].re;
        temp2 = data[array[kk].bits].im;

        data[array[kk].bits].re = data[kk].re;
        data[array[kk].bits].im = data[kk].im;

        data[kk].re = temp1;
        data[kk].im = temp2;
    }

#ifdef _DEBUG_ 
    // print out the FFT
    for (int qq = 0; qq < NN; ++qq)
    {
        cout << right << "dataRe[" << setw(4) << qq << "] = " 
             << setw(8) << fixed << setprecision(4) << data[qq].re << "\t\t";

        cout << right << "dataIm[" << setw(4) << qq << "] = " 
             << setw(8) << fixed << setprecision(4) << data[qq].im << endl;
    }
#endif

}   // end fourier

// function main
int main()
{
    COMPLEX_NUMBER arreglo[16];

    arreglo[0].re = 0.0;
    arreglo[0].im = 0.0;
    arreglo[1].re = 1.0;
    arreglo[1].im = 0.0;
    arreglo[2].re = 2.0;
    arreglo[2].im = 0.0;
    arreglo[3].re = 3.0;
    arreglo[3].im = 0.0;
    arreglo[4].re = 4.0;
    arreglo[4].im = 0.0;
    arreglo[5].re = 5.0;
    arreglo[5].im = 0.0;
    arreglo[6].re = 6.0;
    arreglo[6].im = 0.0;
    arreglo[7].re = 7.0;
    arreglo[7].im = 0.0;
    arreglo[8].re = 8.0;
    arreglo[8].im = 0.0;
    arreglo[9].re = 9.0;
    arreglo[9].im = 0.0;
    arreglo[10].re = 0.0;
    arreglo[10].im = 0.0;
    arreglo[11].re = 0.0;
    arreglo[11].im = 0.0;
    arreglo[12].re = 0.0;
    arreglo[12].im = 0.0;
    arreglo[13].re = 0.0;
    arreglo[13].im = 0.0;
    arreglo[14].re = 0.0;
    arreglo[14].im = 0.0;
    arreglo[15].re = 0.0;
    arreglo[15].im = 0.0;

    // calculate FFT
	fourier(arreglo, 16, 1);

    arreglo[0].re = 45.0;
    arreglo[0].im = 0.0;
    arreglo[1].re = -25.4520;
    arreglo[1].im = 16.6652;
    arreglo[2].re = 10.3640;
    arreglo[2].im = -3.2929;
    arreglo[3].re = -9.0641;
    arreglo[3].im = -2.3285;
    arreglo[4].re = 4.0;
    arreglo[4].im = 5.0;
    arreglo[5].re = -1.2791;
    arreglo[5].im = -5.6422;
    arreglo[6].re = -2.3640;
    arreglo[6].im = 4.7071;
    arreglo[7].re = 3.7951;
    arreglo[7].im = -2.6485;
    arreglo[8].re = -5.0;
    arreglo[8].im = 0;
    arreglo[9].re = 3.7951;
    arreglo[9].im = 2.6485;
    arreglo[10].re = -2.3640;
    arreglo[10].im = -4.7071;
    arreglo[11].re = -1.2791;
    arreglo[11].im = 5.6422;
    arreglo[12].re = 4.0;
    arreglo[12].im = -5.0;
    arreglo[13].re = -9.0641;
    arreglo[13].im = 2.3285;
    arreglo[14].re = 10.3640;
    arreglo[14].im = 3.2929;
    arreglo[15].re = -25.4520;
    arreglo[15].im = -16.6652;

    // calculate IFFT
	fourier(arreglo, 16, -1);
    // need to normalize    

    return 0;
}   // end main
