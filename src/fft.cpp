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
//       * make generic the swapping loops
//       * remove the hard-coded values.
//       * use c++ 14 initializers, constructors, etc.
//       * your contribution is welcome

#include <iostream>
#include <iomanip>
#include "fft.h"

using namespace std;

unsigned int MASK1;
unsigned int MASK2;

int NSTAGES  = log2(NN);     // Number of stages to complete the FFT

void fourier(COMPLEX_NUMBER data[], int nn, int isign)
{
    // variables to perform the bit reversal algorithm
    BITS u, _au, _bu, _cu, _xu, _yu, _u;
    BITS v, _av, _bv, _cv, _xv, _yv, _v;

    // array to hold the values nedeed to swap when
    // butterfly operation are done
    BITS array[NN];

    COMPLEX_NUMBER array_swap[NN];

    BITS _ar, _br, _cr, _xr, _yr, _r;
    BITS _as, _bs, _cs, _xs, _ys, _s;

    // bits for twiddle factor counter
    BITS_TWIDDLE _w;

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

    // masks
    masks(NBITS);

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
            wr = cos(theta);                            // real
            wi = sin(theta);                            // imaginary
                                    
            // bit reversal algorithm
            u.bits = j;
            v.bits = j + 1;

            // bit reversal algorithm to swap values 
            // when butterfly operations are done
            if (0 == i)
            {
                array[j].bits   = u.bits;
                array[j+1].bits = v.bits;
            }
            
            // hold the actual value
            array[j].bits   = array[j].bits;
            array[j+1].bits = array[j+1].bits;

            _ar.bits = array[j].bits & MASK1;
            _br.bits = array[j].bits & (MASK2 >> i);
            _cr.bits = array[j].bits & ((~MASK2) >> i);
            _xr.bits = _ar.bits << (NSTAGES - 1 - i);
            _yr.bits = _br.bits >> 1;

            _as.bits = array[j+1].bits & MASK1;
            _bs.bits = array[j+1].bits & (MASK2 >> i);
            _cs.bits = array[j+1].bits & ((~MASK2) >> i);
            _xs.bits = _as.bits << (NSTAGES - 1 - i);
            _ys.bits = _bs.bits >> 1;

            _r.bits = _xr.bits | _yr.bits | _cr.bits;
            _s.bits = _xs.bits | _ys.bits | _cs.bits;

            array[j].bits   = _r.bits;
            array[j+1].bits = _s.bits;

            // this first approach was intended to do what is actually
            // perform in the code above (see the README file), but, 
            // is actually doing what the butterfly architecture
            // is intended to be, this initial bug is in fact a feature,
            // is a kind of penicillin :)
            _au.bits = u.bits & MASK1;
            _bu.bits = u.bits & (MASK2 >> i);
            _cu.bits = u.bits & ((~MASK2) >> i);
            _xu.bits = _au.bits << (NSTAGES - 1 - i);
            _yu.bits = _bu.bits >> 1;

            _av.bits = v.bits & MASK1;
            _bv.bits = v.bits & (MASK2 >> i);
            _cv.bits = v.bits & ((~MASK2) >> i);
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

    // swapping
    // counter variables for below loops
    int h = 0;
    int q = 0;

    // first half
    //for (h = 0; h < (NN / 8); ++h)
    //{
    //    for (q = h * NN / (NN / 4) + 1; q < (1 + h) * NN / (NN / 4); ++q)
    //    {
    //        swap(&data[array[q].bits].re, &data[q].re);
    //        swap(&data[array[q].bits].im, &data[q].im);
    //    }
    //}

    //// second half
    //for (h = 0; h < (NN / 8) / 2; ++h)
    //{
    //    for (q = NN / 2 + (1 + h) * NN / (NN / 4) - 1; q < NN / 2 + (1 + h) * NN / (NN / 4); ++q)
    //    {
    //        swap(&data[array[q].bits].re, &data[q].re);
    //        swap(&data[array[q].bits].im, &data[q].im);
    //    }
    //}

    // WA - no optimum
    for (h = 0; h < NN; ++h)
    {
        array_swap[h] = data[array[h].bits];
    }
    
    for (q = 0; q < NN; ++q)
    {
        data[q] = array_swap[q];
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

// function to return the mask accordind to the number of bits
int masks(int nbits)
{
    switch (nbits)
    {
        case 4:
            MASK1 = 0x1;
            MASK2 = 0xF;
            break;
        case 5:
            MASK1 = 0x01;
            MASK2 = 0x1F;
            break;
        case 6:
            MASK1 = 0x01;
            MASK2 = 0x3F;
            break;
        case 7:
            MASK1 = 0x01;
            MASK2 = 0x7F;
            break;
        case 8:
            MASK1 = 0x01;
            MASK2 = 0xFF;
            break;
        case 9:
            MASK1 = 0x001;
            MASK2 = 0x1FF;
            break;
        case 10:
            MASK1 = 0x001;
            MASK2 = 0x3FF;
            break;
        case 11:
            MASK1 = 0x001;
            MASK2 = 0x7FF;
            break;
        case 12:
            MASK1 = 0x001;
            MASK2 = 0xFFF;
            break;
        case 13:
            MASK1 = 0x0001;
            MASK2 = 0x1FFF;
            break;
        case 14:
            MASK1 = 0x0001;
            MASK2 = 0x3FFF;
            break;
        case 15:
            MASK1 = 0x0001;
            MASK2 = 0x7FFF;
            break;
        case 16:
            MASK1 = 0x0001;
            MASK2 = 0xFFFF;
            break;
        default:
            // TODO: find a restriction
            MASK1 = 0x0001;
            MASK2 = 0xFFFF;
    }   // end switch

    return EXIT_SUCCESS;
}   // end masks

// function to swap numbers
inline int swap(double *a, double *b)
{
    double temp;    // temporary

    temp = *a;
    *a = *b;
    *b = temp;

    return EXIT_SUCCESS;
}   // end swap
