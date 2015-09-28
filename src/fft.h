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
//       * 
//       * use templates
//       * use c++ 14 initializers, constructors, etc.
//       * your contribution is welcome

#include <math.h>

#define _USE_MATH_DEFINES
#define PI          M_PI            // pi to machine precision, defined in math.h
#define TWOPI       (2.0 * PI)
#define NBITS       10               // WA: hard-coded value to number of bits
#define NBITS_TW    (NBITS - 1)     // 
#define _DEBUG_

const int NN = 1024;  // Number of points of the FFT

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
    unsigned int bits : NBITS_TW;
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

// prototypes

// function to return the mask accordind to the number of bits
int masks(int nbits);
// function to swap numbers
inline int swap(double *a, double *b);
// fast fourier transform
void fourier(COMPLEX_NUMBER data[], int nn, int isign);
