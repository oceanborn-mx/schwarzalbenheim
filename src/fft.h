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
//       * use c++ 14 initializers, constructors, etc.
//       * your contribution is welcome

#include <math.h>

#define _USE_MATH_DEFINES
#define PI          M_PI            // pi to machine precision, defined in math.h
#define TWOPI       (2.0 * PI)
#define NBITS       4               // WA: hard-coded value to number of bits
#define NBITS_TW    (NBITS - 1)     // 
//#define NBITS(X)    (log2(X))
#define _DEBUG_

constexpr int NN       = 16;           // Number of points of the FFT
//constexpr int aux = log2(NN);
//const int nbits = 4;
//constexpr log2
//static constexpr int nbits_n = static_cast<int>(log2(NN));
//static constexpr int nbits = nbits_n;
//constexpr int NBITS    =   nbits;   // Number of bits to fit the number of FFT points
//constexpr int mylog2(int base)
//{
//    int aux = static_cast<int>(log2(base));
//    return aux;
//}
const int aux = static_cast<int>(log2(NN));
//enum nbits
//{
//    NBITS = aux
//};
int NSTAGES  = log2(NN);     // Number of stages to complete the FFT

//constexpr int NBITS = noexcept(mylog2(NN));
//constexpr void set_constants()
//{
//    int NBITS = log2(NN);
//}   // end set_constants

unsigned int MASK1;
unsigned int MASK2;

// needed to eliminate the hard-code value
//enum bits
//{
//    NBITS = static_cast<int>(log2(NN))
//};

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

// function to return the mask accordind to the number of bits
void masks(int nbits)
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
}
