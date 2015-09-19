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

#include <iostream>
#include <iomanip>
#include "fft.h"

// function main
int main()
{
    //masks(4);

    COMPLEX_NUMBER *arreglo; //[NN];

    // allocate memory for FFT complex numbers
    arreglo = (COMPLEX_NUMBER *) malloc((NN) * sizeof(COMPLEX_NUMBER));

    for (int i = 0; i < NN; ++i)
    {
        arreglo[i].re = static_cast<double>(i);
        arreglo[i].im = static_cast<double>(0);
    }

    //arreglo[0].re = 0.0;
    //arreglo[0].im = 0.0;
    //arreglo[1].re = 1.0;
    //arreglo[1].im = 0.0;
    //arreglo[2].re = 2.0;
    //arreglo[2].im = 0.0;
    //arreglo[3].re = 3.0;
    //arreglo[3].im = 0.0;
    //arreglo[4].re = 4.0;
    //arreglo[4].im = 0.0;
    //arreglo[5].re = 5.0;
    //arreglo[5].im = 0.0;
    //arreglo[6].re = 6.0;
    //arreglo[6].im = 0.0;
    //arreglo[7].re = 7.0;
    //arreglo[7].im = 0.0;
    //arreglo[8].re = 8.0;
    //arreglo[8].im = 0.0;
    //arreglo[9].re = 9.0;
    //arreglo[9].im = 0.0;
    //arreglo[10].re = 0.0;
    //arreglo[10].im = 0.0;
    //arreglo[11].re = 0.0;
    //arreglo[11].im = 0.0;
    //arreglo[12].re = 0.0;
    //arreglo[12].im = 0.0;
    //arreglo[13].re = 0.0;
    //arreglo[13].im = 0.0;
    //arreglo[14].re = 0.0;
    //arreglo[14].im = 0.0;
    //arreglo[15].re = 0.0;
    //arreglo[15].im = 0.0;

    // calculate FFT
    fourier(arreglo, NN, 1);

    //arreglo[0].re = 45.0;
    //arreglo[0].im = 0.0;
    //arreglo[1].re = -25.4520;
    //arreglo[1].im = 16.6652;
    //arreglo[2].re = 10.3640;
    //arreglo[2].im = -3.2929;
    //arreglo[3].re = -9.0641;
    //arreglo[3].im = -2.3285;
    //arreglo[4].re = 4.0;
    //arreglo[4].im = 5.0;
    //arreglo[5].re = -1.2791;
    //arreglo[5].im = -5.6422;
    //arreglo[6].re = -2.3640;
    //arreglo[6].im = 4.7071;
    //arreglo[7].re = 3.7951;
    //arreglo[7].im = -2.6485;
    //arreglo[8].re = -5.0;
    //arreglo[8].im = 0;
    //arreglo[9].re = 3.7951;
    //arreglo[9].im = 2.6485;
    //arreglo[10].re = -2.3640;
    //arreglo[10].im = -4.7071;
    //arreglo[11].re = -1.2791;
    //arreglo[11].im = 5.6422;
    //arreglo[12].re = 4.0;
    //arreglo[12].im = -5.0;
    //arreglo[13].re = -9.0641;
    //arreglo[13].im = 2.3285;
    //arreglo[14].re = 10.3640;
    //arreglo[14].im = 3.2929;
    //arreglo[15].re = -25.4520;
    //arreglo[15].im = -16.6652;

    //// calculate IFFT
    fourier(arreglo, NN, -1);
    //// need to normalize
    
    // release allcated memory
    free(arreglo);

    // avoid reuse
    arreglo = NULL;

    return 0;
}   // end main

