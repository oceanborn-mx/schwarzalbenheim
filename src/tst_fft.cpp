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

    // calculate FFT
    fourier(arreglo, NN, 1);

    //// calculate IFFT
    fourier(arreglo, NN, -1);
    //// need to normalize
    
    // release allcated memory
    free(arreglo);

    // avoid reuse
    arreglo = NULL;

    return 0;
}   // end main

