/* dopri8 -- header file

    Copyright (C) 1995-2014, IMCCE,CNRS

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
    Authors : J. Laskar, F. Joutel, M. Gastineau
    e-mail : laskar@imcce.fr
*/

int gdopri8( int *n,void (*fcn)(int n, double x, double *y, double *f),double *y0,
                   double *x, double *xend, double *eps, double *hmax, double *h);
