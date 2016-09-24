/* information
Copyright (C) Wed Sep 21 18:08:08 2016  Jianshan Zhou
Contact: zhoujianshan@buaa.edu.cn	jianshanzhou@foxmail.com
Website: <https://github.com/JianshanZhou>

This program is free software: you can redistribute
 it and/or modify it under the terms of
 the GNU General Public License as published
 by the Free Software Foundation,
 either version 3 of the License,
 or (at your option) any later version.

This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program.
 If not, see <http://www.gnu.org/licenses/>.

In this module, I developed a series of numerical algorithms for solving
 some typical numerical problems such as linear equations, nonlinear function
 approximation, etc.
*/

#ifndef LINEAR_EQUATION_SOLUTIONS_H_INCLUDED
#define LINEAR_EQUATION_SOLUTIONS_H_INCLUDED


int sequential_Gauss_elimination(double **array, double *solution, int rn, int cn);
int Gauss_elimination_with_master_element(double **array, double *solution, int rn, int cn);
int Doolittle_LU_decomposition(double **array, double *solution, int rn, int cn);
int modified_Doolittle_LU_decomposition(double **array, double *solution, int *M, double *s, int rn, int cn);
int Doolittle_LU_on_banded_matrix(double **array, double *solution, int rn, int cn, int s, int r);
int Doolittle_LU_on_compressed_matrix(double **C, double *solution, double *b, int rn, int cn, int s, int r);
int Crout_LU_decomposition(double *array, double *solution, double *b, int rn, int cn, int r, int s);
int modified_Crout_LU_decomposition(double *array, double * array_sr, double *solution, double *b, int rn, int cn, int r, int s);

#endif // LINEAR_EQUATION_SOLUTIONS_H_INCLUDED
