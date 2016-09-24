/* information
Copyright (C) Thu Sep 20 10:11:08 2016  Jianshan Zhou
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

In this module, the function that is used to read double-type data from external
 txt file is developed. Additionally, the data recorded in the txt file should be
 separated by using a space and arranged into the form of ROW_NUM rows and COL_NUM
 columns.
*/


#ifndef IMPORT_EXPORT_DATA_H_INCLUDED
#define IMPORT_EXPORT_DATA_H_INCLUDED

int read_data_from_txt(FILE *filepoint, double **array, int rn, int cn);
int write_data_to_txt(FILE *filepointer, double *array, int n);
int establish_compressed_matrix(FILE *filepoint, double ** array, int rn, int cn, double *b, int r, int s);
int establish_one_array(FILE *filepoint, double ** array, int rn, int cn, double *b, int r, int s);
int establish_one_quasi_array(FILE *filepoint, double *array, int rn, int cn, double *b, int r, int s);

#endif // IMPORT_EXPORT_DATA_H_INCLUDED
