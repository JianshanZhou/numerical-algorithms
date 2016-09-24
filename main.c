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

This is the main file that calls other functions and does run the program.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "import_export_data.h"
#include "linear_equation_solutions.h"

/*Define the numbers of the rows and the columns of a specified array.
These two parameters should be given before calling this module file.
The array corresponds to an augmented matrix of a linear equation system, so that
the row number of the array is specified to satisfy the equation:
row number + 1 = column number.*/
#define ROW_NUM 5
#define COL_NUM 6
#define S 1
#define R 1

/*Note that three cases with data.txt, data1.txt, and data2.txt are specified for
testing the different algorithms developed in "linear_equation_solutions.c" for
solving different types of linear equations.*/
/*the data.txt is used for test_SGE(), test_SEwME(), test_DLU(), test_mDLU().*/
/*the data1.txt is used for test_DLU_on_banded_matrix(), test_DLU_on_compressed_matrix(),
 and test_Crout_LU_decomposition().*/
/*the data2.txt is used for test_modified_Crout_DLU().*/
const char data_file_name[] = "C:\\Users\\zhoujianshan\\Documents\\My_C_Projects\\test1\\data2.txt";
const char data_file_name2[] = "C:\\Users\\zhoujianshan\\Documents\\My_C_Projects\\test1\\solution.txt";


/*statements on some function APIs*/
FILE *get_data_file_pointer(void);
int test1(void);
int test_SGE(void);
int test_SEwME(void);
int test_DLU(void);
int test_mDLU(void);
int test_DLU_on_banded_matrix(void);
int test2(void);
int test_DLU_on_compressed_matrix(void);
int test_Crout_LU_decomposition(void);
int test3(void);
int test_modified_Crout_DLU(void);


int main()
{
    /*test the read and write functions*/
    //test1();
    //test_SGE();
    //test_SEwME();
    //test_DLU();
    //test_mDLU();
    //test_DLU_on_banded_matrix();
    //test2();
    //test_DLU_on_compressed_matrix();
    //test_Crout_LU_decomposition();
    //test3();
    test_modified_Crout_DLU();

    return 0;
}

/*this function is developed to test the solution of a quasi-triple-diagonal
linear equations based on the modified Crout LU decomposition algorithm.*/
/*Note that the test data is from data2.txt.
The corresponding quasi_triple_diagonal system has a solution of
   0.162646675358540
   0.210397653194263
  -0.504237288135593
   0.806551499348110
   0.278031290743155*/
int test_modified_Crout_DLU(void)
{
    /*get the file pointer of the txt file*/
    FILE *filepoint = get_data_file_pointer();
    /*read the data from the txt file*/
    double array[3*ROW_NUM] = {0.0};/*state and initialize an array*/
    if(ROW_NUM<=3){
        printf("The row number of the coefficient matrix should be larger than 3!\n");
        exit(1);
    }
    double array_sr[2*(ROW_NUM-3)] = {0.0};/*used to record the intermediate results r, s*/
    double b[ROW_NUM] = {0.0};
    establish_one_quasi_array(filepoint,(double*)array,ROW_NUM,COL_NUM, (double*)b, R, S);
    printf("The compressed augmented matrix is:\n");
    int i,j;
    for(i=0;i<(3*ROW_NUM);i++){
        printf(" %lf\n",array[i]);
    }
    printf("The constant vector is:\n");
    for(i=0;i<ROW_NUM;i++)
    {
        printf("%lf\n",*(b+i));
    }
    /*initialize the solution array*/
    double solution[ROW_NUM] = {0.0};


    modified_Crout_LU_decomposition(array, array_sr, solution, b, ROW_NUM, COL_NUM, R, S);
    return 0;
}

int test_Crout_LU_decomposition(void)
{
    /*get the file pointer of the txt file*/
    FILE *filepoint = get_data_file_pointer();
    /*read the data from the txt file*/
    double array[3*ROW_NUM-2] = {0.0};/*state and initialize an array*/
    double b[ROW_NUM] = {0.0};
    establish_one_array(filepoint,(double*)array,ROW_NUM,COL_NUM, (double*)b, R, S);
    printf("The compressed augmented matrix is:\n");
    int i,j;
    for(i=0;i<(3*ROW_NUM-2);i++){
        printf(" %lf\n",array[i]);
    }
    printf("The constant vector is:\n");
    for(i=0;i<ROW_NUM;i++)
    {
        printf("%lf\n",*(b+i));
    }
    /*initialize the solution array*/
    double solution[ROW_NUM] = {0.0};
    Crout_LU_decomposition(array, solution, b, ROW_NUM, COL_NUM, R, S);

    return 0;
}


int test3(void)
{
    /*get the file pointer of the txt file*/
    FILE *filepoint = get_data_file_pointer();
    /*read the data from the txt file*/
    double array[3*ROW_NUM] = {0.0};/*state and initialize an array*/
    if(ROW_NUM<=3){
        printf("The row number of the coefficient matrix should be larger than 3!\n");
        exit(1);
    }
    double array_rs[2*(ROW_NUM-3)] = {0.0};/*used to record the intermediate results r, s*/
    double b[ROW_NUM] = {0.0};
    establish_one_quasi_array(filepoint,(double*)array,ROW_NUM,COL_NUM, (double*)b, R, S);
    printf("The compressed augmented matrix is:\n");
    int i,j;
    for(i=0;i<(3*ROW_NUM);i++){
        printf(" %lf\n",array[i]);
    }
    printf("The constant vector is:\n");
    for(i=0;i<ROW_NUM;i++)
    {
        printf("%lf\n",*(b+i));
    }
    /*initialize the solution array*/
    double solution[ROW_NUM] = {0.0};

    return 0;
}

int test_DLU_on_compressed_matrix(void)
{
    /*get the file pointer of the txt file*/
    FILE *filepoint = get_data_file_pointer();
    /*read the data from the txt file*/
    double array[R+S+1][ROW_NUM] = {0.0};/*state and initialize an array*/
    double b[ROW_NUM] = {0.0};
    establish_compressed_matrix(filepoint,(double**)array,ROW_NUM,COL_NUM, (double*)b, R, S);
    printf("The compressed augmented matrix is:\n");
    int i,j;
    for(i=0;i<(R+S+1);i++){
        for(j=0;j<ROW_NUM;j++){
            printf(" %lf",array[i][j]);
        }
        printf("\n");
    }
    printf("The constant vector is:\n");
    for(i=0;i<ROW_NUM;i++)
    {
        printf("%lf\n",*(b+i));
    }
    /*initialize the solution array*/
    double solution[ROW_NUM] = {0.0};

    Doolittle_LU_on_compressed_matrix((double**)array, solution, b, ROW_NUM, COL_NUM, S, R);
    return 0;
}


int test2(void)
{
    /*get the file pointer of the txt file*/
    FILE *filepoint = get_data_file_pointer();
    /*read the data from the txt file*/
    double array[R+S+1][ROW_NUM] = {0.0};/*state and initialize an array*/
    double b[ROW_NUM] = {0.0};
    establish_compressed_matrix(filepoint,(double**)array,ROW_NUM,COL_NUM, (double*)b, R, S);
    printf("The compressed augmented matrix is:\n");
    int i,j;
    for(i=0;i<(R+S+1);i++){
        for(j=0;j<ROW_NUM;j++){
            printf(" %lf",array[i][j]);
        }
        printf("\n");
    }
    printf("The constant vector is:\n");
    for(i=0;i<ROW_NUM;i++)
    {
        printf("%lf\n",*(b+i));
    }
    /*initialize the solution array*/
    double solution[ROW_NUM] = {0.0};

    return 0;
}

int test_DLU_on_banded_matrix(void)
{
    /*get the file pointer of the txt file*/
    FILE *filepoint = get_data_file_pointer();
    /*read the data from the txt file*/
    double array[ROW_NUM][COL_NUM] = {0.0};/*state and initialize an array*/
    read_data_from_txt(filepoint,(double**)array,ROW_NUM,COL_NUM);
    printf("The augmented matrix is:\n");
    int i,j;
    for(i=0;i<ROW_NUM;i++){
        for(j=0;j<COL_NUM;j++){
            printf(" %lf",array[i][j]);
        }
        printf("\n");
    }
    /*initialize the solution array*/
    double solution[ROW_NUM] = {0.0};

    Doolittle_LU_on_banded_matrix((double**)array,solution,ROW_NUM,COL_NUM, 1, 1);

    return 0;
}

int test_mDLU(void){
    /*get the file pointer of the txt file*/
    FILE *filepoint = get_data_file_pointer();
    /*read the data from the txt file*/
    double array[ROW_NUM][COL_NUM] = {0.0};/*state and initialize an array*/
    read_data_from_txt(filepoint,(double**)array,ROW_NUM,COL_NUM);
    printf("The augmented matrix is:\n");
    int i,j;
    for(i=0;i<ROW_NUM;i++){
        for(j=0;j<COL_NUM;j++){
            printf(" %lf",array[i][j]);
        }
        printf("\n");
    }
    /*initialize the solution array*/
    double solution[ROW_NUM] = {0.0}, s[ROW_NUM]={0.0};
    int M[ROW_NUM]={0};
    int ai;
    for(i=1;i<=ROW_NUM;i++)
    {
        ai = i-1;
        M[ai] = i;
    }
    modified_Doolittle_LU_decomposition((double**)array, solution, M, s, ROW_NUM, COL_NUM);

    return 0;
}

int test_DLU(void){
    /*get the file pointer of the txt file*/
    FILE *filepoint = get_data_file_pointer();
    /*read the data from the txt file*/
    double array[ROW_NUM][COL_NUM] = {0.0};/*state and initialize an array*/
    read_data_from_txt(filepoint,(double**)array,ROW_NUM,COL_NUM);
    printf("The augmented matrix is:\n");
    int i,j;
    for(i=0;i<ROW_NUM;i++){
        for(j=0;j<COL_NUM;j++){
            printf(" %lf",array[i][j]);
        }
        printf("\n");
    }
    /*initialize the solution array*/
    double solution[ROW_NUM] = {0.0};
    Doolittle_LU_decomposition((double**)array,solution,ROW_NUM,COL_NUM);

    return 0;
}

int test_SEwME(void){
    /*get the file pointer of the txt file*/
    FILE *filepoint = get_data_file_pointer();
    /*read the data from the txt file*/
    double array[ROW_NUM][COL_NUM] = {0.0};/*state and initialize an array*/
    read_data_from_txt(filepoint,(double**)array,ROW_NUM,COL_NUM);
    printf("The augmented matrix is:\n");
    int i,j;
    for(i=0;i<ROW_NUM;i++){
        for(j=0;j<COL_NUM;j++){
            printf(" %lf",array[i][j]);
        }
        printf("\n");
    }
    /*initialize the solution array*/
    double solution[ROW_NUM] = {0.0};
    Gauss_elimination_with_master_element((double**)array,solution,ROW_NUM,COL_NUM);

    return 0;
}

int test_SGE(void){
    /*get the file pointer of the txt file*/
    FILE *filepoint = get_data_file_pointer();
    /*read the data from the txt file*/
    double array[ROW_NUM][COL_NUM] = {0.0};/*state and initialize an array*/
    read_data_from_txt(filepoint,(double**)array,ROW_NUM,COL_NUM);
    printf("The augmented matrix is:\n");
    int i,j;
    for(i=0;i<ROW_NUM;i++){
        for(j=0;j<COL_NUM;j++){
            printf(" %lf",array[i][j]);
        }
        printf("\n");
    }
    /*initialize the solution array*/
    double solution[ROW_NUM] = {0.0};
    sequential_Gauss_elimination((double**)array,solution,ROW_NUM,COL_NUM);

    return 0;
}

int test1(void){
    /*get the file pointer of the txt file*/
    FILE *filepoint = get_data_file_pointer();
    /*read the data from the txt file*/
    double array[ROW_NUM][COL_NUM] = {0.0};/*state and initialize an array*/
    read_data_from_txt(filepoint,(double**)array,ROW_NUM,COL_NUM);
    printf("The augmented matrix is:\n");
    int i,j;
    for(i=0;i<ROW_NUM;i++){
        for(j=0;j<COL_NUM;j++){
            printf(" %lf",array[i][j]);
        }
        printf("\n");
    }

    /*get a pointer to a txt for writing data*/
    double solution[] = {1.11,2.22,3.33};
    int solution_num = sizeof(solution)/sizeof(solution[0]);
    FILE *pWrite;
    if((pWrite=fopen(data_file_name2,"w"))==NULL){
        printf("Errors occur in creating a txt file pointer!\n");
        exit(1);
    }
    /*write the data to a txt file*/
    write_data_to_txt(pWrite, solution, solution_num);
    return 0;
}

FILE *get_data_file_pointer(void){
    /*get the file pointer of the txt file*/
    FILE *filepoint;
    if((filepoint=fopen(data_file_name,"rb"))==NULL){
        printf("Cannot open the data file!\n");
        printf("Check the existence of the data file on: %s\n",data_file_name);
        exit(1);
    }else{
    return filepoint;
    }
}
