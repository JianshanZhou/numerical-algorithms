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


#include <stdio.h>
#include <stdlib.h>

int establish_one_quasi_array(FILE *filepoint, double *array, int rn, int cn, double *b, int r, int s)
{
    if((r!= 1)||(s!=1))
    {
        printf("the bandwidths r and s should be equal to 1!\n");
        exit(1);
    }

    double temp_value;
    int i=0, j=0, ai, aj, dataCount=0, flag = 0;
    //int flag1=0, flag2=0;
    int n=rn;

    while(1==fscanf(filepoint,"%lf",&temp_value)){

            i = dataCount/cn;
            j = dataCount%cn;

            //read A
            //printf("a[%d,%d]=%lf\n",i+1,j+1,temp_value);
            ai = i+1;
            aj = j+1;

            if(aj<=n)
            {
                if((ai+s)==aj)
                {
                    *(array+i)=temp_value;
                }
                else if(ai == aj)
                {
                    *(array+i+n)=temp_value;
                }
                else if((ai-r)==aj)
                {
                    *(array+i+2*n)=temp_value;
                }
                else if((ai==n)&&(aj==1))
                {
                    *(array+i)=temp_value;
                }
                else if((ai==1)&&(aj==n))
                {
                    *(array+n+n)=temp_value;
                }
            }

            //read b
            if(aj == cn)
            {
                *(b+i) = temp_value;
            }

            dataCount++;
    }
    //printf("total read %d data.",dataCount);
    fclose(filepoint);/*close the file pointer*/
    return 0;
}

int establish_one_array(FILE *filepoint, double *array, int rn, int cn, double *b, int r, int s)
{
    if((r!= 1)||(s!=1))
    {
        printf("the bandwidths r and s should be equal to 1!\n");
        exit(1);
    }

    double temp_value;
    int i=0, j=0, ai, aj, dataCount=0, flag = 0;
    //int flag1=0, flag2=0;
    int n=rn;

    while(1==fscanf(filepoint,"%lf",&temp_value)){

            i = dataCount/cn;
            j = dataCount%cn;

            //read A
            //printf("a[%d,%d]=%lf\n",i+1,j+1,temp_value);
            ai = i+1;
            aj = j+1;

            if(aj<=n)
            {
                if((ai+s)==aj)
                {
                    *(array+i)=temp_value;
                }
                else if(ai == aj)
                {
                    *(array+i+n-s)=temp_value;
                }
                else if((ai-r)==aj)
                {
                    *(array+i+n-s+n-r)=temp_value;
                }
            }

            //read b
            if(aj == cn)
            {
                *(b+i) = temp_value;
            }

            dataCount++;
    }
    //printf("total read %d data.",dataCount);
    fclose(filepoint);/*close the file pointer*/
    return 0;
}

int establish_compressed_matrix(FILE *filepoint, double ** array, int rn, int cn, double *b, int r, int s)
{
    double temp_value;
    int i=0, j=0, ai, aj, dataCount=0, x=0, y=0;
    //int flag1=0, flag2=0;
    int n=rn;

    while(1==fscanf(filepoint,"%lf",&temp_value)){

            i = dataCount/cn;
            j = dataCount%cn;

            //read A
            //printf("a[%d,%d]=%lf\n",i+1,j+1,temp_value);
            ai = i+1;
            aj = j+1;

            if(aj<=n)
            {
                if((aj>(s+1))&&(aj<=n)&&(ai>=1)&&(ai<(aj-s)))
                {
                    ;
                }
                else if((ai>(r+1))&&(ai<=n)&&(aj>=1)&&(aj<(ai-r)))
                {
                    ;
                }
                else
                {
                    x = (i-j)+s;
                    y = j;
                    //printf("a[%d,%d]=%lf\n",ai,aj,temp_value);
                    *((double*)array+n*x+y) = temp_value;
                }
            }
            //read b
            if(aj == cn)
            {
                *(b+i) = temp_value;
            }

            dataCount++;
    }
    //printf("total read %d data.",dataCount);
    fclose(filepoint);/*close the file pointer*/
    return 0;
}

int read_data_from_txt(FILE *filepoint, double **array, int rn, int cn)
{
    double temp_value;
    int i=0, j=0, dataCount=0;

    while(1==fscanf(filepoint,"%lf",&temp_value)){
            i = dataCount/cn;
            j = dataCount%cn;
            *((double*)array + cn*i +j) = temp_value;
            dataCount++;
    }
    fclose(filepoint);/*close the file pointer*/
    return 0;
}

int write_data_to_txt(FILE *filepointer, double *array, int n)
{
    int i=0;
    for(i=0;i<n;i++){
        fprintf(filepointer,"%lf\n",*(array+i));
        printf("Write %lf to the txt file!\n",*(array+i));
    }
    fclose(filepointer);/*close the file pointer*/
    printf("Successfully write the whole data to the txt file!\n");
    return 0;
}
