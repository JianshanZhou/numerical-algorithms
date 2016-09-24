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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int max_index(double **array, int rn, int cn, int k);
int exchange_two_row(double ** array, int rn, int cn, int k, int ik);
double sum_LU(double **array, int rn, int cn, int x, int y, char flag);
double sum_LU_XY(double **array, int rn, int cn, double *solution, int i, char flag);
int max_index_in_s(double *s, int k, int n);
int max_of_two(int x, int y);
int max_of_three(int x, int y, int z);
int min_of_two(int x, int y);
int is_in_band(int x, int y, int m, int n);

const double epsilon = 1.0e-6;


/*the modified Crout LU decomposition algorithm is developed to solve a kind of
special linear equations in which the coefficient matrix is a quasi diagonal matrix.*/
int modified_Crout_LU_decomposition(double *array, double * array_sr, double *solution, double *b, int rn, int cn, int r, int s)
{
    if(rn!=(cn-1))
    {
        printf("The row number of the coefficient matrix is not equal to the column number!\n");
        exit(1);
    }
    if(r != s)
    {
        printf("It should be noted that the Crout LU decomposition algorithm is specified for the case s=r!\n");
        exit(1);
    }
    int n=rn;
    /*Note that:
    the 1 to n elements in the given array correspond to the c_1 to c_n;
    the n+1 to 2n ---> a_1 to a_n;
    the 2n+1 to 3n ---> d_1 to d_n;
    the 1 to n-3 in array_sr are used to record ---> r_2 to r_(n-2);
    the n-3+1 to 2*(n-3) ---> s_2 to s_(n-2)*/
    int i, ai, p_index, q_index, d_index, s_index, r_index;
    //the i and the ai indexes are used for original matrix A
    //p_index --> p and a
    //q_index --> q and c
    //d_index --> d

    /*do LU dedecomposition*/
    i = 1;
    ai = i-1;
    p_index = ai + n;
    *(array + p_index) = *(array + p_index);

    for(i=1;i<=(n-2);i++)
    {
        ai = i-1;
        q_index = ai;
        p_index = ai + n;
        d_index = ai + n+ n;

        *(array+q_index) = 1.0*(*(array+q_index))/(*(array+p_index));
        *(array+p_index+1) = *(array+p_index+1) - (*(array+d_index+1))*(*(array+q_index));
    }
    /*s_1 corresponds to d_1 at (1,n)*/
    i = 1;
    ai = i-1;
    d_index = ai+2*n;
    p_index = ai + n;
    *(array+d_index) = 1.0*(*(array+d_index))/(*(array+p_index));
    for(i=2;i<=n-2;i++)
    {
        ai = i-1;
        s_index = i-2;
        d_index = ai + 2*n;
        p_index = ai + n;
        if((s_index+1)==1)
        {
            *(array_sr + s_index) = -1.0*(*(array+d_index))*(*(array+2*n))/(*(array+p_index));
        }
        else
        {
            *(array_sr + s_index) = -1.0*(*(array+d_index))*(*(array_sr+s_index-1))/(*(array+p_index));
        }
    }
    i = n-1;
    ai = i-1;
    q_index = ai;
    d_index = ai + 2*n;
    s_index = (n-2) - 2;
    p_index = ai + n;
    *(array+q_index) = 1.0*(*(array+q_index)-(*(array+d_index))*(*(array_sr+s_index)))/(*(array+p_index));

    i = 1;
    ai = i-1;
    q_index = n-1;
    *(array+q_index) = *(array+q_index);
    for(i=2;i<=n-2;i++)
    {
        ai = i-1;
        r_index = i-2+(n-3);
        //q_index = n-1;
        if(i==2)
        {
            q_index = n-1;
            *(array_sr+r_index)=-1.0*(*(array+q_index))*(*(array+ai-1));
        }
        else
        {
            *(array_sr+r_index)=-1.0*(*(array_sr+r_index-1))*(*(array+ai-1));
        }
    }
    *(array+n-1+2*n) = *(array+n-1+2*n) - (*(array_sr+2*(n-3)-1))*(*(array+n-2-1));

    double temp_value=0.0;
    for(i=1;i<=n-1;i++)
    {
        if(i==1)
        {
            temp_value += (*(array+n-1))*(*(array+2*n));
        }
        else if(i==(n-1))
        {
            temp_value += (*(array+n-1+2*n))*(*(array+n-1-1));
        }
        else
        {
            temp_value += (*(array_sr+i-2))*(*(array_sr+i-2+(n-3)));
        }
    }
    *(array+n-1+n) = *(array+n-1+n) - temp_value;

    /*solve Ly=b and Ux=y*/
    *(solution+1-1) = 1.0*(*(b+1-1))/(*(array+1-1+n));
    for(i=2;i<=n-1;i++)
    {
        ai = i-1;
        d_index = ai + 2*n;
        p_index = ai + n;
        *(solution+ai) = 1.0*(*(b+ai)-(*(array+d_index))*(*(solution+ai-1)))/(*(array+p_index));
    }
    temp_value = 0.0;
    for(i=1;i<=n-1;i++)
    {
        if(i==1)
        {
            temp_value += (*(array+n-1))*(*(solution+i-1));
        }
        else if(i==(n-1))
        {
            temp_value += (*(array+n-1+2*n))*(*(solution+i-1));
        }
        else
        {
            temp_value += (*(array_sr+n-3+i-2))*(*(solution+i-1));
        }
    }
    *(solution+n-1) = 1.0*(*(b+n-1)-temp_value)/(*(array+n-1+n));

    *(solution+n-1) = *(solution+n-1);
    *(solution+n-1-1) = *(solution+n-1-1) - *(array+n-1-1)*(*(solution+n-1));
    for(i=n-2;i>=1;i--)
    {
        if(i==1)
        {
            *(solution+i-1) = *(solution+i-1)-(*(array+i-1))*(*(solution+i-1+1))-(*(array+i-1+2*n))*(*(solution+n-1));
        }
        else
        {
            *(solution+i-1) = *(solution+i-1)-(*(array+i-1))*(*(solution+i-1+1))-(*(array_sr+i-2))*(*(solution+n-1));
        }
    }


    printf("Successfully achieve the algorithm! The solution is:\n");
    for(i=1; i<=n; i++)
    {
        ai = i-1;
        printf("x[%d]=%lf\n",i,*(solution+ai));
    }

    return 0;
}

/*the Crout LU decomposition algorithm is developed to solve a kind of specified
linear equations, i.e., the so-called tri-diagonal linear equations.*/
int Crout_LU_decomposition(double *array, double *solution, double *b, int rn, int cn, int r, int s)
{
    if(rn!=(cn-1))
    {
        printf("The row number of the coefficient matrix is not equal to the column number!\n");
        exit(1);
    }
    if(r != s)
    {
        printf("It should be noted that the Crout LU decomposition algorithm is specified for the case s=r!\n");
        exit(1);
    }
    int n=rn;
    /*the 1~(n-1) elements of the given array correspond to c_1~c_(n-1) (q_1~q_(n-1))*/
    /*the n~(2n-1) correspond to a_1~a_n (p_1~p_n)*/
    /*the 2n~(3n-2) correspond to d_2~d_n (r_2~r_n)*/

    /*LU decomposition*/
    int i, ai;//i ranges from 1 to 3n-2
    int p_index, d_index, q_index;
    p_index = n-1;
    *(array+p_index) = *(array+p_index);

    for(i=1;i<=n-1;i++)
    {
        ai = i-1;
        p_index = ai + (n-1);
        d_index = (ai) + (n-1+n);
        *(array+ai) = 1.0*(*(array+ai))/(*(array + p_index));
        *(array+p_index+1)=*(array+p_index+1)-(*(array+d_index))*(*(array+ai));
    }

    i=1;
    ai = i-1;
    p_index = ai + n-1;
    *(solution+ai) = 1.0*(*(b+ai))/(*(array+p_index));
    for(i=2;i<=n;i++)
    {
        ai = i-1;
        p_index = ai + n-1;
        d_index = ai + (n-1+n) - 1;
        *(solution+ai) = 1.0*(*(b+ai)-(*(array+d_index))*(*(solution+ai-1)))/(*(array+p_index));
    }
    i=n;
    ai=i-1;
    *(solution+ai) = *(solution+ai);
    for(i=n-1;i>=1;i--)
    {
        ai=i-1;
        *(solution+ai) = *(solution+ai) - (*(array+ai))*(*(solution+ai+1));
    }

    printf("Successfully achieve the algorithm! The solution is:\n");
    for(i=1; i<=n; i++)
    {
        ai = i-1;
        printf("x[%d]=%lf\n",i,*(solution+ai));
    }

    return 0;
}


/*the Doolittle LU decomposition algorithm with the compressed coefficient matrix:
the input array denotes the compressed matrix C related to the original A,
and an additional constant vector b is also provided. It should be noted that
the array is not the augmented matrix here any longer.*/
int Doolittle_LU_on_compressed_matrix(double **C, double *solution, double *b, int rn, int cn, int s, int r)
{
    if(rn!=(cn-1))
    {
        printf("The row number of the coefficient matrix is not equal to the column number!\n");
        exit(1);
    }
    int n=rn, m=s+r+1;

    int i, ai, i_upper, j, aj, j_upper, k, ak, t, at, t_lower;
    int x, y;
    double temp_sum = 0.0;

    /*do LU decomposition*/
    for(k=1;k<=n;k++)
    {
        ak = k-1;
        j_upper = min_of_two(k+s,n);
        for(j=k;j<=j_upper;j++)
        {
            aj = j-1;
            if(k==1){
                ;
            }else{
                t_lower = max_of_three(1,k-r,j-s);
                temp_sum = 0.0;
                for(t=t_lower;t<=(k-1);t++)
                {
                    at = t-1;
                    x = (k-t+s+1)-1;
                    y = (t-j+s+1)-1;

                    temp_sum += (*((double*)C+n*x+at))*(*((double*)C+n*y+aj));
                }
                x = (k-j+s+1)-1;

                *((double*)C+n*x+aj) = *((double*)C+n*x+aj) - temp_sum;
            }
        }

        if(k<n)
        {
            i_upper = min_of_two(k+r,n);
            for(i=k+1;i<=i_upper;i++)
            {
                ai = i-1;
                if(k==1)
                {
                    x=(i-k+s+1)-1;
                    y=s+1-1;

                    *((double*)C+n*x+ak) = 1.0*(*((double*)C+n*x+ak))/(*((double*)C+n*y+ak));
                }
                else
                {
                    t_lower = max_of_three(1,i-r,k-s);
                    temp_sum = 0.0;
                    for(t=t_lower;t<=(k-1);t++)
                    {
                        at = t-1;
                        x = (i-t+s+1)-1;
                        y = (t-k+s+1)-1;

                        temp_sum += (*((double*)C+n*x+at))*(*((double*)C+n*y+ak));
                    }
                    x=(i-k+s+1)-1;
                    y=(s+1)-1;

                    *((double*)C+n*x+ak) = 1.0*(*((double*)C+n*x+ak)-temp_sum)/(*((double*)C+n*y+ak));
                }
            }
        }
    }

    /*solve Ly=b and Ux=y*/
    for(i=2;i<=n;i++)
    {
        ai = i-1;
        t_lower = max_of_two(1,i-r);
        temp_sum = 0.0;
        for(t=t_lower;t<=(i-1);t++)
        {
            at = t-1;
            x = (i-t+s+1)-1;

            temp_sum += (*((double*)C+n*x+at))*(*(b+at));
        }
        *(b+ai) = *(b+ai)-temp_sum;
    }

    i = n;
    ai = i-1;
    y = (s+1)-1;

    *(solution+ai) = (*(b+ai))/(*((double*)C+n*y+ai));
    for(i=n-1;i>=1;i--)
    {
        ai = i-1;
        t_lower = min_of_two(i+s,n);
        temp_sum = 0.0;

        for(t=i+1;t<=t_lower;t++)
        {
            at = t-1;
            x = (i-t+s+1)-1;

            temp_sum += (*((double*)C+n*x+at))*(*(solution+at));
        }
        y = (s+1)-1;
        *(solution+ai) = 1.0*(*(b+ai)-temp_sum)/(*((double*)C+n*y+ai));
    }

    printf("Successfully achieve the algorithm! The solution is:\n");
    for(i=1; i<=n; i++)
    {
        ai = i-1;
        printf("x[%d]=%lf\n",i,*(solution+ai));
    }

    return 0;
}

/*to see whether the (x,y) indexes is in the C*/
int is_in_band(int x, int y, int m, int n)
{
    int ax=x+1,ay=y+1;
    return (ax>=1)&&(ax<=m)&&(ay>=1)&&(ay<=n);
}

/*the Doolittle LU decomposition algorithm for specified
linear equations whose coefficient matrix is banded.
This algorithm requires that the elements located at the banded
region of the coefficient matrix cannot be zero.*/
int Doolittle_LU_on_banded_matrix(double **array, double *solution, int rn, int cn, int s, int r)
{
    if(rn!=(cn-1))
    {
        printf("The row number of the coefficient matrix is not equal to the column number!\n");
        exit(1);
    }
    int n=rn;

    int k,ak, i, ai, j, aj, t, at, j_upper, i_upper, t_lower;
    double temp_LU_sum = 0.0;

    /*do Lu decomposition*/
    for(k=1;k<=n;k++)
    {
        ak = k-1;
        j_upper = min_of_two(k+s,n);
        for(j=k;j<=j_upper;j++)
        {
            aj = j-1;
            temp_LU_sum = 0.0;
            t_lower = max_of_three(1,k-r,j-s);
            for(t=t_lower;t<=k-1;t++)
            {
                at = t-1;
                temp_LU_sum += (*((double*)array+cn*ak+at))*(*((double*)array+cn*at+aj));
            }
            *((double*)array+cn*ak+aj) = *((double*)array+cn*ak+aj) - temp_LU_sum;
        }
        if(k<n)
        {
            i_upper = min_of_two(k+r,n);
            for(i=k+1;i<=i_upper;i++)
            {
                ai = i-1;
                temp_LU_sum = 0.0;
                t_lower = max_of_three(1,i-r,k-s);
                for(t=t_lower;t<=k-1;t++)
                {
                    at = t-1;
                    temp_LU_sum += (*((double*)array+cn*ai+at))*(*((double*)array+cn*at+ak));
                }
                *((double*)array+cn*ai+ak) = 1.0*(*((double*)array+cn*ai+ak)-temp_LU_sum)/(*((double*)array+cn*ak+ak));
            }
        }
    }

    /*solve Ly=b and Ux=y*/
    i = 1;
    ai = i-1;
    *(solution+ai) = *((double*)array+cn*ai+n);
    for(i=2;i<=n;i++){
        ai = i-1;
        temp_LU_sum = 0.0;
        t_lower = max_of_two(1,i-r);
        for(t=t_lower;t<=i-1;t++)
        {
            at = t-1;
            temp_LU_sum += (*((double*)array+cn*ai+at))*(*(solution+at));
        }
        *(solution+ai) = *((double*)array+cn*ai+n) - temp_LU_sum;
    }

    i = n;
    ai = i-1;
    *(solution+ai) = 1.0*(*(solution+ai))/(*((double*)array+cn*ai+ai));
    for(i=n-1;i>=1;i--)
    {
        ai = i-1;
        temp_LU_sum = 0.0;
        t_lower = min_of_two(i+s,n);
        for(t=i+1;t<=t_lower;t++)
        {
            at = t-1;
            temp_LU_sum += 1.0*(*((double*)array+cn*ai+at))*(*(solution+at));
        }
        *(solution+ai) = 1.0*(*(solution+ai)-temp_LU_sum)/(*((double*)array+cn*ai+ai));
    }


    printf("Successfully achieve the algorithm! The solution is:\n");
    for(i=1; i<=n; i++)
    {
        ai = i-1;
        printf("x[%d]=%lf\n",i,*(solution+ai));
    }

    return 0;
}

/*this function outputs the minimum one among two given elements.*/
int min_of_two(int x, int y)
{
    return ((x>y)?y:x);
}

/*this function outputs the maximum one among two given elements.*/
int max_of_two(int x, int y)
{
    return ((x>y)?x:y);
}

/*this function outputs the maximum one among three given elements.*/
int max_of_three(int x, int y, int z)
{
    return max_of_two(max_of_two(x,y), z);
}

/*the modified Doolittle LU decomposition algorithm:
in this modification, the diagonal elements are allowed
not to follow the basic assumption. This version can handle
a more general computation situation based on selection of the
maximum positive element in each column as the master element.
Additional input parameters, i.e., M and s, are needed, one of
which is used to denote the index of each master element and the
other is used to record the intermediate computing results.
M is initialized containing the indexes of the initial master elements.*/
int modified_Doolittle_LU_decomposition(double **array, double *solution, int *M, double *s, int rn, int cn)
{
    if(rn!=(cn-1))
    {
        printf("The row number of the coefficient matrix is not equal to the column number!\n");
        exit(1);
    }
    int n=rn;

    /*step-1: do LU decomposition*/
    int i, j, k, ai, aj, ak, ik, a_ik, t, at;
    double temp_value = 0.0;
    for(k=1;k<=n;k++)
    {
        ak = k-1;
        /*substep-1: calculate s_i for i=k...n*/
        for(i=k;i<=n;i++)
        {
            ai = i-1;
            /*update the intermediate results*/
            if(k==1)
            {
                *(s+ai) = *((double*)array+cn*ai+ak);
            }
            else
            {
                *(s+ai) = *((double*)array+cn*ai+ak) - sum_LU(array, rn, cn, i, k, 'l');
            }
        }
        /*substep-2: select the maximum positive master element*/
        ik = max_index_in_s(s, k, n);
        a_ik = ik-1;
        *(M+ak) = ik;
        /*substep-3: exchange the rows related to the previous and the current master elements*/
        if(ik != k)
        {
            for(t=1;t<=k-1;t++)
            {
                at = t-1;
                temp_value = *((double*)array+cn*ak+at);
                *((double*)array+cn*ak+at) = *((double*)array+cn*a_ik+at);
                *((double*)array+cn*a_ik+at) = temp_value;
            }
            for(t=k;t<=n;t++)
            {
                at = t-1;
                temp_value = *((double*)array+cn*ak+at);
                *((double*)array+cn*ak+at) =*((double*)array+cn*a_ik+at);
                *((double*)array+cn*a_ik+at) = temp_value;
            }
            temp_value = *(s+ak);
            *(s+ak) = *(s+a_ik);
            *(s+a_ik) = temp_value;
        }
        /*substep-4: calculate u_kk, u_kj, l_ik*/
        *((double*)array+cn*ak+ak) = *(s+ak);
        if(k<n)
        {
            for(j=k+1;j<=n;j++)
            {
                aj = j-1;
                if(k==1)
                {
                    *((double*)array+cn*ak+aj) = *((double*)array+cn*ak+aj);
                }
                else
                {
                    *((double*)array+cn*ak+aj) = *((double*)array+cn*ak+aj) - sum_LU(array, rn, cn, k, j, 'u');
                }
            }
            for(i=k+1;i<=n;i++)
            {
                ai = i-1;
                *((double*)array+cn*ai+ak) = 1.0*(*(s+ai))/(*((double*)array+cn*ak+ak));
            }
        }
    }

    /*step-2: do solve Qb*/
    for(k=1;k<=n-1;k++)
    {
        ak = k-1;
        t = *(M+ak);
        at = t-1;
        temp_value = *((double*)array+cn*ak+(cn-1));
        *((double*)array+cn*ak+(cn-1)) = *((double*)array+cn*at+(cn-1));
        *((double*)array+cn*at+(cn-1)) = temp_value;
    }
    i=1;
    ai = i-1;
    *(solution+ai) = *((double*)array+cn*ai+(cn-1));
    for(i=2;i<=n;i++)
    {
        ai = i-1;
        *(solution+ai) = *((double*)array+cn*ai+(cn-1)) - sum_LU_XY(array, rn, cn, solution, i, 'y');
    }
    i = n;
    ai = i-1;
    *(solution+ai) = 1.0*(*(solution+ai))/(*((double*)array+cn*ai+ai));
    for(i=n-1;i>=1;i--)
    {
        ai = i-1;
        *(solution+ai) = 1.0*(*(solution+ai)-sum_LU_XY(array, rn, cn, solution, i, 'x'))/(*((double*)array+cn*ai+ai));
    }

    printf("Successfully achieve the algorithm! The solution is:\n");
    for(i=1; i<=n; i++)
    {
        ai = i-1;
        printf("x[%d]=%lf\n",i,*(solution+ai));
    }

    return 0;
}

/*this function calculates the index of the maximum positive
mater element given k and n.*/
int max_index_in_s(double *s, int k, int n)
{
    int ak = k-1;
    double temp_max_value=fabs(*(s+ak));
    int ik = k, i, ai;
    for(i=k;i<=n;i++)
    {
        ai = i-1;
        if(fabs(*(s+ai))>temp_max_value)
        {
            temp_max_value = fabs(*(s+ai));
            ik = i;
        }
    }
    return ik;
}

/*the Doolittle LU decomposition algorithm:
the following function is developed based on
the well-know Doolittle LU decomposition algorithm
to solve a linear equations system. It should be
noted that this is a basic Doolittle algorithm, and it
requires that the diagonal elements cannot be zero.*/
int Doolittle_LU_decomposition(double **array, double *solution, int rn, int cn)
{
    if(rn!=(cn-1))
    {
        printf("The row number of the coefficient matrix is not equal to the column number!\n");
        exit(1);
    }
    int n=rn;
    /*step-1: do Doolittle decomposition to obtain the L and the U matrices*/
    int k, ak, i, ai, j, aj;
    for(k=1;k<=n;k++){
        ak = k-1;
        if(k==1){
            if(k<n)
            {
                for(i=k+1;i<=n;i++){
                    ai = i-1;
                    if(fabs((*((double*)array+cn*ak+ak)))<epsilon){
                        printf("u_kk is approximatively zero!\n");
                        printf("The basic Doolittle algorithm cannot continue!\n");
                        exit(1);
                    }
                    *((double*)array+cn*ai+ak) = 1.0*(*((double*)array+cn*ai+ak))/(*((double*)array+cn*ak+ak));
                }
            }
        }else{
            for(j=k;j<=n;j++){
                aj = j-1;
                *((double*)array+cn*ak+aj) = *((double*)array+cn*ak+aj) - sum_LU(array, rn, cn, k, j, 'u');
            }
            if(k<n)
            {
                for(i=k+1;i<=n;i++)
                {
                    ai = i-1;
                    *((double*)array+cn*ai+ak) = 1.0*(*((double*)array+cn*ai+ak)-sum_LU(array, rn, cn, i, k, 'l'))/(*((double*)array+cn*ak+ak));
                }
            }
        }
    }
    /*step-2: derive the solution x based on Ly=b and Ux=y*/
    /*solve the y*/
    i=1;
    ai = i-1;
    *(solution + ai) = *((double*)array+cn*ai+(cn-1));
    for(i=2;i<=n;i++)
    {
        ai = i-1;
        *(solution + ai) = *((double*)array+cn*ai+(cn-1)) - sum_LU_XY(array, rn, cn, solution, i, 'y');
    }
    /*solve the x*/
    i = n;
    ai = i-1;
    *(solution + ai) = 1.0*(*(solution + ai))/(*((double*)array+cn*ai+(ai)));
    for(i=n-1;i>=1;i--)
    {
        ai = i-1;
        *(solution + ai) = 1.0*(*(solution + ai)-sum_LU_XY(array, rn, cn, solution, i, 'x'))/(*((double*)array+cn*ai+(ai)));
    }

    printf("Successfully achieve the algorithm! The solution is:\n");
    for(i=1; i<=n; i++)
    {
        ai = i-1;
        printf("x[%d]=%lf\n",i,*(solution+ai));
    }
    return 0;

}

/*this function is used to calculate the sum of a series of
products between l_xt and u_ty.*/
double sum_LU(double **array, int rn, int cn, int x, int y, char flag)
{
    double temp_sum = 0.0;
    int ax = x-1, ay = y-1, t, at, k;
    if(flag == 'u')
    {
        k = x;
    }
    else if(flag == 'l')
    {
        k = y;
    }
    else{
        printf("An error occurs in flag!\n");
        exit(1);
    }
    for(t=1;t<=(k-1);t++){
        at = t-1;
        temp_sum += (*((double*)array+cn*ax+at))*(*((double*)array+cn*at+ay));
    }
    return temp_sum;
}

/*this function calculates the sum of the products between the elements
l_it and yt or between u_it and xt.*/
double sum_LU_XY(double **array, int rn, int cn, double *solution, int i, char flag)
{
    double temp_sum = 0.0;
    int ai=i-1, t, at;
    if(flag == 'y')
    {
        for(t=1;t<=(i-1);t++){
            at = t-1;
            temp_sum += (*((double*)array+cn*ai+at))*(*(solution+at));
        }
    }
    else if(flag == 'x')
    {
        for(t=i+1;t<=rn;t++)
        {
            at = t-1;
            temp_sum += (*((double*)array+cn*ai+at))*(*(solution+at));
        }
    }
    else
    {
        printf("An error occurs in the flag setting!\n");
        exit(1);
    }
    return temp_sum;
}



/*the sequential Gauss elimination algorithm:
this algorithm requires that each diagonal element at
each iteration cannot be zero.*/
int sequential_Gauss_elimination(double **array, double *solution, int rn, int cn)
{
    if(rn!=(cn-1))
    {
        printf("The row number of the coefficient matrix is not equal to the column number!\n");
        exit(1);
    }
    int n=rn;

    /*step-1: eliminate elements of the given augmented matrix array*/
    int k, i, j, ak, ai, aj;
    double m_ik=0.0;
    for(k=1; k<=(n-1); k++)
    {
        /*calculate the actual index k in C*/
        ak = k-1;

        if(fabs(*((double*)array+cn*ak+ak))<=epsilon)
        {
            printf("The diagonal element(%d,%d) is approximately equal to zero!\n The algorithm cannot continue!\n",k,k);
            exit(1);
        }
        else
        {
            for(i=(k+1); i<=n; i++)
            {
                /*calculate the actual index i in C*/
                ai = i-1;

                /*update the element values*/
                m_ik=(1.0*(*((double*)array+cn*ai+ak)))/(*((double*)array+cn*(ak)+ak));
                for(j=k+1; j<=n; j++)
                {
                    /*calculate the actual index j in C*/
                    aj = j-1;

                    *((double*)array+cn*ai+aj) = *((double*)array+cn*ai+aj) - m_ik*(*((double*)array+cn*ak+aj));
                }
                *((double*)array+cn*ai+(n)) = *((double*)array+cn*ai+(n)) - m_ik*(*((double*)array+cn*ak+(n)));
            }
        }
    }

    /*step-2: re-iteration*/
    double sum_temp=0.0;
    for(k=n; k>=1; k--)
    {
        ak = k-1;
        if(k==n)
        {
            *(solution + ak) = 1.0*(*((double*)array+cn*ak+n))/(*((double*)array+cn*ak+ak));
        }
        else
        {
            sum_temp = 0.0;
            for(j=k+1; j<=n; j++)
            {
                aj = j-1;
                sum_temp+=(*((double*)array+cn*ak+aj))*(*(solution + aj));
            }
            *(solution + ak) = 1.0*(*((double*)array+cn*ak+n)-sum_temp)/(*((double*)array+cn*ak+ak));
        }
    }

    printf("Successfully achieve the algorithm! The solution is:\n");
    for(i=1; i<=n; i++)
    {
        ai = i-1;
        printf("x[%d]=%lf\n",i,*(solution+ai));
    }
    return 0;
}


/*the sequential Gauss elimination algorithm based on the master nonnegative element:
this algorithm needs not to assume that each diagonal element at
each iteration cannot be zero.*/
int Gauss_elimination_with_master_element(double **array, double *solution, int rn, int cn)
{

    if(rn!=(cn-1))
    {
        printf("The row number of the coefficient matrix is not equal to the column number!\n");
        exit(1);
    }

    int n=rn;

    /*step-1: eliminate elements of the given augmented matrix array*/
    int k, i, j, ak, ai, aj, ik;
    double m_ik=0.0;
    for(k=1; k<=(n-1); k++)
    {
        /*calculate the actual index k in C*/
        ak = k-1;
        /*find the index of the maximum element in the k-th column*/
        ik = max_index(array, rn, cn, k);
        /*exchange the k-th and the ik-th rows*/
        exchange_two_row(array, rn, cn, k, ik);

        for(i=(k+1); i<=n; i++)
        {
            /*calculate the actual index i in C*/
            ai = i-1;

            /*update the element values*/
            m_ik=(1.0*(*((double*)array+cn*ai+ak)))/(*((double*)array+cn*(ak)+ak));
            for(j=k+1; j<=n; j++)
            {
                /*calculate the actual index j in C*/
                aj = j-1;

                *((double*)array+cn*ai+aj) = *((double*)array+cn*ai+aj) - m_ik*(*((double*)array+cn*ak+aj));
            }
            *((double*)array+cn*ai+(n)) = *((double*)array+cn*ai+(n)) - m_ik*(*((double*)array+cn*ak+(n)));
        }
    }

    /*step-2: re-iteration*/
    double sum_temp=0.0;
    for(k=n; k>=1; k--)
    {
        ak = k-1;
        if(k==n)
        {
            *(solution + ak) = 1.0*(*((double*)array+cn*ak+n))/(*((double*)array+cn*ak+ak));
        }
        else
        {
            sum_temp = 0.0;
            for(j=k+1; j<=n; j++)
            {
                aj = j-1;
                sum_temp+=(*((double*)array+cn*ak+aj))*(*(solution + aj));
            }
            *(solution + ak) = 1.0*(*((double*)array+cn*ak+n)-sum_temp)/(*((double*)array+cn*ak+ak));
        }
    }

    printf("Successfully achieve the algorithm! The solution is:\n");
    for(i=1; i<=n; i++)
    {
        ai = i-1;
        printf("x[%d]=%lf\n",i,*(solution+ai));
    }
    return 0;
}


/*this function is used to find the index of the maximum element in a column vector:
k is the given index of the column of the given array*/
int max_index(double **array, int rn, int cn, int k)
{
    int i, ai, ak, ik;
    i = k;
    ai = i-1;
    ak = k-1;
    ik = i;

    double temp_max_element = fabs(*((double*)array+cn*ai+ak)), temp_value;
    for(i=k;i<=rn;i++){
        ai = i-1;
        temp_value = fabs(*((double*)array+cn*ai+ak));
        if(temp_max_element<temp_value){
            temp_max_element = temp_value;
            ik = i;
        }
    }
    return ik;/*output the index of the maximum element in the k-th column*/
}


/*this function is used to respectively exchange the corresponding elements in
the given k-th and the ik-th rows of the array.*/
int exchange_two_row(double ** array, int rn, int cn, int k, int ik)
{
    int ak = k-1, a_ik = ik-1, j, aj;
    int n=cn-1;
    double temp_value = 0.0;
    for(j=k;j<=rn;j++){
        aj = j-1;
        temp_value = *((double*)array + cn*ak + aj);
        *((double*)array + cn*ak + aj) =*((double*)array + cn*a_ik + aj);
        *((double*)array + cn*a_ik + aj) = temp_value;
    }
    temp_value = *((double*)array + cn*ak + n);
    *((double*)array + cn*ak + n) = *((double*)array + cn*a_ik + n);
    *((double*)array + cn*a_ik + n) = temp_value;
    return 0;
}
