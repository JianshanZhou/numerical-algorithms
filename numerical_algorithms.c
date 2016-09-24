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

const double epsilon = 1.0e-6;


/*the modified Doolittle LU composition algorithm:
in this modification, the diagonal elements are allowed
not to follow the basic assumption. This version can handle
a more general computation situation based on selection of the
maximum positive element in each column as the master element.
Additional input parameters, i.e., M and s, are needed, one of
which is used to denote the index of each master element and the
other is used to record the intermediate computing results.
M is initialized containing the indexes of the initial master elements.*/
int modified_Doolittle_LU_composition(double **array, double *solution, int *M, double *s, int rn, int cn)
{
    if(rn!=(cn-1))
    {
        printf("The row number of the coefficient matrix is not equal to the column number!\n");
        exit(1);
    }
    int n=rn;

    /*step-1: do LU composition*/
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

/*the Doolittle LU composition algorithm:
the following function is developed based on
the well-know Doolittle LU composition algorithm
to solve a linear equations system. It should be
noted that this is a basic Doolittle algorithm, and it
requires that the diagonal elements cannot be zero.*/
int Doolittle_LU_composition(double **array, double *solution, int rn, int cn)
{
    if(rn!=(cn-1))
    {
        printf("The row number of the coefficient matrix is not equal to the column number!\n");
        exit(1);
    }
    int n=rn;
    /*step-1: do Doolittle composition to obtain the L and the U matrices*/
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
