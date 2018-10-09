#include <iostream>
#include <math.h>
#include <boost/math/special_functions/digamma.hpp>
using namespace std;
long double l1(int x, int y, double a0, double a1, double a2, int k, int m, double adjj = 0)
{
	return(exp(log(tgamma(a1 + k))- log(tgamma(k+1))- log(tgamma(a1)) + log(tgamma(x + y + a0 -m -k)) - log(tgamma(x -k +1)) - log(tgamma(a0 + y - m)) 
                                   + log(tgamma(m + a2)) - log(tgamma(m+1)) - log(tgamma(a2)) + log(tgamma(y +a0 -m)) - log(tgamma(y -m +1)) - log(tgamma(a0)) - adjj));
}
long double l1_c (double t1, double t2, int k, int m, double adjj = 0)
{
	return(exp(k *log(t1) + m *log(t2) - adjj));
}
long double l1_AC (double t1, double t2, double x, double y, double a0, double a1, double a2, int k, int m, double adjj = 0) 
{
	return(exp(log(tgamma(a1 + k))- log(tgamma(k+1))- log(tgamma(a1)) + log(tgamma(x + y + a0 -m -k)) - log(tgamma(x -k +1)) - log(tgamma(a0 + y - m)) 
                                   + log(tgamma(m + a2)) - log(tgamma(m+1)) - log(tgamma(a2)) + log(tgamma(y +a0 -m)) - log(tgamma(y -m +1)) - log(tgamma(a0)) + k *log(t1) + m *log(t2) - adjj));
}
long double R0_E1(int x, int y, int k, int m, double a0)
{
	return(x - k + y - m + a0);
} 
long double log_R0_E1(int x, int y, int k, int m, double a0)
{
	return(digamma(x - k + y - m + a0));
} 
long double R1_E1(int k, double a1) 
{
	return(k + a1);
}
long double log_R1_E1(int k, double a1) 
{
	return(digamma(k + a1));
}
long double R2_E1(int m, double a2) 
{
	return(m + a2);
}
long double log_R2_E1(int m, double a2) 
{
	return(digamma(m + a2));
}
double dBvZINB4_Expt(int x, int y, double a0, double a1, double a2, double b1, double b2, double p1, double p2, double p3, double p4) 
{
	double t1 = (b1 + b2 + 1) /(b1 + 1);
	double t2 = (b1 + b2 + 1) /(b2 + 1);
	double adj_A = 0;
	double adj_B1 = 0;
	double adj_C = 0;
	double adj_sum = 0;
	long double l1_B;
	long double l1_b = - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1 + b1) - a2 * log(1 + b2);
	long double l2_B = exp(- (x + a0 + a1)*log(1 + b1) + x * log(b1) + adj_B1) * p2 * y==0?1:0;
	long double l3_B = exp(- (y + a0 + a2)*log(1 + b2) + y * log(b2) + adj_B1) * p3 * x==0?1:0;
	long double l4_B = p4 * (x+y==0)?1:0 * exp(adj_B1);
	long double l_A_mat[x+1][y+1];
	long double l_C_mat[x+1][y+1];
	int sum_AC,sum_A_mat,sum_C_mat;
	for(int i = 0;i <= x;i++)
	{
		for(int j = 0;j <= y;j++)
		{
			l_A_mat[i][j] = l1(x,y,a0,a1,a2,i,j,adj_A);
			l_C_mat[i][j] = l1_c(t1,t2,i,j,adj_C);
			sum_A_mat += l_A_mat[i][j];
			sum_C_mat += l_C_mat[i][j];
			sum_AC += l_A_mat[i][j]*l_C_mat[i][j];
		}
	}
	if (l1_B < -200 && log(l2_B + l3_B + l4_B) < 0)
	{
		adj_B1 = ((-l1_B - 200) / 100) * 100; // prevent exp(l1_B) from being 0
	    l1_B = l1_B + adj_B1;
	}
	l1_B = exp(l1_B) * p1;
	while (log(sum_A_mat) > 250) 
	{
		adj_A = adj_A + 200;
		for(int i = 0;i <= x;i++)
		{
			for(int j = 0;j <= y;j++)
			{
				l_A_mat[i][j] = l1(x,y,a0,a1,a2,i,j,adj_A);  	
			}
	    }
	}
	while (log(sum_C_mat) > 250) 
	{
		adj_C = adj_C + 200;
		for(int i = 0;i <= x;i++)
		{
			for(int j = 0;j <= y;j++)
			{
				l_C_mat[i][j] = l1_c(t1,t2,i,j,adj_C);  	
			}
	    }
	}
	if(log(sum_AC) > 200)
	{
		adj_A = adj_A + 100;
    	adj_C = adj_C + 100;
    	for(int i = 0;i <= x;i++)
		{
			for(int j = 0;j <= y;j++)
			{
				l_A_mat[i][j] = l1(x,y,a0,a1,a2,i,j,adj_A);  
				l_C_mat[i][j] = l1_c(t1,t2,i,j,adj_C);  
				sum_AC += l_A_mat[i][j]*l_C_mat[i][j];	
			}
	    }
	}
	else if(log(sum_AC) < - 100)
	{
		adj_A = adj_A - 200;
    	adj_C = adj_C - 200;
    	for(int i = 0;i <= x;i++)
		{
			for(int j = 0;j <= y;j++)
			{
				l_A_mat[i][j] = l1(x,y,a0,a1,a2,i,j,adj_A);  
				l_C_mat[i][j] = l1_c(t1,t2,i,j,adj_C);  
				l_AC_mat[i][j] = l1_AC(t1,t2,x,y,a0,a1,a2,i,j,adj_A+adj_C);
				sum_AC += l_AC_mat[i][j];
				sum_A +=  l_A_mat[i][j]; 
			}
	    }
	}
	l_sum <- sum_AC * l1_B + sum_A * (l2_B +  l3_B +  l4_B) * exp(-adj_C);
	if (l_sum == 0) 
	{
	    adj_sum = -floor(log(sum_AC)*2/3 + log(l1_B)*2/3);
	    l_sum = sum_AC * exp(adj_sum) * l1_B + sum_A * (exp(adj_sum) * (l2_B +  l3_B +  l4_B)) * exp(-adj_C);
  	}
	// expectation components
	  R0_E1_B = b1/(1 + b1 + b2);
	  R0_E2_B = b1/(1 + b1);
	  R0_E3_B = b1/(1 + b2);
	  R0_E4_B = b1;
	  
	  R1_E1_B = b1/(1 + b1);
	  R1_E2_B = b1/(1 + b1);
	  R1_E3_B = b1;
	  R1_E4_B = b1;
	  
	  R2_E1_B = b1/(1 + b2); 
	  R2_E2_B = b1;
	  R2_E3_B = b1/(1 + b2);
	  R2_E4_B = b1;
	  for(i = 0;i <= x;i++)
	  {
	  	for(j = 0;j <= y;j++)
	  	{
	  		R0_mat[i][j] = R0_E1(x,y,i,j,a0);
	  		R0_mat[i][j] = R0_mat[i][j]*l_A_mat[i][j];
			R0_E += R0_mat[i][j] * l_C_mat[i][j] * exp(adj_sum) * l1_B * R0_E1_B + R0_mat[i][j] * (l2_B * R0_E2_B + l3_B * R0_E3_B + l4_B * R0_E4_B)*exp(-adj_C + adj_sum));
			R1_mat[i][j] = R1_E1(i,a1);
			R2_mat[i][j] = R2_E1(j,a2);
			R1_mat[i][j] = R1_mat[i][j] * l_A_mat[i][j];
			R2_mat[i][j] = R2_mat[i][j] * l_A_mat[i][j];
			log_R0_mat[i][j] = log_R0_E1(x,y,i,j,a0) * l_A_mat[i][j];
			log_R1_mat[i][j] = log_R1_E1(i,a1) * l_A_mat[i][j];
			log_R2_mat[i][j] = log_R2_E1(j,a2) * l_A_mat[i][j];
			log_R0_E += log_R0_mat[i][j] * l_C_mat[i][j] * exp(adj_sum) * l1_B + log_R0_mat[i][j] * (l2_B +  l3_B +  l4_B)*exp(-adj_C + adj_sum);
			log_R1_E += log_R1_mat[i][j] * l_C_mat[i][j] * exp(adj_sum) * l1_B + log_R1_mat[i][j] * (l2_B +  l3_B +  l4_B)*exp(-adj_C + adj_sum);
			log_R2_E += log_R2_mat[i][j] * l_C_mat[i][j] * exp(adj_sum) * l1_B + log_R2_mat[i][j] * (l2_B +  l3_B +  l4_B)*exp(-adj_C + adj_sum);
			R1_E += R1_mat[i][j] * l_C_mat[i][j] * exp(adj_sum) * l1_B * R1_E1_B + R1_mat[i][j] * (l2_B * R1_E2_B + l3_B * R1_E3_B + l4_B * R1_E4_B)*exp(-adj_C + adj_sum));
			R2_E += R2_mat[i][j] * l_C_mat[i][j] * exp(adj_sum) * l1_B * R2_E1_B + R2_mat[i][j] * (l2_B * R2_E2_B + l3_B * R2_E3_B + l4_B * R2_E4_B)*exp(-adj_C + adj_sum));
		}
	  }
	  R0_E = R0_E / l_sum;
	  R1_E = R1_E / l_sum;
	  R2_E = R2_E / l_sum;
	  log_R0_E = log_R0_E + sum_AC * exp(adj_sum) * l1_B * log (R0_E1_B) + sum_A * (l2_B * log(R0_E2_B) + l3_B * log(R0_E3_B) + l4_B * log(R0_E4_B)) *exp(-adj_C + adj_sum);
	  log_R0_E = log_R0_E / l_sum; 
	  log_R1_E = log_R1_E + sum_AC * exp(adj_sum) * l1_B * log (R1_E1_B) + sum_A * (l2_B * log(R1_E2_B) + l3_B * log(R1_E3_B) + l4_B * log(R1_E4_B)) *exp(-adj_C + adj_sum);
	  log_R1_E = log_R1_E / l_sum; 
	  log_R2_E = log_R2_E + sum_AC * exp(adj_sum) * l1_B * log (R2_E1_B) + sum_A * (l2_B * log(R2_E2_B) + l3_B * log(R2_E3_B) + l4_B * log(R2_E4_B)) *exp(-adj_C + adj_sum);
	  log_R2_E = log_R2_E / l_sum; 
	  E_E1 = sum_AC * exp(adj_sum) * l1_B;
	  E_E2 = sum_A * l2_B *exp(-adj_C + adj_sum);
	  E_E3 = sum_A * l3_B *exp(-adj_C + adj_sum);
	  E_E4 = sum_A * l4_B *exp(-adj_C + adj_sum);
	  v_E <- (sum_AC * exp(adj_sum) * l1_B * y + 
            sum_A * l2_B * a2 * b2*exp(-adj_C + adj_sum) +
            dnbinom(x, a0 + a1 + 1, b1/(1+b1)) * exp(-adj_A - adj_C + adj_sum) * a0 * b2 * p2 * y==0?1:0 +
            sum_A * l3_B * y *exp(-adj_C + adj_sum) +
            sum_A * l4_B * (a0 + a2) * b2 *exp(-adj_C + adj_sum)) / l_sum;
}
int main()
{
    long double a,b;
	b = 0_23;
	a=tgamma(b);
	cout << boost::math::digamma(3.14)<< endl;
	return 0;
 } 
