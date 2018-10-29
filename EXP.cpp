#include <iostream>
#include <math.h>
#include <string>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
using namespace std;
double R0_mat[1001][1001],R1_mat[1001][1001],R2_mat[1001][1001],log_R0_mat[1001][1001],log_R1_mat[1001][1001],log_R2_mat[1001][1001];
long double l1(int x, int y, int a0, int a1, int a2, int k, int m, double adjj = 0)
{
	return(exp(log(tgamma(a1 + k))- log(tgamma(k+1))- log(tgamma(a1)) + log(tgamma(x + y + a0 -m -k)) - log(tgamma(x -k +1)) - log(tgamma(a0 + y - m)) 
                                   + log(tgamma(m + a2)) - log(tgamma(m+1)) - log(tgamma(a2)) + log(tgamma(y +a0 -m)) - log(tgamma(y -m +1)) - log(tgamma(a0)) - adjj));
}
long double l1_c (double t1, double t2, int k, int m, double adjj = 0)
{
	return(exp(k *log(t1) + m *log(t2) - adjj));
}
long double l1_AC (double t1, double t2, int x, int y, int a0, int a1, int a2, int k, int m, double adjj = 0) 
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
	return(boost::math::digamma(x - k + y - m + a0));
} 
long double R1_E1(int k, double a1) 
{
	return(k + a1);
}
long double log_R1_E1(int k, double a1) 
{
	return(boost::math::digamma(k + a1));
}
long double R2_E1(int m, double a2) 
{
	return(m + a2);
}
long double log_R2_E1(int m, double a2) 
{
	return(boost::math::digamma(m + a2));
}
double po(double X,int Y)
{
	double Q = 1;
	for(int i = 1;i <= Y;i++)
	{
		Q = Q * X;
	}
	return(Q);
}
double dnbinom(int k,int r, double p)
{
	return(exp(log(tgamma(r + k)) - log(tgamma(r)) - log(tgamma(k+1))) * po(p, r) * po((1-p), k));
}
double* dBvZINB4_Expt(int x, int y, int a0, int a1, int a2, int b1, int b2, double p1, double p2, double p3, double p4) 
{
		double t1 = (float)(b1 + b2 + 1) /(b1 + 1);
	double t2 = (float)(b1 + b2 + 1) /(b2 + 1);
	double adj_A = 0;
	double adj_B1 = 0;
	double adj_C = 0;
	double adj_sum = 0;
	double l1_B = - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1 + b1) - a2 * log(1 + b2);
	double l2_B = exp(- (x + a0 + a1)*log(1 + b1) + x * log(b1) + adj_B1) * p2 * y==0?1:0;
	double l3_B = exp(- (y + a0 + a2)*log(1 + b2) + y * log(b2) + adj_B1) * p3 * x==0?1:0;
	double l4_B = p4 * (x+y==0)?1:0 * exp(adj_B1);
	double l_A_mat[x+1][y+1];
	double l_C_mat[x+1][y+1];
	double l_AC_mat[x+1][y+1];
	//cout << "@@"<< l2_B << endl;
	double sum_A = 0,sum_AC = 0,sum_A_mat = 0,sum_C_mat = 0;
	for(int i = 0;i <= x;i++)
	{
		for(int j = 0;j <= y;j++)
		{
			l_A_mat[i][j] = l1(x,y,a0,a1,a2,i,j,adj_A);
			l_C_mat[i][j] = l1_c(t1,t2,i,j,adj_C);
			sum_A_mat = sum_A_mat + l_A_mat[i][j];
			sum_C_mat += l_C_mat[i][j];
			sum_AC = sum_AC + l_A_mat[i][j]*l_C_mat[i][j];
		}
	}
//	cout << sum_AC << endl;
	
	if (l1_B < -200 && log(l2_B + l3_B + l4_B) < 0)
	{
		adj_B1 = ((-l1_B - 200)*1.0 / 100) * 100; // prevent exp(l1_B) from being 0
	    l1_B = l1_B + adj_B1;
	}
	l1_B = exp(l1_B) * p1;
	//cout << l1_B <<endl;
	
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
	double l_sum =0;
	l_sum = sum_AC * l1_B + sum_A * (l2_B +  l3_B +  l4_B) * exp(-adj_C);

	if (l_sum == 0) 
	{
	    adj_sum = -floor(log(sum_AC)*2*1.0/3 + log(l1_B)*2*1.0/3);
	    //cout << adj_sum<<"adj";
	    //cout << sum_AC <<"AC"<< adj_sum <<"@"<< l1_B<<"c"<<adj_C<<endl;
	    l_sum = sum_AC * exp(adj_sum) * l1_B + sum_A * (exp(adj_sum) * (l2_B +  l3_B +  l4_B)) * exp(-adj_C);
  	}
    // cout << l_sum << "$$"<<endl;
  	long double R0_E1_B,R0_E2_B,R0_E3_B,R0_E4_B,R1_E1_B,R1_E2_B,R1_E3_B,R1_E4_B,R2_E1_B,R2_E2_B,R2_E3_B,R2_E4_B;
	// expectation components
	  R0_E1_B = (float)b1/(1 + b1 + b2);
	  R0_E2_B = (float)b1/(1 + b1);
	  R0_E3_B = (float)b1/(1 + b2);
	  R0_E4_B = (float)b1;
	  
	  R1_E1_B = (float)b1/(1 + b1);
	  R1_E2_B = (float)b1/(1 + b1);
	  R1_E3_B = (float)b1;
	  R1_E4_B = (float)b1;
	  
	  R2_E1_B = (float)b1/(1 + b2); 
	  R2_E2_B = (float)b1;
	  R2_E3_B = (float)b1/(1 + b2);
	  R2_E4_B = (float)b1;
	  
	  double log_R0_E = 0,log_R1_E = 0,log_R2_E = 0,R0_E = 0,R1_E = 0,R2_E = 0;
	  for(int i = 0;i <= x;i++)
	  {
	  	for(int j = 0;j <= y;j++)
	  	{
	  		R0_mat[i][j] = R0_E1(x,y,i,j,a0);
	  		R0_mat[i][j] = R0_mat[i][j]*l_A_mat[i][j];
			R0_E += R0_mat[i][j] * l_C_mat[i][j] * exp(adj_sum) * l1_B * R0_E1_B + R0_mat[i][j] * (l2_B * R0_E2_B + l3_B * R0_E3_B + l4_B * R0_E4_B)*exp(-adj_C + adj_sum);
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
			R1_E += R1_mat[i][j] * l_C_mat[i][j] * exp(adj_sum) * l1_B * R1_E1_B + R1_mat[i][j] * (l2_B * R1_E2_B + l3_B * R1_E3_B + l4_B * R1_E4_B)*exp(-adj_C + adj_sum);
			R2_E += R2_mat[i][j] * l_C_mat[i][j] * exp(adj_sum) * l1_B * R2_E1_B + R2_mat[i][j] * (l2_B * R2_E2_B + l3_B * R2_E3_B + l4_B * R2_E4_B)*exp(-adj_C + adj_sum);
		}
	  }
	  
	  
	  R0_E = R0_E*1.0 / l_sum;
	  R1_E = R1_E*1.0 / l_sum;
	  R2_E = R2_E*1.0 / l_sum;
	  log_R0_E = log_R0_E + sum_AC * exp(adj_sum) * l1_B * log (R0_E1_B) + sum_A * (l2_B * log(R0_E2_B) + l3_B * log(R0_E3_B) + l4_B * log(R0_E4_B)) *exp(-adj_C + adj_sum);
	  log_R0_E = log_R0_E * 1.0 / l_sum; 
	  log_R1_E = log_R1_E + sum_AC * exp(adj_sum) * l1_B * log (R1_E1_B) + sum_A * (l2_B * log(R1_E2_B) + l3_B * log(R1_E3_B) + l4_B * log(R1_E4_B)) *exp(-adj_C + adj_sum);
	  log_R1_E = log_R1_E * 1.0 / l_sum; 
	  log_R2_E = log_R2_E + sum_AC * exp(adj_sum) * l1_B * log (R2_E1_B) + sum_A * (l2_B * log(R2_E2_B) + l3_B * log(R2_E3_B) + l4_B * log(R2_E4_B)) *exp(-adj_C + adj_sum);
	  log_R2_E = log_R2_E * 1.0 / l_sum; 
	
	  double E_E1 = sum_AC * exp(adj_sum) * l1_B;
	  double E_E2 = sum_A * l2_B *exp(-adj_C + adj_sum);
	  double E_E3 = sum_A * l3_B *exp(-adj_C + adj_sum);
	  double E_E4 = sum_A * l4_B *exp(-adj_C + adj_sum);
	  
	  double su = E_E1+E_E2+E_E3+E_E4;
	  E_E1 = E_E1/su;
	  E_E2 = E_E2/su;
	  E_E3 = E_E3/su;
	  E_E4 = E_E4/su;
      //E_E{2} <- E.E/sum(E.E)
	  
            double xx = sum_AC * exp(adj_sum) * l1_B * y ; 
            double yy = sum_A * l2_B * a2 * b2*exp(-adj_C + adj_sum) +
            dnbinom(x, a0 + a1 + 1, b1*1.0/(1+b1)) * exp(-adj_A - adj_C + adj_sum) * a0 * b2 * p2 * y==0?1:0 +
            sum_A * l3_B * y *exp(-adj_C + adj_sum) +
            sum_A * l4_B * (a0 + a2) * b2 *exp(-adj_C + adj_sum);
            double v_E = xx + yy;
            v_E= v_E/l_sum;
      double *result = new double[12];
      result[0] = log(l_sum) + adj_A -adj_B1 + adj_C - adj_sum;
      result[1] = R0_E;
	  result[2] = R1_E;
	  result[3] = R2_E;
	  result[4] = log_R0_E;
	  result[5] = log_R1_E;
	  result[6] = log_R2_E;
	  result[7] = v_E;
	  result[8] = E_E1;
	  result[9] = E_E2;
	  result[10] = E_E3;
	  result[11] = E_E4;
	  return(result);
}
//void ML_BvZINB4(int xvec[], int yvec[],double tol = 1e-8,int maxiter = 200)

int main()
{
	int x,y,a0,a1,a2,b1,b2;
	double p1,p2,p3,p4;
	cin >> x >> y >> a0 >> a1 >> a2>> b1 >> b2 >> p1 >> p2>> p3>>p4;
	double *str = new double[12];
	string name[12] = {"logdensity","R0_E","R1_E","R2_E", "log_R0_E","log_R1_E","log_R2_E","v_E","E_E1","E_E2","E_E3","E_E4"};
	str = dBvZINB4_Expt(x,y,a0,a1,a2,b1,b2,p1,p2,p3,p4);
	for(int i = 0;i <= 11;i++)
	{
		cout << name[i] << " " << str[i] << endl;
	}
	return 0;
    /*double a,b;
	b = 0.23;
	a=tgamma(b);
	cout << boost::math::digamma(3.14)<< endl;
	
	int x=20, y=5, a0=1,  a1=2,  a2=1, b1=2, b2=3;
	double p1=0.2, p2=0.2,  p3=0.3,  p4=0.3;
	double t1 = (float)(b1 + b2 + 1) /(b1 + 1);
	double t2 = (float)(b1 + b2 + 1) /(b2 + 1);
	double adj_A = 0;
	double adj_B1 = 0;
	double adj_C = 0;
	double adj_sum = 0;
	double l1_B = - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1 + b1) - a2 * log(1 + b2);
	double l2_B = exp(- (x + a0 + a1)*log(1 + b1) + x * log(b1) + adj_B1) * p2 * y==0?1:0;
	double l3_B = exp(- (y + a0 + a2)*log(1 + b2) + y * log(b2) + adj_B1) * p3 * x==0?1:0;
	double l4_B = p4 * (x+y==0)?1:0 * exp(adj_B1);
	double l_A_mat[x+1][y+1];
	double l_C_mat[x+1][y+1];
	double l_AC_mat[x+1][y+1];
	cout << "@@"<< l2_B << endl;
	double sum_A = 0,sum_AC = 0,sum_A_mat = 0,sum_C_mat = 0;
	for(int i = 0;i <= x;i++)
	{
		for(int j = 0;j <= y;j++)
		{
			l_A_mat[i][j] = l1(x,y,a0,a1,a2,i,j,adj_A);
			l_C_mat[i][j] = l1_c(t1,t2,i,j,adj_C);
			sum_A_mat = sum_A_mat + l_A_mat[i][j];
			sum_C_mat += l_C_mat[i][j];
			sum_AC = sum_AC + l_A_mat[i][j]*l_C_mat[i][j];
		}
	}
	cout << sum_AC << endl;
	
	if (l1_B < -200 && log(l2_B + l3_B + l4_B) < 0)
	{
		adj_B1 = ((-l1_B - 200)*1.0 / 100) * 100; // prevent exp(l1_B) from being 0
	    l1_B = l1_B + adj_B1;
	}
	l1_B = exp(l1_B) * p1;
	cout << l1_B <<endl;
	
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
	double l_sum =0;
	l_sum = sum_AC * l1_B + sum_A * (l2_B +  l3_B +  l4_B) * exp(-adj_C);

	if (l_sum == 0) 
	{
	    adj_sum = -floor(log(sum_AC)*2*1.0/3 + log(l1_B)*2*1.0/3);
	    //cout << adj_sum<<"adj";
	    //cout << sum_AC <<"AC"<< adj_sum <<"@"<< l1_B<<"c"<<adj_C<<endl;
	    l_sum = sum_AC * exp(adj_sum) * l1_B + sum_A * (exp(adj_sum) * (l2_B +  l3_B +  l4_B)) * exp(-adj_C);
  	}
    // cout << l_sum << "$$"<<endl;
  	long double R0_E1_B,R0_E2_B,R0_E3_B,R0_E4_B,R1_E1_B,R1_E2_B,R1_E3_B,R1_E4_B,R2_E1_B,R2_E2_B,R2_E3_B,R2_E4_B;
	// expectation components
	  R0_E1_B = (float)b1/(1 + b1 + b2);
	  R0_E2_B = (float)b1/(1 + b1);
	  R0_E3_B = (float)b1/(1 + b2);
	  R0_E4_B = (float)b1;
	  
	  R1_E1_B = (float)b1/(1 + b1);
	  R1_E2_B = (float)b1/(1 + b1);
	  R1_E3_B = (float)b1;
	  R1_E4_B = (float)b1;
	  
	  R2_E1_B = (float)b1/(1 + b2); 
	  R2_E2_B = (float)b1;
	  R2_E3_B = (float)b1/(1 + b2);
	  R2_E4_B = (float)b1;
	  
	  double log_R0_E = 0,log_R1_E = 0,log_R2_E = 0,R0_E = 0,R1_E = 0,R2_E = 0;
	  for(int i = 0;i <= x;i++)
	  {
	  	for(int j = 0;j <= y;j++)
	  	{
	  		R0_mat[i][j] = R0_E1(x,y,i,j,a0);
	  		R0_mat[i][j] = R0_mat[i][j]*l_A_mat[i][j];
			R0_E += R0_mat[i][j] * l_C_mat[i][j] * exp(adj_sum) * l1_B * R0_E1_B + R0_mat[i][j] * (l2_B * R0_E2_B + l3_B * R0_E3_B + l4_B * R0_E4_B)*exp(-adj_C + adj_sum);
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
			R1_E += R1_mat[i][j] * l_C_mat[i][j] * exp(adj_sum) * l1_B * R1_E1_B + R1_mat[i][j] * (l2_B * R1_E2_B + l3_B * R1_E3_B + l4_B * R1_E4_B)*exp(-adj_C + adj_sum);
			R2_E += R2_mat[i][j] * l_C_mat[i][j] * exp(adj_sum) * l1_B * R2_E1_B + R2_mat[i][j] * (l2_B * R2_E2_B + l3_B * R2_E3_B + l4_B * R2_E4_B)*exp(-adj_C + adj_sum);
		}
	  }
	  
	  
	  R0_E = R0_E*1.0 / l_sum;
	  R1_E = R1_E*1.0 / l_sum;
	  R2_E = R2_E*1.0 / l_sum;
	  log_R0_E = log_R0_E + sum_AC * exp(adj_sum) * l1_B * log (R0_E1_B) + sum_A * (l2_B * log(R0_E2_B) + l3_B * log(R0_E3_B) + l4_B * log(R0_E4_B)) *exp(-adj_C + adj_sum);
	  log_R0_E = log_R0_E * 1.0 / l_sum; 
	  log_R1_E = log_R1_E + sum_AC * exp(adj_sum) * l1_B * log (R1_E1_B) + sum_A * (l2_B * log(R1_E2_B) + l3_B * log(R1_E3_B) + l4_B * log(R1_E4_B)) *exp(-adj_C + adj_sum);
	  log_R1_E = log_R1_E * 1.0 / l_sum; 
	  log_R2_E = log_R2_E + sum_AC * exp(adj_sum) * l1_B * log (R2_E1_B) + sum_A * (l2_B * log(R2_E2_B) + l3_B * log(R2_E3_B) + l4_B * log(R2_E4_B)) *exp(-adj_C + adj_sum);
	  log_R2_E = log_R2_E * 1.0 / l_sum; 
	
	  double E_E1 = sum_AC * exp(adj_sum) * l1_B;
	  double E_E2 = sum_A * l2_B *exp(-adj_C + adj_sum);
	  double E_E3 = sum_A * l3_B *exp(-adj_C + adj_sum);
	  double E_E4 = sum_A * l4_B *exp(-adj_C + adj_sum);
	  E_E[1] = sum_AC * exp(adj_sum) * l1_B;
	  E_E[2] = sum.A * l2_B *exp(-adj_C + adj_sum);
	  E_E[3] = sum.A * l3_B *exp(-adj_C + adj_sum);
      //E_E{2} <- E.E/sum(E.E)
	  
            double xx = sum_AC * exp(adj_sum) * l1_B * y ; 
            double yy = sum_A * l2_B * a2 * b2*exp(-adj_C + adj_sum) +
            dnbinom(x, a0 + a1 + 1, b1*1.0/(1+b1)) * exp(-adj_A - adj_C + adj_sum) * a0 * b2 * p2 * y==0?1:0 +
            sum_A * l3_B * y *exp(-adj_C + adj_sum) +
            sum_A * l4_B * (a0 + a2) * b2 *exp(-adj_C + adj_sum);
            double v_E = xx + yy;
            v_E= v_E/l_sum;
      double result[8];
      result[0] = log(l_sum) + adj_A -adj_B1 + adj_C - adj_sum;
      result[1] = R0_E;
	  result[2] = R1_E;
	  result[3] = R2_E;
	  result[4] = log_R0_E;
	  result[5] = log_R1_E;
	  result[6] = log_R2_E;
	  //result[7] = E_E;
	  result[7] = v_E;
	  for(int i =0;i<=7;i++)
	  {
	  	cout << result[i] << endl;
	  }
	return 0;*/
 } 
