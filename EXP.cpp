#include <iostream>
#include <math.h>
#include <string>
#include <time.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
using namespace std;
long double R0_mat[1001][1001],R1_mat[1001][1001],R2_mat[1001][1001],log_R0_mat[1001][1001],log_R1_mat[1001][1001],log_R2_mat[1001][1001];
long double l_A_mat[1001][1001];
long double l_C_mat[1001][1001];
long double l_AC_mat[1001][1001];
long double l1(int x, int y, double a0, double a1, double a2, int k, int m, double adjj = 0)
{
	return(exp(lgamma(a1 + k)- lgamma(k+1)- lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m)
                                   + lgamma(m + a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) - adjj));
}
long double l1_c (double t1, double t2, int k, int m, double adjj = 0)
{
	return(exp(k *log(t1) + m *log(t2) - adjj));
}
long double l1_AC (double t1, double t2, int x, int y, int a0, int a1, int a2, int k, int m, double adjj = 0) 
{
	return(exp(lgamma(a1 + k)- lgamma(k+1)- lgamma(a1) + lgamma(x + y + a0 -m -k) - lgamma(x -k +1) - lgamma(a0 + y - m) 
                                   + lgamma(m + a2) - lgamma(m+1) - lgamma(a2) + lgamma(y +a0 -m) - lgamma(y -m +1) - lgamma(a0) + k *log(t1) + m *log(t2) - adjj));
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
	return(exp(lgamma(r + k) - lgamma(r) - lgamma(k+1)) * po(p, r) * po((1-p), k));
}
long double* dBvZINB4_Expt(int x, int y, int a0, int a1, int a2, int b1, int b2, double p1, double p2, double p3, double p4) 
{
	double t1 = (float)(b1 + b2 + 1) /(b1 + 1);
	double t2 = (float)(b1 + b2 + 1) /(b2 + 1);
	double adj_A = 0;
	double adj_B1 = 0;
	double adj_C = 0;
	double adj_sum = 0;
	long double l1_B = - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1 + b1) - a2 * log(1 + b2);
	long double l2_B = y==0?1:0 * exp(- (x + a0 + a1)*log(1 + b1) + x * log(b1) + adj_B1) * p2 ;
	long double l3_B = x==0?1:0 * exp(- (y + a0 + a2)*log(1 + b2) + y * log(b2) + adj_B1) * p3 ;
	long double l4_B = p4 * (x+y==0)?1:0 * exp(adj_B1);

	long double sum_A = 0,sum_AC = 0,sum_A_mat = 0,sum_C_mat = 0;
	for(int i = 0;i <= x;i++)
	{
		for(int j = 0;j <= y;j++)
		{
			l_A_mat[i][j] = l1(x,y,a0,a1,a2,i,j,adj_A);
			l_C_mat[i][j] = l1_c(t1,t2,i,j,adj_C);
			sum_A_mat = sum_A_mat + l_A_mat[i][j];
			sum_C_mat += l_C_mat[i][j];
			//cout << i <<" "<<j<<" "<< l_A_mat[i][j] << endl;
			//system("pause");
		}
	}
	
	if (l1_B < -200 && log(l2_B + l3_B + l4_B) < 0)
	{
		adj_B1 = ((-l1_B - 200)*1.0 / 100) * 100; // prevent exp(l1_B) from being 0
	    l1_B = l1_B + adj_B1;
	}
	l1_B = exp(l1_B) * p1;
	cout << l1_B <<endl;
	
	while (log(sum_A_mat) > 250) 
	{
		sum_A_mat = 0;
		adj_A = adj_A + 200;
		for(int i = 0;i <= x;i++)
		{
			for(int j = 0;j <= y;j++)
			{
				l_A_mat[i][j] = l1(x,y,a0,a1,a2,i,j,adj_A); 
				sum_A_mat += l_A_mat[i][j];	
			}
	    }
	}
	while (log(sum_C_mat) > 250) 
	{
		sum_C_mat = 0;
		adj_C = adj_C + 200;
		for(int i = 0;i <= x;i++)
		{
			for(int j = 0;j <= y;j++)
			{
				l_C_mat[i][j] = l1_c(t1,t2,i,j,adj_C);  
				sum_C_mat += l_C_mat[i][j];	
			}
	    }
	}
	for(int i = 0;i <= x;i++)
	{
		for(int j = 0;j <= y;j++)
		{
			sum_AC = sum_AC + l_A_mat[i][j]*l_C_mat[i][j];
		}
	}
	if(log(sum_AC) > 200)
	{
		adj_A = adj_A + 100;
    	adj_C = adj_C + 100;
    	sum_AC = 0;
		sum_A = 0; 
    	for(int i = 0;i <= x;i++)
		{
			for(int j = 0;j <= y;j++)
			{
				l_A_mat[i][j] = l1(x,y,a0,a1,a2,i,j,adj_A);  
				l_C_mat[i][j] = l1_c(t1,t2,i,j,adj_C);  
				sum_AC += l_A_mat[i][j]*l_C_mat[i][j];	
				sum_A +=  l_A_mat[i][j]; 
			}
	    }
	}
	else if(log(sum_AC) < - 100)
	{
		adj_A = adj_A - 200;
    	adj_C = adj_C - 200;
    	sum_AC=0;
    	sum_A= 0;
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
	long double l_sum =0;
	cout << sum_AC << " "<< l2_B << endl;
	l_sum = sum_AC * l1_B + sum_A * (l2_B +  l3_B +  l4_B) * exp(-adj_C);

	if (l_sum == 0) 
	{
	    adj_sum = -floor(log(sum_AC)*2*1.0/3 + log(l1_B)*2*1.0/3);
	    //cout << adj_sum<<"adj";
	    //cout << sum_AC <<"AC"<< adj_sum <<"@"<< l1_B<<"c"<<adj_C<<endl;
	    l_sum = sum_AC * exp(adj_sum) * l1_B + sum_A * (exp(adj_sum) * (l2_B +  l3_B +  l4_B)) * exp(-adj_C);
  	}
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
	  
            long double xx = sum_AC * exp(adj_sum) * l1_B * y ; 
            cout<<xx<<"&"<<endl;
            long double yy = sum_A * l2_B * a2 * b2*exp(-adj_C + adj_sum) +
            dnbinom(x, a0 + a1 + 1, b1*1.0/(1+b1)) * exp(-adj_A - adj_C + adj_sum) * a0 * b2 * p2 * y==0?1:0 +
            sum_A * l3_B * y *exp(-adj_C + adj_sum) +
            sum_A * l4_B * (a0 + a2) * b2 *exp(-adj_C + adj_sum);
            cout << yy <<"%"<<endl;
            long double v_E = xx + yy;
            v_E= v_E*1.0/l_sum;
      long double *result = new long double[12];
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
	long double *str = new long double[12];
	string name[12] = {"logdensity","R0_E","R1_E","R2_E", "log_R0_E","log_R1_E","log_R2_E","v_E","E_E1","E_E2","E_E3","E_E4"};
	clock_t start = clock();
	str = dBvZINB4_Expt(x,y,a0,a1,a2,b1,b2,p1,p2,p3,p4);
	    clock_t finish = clock();
	for(int i = 0;i <= 11;i++)
	{
		cout << name[i] << " " << str[i] << endl;
	}
	cout << "start :\t\t" << (double)start/CLOCKS_PER_SEC << 's' << endl;
    cout << "finish:\t\t" << (double)finish/CLOCKS_PER_SEC << 's' << endl;

    cout << "finish-start:\t" << (double)(finish - start) / CLOCKS_PER_SEC << 's' << endl;
	return 0;
    /*double a,b;
	b = 0.23;
	a=tgamma(b);
	cout << boost::math::digamma(3.14)<< endl;
	
	int x=150, y=150, a0=1,  a1=2,  a2=1, b1=2, b2=3;
    double p1=0.2, p2=0.2,  p3=0.3,  p4=0.3;
	double t1 = (float)(b1 + b2 + 1) /(b1 + 1);
	double t2 = (float)(b1 + b2 + 1) /(b2 + 1);
	double adj_A = 0;
	double adj_B1 = 0;
	double adj_C = 0;
	double adj_sum = 0;
	long double l1_B = - (+x+y+a0)*log(1 + b1 + b2) + x * log(b1) + y * log(b2) - a1 * log(1 + b1) - a2 * log(1 + b2);
	long double l2_B = y==0?1:0 * exp(- (x + a0 + a1)*log(1 + b1) + x * log(b1) + adj_B1) * p2 ;
	long double l3_B = x==0?1:0 * exp(- (y + a0 + a2)*log(1 + b2) + y * log(b2) + adj_B1) * p3 ;
	long double l4_B = p4 * (x+y==0)?1:0 * exp(adj_B1);

	long double sum_A = 0,sum_AC = 0,sum_A_mat = 0,sum_C_mat = 0;
	for(int i = 0;i <= x;i++)
	{
		for(int j = 0;j <= y;j++)
		{
			l_A_mat[i][j] = l1(x,y,a0,a1,a2,i,j,adj_A);
			l_C_mat[i][j] = l1_c(t1,t2,i,j,adj_C);
			sum_A_mat = sum_A_mat + l_A_mat[i][j];
			sum_C_mat += l_C_mat[i][j];
			//cout << i <<" "<<j<<" "<< l_A_mat[i][j] << endl;
			//system("pause");
		}
	}
	
	if (l1_B < -200 && log(l2_B + l3_B + l4_B) < 0)
	{
		adj_B1 = ((-l1_B - 200)*1.0 / 100) * 100; // prevent exp(l1_B) from being 0
	    l1_B = l1_B + adj_B1;
	}
	l1_B = exp(l1_B) * p1;
	cout << l1_B <<endl;
	
	while (log(sum_A_mat) > 250) 
	{
		sum_A_mat = 0;
		adj_A = adj_A + 200;
		for(int i = 0;i <= x;i++)
		{
			for(int j = 0;j <= y;j++)
			{
				l_A_mat[i][j] = l1(x,y,a0,a1,a2,i,j,adj_A); 
				sum_A_mat += l_A_mat[i][j];	
			}
	    }
	}
	while (log(sum_C_mat) > 250) 
	{
		sum_C_mat = 0;
		adj_C = adj_C + 200;
		for(int i = 0;i <= x;i++)
		{
			for(int j = 0;j <= y;j++)
			{
				l_C_mat[i][j] = l1_c(t1,t2,i,j,adj_C);  
				sum_C_mat += l_C_mat[i][j];	
			}
	    }
	}
	for(int i = 0;i <= x;i++)
	{
		for(int j = 0;j <= y;j++)
		{
			sum_AC = sum_AC + l_A_mat[i][j]*l_C_mat[i][j];
		}
	}
	if(log(sum_AC) > 200)
	{
		adj_A = adj_A + 100;
    	adj_C = adj_C + 100;
    	sum_AC = 0;
		sum_A = 0; 
    	for(int i = 0;i <= x;i++)
		{
			for(int j = 0;j <= y;j++)
			{
				l_A_mat[i][j] = l1(x,y,a0,a1,a2,i,j,adj_A);  
				l_C_mat[i][j] = l1_c(t1,t2,i,j,adj_C);  
				sum_AC += l_A_mat[i][j]*l_C_mat[i][j];	
				sum_A +=  l_A_mat[i][j]; 
			}
	    }
	}
	else if(log(sum_AC) < - 100)
	{
		adj_A = adj_A - 200;
    	adj_C = adj_C - 200;
    	sum_AC=0;
    	sum_A= 0;
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
	long double l_sum =0;
	cout << sum_AC << " "<< l2_B << endl;
	l_sum = sum_AC * l1_B + sum_A * (l2_B +  l3_B +  l4_B) * exp(-adj_C);

	if (l_sum == 0) 
	{
	    adj_sum = -floor(log(sum_AC)*2*1.0/3 + log(l1_B)*2*1.0/3);
	    //cout << adj_sum<<"adj";
	    //cout << sum_AC <<"AC"<< adj_sum <<"@"<< l1_B<<"c"<<adj_C<<endl;
	    l_sum = sum_AC * exp(adj_sum) * l1_B + sum_A * (exp(adj_sum) * (l2_B +  l3_B +  l4_B)) * exp(-adj_C);
  	}
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
	  
            long double xx = sum_AC * exp(adj_sum) * l1_B * y ; 
            cout<<xx<<"&"<<endl;
            long double yy = sum_A * l2_B * a2 * b2*exp(-adj_C + adj_sum) +
            dnbinom(x, a0 + a1 + 1, b1*1.0/(1+b1)) * exp(-adj_A - adj_C + adj_sum) * a0 * b2 * p2 * y==0?1:0 +
            sum_A * l3_B * y *exp(-adj_C + adj_sum) +
            sum_A * l4_B * (a0 + a2) * b2 *exp(-adj_C + adj_sum);
            cout << yy <<"%"<<endl;
            long double v_E = xx + yy;
            v_E= v_E*1.0/l_sum;
      long double result[12];
      result[0] = log(l_sum) + adj_A -adj_B1 + adj_C - adj_sum;
      result[1] = R0_E;
	  result[2] = R1_E;
	  result[3] = R2_E;
	  result[4] = log_R0_E;
	  result[5] = log_R1_E;
	  result[6] = log_R2_E;
	  //result[7] = E_E;
	  result[7] = v_E;
	  result[8] = E_E1;
	  result[9] = E_E2;
	  result[10] = E_E3;
	  result[11] = E_E4;
	  for(int i =0;i<=11;i++)
	  {
	  	cout << result[i] << endl;
	  }
	return 0;*/
 } 
