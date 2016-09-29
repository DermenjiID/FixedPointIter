
#include "stdafx.h"
#include <vector>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <locale.h> 
#include <stdio.h>
#include <fstream>
#include <Windows.h>

using namespace std;

int main()
{
	
	int k=1; 
	
	char *str = new char [256];
    int n=0;
    ifstream base("smile.txt");
    while (!base.eof())
    {
        base.getline(str, 256, '\n');
        n++;
    }
    base.close();
    delete str;
	
	fstream fi;
	fi.open("smile.txt", ios::in);
	long double ** matr=new long double *[n];
	for (int i=0; i<n; i++) 
	{
        matr[i]=new long double [n];
	}
	int i, j;
	if (!fi)
	{ 
		cerr << "error file doesn't exist\n"; 
        exit(1);
    }
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			fi>>matr[i][j];
		}
	}
	for (i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
		{
			cout<<setw(10)<<matr[i][j]<<"    ";
		}
		cout<<"\n";
	}
	fi.close();
	long double eps=0.001;
	long double del=0;
	for (i=0;i<n;i++)
		{
			for (j=0;j<n;j++)
			{
				if(i!=j)
				{
				   del=del+matr[i][j]*matr[i][j];
				}
			}
		}
	long double ** multVrash=new long double *[n];
	long double ** prevmatr=new long double *[n];
	long double ** mult=new long double *[n];
	for (int i=0; i<n; i++) 
	{
		multVrash[i]=new long double [n];
		prevmatr[i]=new long double [n];
		mult[i]=new long double [n];
	}
	
	for (i=0;i<n;i++)
	{
			for (j=0;j<n;j++)
			{
				prevmatr[i][j]=matr[i][j];
				multVrash[i][j]=0;
				multVrash[i][i]=1;
			}
	 }

	while (true)
	{
		long double ** U=new long double *[n];
		long double ** currmatr=new long double *[n];
		long double ** tU=new long double *[n];
		for (int i=0; i<n; i++) 
	    {
			U[i]=new long double [n];
			currmatr[i]=new long double [n];
			tU[i]=new long double [n];
		}
	    for (i=0;i<n;i++)
		{
			for (j=0;j<n;j++)
			{
				currmatr[i][j]=0;
				U[i][j]=0;
				U[i][i]=1;
			}
		}
	
		long double phi=0;
		int a=0;
		int b=1;
		long double maxd=prevmatr[0][1];
		
		for(i=0;i<n;i++)
		{
			for(j=1;j<n;j++)
				if(i<j)
				{
					if (abs(maxd)<abs(prevmatr[i][j]))
					{
						maxd=prevmatr[i][j];
						a=i;
						b=j;
					}
				}
		}
		
		del=del-2*maxd*maxd;
		
		if(maxd==0)
		{
			break;
		}
		if(del < eps)
		{
			break;
		}
		
		phi=(atan((2*prevmatr[a][b])/(prevmatr[a][a]-prevmatr[b][b])))/2;
		U[a][a]=cos(phi);
		U[b][b]=U[a][a];
		U[a][b]=-sin(phi);
		U[b][a]=sin(phi);

	    for (i=0;i<n;i++)
		{
			for (j=0;j<n;j++)
			{
				tU[i][j]=U[j][i];
			}
		}

		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				mult[i][j] = 0;
				for(int z = 0; z < n; z++)
				{
					mult[i][j] += tU[i][z] *prevmatr[z][j];
				}
			}
		}
		
		for( i = 0; i < n; i++)
		{
			for( j = 0; j < n; j++)
			{
				for(int z = 0; z < n; z++)
				{
					currmatr[i][j] += mult[i][z] *U[z][j];
				}
			}
		}
		
		for (i=0;i<n;i++)
		{
			for (j=0;j<n;j++)
			{
				prevmatr[i][j]=currmatr[i][j];
			}
		}

		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				mult[i][j]=0;
				for(int z = 0; z < n; z++)
				{
					mult[i][j] += multVrash[i][z] *U[z][j];
				}
			}
		}
		for (i=0;i<n;i++)
		{
			for (j=0;j<n;j++)
			{
				multVrash[i][j]=mult[i][j];
			}
		}
	
		cout<<"Resultat"<<' '<<k<<' '<<"iteracii"<<endl;
		cout<<"max naddiagonalnii element"<<' '<<maxd<<endl;
		cout<<"FI:"<<' '<<phi<<endl;
		cout<<"matr vrasheniya:"<<endl;
		for (i=0;i<n;i++)
		{
			for (j=0;j<n;j++)
			{
				cout<<setw(10)<<U[i][j]<<"  ";
			}
			cout<<"\n";
		}
		cout<<"Matrica na"<<' '<<k<<' '<< "iteracii :"<<endl;
		for (i=0;i<n;i++)
		{
			for (j=0;j<n;j++)
			{
				cout<<setw(10)<<prevmatr[i][j]<<"  ";
			}
			cout<<"\n";
		}
		cout<<endl;
		k++;
	}


	vector <long double> sobstvVec (n,0);
	
	for(i=0;i<n;i++)
	{
		cout<<"Sobstvennoe znachenie"<<' '<<i+1<<":"<<' '<<prevmatr[i][i]<<endl;
		cout<<endl;
	}
	
	for ( i = 0; i < n; i++) 
	{
		cout<<"Sobstvennii vektor"<<' '<<i+1<<endl;
		cout<<"{"<<' ';
        for ( j = 0; j < n; j++) 
		{
            sobstvVec[i] = multVrash[j][i]/multVrash[i][i];
			cout << ' ' << sobstvVec[i];
        }
		cout<<' '<<"}";
        cout<<endl;
		cout<<endl;
    }
	system("pause");
	return 0;
}

