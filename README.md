//11.12.15
//MS VS2015

#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iomanip>

#define anisotropy 0.5 
struct element
{
	double value;
	int column;
	int row;
};

struct data
{
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double t_min;
	double t_max;
	double **prev_func_values; //1ts(close) index - X, 2nd - Y, size!!!: (Nx+1)*(Ny+1)
	int Nx;
	int Ny;
	int Nt;
	double curr_t;
	element *matrix;
	int is_triag;
	std::ofstream *debug;
};
int print_matrix(std::ofstream &str, data D);

double func(double x, double y, double t)
{
	return	(1 - x*x)*(1 - y*y) +(y + 1)*(1 - y*y)*sin(2 * M_PI*t);
}

double func_x_min(double y, double t)
{
	return 0;
}

double func_x_max(double y, double t)
{
	return 0;
}

double func_y_min(double x, double t)
{
	return 0;
}

double func_y_max(double x, double t)
{
	return 0;
}

int build(data* D)
{
//	std::cout << "build started" << std::endl;
	int SIZE = (D->Nx - 1)*(D->Ny - 1);
	if (D->matrix)
	{
		delete D->matrix;
	}
	D->matrix = new element[2*D->Nx*SIZE];
	D->is_triag = 0;
	double tau = (D->t_max - D->t_min) /(D->Nt-1);
	double h = (D->x_max - D->x_min) / (D->Nx);
	double g = (D->y_max - D->y_min) / (D->Ny);
	D->curr_t += tau;
	if (D->debug)
	{
		*(D->debug) <<"func values"<< std::endl;
		for (int j = 0; j < SIZE; j++)
		{
			int m = j / (D->Nx - 1);
			int n = j - m*(D->Nx - 1);
			double x = D->x_min + h*(n + 1);
			double y = D->y_min + g*(m + 1);
			*(D->debug) <<std::setprecision(3)<< -func(x, y, D->curr_t) << "\t";
		}
		*(D->debug) << std::endl;
	}
	for (int j = 0; j < SIZE; j++)
	{
		int m = j / (D->Nx - 1);
		int n = j - m*(D->Nx - 1);
		double x = D->x_min + h*(n+1);
		double y = D->y_min + g*(m+1);
		for (int k = 0; k <= 2 * (D->Nx-1); k++)
		{
			D->matrix[(2 * (j)*D->Nx+k)].value = 0;
			D->matrix[(2 * (j)*D->Nx + k)].row = j;
			D->matrix[(2 * (j)*D->Nx + k)].column = (((j + k - (D->Nx - 1)) >= SIZE)||((j + k - (D->Nx - 1)) <0)) ? (-1) : (j + k - (D->Nx - 1));
		}
		D->matrix[(2 * (j + 1)*D->Nx) - 1].value = -*(*(D->prev_func_values + n+1) + m+1) / tau - func(x, y, D->curr_t);//F
		D->matrix[(2 * (j + 1)*D->Nx) - 1].row = j;
		D->matrix[(2 * (j + 1)*D->Nx) - 1].column = SIZE;
		
		if (m== 0)
			D->matrix[(2 * (j + 1)*D->Nx) - 1].value -= anisotropy*func_y_min(x, D->curr_t) / (g*g);
		else
			D->matrix[(2 * (j)*D->Nx)].value = anisotropy / (g*g);

		if (n== 0)
			D->matrix[(2 * (j + 1)*D->Nx) - 1].value -= func_x_min(y, D->curr_t) / (h*h);
		else
			D->matrix[(2 * (j)*D->Nx) + D->Nx - 2].value = 1 / (h*h);

		D->matrix[(2 * j*D->Nx) + D->Nx - 1].value = -2 * anisotropy / (g*g) - 2 / (h*h) - 1 / tau;

		if (n== (D->Nx-2))
			D->matrix[(2 * (j + 1)*D->Nx) - 1].value -= func_x_max(y, D->curr_t) / (h*h);
		else
			D->matrix[(2 * j*D->Nx) + D->Nx].value = 1 / (h*h);

		if (m== (D->Ny-2))
			D->matrix[(2 * (j + 1)*D->Nx) - 1].value -= anisotropy*func_y_max(x, D->curr_t) / (g*g);
		else
			D->matrix[(2 * j*D->Nx) + 2 * (D->Nx - 1)].value = anisotropy / (g*g);
	}
	return 0;
}

int print_func(std::ofstream &str, data D)
{
	if (!D.is_triag) return -2;
	double x, y;
	for (int n = 0; n < D.Nx + 1; n++)
	{
		for (int m = 0; m < D.Ny + 1; m++)
		{
			x = D.x_min + n*(D.x_max - D.x_min) / (D.Nx);
			y = D.y_min + m*(D.y_max - D.y_min) / (D.Ny);
			str << std::setprecision(5)<<D.curr_t <<"\t"<< x << "\t" << y << "\t" << *(*(D.prev_func_values + n) + m) << "\t" << std::endl;
		}
	}
	return 0;
}


int print_center(std::ofstream &str, data D)
{
	if (!D.is_triag) return -2;
	int m, n;

	n = D.Nx / 2;
	m = D.Ny / 2;

	str << std::setprecision(5) << D.curr_t << "\t" << *(*(D.prev_func_values + n) + m) << "\t" << std::endl;

	return 0;
}


element* search_row_col(element* els, int row, int col,int SIZE,int Nx)
{
	for (int k = 0; k <2 * (Nx); k++)
	{
		if ((els[2*row*Nx+k].column == col))
		{
			return els + 2 * row*Nx + k;
		}
	}
	return 0;
}

int triag_alg(data *D)
{
	//std::cout << "triag_alg started" << std::endl;
	if (D->debug)	*(D->debug) << "start triag" << std::endl;
	element* els = D->matrix;
	if (els == 0) return -1;//error
	int SIZE = (D->Nx - 1)*(D->Ny - 1);
	for (int j = 0; j < SIZE; j++)
	{
		//std::cout << "row: "<<j << std::endl;
		if (D->debug)
		{
			*(D->debug) << "row: " << j << std::endl;
			print_matrix(*(D->debug), *D);
		}
		for (int Raw = j + 1; Raw <((j+ D->Nx<SIZE)?(j+D->Nx):(SIZE)); Raw++)
		{
			element* v,*diag;
			v = search_row_col(els, Raw, j,SIZE,D->Nx);
			if (v)
			{
				if (v->value != 0)
				{
					diag = search_row_col(els, j, j, SIZE, D->Nx);
					if (diag == 0) return -1;
					double divider = v->value / diag->value;
					element* low_row, *high_row;
					for (int Col = j; Col < ((SIZE<j+D->Nx)?(SIZE):(j+D->Nx)); Col++)
					{
						low_row = search_row_col(els, Raw, Col, SIZE, D->Nx);
						high_row = search_row_col(els, j, Col, SIZE, D->Nx);
						//if (!((low_row) && (high_row)) && ((low_row) || (high_row))) return -2;//no need
						if ((low_row) && (high_row))
						{
							low_row->value -= divider*high_row->value;
						}
					}
					low_row = search_row_col(els, Raw, SIZE, SIZE, D->Nx);
					high_row = search_row_col(els, j, SIZE, SIZE, D->Nx);
					if ((low_row) && (high_row))
					{
						low_row->value -= divider*high_row->value;
					}
				}
			}
		}
	}
	D->is_triag = 1;
	return 0;
}

int fill_func(data* D)
{
	//std::cout << "fill_func started" << std::endl;
	int SIZE = (D->Nx - 1)*(D->Ny - 1);
	double *temp_f = new double[SIZE];
	for (int j = SIZE - 1; j > -1; j--)
	{
		element* v = search_row_col(D->matrix,j,SIZE,SIZE,D->Nx);
		if (v == 0)
		{
			delete temp_f;
			return -1;
		}
		temp_f[j] = v->value;
		for (int col = ((SIZE - 1)<(j + D->Nx - 1)) ? (SIZE - 1) : (j + D->Nx - 1); col >j; col--)
		{
			v = search_row_col(D->matrix, j, col, SIZE, D->Nx);
			if (v == 0)
			{
				delete temp_f;
				return -2;
			}
			temp_f[j] -= temp_f[col]*v->value;
		}
		v = search_row_col(D->matrix, j, j, SIZE, D->Nx);
		if (v == 0)
		{
			delete temp_f;
			return -3;
		}
		temp_f[j] = temp_f[j]/v->value;
	}
	if (D->debug)
	{
		*(D->debug) << "func: " << std::endl;
		for (int j = 0; j < SIZE; j++)
		{
			*(D->debug) << std::setprecision(3) << temp_f[j] << "\t";
		}
		*(D->debug) << std::endl;
	}
	for (int n = 0; n < D->Nx + 1; n++)
	{
		double x = D->x_min + n*(D->x_max - D->x_min) / (D->Nx);
		**(D->prev_func_values + n) = func_y_min(x, D->curr_t);
		*(*(D->prev_func_values + n)+D->Ny) = func_y_max(x, D->curr_t);
	}
	for (int m = 0; m < D->Ny + 1; m++)
	{
		double y = D->y_min + m*(D->y_max - D->y_min) / (D->Ny);
		*(*(D->prev_func_values)+m) = func_x_min(y, D->curr_t);
		*(*(D->prev_func_values + D->Nx) + m) = func_x_max(y, D->curr_t);
	}
	for (int j = SIZE - 1; j > -1; j--)
	{
		int m = j / (D->Nx - 1);
		int n = j - m*(D->Nx - 1);
		*(*(D->prev_func_values + n+1) + m+1) = temp_f[j];
	}
	delete temp_f;
	return 0;
}

int print_matrix(std::ofstream &str, data D)
{
	int SIZE = (D.Nx - 1)*(D.Ny - 1);
	element *v;
	for (int row = 0; row <SIZE; row++)
	{
		for (int col= 0; col <= SIZE; col++)
		{
			v = search_row_col(D.matrix, row, col, SIZE, D.Nx);
			if (v)
			{
				if (abs(v->value)>1.0e-11) str << std::setprecision(3) << v->value << "\t";
				else
				{
					str << std::setprecision(3) << 0 << "\t";
				}
			}
			else str << std::setprecision(3) << 0 << "\t";
		}
		str << std::endl;
	}
	return 0;
}

double init_func(double x, double y)
{
	return 0;
}

int set_initial_values(data* D)
{
	for (int n = 0; n < D->Nx + 1; n++)
	{
		double x = D->x_min + n*(D->x_max - D->x_min) / (D->Nx);
		**(D->prev_func_values + n) = func_y_min(x, D->curr_t);
		*(*(D->prev_func_values + n) + D->Ny) = func_y_max(x, D->curr_t);
	}
	for (int m = 0; m < D->Ny + 1; m++)
	{
		double y = D->y_min + m*(D->y_max - D->y_min) / (D->Ny);
		*(*(D->prev_func_values) + m) = func_x_min(y, D->curr_t);
		*(*(D->prev_func_values + D->Nx) + m) = func_x_max(y, D->curr_t);
	}
	for (int n =1; n<D->Nx; n++)
	{
		for (int m =1; m < D->Ny; m++)
		{
			double x = D->x_min + n*(D->x_max - D->x_min) / (D->Nx);
			double y = D->y_min + m*(D->y_max - D->y_min) / (D->Ny);
			*(*(D->prev_func_values + n) + m) = init_func(x,y);
		}
	}
	return 0;
}

int main(void)
{
	int Nx, Ny,Nt;
	std::cout << "Enter X discretization number" << std::endl;
	do
	{
		std::cin >> Nx;
	} while ((Nx > 2) ? (0) : (std::cout << "value must be more than 2"<<std::endl, 1));
	std::cout << "Enter Y discretization number" << std::endl;
	do
	{
		std::cin >> Ny;
	} while ((Ny > 2) ? (0) : (std::cout << "value must be more than 2" << std::endl, 1));
	std::cout << "Enter T discretization number" << std::endl;
	do
	{
		std::cin >> Nt;
	} while ((Nt >30) ? (0) : (std::cout << "value must be more than 30" << std::endl, 1));
	data DATA;
	DATA.t_min = 0;
	DATA.x_min = -1;
	DATA.y_min = -1;
	DATA.t_max = 2;
	DATA.x_max = 1;
	DATA.y_max = 1;
	DATA.curr_t = DATA.t_min;
	DATA.Nt = Nt;
	DATA.Nx = Nx - 1;
	DATA.Ny = Ny - 1;
	DATA.matrix = 0;
	DATA.is_triag = 0;
	DATA.prev_func_values = new double*[Nx];
	for (int j = 0; j < Nx; j++)
	{
		*(DATA.prev_func_values + j) = new double[Ny];
	}
	set_initial_values(&DATA);
	std::ofstream str1, str2, str3, test1,test2;
	str1.open("ff.txt", std::ios_base::trunc);
	str2.open("center.txt", std::ios_base::trunc);
	//str3.open("time_3.txt", std::ios_base::trunc);
	test1.open("matrix1.txt", std::ios_base::trunc);
	//test2.open("matrix2.txt", std::ios_base::trunc);
	DATA.debug = &test1;
	for (int tt = 0; tt < Nt; tt++)
	{
		if (tt == 0)
		{
			DATA.debug = 0;
		}
		build(&DATA);
		triag_alg(&DATA);
		fill_func(&DATA);
		print_func(str1, DATA);
		print_center(str2,DATA);
	}
	str1.close();
	//str2.close();
	//str3.close();
	test1.close();
	//test2.close();
	for (int j = 0; j < Nx; j++)
	{
		if (*(DATA.prev_func_values+j)) delete *(DATA.prev_func_values+j);
	}
	if (DATA.prev_func_values) delete DATA.prev_func_values;
	if (DATA.matrix) delete DATA.matrix;
	std::cout << "Enter 0 to quit" << std::endl;
	std::cin >> Nx;
	return 0;
}
