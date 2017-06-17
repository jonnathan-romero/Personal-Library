// Jonnathan Romero, May 2017
// rand0m.cpp

#ifndef JONN_H
#define JONN_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <curl/curl.h> // -lcurl, (terminal)-> sudo apt-get install libcurl4-gnutls-dev
#include <ctime>  

#define ZERO 0.00000001

class Matrix{
public:

	double **data;

	int rows;
	int columns;

	Matrix(const int& r, const int& c);
	Matrix(const Matrix& copyFrom);
	void operator=(const Matrix& copyFrom);

	Matrix(Matrix&& moveFrom);
	void operator=(Matrix&& moveFrom);

	int getRows()const;
	int getColumns()const;

	void fill();
	void show();

	~Matrix();
};

Matrix scalarMultiple(const Matrix& A, const int& num);
Matrix scalarMultiple(const int& num, const Matrix& A);
Matrix transpose(const Matrix& A);
Matrix identity(const int& len);
Matrix invert(const Matrix& A);
Matrix cholesky(Matrix A);
Matrix add(const Matrix& A, const Matrix& B);
Matrix multiply(const Matrix& A, const Matrix& B);

double drop(const Matrix& A);                                                       // Given 1x1 Matrix returns double

double time();                                                                      // Timer

void   printToFile(std::string filename, const std::string &text);                  // Prints string to a file 
void   quandlData(const std::string &tickerName, const std::string &days);          // Download csv file with daily adj close prices
void   stockData(Matrix& returns, Matrix& cov);                                     // Create return and cov Matrix
void   yieldData(Matrix& Yield);                                                    // Download Yield Data

double uniform (unsigned int seed);                                                 // Uniform[0,1] 
double psiInv (double u);                                                           // Standard Normal associated with value [0,1] 
double psi (double x);                                                              // [0,1] associated with Standard Normal

//***********************************************************************************************************************************

double time () {
   static clock_t time = clock();

   return (clock() - time) / CLOCKS_PER_SEC;
}

void printToFile(std::string filename, const std::string &text){
	
	std::ofstream mfile;
	
	mfile.open(filename);
	
	for(int i=0;i<text.length();++i)
		mfile<<text[i];
	
	mfile.close();
}

static size_t writerF(void *ptr, size_t size, size_t nmemb, void *userdata){

	((std::string*)userdata)->append((char*)ptr, size * nmemb);

	return size * nmemb;
}

void quandlData(const std::string &tickerName, const std::string &days){
	
	std::string mainLink = "https://www.quandl.com/api/v3/datasets/";
	std::string database = "WIKI";
	std::string quandl_auth_token = "8N5sVKPAKEqu3FiMm5nz";

	mainLink+=database;
	mainLink+="/"+tickerName;
	mainLink+=".csv";
	mainLink+="?column_index=11";
	mainLink+="&rows="+days;
	mainLink+="&sort_order=asc&auth_token=";
	mainLink+=quandl_auth_token;
	
	CURL *curl;
	std::string quandlData;
	std::string fName=tickerName;
	fName+=".csv";

	curl = curl_easy_init();

	if(curl){
		const char* linkArrayChar=mainLink.c_str();
		curl_easy_setopt(curl, CURLOPT_URL, linkArrayChar);
		curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writerF);
		curl_easy_setopt(curl, CURLOPT_WRITEDATA, &quandlData);
		curl_easy_perform(curl);
		curl_easy_cleanup(curl);
		printToFile(fName,quandlData);
	}
}

// Yield.rows = 11
// Yield.columns = days
void yieldData(Matrix& Yield){

	int i,j;
	int col = Yield.getColumns();
	int row = Yield.getColumns();

	double temp;

	std::ifstream file;
	std::string str;

	str=std::to_string(Yield.getColumns());

	std::cout<<"Your matrix of Yields will contain "<<str<<" points(columns) of\n1 MO, 3 M, 6 MO, 1 YR, 2 YR, 3 YR, 5 YR, 7 YR, 10 YR, 20 YR, & 30 YR (rows)\n";

	std::string mainLink = "https://www.quandl.com/api/v3/datasets/";
	std::string database = "USTREASURY";
	std::string quandl_auth_token = "8N5sVKPAKEqu3FiMm5nz";

	mainLink+=database;
	mainLink+="/YIELD.csv?&rows=";
	mainLink+= str;
	mainLink+="&api_key=ZpSA483KUzTC8oDzPzAq&auth_token=";
	mainLink+=quandl_auth_token;

	CURL *curl;
	std::string quandlData;
	std::string fName="Yield.csv";	

	curl = curl_easy_init();

	if(curl){
		const char* linkArrayChar=mainLink.c_str();
		curl_easy_setopt(curl, CURLOPT_URL, linkArrayChar);
		curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writerF);
		curl_easy_setopt(curl, CURLOPT_WRITEDATA, &quandlData);
		curl_easy_perform(curl);
		curl_easy_cleanup(curl);
		printToFile(fName,quandlData);
	}

	str="Yield.csv";

	file.open(str,std::ifstream::in);

	if(file.good()){

		getline(file,str);

		for(j = 0; j<col; ++j){
			std::cout<<j<<": ";
			getline(file,str);
			str = (str.replace(str.begin(),str.begin()+11,""));
			str +=",";
			std::stringstream ss(str);
			for(i=0;ss >> temp;++i)
    		{
        		std::cout<<temp<<" ";
        		Yield.data[i][j]=temp;
        		if (ss.peek() == ',')
            	ss.ignore();
    		}
    		std::cout<<"\n";
		}
		file.close();

	}else{
		std::cout<<"File "<<str<<" could not be open!\n";
	}

	std::cout<<"\n\"-1\" to delete files and \"1\" to keep files... ";
  	std::cin>>str;
	
	if(str=="-1"){
		str ="Yield.csv";
		if (remove(str.c_str( )) !=0)
        	std::cout<<"Remove operation failed for "<<str.c_str( )<<"!\n";
	}
}
// returns.rows          = numOfStocks
// returns.columns       = days
// cov.rows = cov.colums = numOfStocks
void stockData(Matrix& returns, Matrix& cov){ 

	int i, j, k, len = returns.getRows(), days = returns.getColumns();
	long double average;

	Matrix tempReturns(returns);

	std::string str;

	std::vector <std::string> tickers;
	std::vector <long double> adjClose;

	std::ifstream file;

	std::cout<<"Input "<<len<<" tickers...\n\n";
	std::cout<<"For example...\nHD IBM AAPL NKE MMM UNH UTX PFE DIS JNJ MSFT GE KO WMT MRK PG V DD XOM MCD JPM BA CSCO AXP VZ INTC\n\n";

	for(i=0;i<len;++i){
		std::cin>>str;
		tickers.push_back(str);
	}

	std::cout<<"\nYour portfolio will include:\n";
  	for(i=0; i<len; ++i)
      std::cout<<tickers[i]<<(i<len-1?", ":"\n");

	std::cout<<"\nData for "<<days<<" days will be used... ";
    
	str = std::to_string(days);

	for(i=0; i<len; ++i)
		quandlData(tickers.at(i), str);
	
	
	for(i=0; i<len; ++i){

		average = 0;

		adjClose.clear();

		str = tickers[i]+".csv";

		file.open(str,std::ifstream::in);

		if(file.good()){

			getline(file,str);

			for(j = 0; j<days; ++j){
				getline(file,str);
				str = (str.replace(str.begin(),str.begin()+11,""));
				adjClose.push_back(atof (str.c_str()));
			}
			file.close();
		}else{
			std::cout<<"File "<<str<<" could not be open!\n";
		}

		for(j=0; j < days-1; ++j){
			returns.data[i][j]=log(adjClose[j]/adjClose[j+1]);
			average += returns.data[i][j];
    	}

    	average /= (days-1);

    	for(j=0; j < days-1; ++j){
			tempReturns.data[i][j]=returns.data[i][j]-average;
    	}
	}

	for(i = 0; i < len; ++i)
	    for(j = 0; j < len; ++j){
	    	for(k=0; k<days; ++k)
	    		cov.data[i][j] += tempReturns.data[i][k]*tempReturns.data[j][k];
	    	cov.data[i][j]/=(days-1);
	    	cov.data[i][j]*=252;
	    }
  	std::cout<<"\n\"-1\" to delete files and \"1\" to keep files... ";
  	std::cin>>str;
	
	if(str=="-1"){

		for(i=0;i<len;++i){

  			str = tickers[i]+".csv";

  			if (remove(str.c_str( )) !=0)
           		std::cout<<"Remove operation failed for "<<str.c_str( )<<"!\n";
  		}
	}
}

unsigned int temper (unsigned int N) {

   N ^= (N >> 11);
   N ^= (N << 7) & 0x9d2c5680;
   N ^= (N << 15) & 0xefc60000;
   N ^= (N >> 18);

   return N;
}

// seed = 1 to initialize
// seed = 0 for next integer
double uniform (unsigned int seed) {

	static unsigned int X[1248], m[2], initialized = 0, k;
	unsigned int N, Y;


	if (seed || !initialized) {
	 
		X[0] = (seed ? seed : 1);

		for (k = 1; k < 624; k++) {
			X[k] = 22695477 * X[k-1] + 1;
		}
		m[0] = 0; m[1] = 0x9908b0df;

		initialized = 1;
	}

	Y = (X[k-624] & 0x80000000) | (X[k-623] & 0x7fffffff);

	X[k] = ((Y >> 1) ^ m[Y & 1] ^ X[k-227]);

	X[k-624] = X[k];

	N = temper (X[k]);

	k ++;
	if (k == 1248) k = 624;

	return ( (N + 0.5) / 4294967296.0 );

}

double psi (double x) {

   static double    A =    0.0293868870,
                    B =    0.2161934875,
                    C =    0.6503029656,
                    D =    0.7978845608,
                    E =    0.0594864800,
                    F =    0.4160822924;

   int sign=1;
   double R;

   if (x < 0) {
      sign = -1;
      x *= -1;
   }

   x = fabs(x);
   R = (((A*x+B)*x+C)*x+D)*x / ((E*x+F)*x+1);

   return (0.5 + sign * 0.5 * (1.0 - exp (-R)));
}

double psiInv (double u) {

   static double
    A1 =  -3.969683028665376e+01,
    A2 =   2.209460984245205e+02,
    A3 =  -2.759285104469687e+02,
    A4 =   1.383577518672690e+02,
    A5 =  -3.066479806614716e+01,
    A6 =   2.506628277459239e+00,
    B1 =  -5.447609879822406e+01,
    B2 =   1.615858368580409e+02,
    B3 =  -1.556989798598866e+02,
    B4 =   6.680131188771972e+01,
    B5 =  -1.328068155288572e+01,
    C1 =  -7.784894002430293e-03,
    C2 =  -3.223964580411365e-01,
    C3 =  -2.400758277161838e+00,
    C4 =  -2.549732539343734e+00,
    C5 =   4.374664141464968e+00,
    C6 =   2.938163982698783e+00,
    D1 =   7.784695709041462e-03,
    D2 =   3.224671290700398e-01,
    D3 =   2.445134137142996e+00,
    D4 =   3.754408661907416e+00,
    P0 =  0.02425,
    P1 =  0.97575;

   double N, q, r;

   if (u < P0) {
      q = sqrt(-2*log(u));
      N = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
   }

   else if (u <= P1) {
      q = u - 0.5;
      r = q*q;
      N = (((((A1*r+A2)*r+A3)*r+A4)*r+A5)*r+A6)*q /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1);
   }

   else {
      q = sqrt(-2*log(1.0-u));
   	  N = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
   }

   return (N);

}

Matrix::Matrix(const int& r, const int& c){
	rows = r;
	columns = c;

	int i,j;

	if(getRows()>0 && getColumns()>0){
		data = new double*[getRows()];
	
	for (i = 0; i < getRows(); ++i) 
    	data[i] = new double[getColumns()];
	}else{
		rows=0; columns=0;
		printf("Invalid Dimensions in constructor.");
		
	}
}

Matrix::Matrix(const Matrix& copyFrom){
	rows = copyFrom.getRows();
	columns = copyFrom.getColumns();

	int i,j;

	data = new double*[getRows()];
	
	for (i = 0; i < getRows(); ++i) 
    	data[i] = new double[getColumns()];

    for (i = 0; i < getRows(); ++i) 
    	for (j = 0; j < getColumns(); ++j) 
        	data[i][j] = copyFrom.data[i][j];
}

void Matrix::operator=(const Matrix& copyFrom){

	int i,j;

	if(data!=nullptr){
		for (i = 0; i < getRows(); ++i)
        	delete data[i];
   		delete[] data;
	}

	rows = copyFrom.getRows();
	columns = copyFrom.getColumns();

	data = new double*[getRows()];
	
	for (i = 0; i < getRows(); ++i) 
    	data[i] = new double[getColumns()];

    for (i = 0; i < getRows(); ++i) 
    	for (j = 0; j < getColumns(); ++j) 
        	data[i][j] = copyFrom.data[i][j];
}

Matrix::Matrix(Matrix&& moveFrom){
	rows = moveFrom.getRows();
	columns = moveFrom.getColumns();

	data=moveFrom.data;

	moveFrom.data = nullptr;
}

void Matrix::operator=(Matrix&& moveFrom){
	rows = moveFrom.getRows();
	columns = moveFrom.getColumns();

	data=moveFrom.data;

	moveFrom.data = nullptr;
}

int Matrix::getRows()const{
	return rows;
}
int Matrix::getColumns()const{
	return columns;
}

void Matrix::fill(){
	int i, j;

	std::cout<<"("<<getRows()<<"x"<<getColumns()<<") fill it up\n";

	for(i=0;i<getRows();++i){
		for(j=0;j<getColumns();++j)
			std::cin>>data[i][j];
	}
}
void Matrix::show(){
	int i, j;

	printf("\n     ");
	for(j=0;j<getColumns();++j)
		printf("%i          ",(j+1));

	for(i=0;i<getRows();++i){
		std::cout<<"\n"<<(i+1)<<": ";

		for(j=0;j<getColumns();++j)
			printf("%8.4f   ", data[i][j]);;
	}
		
	std::cout<<"\n";
}

Matrix::~Matrix(){

	if(data!=nullptr){
		int i;
		for (i = 0; i < getRows(); ++i)
        	delete (data[i]);
   		delete[] (data);
	}
}

Matrix scalarMultiple(const Matrix& A, const int& c){

	int i, j;

	Matrix cA(A.getRows(),A.getColumns());

	for (i = 0; i < cA.getRows(); ++i) 
		for (j = 0; j < cA.getColumns(); ++j) 
			cA.data[i][j] = c * A.data[i][j];

	return cA;
}

Matrix scalarMultiple(const int& c, const Matrix& A){

	int i, j;

	Matrix cA(A.getRows(),A.getColumns());

	for (i = 0; i < cA.getRows(); ++i) 
		for (j = 0; j < cA.getColumns(); ++j) 
			cA.data[i][j] = c * A.data[i][j];

	return cA;
}
Matrix transpose(const Matrix& A){

	int i, j;
	Matrix At(A.getColumns(),A.getRows());

	for (i = 0; i < A.getColumns(); ++i) 
		for (j = 0; j < A.getRows(); ++j) 
			At.data[i][j] = A.data[j][i];

	return At;
}
Matrix identity(const int& len){

	Matrix I(len,len);
	int i, j;

	for(int i=0; i<I.getRows(); ++i)
		for(int j=0; j<I.getColumns(); ++j)
			I.data[i][j]=(i==j);

	return I;
}
Matrix invert(const Matrix& A0){

	int i, j, k, rmax;
	double cc;

	Matrix A = A0;
	Matrix Ainv(A.getRows(),A.getColumns());
	Ainv = identity(A.getRows());

	if(A.getRows()==A.getColumns()){

		for (j = 0; j < Ainv.getRows(); ++j) {

			cc = 0;
			rmax = 0;

			for (i = j; i < Ainv.getRows(); ++i)
				if (fabs(A.data[i][j]) > cc) {
					cc = fabs (A.data[i][j]);
					rmax = i;
				}
			
			if (cc >= ZERO) {

				i = j;
				for (k = 0; k < Ainv.getRows(); ++k) {
					cc = A.data[i][k];
					A.data[i][k] = A.data[rmax][k];
					A.data[rmax][k] = cc;
					cc = Ainv.data[i][k];
					Ainv.data[i][k] = Ainv.data[rmax][k];
					Ainv.data[rmax][k] = cc;
				}

				cc = A.data[i][i];
				for (k = 0; k < Ainv.getRows(); ++k) {
					A.data[i][k] /= cc;
					Ainv.data[i][k] /= cc;
				}

				for (i = 0; i < Ainv.getRows(); ++i) if (i != j) {
					cc = A.data[i][j];
					for (k = 0; k < Ainv.getRows(); ++k) {
						A.data[i][k] += -cc * A.data[j][k];
						Ainv.data[i][k] += -cc * Ainv.data[j][k];
					}
				}
				
			}else{
				printf ("Trying to invert a singular matrix.\n");
				Matrix Z(A.getRows(),A.getColumns());
				return Z;
			}
		}

		return Ainv;
	}else{
		printf ("Trying to invert a non-square matrix.\n");
		
		Matrix Z(A.getRows(),A.getColumns());
		return Z;
	}
}
Matrix cholesky(Matrix V){
   int n, i, j, k;
   Matrix L(V.getRows(),V.getColumns()); 
   double sum;

   // Check that the matrix V is square and symmetric.
   n = V.getRows();
    if (n != V.getColumns()) {
      printf ("In Cholesky function V is not square.\n");
      return L;
   }

   for (i = 0; i < n; i++) 
      for (j = 0; j < n; j++) 
         if (V.data[i][j] != V.data[j][i]) {
            printf ("In Cholesky function V is not symmetric.\n");
            return L;
         }

   for (k = 0; k < n; k++) {

      for (j = 0; j < k; j++) {
         sum = 0;
         for (i = 0; i < j; i++) {
            sum += L.data[k][i] * L.data[j][i];
         }
         L.data[k][j] = (V.data[k][j] - sum) / L.data[j][j];
     }

      sum = V.data[k][k];
      for (i = 0; i < k; i++) 
         sum -= L.data[k][i] * L.data[k][i];

      if (sum <= ZERO){
         printf ("In Cholesky function V is not positive definite.\n");
         Matrix Z(V.getRows(),V.getColumns());
         return Z;
      }

      L.data[k][k] = sqrt (sum);

   }

   return L;
}
Matrix add(const Matrix& A, const Matrix& B){

	int i, j;

	Matrix S(A.getRows(),B.getColumns());

	if(A.getRows()==B.getRows()&&A.getColumns()==B.getColumns()){

		for (i = 0; i < S.getRows(); ++i)
			for (j = 0; j < S.getColumns(); ++j)
		    	S.data[i][j]= A.data[i][j] + B.data[i][j];
	}else{
		printf("Dimensions don't match in matrix additions...\n");
	}

	return S;
}
Matrix multiply(const Matrix& A, const Matrix& B){

	Matrix AB(A.getRows(),B.getColumns());
	int i, j, k;

	if(A.getColumns()==B.getRows()){

		for (i = 0; i < A.getRows(); ++i) 
			for (j = 0; j < B.getColumns(); ++j) 
				for (k = 0; k < A.getColumns(); ++k) 
					AB.data[i][j] += A.data[i][k] * B.data[k][j];
	}else{
		printf ("Dimensions don't match in matrix multiplication.\n");
		
	}

	return AB;	
}
double drop(const Matrix& A){
	if(A.getRows()==1 && A.getColumns()==1)
		return A.data[0][0];
	printf ("Matrix is not 1x1...\n");
	
	return -999;
}
#endif // JONN_H
