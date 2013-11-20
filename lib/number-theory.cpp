#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <cstring>
#include <algorithm>
#define pb push_back
#define mp make_pair

//Fatorizacao - Funcao phi - Funcao sigma
 
using namespace std;

typedef pair<int, int> prime_pot;
typedef long long int64;

//aqui que dtermina o limite dos primos
const unsigned MAX = 20001000/60, MAX_S = sqrt(MAX/60);

unsigned w[16] = {1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59};
unsigned short composite[MAX];
vector<int> primes;

void sieve() {
   unsigned mod[16][16], di[16][16], num;
   for(int i = 0; i < 16; i++)
       for(int j = 0; j < 16; j++) {
           di[i][j] = (w[i]*w[j])/60;
           mod[i][j] = lower_bound(w, w + 16, (w[i]*w[j])%60) - w;
       }

   primes.push_back(2); primes.push_back(3); primes.push_back(5);

   for(unsigned i = 0; i < MAX; i++)
       for(int j = (i==0); j < 16; j++) {
           if(composite[i] & (1<<j)) continue;
           primes.push_back(num = 60*i + w[j]);

           if(i > MAX_S) continue;
           for(unsigned k = i, done = false; !done; k++)
               for(int l = (k==0); l < 16 && !done; l++) {
                   unsigned mult = k*num + i*w[l] + di[j][l];
                   if(mult >= MAX) done = true;
                   else composite[mult] |= 1<<mod[j][l];
               }
       }
}

//Fatorizacao do numero em fatores primos
void fact(int n, vector<prime_pot>& fat){
    if(n == 1){
        fat.pb(mp(1, 0));
        return;
    }
    int limit = (int)sqrt(n) + 1;
    prime_pot p;
    for(int i = 0; primes[i] <= limit; i++){
        if((n % primes[i]) == 0){
            p = mp(primes[i], 0);
            while((n % primes[i]) == 0){
                p.second++;
                n /= primes[i];
            }
            fat.pb(p);
            if(n == 1)
                return;
        }
    }

    if(n != 1)
        fat.pb(mp(n, 1));
}

//exponenciacao logaritmica
int intpow(int a, int n){
    int res = 1, pot = a;
    
    while(n){
        if(n & 1)
            res *= pot;
        pot *= pot;
        n >>= 1;
    }
    return res;  
}

//Funcao phi de Euler
//Quantidade de numeros menores que n tais que mdc(x, n) = 1
int phi(vector<prime_pot>& fat){
    int res = 1;
    vector<prime_pot>::iterator it;
    for(it = fat.begin(); it != fat.end(); it++){
        res *=  intpow(it->first, it->second - 1)*(it->first - 1);
    }
    return res;
}

//Quantidade de divisores de n
int div_cont(vector<prime_pot>& fat){
    int res = 1;
    vector<prime_pot>::iterator it;
    for(it = fat.begin(); it != fat.end(); it++){
        res *= it->second + 1;
    }
    return res;
}

//Funcao sigma
//Soma dos divisores de n
int sigma(vector<prime_pot>& fat){
    if(fat[0].first == 1) return 1;
    
    int res = 0;
    vector<prime_pot>::iterator it;
    for(it = fat.begin(); it != fat.end(); it++){
        res *= (intpow(it->first, it->second + 1) - 1)/(it->first - 1);
    }
    return res;
}

/////////////////////////////////////////////////
// Cycle finding - Algoritmo de Brent

//Helper function
int f(int x){
	//colocar a funcao aqui
	return 0;
}

//Brend Cycle Finding algorithm - (Sempre garanta que o algoritmo pare)
//Based in the tortoise and hare principle
//retorna o tamanho do ciclo
int cycle_find(int x0){
    int power = 1, len = 1;
        
    //tortoise = x0
    int tortoise = x0, hare = f(x0);
  
	while(tortoise != hare){
		if(power == len){
			tortoise = hare;
			power <<= 1;
			len = 0;
		}
		hare = f(hare);
		len += 1;
	}
	return len;
}

////////////////////////////////////////////////
// Euclides extendido

int gcd(int a, int b){
	while(b){
		a = a%b;
		swap(a, b);
	}
	return a;
}

//return d, x, y tais que d = a*x + b*y
//cobre o caso de numeros negativos
int64 gcd_extended(int64 a, int64 b, int64& x, int64& y){
	if(a < 0){
		int64 d = gcd_extended(-a, b, x, y);
		x = -x; return d;
	}
	if(b < 0){
		int64 d = gcd_extended(a, -b, x, y);
		y = -y; return d;
	}
	
	x = 1; y = 0;
	int64 nx = 0, ny = 1, q;
	
	while(b){
		q = a/b;
		x -= q*nx; swap(x, nx);
		y -= q*ny; swap(y, ny);
		a -= q*b; swap(a, b);
	}
	return a;
}

///////////////////////////////////////////////////
// Teorema do resto chines

//dados n equacoes modulares do tipo
//x = ai mod mi onde:
//ai eh relativamente primo com mi
//mi e mj sao relativamente primos se i != j
//
//Retorna a unica solucao modulo prod(mi)
//x = a[i] mod m[i]

/*int chinese_remainder(vector<int>& a, vector<int>& m){
	vector<int> y(a.size());
	
	int prod = 1, aux, res = 0, M;
	for(int i = 0; i < a.size(); i++){
		prod *= m[i];
	}
	
	for(int i = 0; i < a.size(); i++){
		gcd_extended((prod/m[i])%m[i], m[i], y[i], aux);
		res = (res + a[i]*y[i]*(prod/m[i])) % prod;
	}
	
	return res;
}*/

//vesao mais segura
int64 chinese_remainder(vector<int64>& a, vector<int64>& m){
	vector<int64> y(a.size()), M(a.size(), 1);
	
	int64 aux, res = 0, MM = 1;
	
	for(int i = 0; i < a.size(); i++){
		for(int j = 0; j < a.size(); j++){
			if(i != j)
			M[j] = (M[j]*m[i])%m[j];
		}
		//sem o mod para modulo prod(mi)
		MM *= m[i];
	}
	
	for(int i = 0; i < a.size(); i++){
		gcd_extended(M[i], m[i], y[i], aux);
		aux = a[i]*y[i] % MM;
		aux = aux * MM/m[i] % MM;
		res = (res + aux) % MM;
	}
	
	return res;
	
}

/////////////////////////////////////
// Equacoes Diofantinas

//Resolve uma equacao diofantina com duas variaveis ax + by = c
//
//Retorna a solucao se existir
//
//x = x0 + n*dx y = y0 + n*dy
bool diophantus(int64 a, int64 b, int64 c, int64& x0, int64& dx, int64& y0, int64& dy){
	int64 d = gcd_extended(a, b, x0, y0);
	if(c % d) return false;
	int64 q = c/d;
	x0 *= q;  y0 *= q;
	dx = b/d; dy = -a/d;
	return true;
}

///////////////////////////////////////
//Equacoes diofantinas gerais

//dada a equacao a1*x1 + a1*x1 + .. + an*xn = c
//encontrar uma solucao se possivel
//
//a solucao sera dada na forma de sol + t1*bas1 + t2*bas2 + .. + tn-1*basn-1
//WARNING: nenhum a pode ser zero!!!

//resolve a1*x1 + .. + an*xn = gcd(x1, .., xn)
int64 gcd_solve(vector<int64>& a, int n, vector<int64>& sol, vector<vector<int64> >& bases){
	int64 ret;
	
	if(n == 2){
		ret = gcd_extended(a[0], a[1], sol[0], sol[1]);
		bases[0][0] = a[1]/ret;
		bases[0][1] = -a[0]/ret;
		return ret;
	}
	
	int64 v, d = gcd_solve(a, n-1, sol, bases);
	ret = gcd_extended(d, a[n-1], v, sol[n-1]);
	
	for(int i = 0; i < n-1; i++){
		bases[n-2][i] = a[n-1]/ret * sol[i];
		sol[i] *= v;
	}
	bases[n-2][n-1] = -d/ret;
	
	return ret;
}

//resolve a equacao diofantina propriamente dita
bool diophantus(vector<int64> a, int64 c, vector<int64>& sol, vector<vector<int64> >& bases){
	bases.resize(a.size());
	sol.resize(a.size());
	for(int i = 0; i < a.size()-1; i++){
		sol[i] = 0;
		bases[i].resize(a.size());
		for(int j = 0; j < a.size(); j++)
			bases[i][j] = 0;
	}
	
	int64 d = gcd_solve(a, a.size(), sol, bases);
	
	if( c % d ) return false;
	
	int64 q = c/d;
	for(int i = 0; i < a.size(); i++)
		sol[i] *= q;
	return true;
}
//////////////////////////////////////
// Integrais

double f(double x){
	//colocar funcao aqui
	return pow(x, 4) + pow(x, 3) + pow(x, 2) + x;
}

//metodo de integracao de f(x) no intervalo [a,b] O(k)
double simpson(double a, double b, int k){
	double dx, x, t;
	dx = (b-a)/(2.0*k);
	t = 0;

	for( int i=0; i<k; i++ ) {
		t += (i==0 ? 1.0 : 2.0) * f( a+ 2.0*i*dx );
		t += 4.0 * f(a + (2.0*i+1.0)*dx);
	}
	t += f(b);

	return t * (b-a)/6.0/k;
}

int main(){
	
}
    
    