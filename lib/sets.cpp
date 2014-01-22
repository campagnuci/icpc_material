#include <set>
#include <list>
#include <stack>
#include <cmath>
#include <queue>
#include <ctime>
#include <cfloat>
#include <vector>
#include <string>
#include <cstdio>
#include <bitset>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <utility>
#include <iostream>
#include <algorithm>

using namespace std;

#define FOR(i, a, b) for(int i = a; i <= b; ++i)
#define RFOR(i, b, a) for(int i = b; i >= a; --i)
#define REP(i, N) for(int i = 0; i < N; ++i)
#define RREP(i, N) for(int i = N-1; i >= 0; --i)
#define FORIT(i, a) for( TI(a) i = a.begin(); i != a.end(); i++ )
#define MAXN 10000
#define INF 0x3F3F3F3F
#define LINF 0x3F3F3F3FFFFFFFFFLL
#define FILL(X, V) memset( X, V, sizeof(X) )
#define TI(X) __typeof((X).begin())
#define ALL(V) V.begin(), V.end()
#define SIZE(V) int((V).size())
#define pb push_back
#define mp make_pair
#define pair < int, int > pii;


int n, grupo;

int pset[40000], vet[40000];

void initSet(){
	for(int i = 0; i < n; i++){ pset[i] = i; vet[i] = 1; }
}

int findSet(int i){
	int x = i;
	while(x != pset[x]) x = pset[x];
	int aux;
	while(i != x){
		aux = pset[i];
		pset[i] = x;
		i = aux;
	}
	return x ;
}

bool isSameSet (int i, int j){
	return ( findSet(i) == findSet(j) );
}

void unionSet(int i, int j){
	i = findSet(i);
	j = findSet(j);
	if(i == j) return ;
	pset[i] = j;
	vet[j] += vet[i];
}

int main(){
	int qt, no, qtNo;
	ios::sync_with_stdio(false);
	cin >> n >> grupo;
	while(n || grupo){
		initSet();
		for(int i = 0; i < grupo; i++) {
			cin >> qt;
			if(qt){
				cin >> qtNo;
				qt--;
				for(int j = 0; j < qt; j++){
					cin >> no;
					unionSet(qtNo,no);
				}
			}
		}
		int zero = findSet(0); 
		cout << vet[zero] << endl;//printf("%d\n", ans);
		cin >> n >> grupo;
	}
	return 0;
}
