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

struct tri{
	int atual, tipo;
	tri ( int atual = 0, int tipo = 0) : atual(atual), tipo(tipo) { }
};

int n, grupo, g = 1;
typedef vector < int > vi;
typedef vector < vi > vii;
typedef vector < tri > vtri;
typedef vector < vtri > vvtri;
typedef long long int64;
typedef unsigned long long uint64;

int x[] = {-1,-1,-1, 0, 1, 1, 1, 0 };
int y[] = {-1, 0, 1, 1, 1, 0,-1,-1 };

vii grafo;

int visitado[40000], vis = 1, ans = 1;

void dfs_visit( int v ){
	REP(i, SIZE(grafo[v])){
		if(visitado[grafo[v][i]] != vis){
			visitado[grafo[v][i]] = vis;
			if(grafo[v][i] < n ) ans++; 
			dfs_visit( grafo[v][i] );
		}
	}
}

int main(){
	int qt, no, qtNo;
	ios::sync_with_stdio(false);
	while(cin >> n >> grupo && (n + grupo)){
		grafo.resize(n + 2 * grupo);
		qtNo = n;
		REP(i,grupo){
			cin >> qt;
			REP(j, qt){
				cin >> no;
				grafo[no].pb(qtNo);
				grafo[qtNo].pb(no);
			}
			qtNo++;
		}
		visitado[0] = vis;
		dfs_visit(0);
		cout << ans << endl;
		vis++;
		ans = 1;
		grafo.clear();
	}
	return 0;
}