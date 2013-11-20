#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>
#include <vector>
#include <algorithm>
#include <ctime>

using namespace std;

const double pi = acos(-1.0);
const double EPS = 1e-9;
const double INF = 1e50;

struct pt;
typedef pair<pt, pt> line;
typedef vector<pt> polygon;

//funcao de comparacao
int cmp(double a, double b = 0.0){
	if(fabs(a-b) < EPS) return 0;
	return a > b ? 1 : -1;
}

//estrutura que representa um ponto
// ou um vetor dependendo da necessidade
struct pt{
	double x, y;
	
	pt(double x = 0.0, double y = 0.0) : x(x), y(y) {}
	
	double length() { return sqrt(x*x + y*y); }
	double length2() { return x*x + y*y; }
	
	pt normalize(){ return (*this)/length(); }
	
	pt operator - (pt p){ return pt(x - p.x, y - p.y); }
	pt operator + (pt p){ return pt(x + p.x, y + p.y); }
	pt operator * (double k){ return pt(x * k, y * k); }
	pt operator / (double k){ return pt(x / k, y / k); }
	bool operator < (const pt& p) const{
		if(fabs( x - p.x ) >= EPS) return x < p.x;
		return y < p.y;
	}
	bool operator == (const pt p){
		return ( fabs(x - p.x) < EPS && fabs(y - p.y) < EPS );
	}
};

////////////////////
// Funcoes basicas

double dist(pt a, pt b){ return (a - b).length(); }
double dot(pt a, pt b){ return a.x*b.x + a.y*b.y; }
double cross(pt a, pt b){ return a.x*b.y - a.y*b.x; }
double signed_area(pt a, pt b, pt c){ return cross((a-c),(b-c))/2.0; }

/////////////////////
// Orientacao

//Determina o lado que c esta em relacao ao vetor a->b
//1  -> left
//0  -> on
//-1 -> right
int side_sign(pt a, pt b, pt c){ return cmp(signed_area(a, b, c)); }
bool cw(pt a, pt b, pt c){ return cmp(signed_area(a, b, c)) > 0; }
bool ccw(pt a, pt b, pt c){	return cmp(signed_area(a, b, c)) < 0; }
bool colinear(pt a, pt b, pt c){ return !cmp(signed_area(a, b, c)); }

//testa se c esta na caixa limitada por a e b
bool in_box(pt a, pt b, pt c){
	return ( cmp(c.x, min(a.x, b.x)) >= 0 && cmp(c.x, max(a.x, b.x)) <= 0
		&&   cmp(c.y, min(a.y, b.y)) >= 0 && cmp(c.y, max(a.y, b.y)) <= 0 );
}

//////////////////////////////
// Ponto mais proximo

//determina o ponto mais proximo de c
//na reta ab
pt closest_point(pt a, pt b, pt c){
	//se colinares
	if( fabs( signed_area(a, b, c) ) < EPS) return c;
	pt dir = b-a;
	return a + ( dir*dot(dir, c-a)/dir.length2() );
}

//determina o ponto mais proximo de c
//no segmento ab
pt closest_point_seg(pt a, pt b, pt c){
	if( fabs( signed_area(a, b, c) ) < EPS) return c;
	pt close = closest_point(a, b, c);
	if( in_box(a, b, close) ) return close;
	return dist(a, c) < dist(b, c) ? a : b;
}

////////////////////////////////
//Angulos

//retorna o angulo entre os vetores a->b e a->c
//angulo com sinal
double angle(pt a, pt b, pt c){
	pt va = b-a, v = c-a;
	//va.normalize(); **nao parece ser necessario
	pt vb(-va.y, va.x);
	return atan2( dot(v, vb), dot(v, va) );
}

//retorna o angulo entre os vetorre a-> e a->c
//angulo sem sinal
double angle_2(pt a, pt b, pt c){
	pt va = b-a, vb = c-a;
	return acos(dot(va, vb) / (va.length() * vb.length()));
}

//Lei dos cossenos determina - o angulo A
double angle(double a, double b, double c){
	return acos( (b*b + c*c - a*a)/(2.0*b*c) );
}

//gira a em torno da origem por theta radianos
pt rotate(pt a, double theta){ return pt(cos(theta)*a.x - sin(theta)*a.y, sin(theta)*a.x + cos(theta)*a.y); }

//gira b em torno de a por theta radianos
pt rotate(pt a, pt b, double theta){ return a + rotate(b-a, theta); }

//reflete c em torno de a->b
pt reflect(pt a, pt b, pt c){
	double theta = angle(a, b, c);
	return rotate(c-a, -2.0*theta) + a;
}


///////////////////////////
// Intersecoes

//verifica se o ponto c se encontra no segmento ab
bool point_and_seg(pt a, pt b, pt c){
	if( !colinear(a, b, c) ) return false;
	return in_box(a, b, c); 
}

//interseccao entre as retas a->b e c->d e guarda em inter
bool intersect(pt a, pt b, pt c, pt d, pt& inter){
	double det = cross(b-a, d-c);
	if(fabs(det) < EPS){
		if( fabs(signed_area(a, b, c) ) < EPS){
			inter = pt(INF, INF); 
			return true; //retas coincidentes
		}
		return false; //retas paralelas
	}
	//retas concorrentes
	double ua = ( cross(d-c, a-c) )/det; //parametros (pode ser importante)
//  double ub = ( cross(b-a, a-c) )/det;
	inter = a + (b-a)*ua;
	return true;
}

//testa se exsite interseccao entre os dois segmentos
bool intersect_seg(pt a, pt b, pt c, pt d, pt& inter){
	if( !intersect(a, b, c, d, inter) ) return false; //segmentos paralelos
	
	if( inter == pt(INF, INF) )
		return in_box(a, b, min(c, d)) || in_box(c, d, min(a, b)); //segmentos sobrepostos
		
	return in_box(a, b, inter) && in_box(c, d, inter); //segmentos concorrentes
}

//intersecao entre reta e circunferencia
bool intersect(pt center, double r, pt a, pt b, pt& r1, pt& r2){
	pt close = closest_point(a, b, center);
	double x = dist(center, close);
	
	if( cmp(x, r) > 0 ) return false;
	
	double d = sqrt(r*r - x*x);
	
	pt v = (b - a);
	r1 = close - (v*d/v.length());
	r2 = close + (v*d/v.length());
	return true;
}

//intersecao de ciucunferencias
bool intersect(pt c1, double r1, pt c2, double r2, pt& p1, pt& p2){
	double d = dist(c1, c2);
	if(r1 < r2) { swap(c1, c2); swap(r1, r2); }
	if( cmp(d, r1+r2) > 0 || cmp(d, r1-r2) < 0) return false; //nao tem intersecao
	if( !cmp(d) && !cmp(r1, r2) ) return true; //sao iguais
	
	pt v = c2 - c1;
	p1 = c1 + v*r1/v.length();
	if( !cmp(d-r1-r2) || !cmp(d+r2-r1) ){ p2 = p1; return true; } //tengencia interna/externa
	
	double theta = angle(r2, d, r1);
	p2 = rotate(c1, p1, theta);
	p1 = rotate(c1, p1, -theta);
	return true;
}

///////////////////////////////////
// Areas, perimetros e convexidade

double trap(pt& a, pt& b){ return 0.5*(b.x - a.x)*(b.y + a.y); }

double area(polygon& poly){
	double ret = 0.0;
	for(int i = 0; i < poly.size(); i++){
		ret += trap(poly[i], poly[(i+1)%poly.size()]);
	}
	return fabs(ret);
}

//Determina se o poligono simples eh convexo
bool is_convex(polygon& p){
	int left = 0, right = 0, side;
	for(int i = 0; i < p.size(); i++){
		side = side_sign(p[i], p[(i+1)%p.size()], p[(i+2)%p.size()]);
		if(side < 0) right++;
		if(side > 0) left++;
	}
	return !(left && right);
}

double perimeter(polygon p){
	double per = 0.0;
	for(int i = 0; i < p.size(); i++){
		per += dist(p[i], p[(i+1)%p.size()]);
	}
	return per;
}

//testa se o ponto esta no poligono convexo 
bool inside_convex_poly(pt p, polygon& poly){
	int left = 0, right = 0, side;
	for(int i = 0; i < poly.size(); i++){
		side = side_sign(p, poly[i], poly[(i+1)%poly.size()]);
		if(side < 0) right++;
		if(side > 0) left++;
	}
	return !(left && right);
}

//testa se o ponto esta dentro de um poligono (nao necessariamente convexo)
bool inside_poly(pt p, polygon poly){
	poly.push_back(poly[0]);
	
	for(int i = 0; i < poly.size()-1; i++)
		if(point_and_seg(poly[i], poly[i+1], p)) return true; //na borda
	
	for(int i = 0; i < poly.size()-1; i++) poly[i] = poly[i] - p;
	p = pt(0, 0);
	
	double theta, y;
	
	while(true){
		theta = (double)rand()/10000.0;
		
		bool inter = false;
		//evita que um ponto fique no eixo x
		for(int i = 0; i < poly.size()-1; i++){
			poly[i] = rotate(poly[i], theta);
			if( !cmp(poly[i].x) ) inter = true;
		}
		
		if( !inter ){
			poly[poly.size()-1] = poly[0];
			//testa as possiveis intersecoes
			for(int i = 0; i < poly.size()-1; i++){
				if( cmp( poly[i].x * poly[i+1].x ) < 0 ){
					y = poly[i+1].y - poly[i+1].x * (poly[i].y - poly[i+1].y) / (poly[i].x - poly[i+1].x);
					if( cmp(y) > 0 ) inter = !inter; //se interecao valida
				}
			}
			return inter; //testa a paridade da semi-reta vertical partindo de p
		}
	}
	return true;
}

/////////////////////////////////////////
//Triangulos

double triangle_area(pt a, pt b, pt c){
	return fabs(signed_area(a, b, c));
}

//altura de a para a linha bc
double height(pt a, pt b, pt c){
	return 2.0*triangle_area(a, b, c)/dist(b, c);
}

//Encontro das alturas
bool hcenter(pt a, pt b, pt c, pt& hc){
	if( !cmp(triangle_area(a, b, c)) ) return false;
	pt p1 = closest_point(b, c, a);
	pt p2 = closest_point(a, c, b);
	intersect(a, p1, b, p2, hc);
	return true;
}

//Centro do circulo circunscrito
bool ccenter(pt a, pt b, pt c, pt& cc){
	if( !cmp(triangle_area(a, b, c)) ) return false;
	pt r1 = (b + c)*0.5;
	pt r2 = (a + c)*0.5;
	
	pt s1(r1.x - (c.y - b.y), r1.y + (c.x - b.x));
	pt s2(r2.x - (c.y - a.y), r2.y + (c.x - a.x));
	
	return intersect(r1, s1, r2, s2, cc);
}

//Encontro das bissetrizes
bool bcenter(pt a, pt b, pt c, pt& bc){
	if( !cmp(triangle_area(a, b, c)) ) return false;
	double s1 = dist(b, c), s2 = dist(a, c), s3 = dist(a, b);
	
	double rt1 = s2/(s2 + s3), rt2 = s1/(s1 + s3);
	
	pt p1 = b*rt1 + c*(1.0 - rt1);
	pt p2 = a*rt2 + c*(1.0 - rt2);
	
	return intersect(a, p1, b, p2, bc);
}

///////////////////////////////////////////
// Centroide (centro de massa do poligono)

pt centroid(polygon p){
	double a = area(p);
	double xc = 0.0, yc = 0.0;
	
	for(int i = 0; i < p.size(); i++){
		int next = (i+1)%p.size();
		xc += (p[i].x + p[next].x)*(p[i].x*p[next].y - p[next].x*p[i].y);
		yc += (p[i].y + p[next].y)*(p[i].x*p[next].y - p[next].x*p[i].y);
	}
	
	return pt(xc/(6.0*a), yc/(6.0*a));
}

///////////////////////////////////////////////////////////
//Fecho convexo (Graham Scan - Nao pega pontos colineares)

pt refer;
bool cmp_angle(pt p1, pt p2){
	double det = signed_area(refer, p1, p2); 
	if(fabs(det) < EPS){
		return (dist(refer, p1) >= dist(refer, p2));
	}
	return (det > EPS);
}

void convex_hull(polygon in, polygon& hull){
	hull.clear();
	
	if(in.size() < 3){ hull = in; return; }
	
	int pos = 0;
	for(int i = 1; i < in.size(); i++) if(in[i] < in[pos]) pos = i;
	swap(in[0], in[pos]);
	refer = in[0];
	
	sort(in.begin() + 1, in.end(), cmp_angle);
	in.resize(unique(in.begin(), in.end()) - in.begin());
	
	hull.push_back(in[0]); hull.push_back(in[1]);
	
	in.push_back(in[0]); //isto evita pontos colineares no final do poligono
	for(int i = 2; i < in.size(); ){
		if(hull.size() > 2 && side_sign(hull[hull.size() - 2], hull[hull.size() - 1], in[i]) <= 0){
			hull.pop_back();
		}
		else hull.push_back(in[i++]);
	}
	//tira a duplicata
	hull.pop_back();
}

//////////////////////////////////////////////////////////////////
//Fecho convexo (Lower hull - upper hull) - Pega pontos colineares
//(começa do ponto maior e em cw)
//para inverter para ccw e ponto menor inverter upper hull e lower hull

// Returns a list of points on the convex hull in counter-clockwise order.
// NOTE: the last point in the returned list is the same as the first one.
void convex_hull_2(polygon P, polygon& hull) {
	hull.clear();
	
	// Sort points lexicographically
	sort(P.begin(), P.end());
	P.resize(unique(P.begin(), P.end()) - P.begin());
	
	// Build lower hull
	for (int i = 0; i < P.size(); i ++) {
		while (hull.size() >= 2 && side_sign(hull[hull.size() - 2], hull[hull.size() - 1], P[i]) <= 0)
			hull.pop_back();
		hull.push_back(P[i]);
	};
 
	// Build upper hull
	for (int i = P.size()-2, t = hull.size() + 1; i >= 0; i --) {
		while (hull.size() >= t && side_sign(hull[hull.size()-2], hull[hull.size()-1], P[i]) <= 0)
			hull.pop_back();
		hull.push_back(P[i]);
	};
}

///////////////////////////////////
// Operacoes com poligonos

//Corta o poligono pela reta ab
//pol1 - poligono do lado esquerdo
//pol2 - poligono do lado direito
void cut_polygon(polygon pol, pt a, pt b, polygon& pol1, polygon& pol2){
	polygon pp, pn;
	pt p1, p2, r;
	
	for(int i = 0; i < pol.size(); i++){
		p1 = pol[i]; p2 = pol[(i+1)%pol.size()];
		int side = side_sign(a, b, p1);
		if(side >= 0) pp.push_back(p1);
		if(side <= 0) pn.push_back(p1);
		
		if(intersect(a, b, p1, p2, r) && (r == pt(INF, INF))){
			if(point_and_seg(p1, p2, r)){
				pp.push_back(r);
				pn.push_back(r);
			}
		}
		if(pp.size() > 2) convex_hull(pp, pol1);
		if(pn.size() > 2) convex_hull(pn, pol2);
	}
}

//intersecao de dois poligonos convexos
bool intersect(polygon a, polygon b, polygon& inter){
	polygon pts; pt p;
	
	//pontos dentro da regiao
	for(int i = 0; i < a.size(); i++)
		if(inside_poly(a[i], b))pts.push_back(a[i]);
	for(int i = 0; i < b.size(); i++)
		if(inside_poly(b[i], a)) pts.push_back(b[i]);
	
	//pontos gerados pela intersecao de segmentos
	for(int i = 0; i < a.size(); i++)
	for(int j = 0; j < b.size(); j++){
		if(intersect_seg(a[i], a[(i+1)%a.size()], b[j], b[(j+1)%b.size()], p) && !(p == pt(INF, INF)) ){
			pts.push_back(p);
		}
	}
	//sem poligono de intersecao
	if(pts.size() <= 1) return false;
	//faz o fecho convexo da solucao
	convex_hull(pts, inter); return true;
}


/////////////////////////////////////////////////////////////////////
// Great Circle Distance - Distância entre dois pontos em uma esfera

//calcula a great circle distance em uma esfera de raio r
//angulos em radianos
double great_circle_distance(double lat1, double lon1, double lat2, double lon2, double r){
	return r*acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2 - lon1));
}


int main(){
	cout << fixed << setprecision(20);
	pt a(354, 287587), b(10000, 2340), c(1435, 2444.45535);
	double theta1 = angle(a, b, c);
	double theta2 = angle_2(a, b, c);
	
	cout << theta1 << endl << theta2 << endl;
}