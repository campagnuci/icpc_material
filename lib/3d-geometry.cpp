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

struct pt3;
struct line3;

//funcao de comparacao
int cmp(double a, double b = 0.0){
	if(fabs(a-b) < EPS) return 0;
	return a > b ? 1 : -1;
}

////////////////////////////////////
/// Ponto 3D

struct pt3{
	double x, y, z;
	pt3(double x = 0.0, double y = 0.0, double z = 0.0) : x(x), y(y), z(z) {}
	
	double length(){ return sqrt(x*x + y*y + z*z); }
	double length2() { return x*x + y*y + z*z; }
	
	pt3 operator + (pt3 p) { return pt3(x + p.x, y + p.y, z + p.z); }
	pt3 operator - (pt3 p) { return pt3(x - p.x, y - p.y, z - p.z); }
	pt3 operator * (double k) { return pt3(x * k, y * k, z * k); }
	pt3 operator / (double k) { return pt3(x / k, y / k, z / k); }
	pt3 normalize() { return (*this)/length(); };
};

//distancia
double dist(pt3 a, pt3 b){ return (b-a).length(); }

//produto escalar
double dot(pt3 a, pt3 b){ return a.x*b.x + a.y*b.y + a.z* b.z; }

//produto vetorial
pt3 cross(pt3 a, pt3 b){ return pt3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }

//////////////////////////////
/// Reta 3D

struct line3{
	pt3 a, b;
	line3(pt3 a = pt3(), pt3 b = pt3()) : a(a), b(b) {}
	//direcao da reta - (nao normalizado)
	pt3 dir() { return (b-a); }
};

//ponto mais proximo de uma reta
//retorna o ponto da reta mais proximo de p
pt3 closest_point_line(pt3 p, line3 l){
	pt3 dir = l.dir();
	return l.a + dir*dot(p - l.a, dir)/dir.length2();
}


//distancia entre retas
//retorna a distancia minima entre duas retas
double dist(line3 r, line3 s){
	pt3 ort = cross(r.dir(), s.dir());
	if(!cmp(ort.length()))
		return dist(closest_point_line(r.a, s), r.a);
	return dot(s.a - r.a, ort)/ort.length();
}

//encontra o ponto mais proximo entre duas retas
//retorna o ponto em r mais proximo da resta s
//assume retas nao paralelas
bool closest_point_line_line(line3 r, line3 s, pt3& close){
	pt3 rdir = r.dir(), sdir = s.dir();
	double rr = dot(rdir, rdir);
	double ss = dot(sdir, sdir);
	double rs = dot(rdir, sdir);
	double t = dot(r.a - s.a, rdir)*ss - dot(r.a - s.a, sdir)*rs;
	
	if(!cmp(rs*rs - rr*ss)) return false; //retas paralelas
	
	t /= (rs*rs - rr*ss);
	close = r.a + rdir*t;
	return true;
}

//ponto mais proximo do segmento
//retorna o ponto do segmento mais proximo de p
pt3 closest_point_seg(pt3 p, line3 l){
	pt3 ldir = (l.b - l.a);
	double s = dot(l.a - p, ldir)/ldir.length2();
	if(s < -1.0) return l.b;
	if(s > 0.0) return l.a;
	return l.a - ldir*s;
}

///////////////////////////////////
/// Plano 3D

//define um plano no 3D
//p eh um ponto no plano
//n eh a normal do plano a partir de p
//(representacao util para calculos algebricos)
struct plane{
	pt3 n, p;
	plane(pt3 n = pt3(), pt3 p = pt3()) : n(n), p(p) {}
	plane(pt3 a, pt3 b, pt3 c) : n(cross(b-a, c-a)), p(a) {};
	//produto misto
	double d() { return -dot(n , p); }
};

//ponto do plano mais proximo de p
pt3 closest_point_plane(pt3 p, plane pl){
	return p - pl.n*(dot(pl.n, p - pl.p))/pl.n.length2();
}

//ponto de intersecao entre uma reta e um plano
//assume que reta nao eh paralela ao plano
bool intersect(line3 l, plane pl, pt3& inter){
	pt3 ldir = l.dir();
	
	if(!cmp(dot(pl.n, ldir))) return false; //reta paralela ao plano
	
	inter =  l.a - ldir*( (dot(pl.n, l.a) + pl.d())/dot(pl.n, ldir) );
	return true;
}

//intersecao de planos
//assume que os planos nao sao paralelos
bool intersect(plane u, plane v, line3& inter){
	pt3 p1 = u.n*(-u.d()), uv = cross(u.n, v.n);
	pt3 uvu = cross(uv, u.n);
	
	if(!cmp(dot(v.n, uvu))) return false; //planos paralelos
	
	pt3 p2 = p1 - uvu*((dot(v.n, p1) + v.d())/(dot(v.n, uvu)));
	inter.a = p2;
	inter.b = p2 + uv;
	return true;
}

//angulo entre dois vetores
double angle(pt3 a, pt3 b){
	return acos(dot(a, b)/(a.length()*b.length()));
}

//determina o volume formado pelo tetraedro delimitado
//pelos pontos a, b, c, d
double signed_volume(pt3 a , pt3 b, pt3 c, pt3 d){
	plane pl(b-a, c-a, d-a);
	return dot(cross(b-a, c-a), (d-a))/6.0;
	return pl.d()/6.0;
}

//area formada pelo triangulo a, b, c
double triangle_area(pt3 a, pt3 b, pt3 c){
	double h = dist(a, closest_point_line(a , line3(b, c)));
	return dist(b, c)*h/2.0;
}

//verifica se o ponto de esta no triangulo a, b, c
bool inside_triangle(pt3 a, pt3 b, pt3 c, pt3 d){
	if(cmp(signed_volume(a, b, c, d))) return false; //nao coplanares
	int left = 0, right = 0;
	int sign;
	sign = cmp(signed_area(a, b, d)); sign > 0 ? left++ : sign < 0 ? right++ : 0;
	sign = cmp(signed_area(b, c, d)); sign > 0 ? left++ : sign < 0 ? right++ : 0;
	sign = cmp(signed_area(c, a, d)); sign > 0 ? left++ : sign < 0 ? right++ : 0;
	return !(left && right);
}
void print(pt3 p){
	cout << "(" << p.x << "," << p.y << "," << p.z << ")";
}

int main(){
	cout << fixed << setprecision(6);
	
	//declara os pontos de dois planos
	pt3 a1(100, 0, 0), b1(0, 100, 0), c1(0, 0 ,0);
	pt3 a(100, 0, 0), b(0, 100, 0), c(0, 0 ,100);

	//cria os dois planos
	plane u(a1, b1, c1), v(a, b, c);
	
	//cria duas retas
	line3 l(pt3(0, 0, 0), pt3(0, 0, 10)), s(pt3(0, 2, 2), pt3(2, 0, 2));
	
	pt3 close;
	
	//calcula o ponto mais proximo entre duas retas
	closest_point_line_line(l, s, close); print(close); cout << endl;
	
	//calcula o ponto mais proximo na reta s
	print(closest_point_line(pt3(0, 2, 0), s)); cout << endl;
	
	//calcula a distancia entre as duas retas
	cout << dist(l, s) << endl;
	
	//calcula a intersecao de dois planos
	cout << intersect(u, v, l) << endl;
	print(l.a); cout << " "; print(l.b); cout << endl;
	
	//verifica se a solucao eh ortogonal as normais dos dois planos
	//dot = 0
	cout << dot(l.dir(), v.n) << endl;
	cout << dot(l.dir(), u.n) << endl;
	
	//calcula os angulos entre vetores coplanares
	//soma deve ser igual a pi neste caso (triangulo)
	cout << angle(a-c, b-c) << " " << angle(b-a, c-a) << " " << angle(c-b, a-b) << endl;
	cout << angle(a-c, b-c) + angle(b-a, c-a) + angle(c-b, a-b) << endl;
	
	//testa se os pontos da solucao estao certos
	//volume = 0 (pontos coplanares)
	cout << signed_volume(l.a, a, b, c) << endl;
	cout << signed_volume(l.b, a, b, c) << endl;
	cout << signed_volume(l.a, a1, b1, c1) << endl;
	cout << signed_volume(l.b, a1, b1, c1) << endl;
	
	//calcula a area do triangulo no espaco
	cout << signed_area(a, b, c) << endl;
	
	cout << inside_triangle(a, b, c, pt3(0, 100, 0));
}