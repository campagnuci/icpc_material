//Crivo otimizado

const unsigned MXX = 300000000/60, MAX_S = sqrt(MXX/60);

unsigned w[16] = {1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59};
unsigned short composite[MXX];
vector<int> primes;

void sieve() {
    unsigned mod[16][16], di[16][16], num;
    for(int i = 0; i < 16; i++)
        for(int j = 0; j < 16; j++) {
            di[i][j] = (w[i]*w[j])/60;
            mod[i][j] = lower_bound(w, w + 16, (w[i]*w[j])%60) - w;
        }

    primes.push_back(2);
    primes.push_back(3);
    primes.push_back(5);

    memset(composite, 0, sizeof composite);
    for(unsigned i = 0; i < MXX; i++)
        for(int j = (i==0); j < 16; j++) {
            if(composite[i] & (1<<j)) continue;
            num = 60*i + w[j];
            primes.push_back( num );

            if(i > MAX_S) continue;
            for(unsigned k = i, done = false; !done; k++)
                for(int l = (k==0); l < 16 && !done; l++) {
                    unsigned mult = k*num + i*w[l] + di[j][l];
                    if(mult >= MXX) done = true;
                    else composite[mult] |= 1<<mod[j][l];
                }
        }
}

//Crivo normal

vector<int> primes;
void sieve(int lim){
	bitset<INSERT_LIM> isp;
	primes.pb(2);
	for(int i = 3; i*i <= lim; i += 2){
		if(!isp[i]){
			primes.pb(i);
			for(int j = i+i; j <= lim; j += i){
				isp[i] = 1;
			}
		}
	}
	
	int val = sqrt(lim) + 1;
	if(!(n & 1)) val++;
	
	for(int i = val; i <= lim; i += 2){
		if(!isp[i]) primes.pb(i);
	}
}

//is prime

bool isPrime(int number){
	if(number < 2) return false;
	if(number == 2) return true;
	
	for(int i = 3; i*i <= number; i++){
		if(number % i == 0) return true;
	}
	
	return false;
}