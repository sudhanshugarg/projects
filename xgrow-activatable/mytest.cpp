#include<iostream>
#include<cstdlib>
#include<cstdio>
using namespace std;

int main(){
    double k00, k01, k10, k11;
    double r, sum;
    bool oops;
    int di=-1,dj=-1;

    r = drand48();
    r = drand48();
    for(int p=0; p<20; p++){
        getchar();
        k00 = drand48()*100;
        k01 = drand48()*100;
        k10 = drand48()*100;
        k11 = drand48()*100;
        sum = k00+k01+k10+k11;

        cout << "Start of Loop:"
            << p
            << ", k00=" << k00
            << ", k01=" << k01
            << ", k10=" << k10
            << ", k11=" << k11
            << ", di=" << di
            << ", dj=" << dj
            << ", r=" << r
            << ", sum=" << sum
            << endl;

        do {
            r = r*sum;  oops=0;
            if ( (r-=k00) < 0) { di=0; dj=0; r=(r+k00)/k00; } else
                if ( (r-=k10) < 0) { di=1; dj=0; r=(r+k10)/k10; } else
                    if ( (r-=k01) < 0) { di=0; dj=1; r=(r+k01)/k01; } else
                        if ( (r-=k11) < 0) { di=1; dj=1; r=(r+k11)/k11; } else
                        { r=drand48(); oops=1; }
        } while (oops);

    }
        cout << "End of For: 7"
            << ", k00=" << k00
            << ", k01=" << k01
            << ", k10=" << k10
            << ", k11=" << k11
            << ", di=" << di
            << ", dj=" << dj
            << ", r=" << r
            << ", sum=" << sum
            << endl;

}
