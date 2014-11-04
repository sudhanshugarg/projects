#include<iostream>
using namespace std;
int main(){
    string s;
    int i;
    cin >> i >> s;
    if(i==9){
        cout << "h" << i << "toehold " << s.substr(47,8) << endl;
        cout << "h" << i << "loop " << s.substr(12,25) << endl;
        cout << "h" << i << "hairpin" << s.substr(0,55) << endl;
    }
    else
    if(i%2){
        cout << "h" << i << "toehold " << s.substr(44,8) << endl;
        cout << "h" << i << "loop " << s.substr(12,22) << endl;
        cout << "h" << i << "hairpin" << s.substr(0,52) << endl;
    }
    else{
        cout << "h" << i << "toehold " << s.substr(24,8) << endl;
        cout << "h" << i << "loop " << s.substr(42,22) << endl;
        cout << "h" << i << "hairpin" << s.substr(24,52) << endl;
    }
}
