#include<iostream>
#include<string>
#include<map>
#include<vector>
#include<string.h>
#include<sstream>
#include<fstream>
#include<unistd.h>
#include<set>
#include<cstdlib>
#include<algorithm>
//##include<ctype>
#include "STRAND.h"
using namespace std;

class BIT;

string gate (const string &s, const int &a){
    ostringstream oss;
    oss << "G";
    oss << s;
    oss << a;
    return oss.str();
}

string intermediate (
        const string &input1, 
        const string &input2, 
        const string &output, 
        const int &i1,
        const int &i2,
        const int &o){
    ostringstream oss;
    oss << input1 << i1;
    oss << DNA::sep;
    oss << input2 << i2;
    oss << DNA::sep;
    oss << "G" << output << o;
    return oss.str();
}

string getFanOutDomainName(const string &domain, const int &fanOutPosn){
    ostringstream oss;
    oss << domain << "_" << fanOutPosn << "_";
    return oss.str();
}

string fanOutOneTime(const string &domain, const int &multiplier){
    ostringstream oss;
    if(multiplier > 1){
        for(int j=0;j<2;j++){
            int count = 1;
            oss << domain << j << " -> {uni_forward} " 
                << getFanOutDomainName(domain, count) << j;
            for(int i=1;i<multiplier;i++){
                count++;
                oss << " + " 
                << getFanOutDomainName(domain, count) << j;
            }
            oss << " |" << endl;
        }
    }
    return oss.str();
}

string fanOutANDReusable(const string &domain, const int &multiplier, const int &conc){
    if(multiplier <= 1) return "";

    ostringstream oss;
    for(int j=0;j<2;j++){

        //Set of forward reactions.
        string dname = getFanOutDomainName(domain, 1);
        oss << "init "
            << gate(dname,j) << DNA::sep 
            << dname << j << DNA::sep
            << dname << DNA::sep2 << j
            << " " << conc << " |"
            << endl;
        oss << domain << j << " + "
            << gate(dname,j) << DNA::sep 
            << dname << j << DNA::sep
            << dname << DNA::sep2 << j
            << " -> "
            << domain << j << DNA::sep
            << gate(dname,j) << " + "
            << dname << j << " + "
            << dname << DNA::sep2 << j
            << " |" << endl;

        string next_dname;
        for(int k=1; k<multiplier-1; k++){
            next_dname = getFanOutDomainName(domain, k+1);

            oss << "init "
                << gate(next_dname,j) << DNA::sep 
                << next_dname << j << DNA::sep
                << next_dname << DNA::sep2 << j
                << " " << conc << " |"
                << endl;
            oss << dname << DNA::sep2 << j << " + "
                << gate(next_dname, j) << DNA::sep
                << next_dname << j << DNA::sep
                << next_dname << DNA::sep2 << j
                << " -> "
                << dname << DNA::sep2 << j << DNA::sep
                << gate(next_dname, j) << " + "
                << next_dname << j << " + "
                << next_dname << DNA::sep2 << j
                << " |" << endl;
            dname = next_dname;
        }
        //Last species.
        dname = getFanOutDomainName(domain, multiplier-1);
        next_dname = getFanOutDomainName(domain, multiplier);
        oss << "init "
            << gate(next_dname,j) << DNA::sep 
            << next_dname << j
            << " " << conc << " |"
            << endl;
        oss << dname << DNA::sep2 << j << " + "
            << gate(next_dname,j) << DNA::sep 
            << next_dname << j
            << " -> "
            << dname << DNA::sep2 << j << DNA::sep
            << gate(next_dname,j) << " + "
            << next_dname << j 
            << " |" << endl;


        //Set of reverse reactions.
        dname = getFanOutDomainName(domain, 1);

        oss << gate(domain,j) << " + "
            << domain << j << DNA::sep
            << gate(dname,j)
            << " -> "
            << gate(domain, j) << DNA::sep
            << domain << j << " + "
            << gate(dname,j)
            << " |" << endl;
        for(int k=1; k<multiplier; k++){
            next_dname = getFanOutDomainName(domain, k+1);
            oss << gate(dname, j) << DNA::sep
                << dname << j << " + "
                << dname << DNA::sep2 << j << DNA::sep
                << gate(next_dname, j)
                << " -> "
                << gate(dname, j) << DNA::sep
                << dname << j << DNA::sep
                << dname << DNA::sep2 << j << " + "
                << gate(next_dname, j)
                << " |" << endl;

            dname = next_dname;
        }
    }
    return oss.str();
}


//Initialize global variables
void initialize(void){
    getGatePrefix[AND_MOTIF] = "AND";
    getGatePrefix[NAND_MOTIF] = "NAND";
    getGatePrefix[OR_MOTIF] = "OR";
    getGatePrefix[NOR_MOTIF] = "NOR";
    getGatePrefix[BIT_MOTIF] = "BIT";
    getGatePrefix[NOT_MOTIF] = "NOT";

    getTruthTable[AND_MOTIF] = TT_AND;
    getTruthTable[NAND_MOTIF] = TT_NAND;
    getTruthTable[OR_MOTIF] = TT_OR;
    getTruthTable[NOR_MOTIF] = TT_NOR;
    getTruthTable[BIT_MOTIF] = NULL;
    getTruthTable[NOT_MOTIF] = TT_NOT;
}

//DNAMotif methods

DNAMotif::DNAMotif(){
    istype = DNA_MOTIF;
    // fan out multiplier - default 0
    fanCount = 0;
    //concentration multiplier - default 1
    concX = 1;
    return;
}
MOTIF_TYPE DNAMotif::getType(void){
    return istype;
}
void DNAMotif::setType(MOTIF_TYPE type){
    istype = type;
}

void DNAMotif::setFanOutMultiplier(int m){
    fanCount = m;
}

int DNAMotif::getFanOutMultiplier(void){
    return fanCount;
}

void DNAMotif::incrementMultiplier(void){
    fanCount++;
}

void DNAMotif::setConcMultiplier(int m){
    concX = m;
}

int DNAMotif::getConcMultiplier(void){
    return concX;
}

string DNAMotif::getID(void){
    return id[0];
    cerr << "in the dnamotif getID" << endl;
}

void DNAMotif::print(void){
    cout << "2inputAndStrand:: ";
    for(int i=0;i<4;i++)
        andStrand[i]->print();
    cout << "gateStrand:: ";
    for(int i=0;i<4;i++)
        gateStrand[i]->print();
}

void DNAMotif::printConcentration(void){
    cout << "Conc: " << DNAMotif::getConcMultiplier()
        << "X" << endl;
}

string DNAMotif::getDSDName(int num){
    string motif = getGatePrefix[istype] + "_" + id[0] + "_0";
    motif[motif.length()-1] += num;
    return motif;
}

void DNAMotif::printForDSD(string motif){

    for(int i=0;i<4;i++){
        motif = getDSDName(i);
        cout << "def " 
            << motif
            << "() = {"
            << gateStrand[i]->getDomain(0,6, DSD) 
            << "}["
            << andStrand[i]->getDomain(0,0,DSD) << " "
            << andStrand[i]->getDomain(0,1,DSD) << " "
            << andStrand[i]->getDomain(0,2,DSD) << " "
            << andStrand[i]->getDomain(0,3,DSD) << " "
            << andStrand[i]->getDomain(0,4,DSD)
            << "]{"
            << gateStrand[i]->getDomain(0,0, DSD) 
            << "}"
            << endl;
    }

    cerr << "in the dnamotif printfordsd" << endl;
}

void DNAMotif::printForDSDWithConc(void){

    string motif;
    for(int i=0;i<4;i++){
        motif = getDSDName(i);
        cout << "| " 
            << DNA::conc * DNAMotif::getConcMultiplier()
            << " * "
            << motif
            << endl;
    }
}

void printGate(const int *t){

    for(int i=0;i<4;i++){
        for(int j=0;j<3;j++){
            cout << t[3*i+j] << ",";
        }
        cout << endl;
    }
}

void DNAMotif::printForCRN(CIRCUIT_TYPE version){
    int conc;
    const int *TT_GATE;
    /*
    TT_GATE = getTruthTable[istype];
    */
    switch(istype){
        case NAND_MOTIF: TT_GATE = TT_NAND;
                         break;
        case AND_MOTIF: TT_GATE = TT_AND;
                        break;
        case OR_MOTIF: TT_GATE = TT_OR;
                       break;
        case NOR_MOTIF: TT_GATE = TT_NOR;
                        break;
        case NOT_MOTIF: TT_GATE = TT_NOT;
                        break;
        default: TT_GATE = TT_AND;
    }
    //printGate(TT_GATE);

    if (version == REUSABLE){
        cout << gate(id[0],0) << " + " << id[0] << "0 -> {bi_forward} " << gate(id[0],0) << DNA::sep << id[0] << "0" << " |" << endl;
        cout << gate(id[0],1) << " + " << id[0] << "1 -> {bi_forward} " << gate(id[0],1) << DNA::sep << id[0] << "1" << " |" << endl;

        if (DNA::fanOutEnabled){
            conc = DNA::conc;
            cout << fanOutANDReusable(id[0], DNAMotif::getFanOutMultiplier(), conc);
        }
        else{
            conc = DNA::conc * DNAMotif::getConcMultiplier();
        }

        //initialize bit 0 by default.
        if(istype == BIT_MOTIF){
            cout << "init " << id[0] << "0 " << conc << " |" << endl;
            return;
        }

        cout << "init " << gate(id[0],0) << DNA::sep << id[0] << "0 " << conc << " |" << endl;
        cout << "init " << gate(id[0],1) << DNA::sep << id[0] << "1 " << conc << " |" << endl;

        if (istype == NOT_MOTIF){
            for(int i=0;i<2;i++){
                cout << id[1] << TT_GATE[2*i] << " + " 
                    << gate(id[0], TT_GATE[2*i+1]) << DNA::sep << id[0] << TT_GATE[2*i+1] 
                    << " -> {bi_forward} " << id[1] << TT_GATE[2*i] << DNA::sep << gate(id[0], TT_GATE[2*i+1])
                    << " + " << id[0] << TT_GATE[2*i+1]
                    << " |" << endl;
                cout << gate(id[1], TT_GATE[2*i]) << " + " 
                    << id[1] << TT_GATE[2*i] << DNA::sep << gate(id[0],TT_GATE[2*i+1])
                    << " -> {bi_forward} " << gate(id[1], TT_GATE[2*i]) << DNA::sep << id[1] << TT_GATE[2*i]
                    << " + " << gate(id[0],TT_GATE[2*i+1])
                    << " |" << endl;
            }
        }
        else{
            for(int i=0;i<4;i++){
                cout << id[1] << TT_GATE[3*i] << " + " << id[2] << TT_GATE[3*i+1] << " + "
                    << gate(id[0], TT_GATE[3*i+2]) << DNA::sep << id[0] << TT_GATE[3*i+2] 
                    << " -> {tri_forward} " << id[0] << TT_GATE[3*i+2] << " + "
                    << intermediate(id[1],id[2],id[0],TT_GATE[3*i], TT_GATE[3*i+1], TT_GATE[3*i+2])
                    << " |" << endl;
                cout << gate(id[1],TT_GATE[3*i]) << " + " 
                    << intermediate(id[1],id[2],id[0],TT_GATE[3*i], TT_GATE[3*i+1], TT_GATE[3*i+2]) 
                    << " -> {bi_forward} " << gate(id[1],TT_GATE[3*i]) << DNA::sep << id[1] << TT_GATE[3*i] 
                    << " + " << id[2] << TT_GATE[3*i+1] << " + " << gate(id[0],TT_GATE[3*i+2]) << " |" << endl;
                cout << gate(id[2],TT_GATE[3*i+1]) << " + " 
                    << intermediate(id[1],id[2],id[0],TT_GATE[3*i], TT_GATE[3*i+1], TT_GATE[3*i+2]) 
                    << " -> {bi_forward} " << gate(id[2],TT_GATE[3*i+1]) << DNA::sep << id[2] << TT_GATE[3*i+1] 
                    << " + " << id[1] << TT_GATE[3*i] << " + " << gate(id[0],TT_GATE[3*i+2]) << " |" << endl;
            }}
    }
    else if (version == ONE_TIME){
        if (DNA::fanOutEnabled){
            cout << fanOutOneTime(id[0], DNAMotif::getFanOutMultiplier());
            conc = DNA::conc;
        }
        else conc = DNA::conc * DNAMotif::getConcMultiplier();

        if(istype == STRAND_MOTIF){
            cout << "init " << id[0] << "0 " << conc << " |" << endl;
            return;
        }
        else if (istype == BIT_MOTIF){
            return;
        }
        else if(istype == NOT_MOTIF){
            for(int i=0;i<2;i++){
                cout << id[1] << TT_GATE[2*i] 
                    << " -> {uni_forward} " 
                    << id[0] << TT_GATE[2*i+1] << " |" << endl;
            }
        }
        else{
            for(int i=0;i<4;i++){
                cout << id[1] << TT_GATE[3*i] << " + " 
                    << id[2] << TT_GATE[3*i+1] << " -> {bi_forward} " 
                    << id[0] << TT_GATE[3*i+2] << " |" << endl;
            }
        }
    }
}

void DNAMotif::constructFromInput(DNAMotif *I1, DNAMotif *I2, int f1, int f2){
    cout << "for gates " << I1->getID() << " and " << I2->getID() 
        << ", the nums are : " << f1 << "," << f2 << endl;
    //andStrand

    string ids[2];
    ids[0] = id[0] + "0";
    ids[1] = id[0] + "1";
    const int *TT_gate;
    TT_gate = getTruthTable[istype];
    printGate(TT_gate);
    I1->print();
    I2->print();

    for(int i=0;i<4;i++){
        //andStrand
        andStrand[i]->push(I1->getDomain(TT_gate[3*i],2));
        cerr << "done till here " << i << endl;
        andStrand[i]->push(I1->getDomain(TT_gate[3*i],3));
        andStrand[i]->push(ids[TT_gate[3*i+2]]);
        andStrand[i]->push(I2->getDomain(TT_gate[3*i+1],1));
        andStrand[i]->push(I2->getDomain(TT_gate[3*i+1],2));

        //gateStrand
        gateStrand[i]->push(I2->getComplementDomain(TT_gate[3*i+1],3));
        gateStrand[i]->push(I2->getComplementDomain(TT_gate[3*i+1],2));
        gateStrand[i]->push(I2->getComplementDomain(TT_gate[3*i+1],1));
        gateStrand[i]->push(andStrand[i]->getComplementDomain(0,2));
        gateStrand[i]->push(I1->getComplementDomain(TT_gate[3*i],3));
        gateStrand[i]->push(I1->getComplementDomain(TT_gate[3*i],2));
        gateStrand[i]->push(I1->getComplementDomain(TT_gate[3*i],1));
    }

    constructCRN(I1->getID(), I2->getID(), f1, f2);
    cerr << "in the dnamotif construct from input " << endl;
}

//The assumption is that only the 1,2,3 domains can be got
//from a 2-input gate.
string DNAMotif::getDomain(int bit, int idx, TOEHOLD_TYPE type){
    const int *tt;
    tt = getTruthTable[istype];
    if(bit == 1){
        for(int i=0;i<4;i++)
            if(tt[3*i+2] == 1)
                return andStrand[i]->getDomain(-1, idx, type);
    }

    for(int i=0;i<4;i++)
        if(tt[3*i+2] == 0)
            return andStrand[i]->getDomain(-1, idx, type);
    cerr << "in the dnamotif getDomain " << endl;
}

string DNAMotif::getComplementDomain(int bit, int idx){
    const int *tt;
    tt = getTruthTable[istype];
    if(bit == 1){
        for(int i=0;i<4;i++)
            if(tt[3*i+2] == 1)
                return andStrand[i]->getComplementDomain(-1, idx);
    }

    for(int i=0;i<4;i++)
        if(tt[3*i+2] == 0)
            return andStrand[i]->getComplementDomain(-1, idx);
    cerr << "in the dnamotif getComplementDomain " << endl;
}

vector<STRAND> DNAMotif::getStrands(int bit){
    vector<STRAND> ret;
    if(bit == -1){
        for(int i=0;i<4;i++){
            ret.push_back(*andStrand[i]);
            ret.push_back(*gateStrand[i]);
        }
    }
    else{
        ret.push_back(*andStrand[bit]);
        ret.push_back(*gateStrand[bit]);
    }
    return ret;
    cerr << "in the dnamotif getStrands " << endl;
}

set<pair<string, string> > DNAMotif::getUniqueDomains(void){
    cerr << "in the dnamotif getUniqueDomains " << endl;
    set<pair<string, string> > ret, current;
    set<pair<string, string> >::iterator it;
    for(int i=0;i<4;i++){
        current = andStrand[i]->getUniqueDomains();
        for(it = current.begin(); it != current.end(); it++)
            ret.insert(*it);
        current = gateStrand[i]->getUniqueDomains();
        for(it = current.begin(); it != current.end(); it++)
            ret.insert(*it);
    }
    return ret;
}

vector<string> DNAMotif::printNupackStructureAndSequence(
        map<int, DNAMotif*> &m, 
        map<int, int> &names,
        vector<vector<int> >&g){
    vector<string> ret;
    return ret;
}

vector<STRAND> DNAMotif::getOutputStrandList(int bit){
    cerr << "or having trouble here it seems" << endl;
    vector<STRAND> ret;
    const int *TT_gate;
    TT_gate = getTruthTable[getType()];
    for(int i=0;i<4;i++)
        if(TT_gate[3*i+2] == bit){
            ret.push_back(*andStrand[i]);
        }
    return ret;
}


void DNAMotif::constructCRN(const string &id1, const string &id2, const int &f1, const int &f2){
    if(f1){
        id[1] = getFanOutDomainName(id1, f1);
    }
    else
        id[1] = id1;

    if(f2){
        id[2] = getFanOutDomainName(id2, f2);
    }
    else
        id[2] = id2;
}

int DNAMotif::parseStrand(const string &name, const string &seq){
    cerr << "in dnamotif parseStrand" << endl;
    return 0;
}

//STRAND methods

STRAND::STRAND(){}
STRAND::STRAND(int num){
    DNAMotif::setType(STRAND_MOTIF);
    idNum = num;
    id[0] = "A";
    id[0][0] += num;
}

STRAND::STRAND(const string &s, int num){
    DNAMotif::setType(STRAND_MOTIF);
    idNum = num;
    id[0] = "A";
    id[0][0] += num;
    init(s);
}

void STRAND::init(const string &s){
    char *tokens = new char[s.length()+1];
    s.copy(tokens,s.length(),0);
    tokens[s.length()] = '\0';

    vector<string> doms;
    char *token = strtok(tokens, " \t");
    while(token != NULL){
        doms.push_back(token);
        token = strtok(NULL, " \t");
    }
    createFromVector(doms);
}

void STRAND::print(void){
    for(int i=0;i<name.size();i++){
        cout << name[i];
        if(complement[i]) cout << "'";
        cout << " ";
    }
    cout << endl;
}
void STRAND::printConcentration(void){
    cout << "Concentration: " << DNAMotif::getConcMultiplier()
        << "X" << endl;
}
void STRAND::printForDSD(string motif){
    cout << "def " 
        << motif
        << "() = < ";
    for(int i=0;i<name.size();i++){
        cout << name[i] << "^";
        if(complement[i]) cout << "*";
        cout << " ";
    }
    cout << ">";
    cout << endl;
}

string STRAND::getDomain(int bit, int idx, TOEHOLD_TYPE type){
    cerr << "in strand getdomain" << endl;
    if(name.size() <= idx){
        throw "REQUEST FOR INDEX NOT PRESENT IN STRAND";
    }

    string suffix="";
    if(type == DSD)
        suffix = "^";

    if(complement[idx] == 0)
        return name[idx]+suffix;
    else return name[idx]+suffix+"*";
}

string STRAND::getComplementDomain (int bit, int idx){
    if(name.size() <= idx){
        throw "REQUEST FOR INDEX NOT PRESENT IN STRAND";
    }

    if(complement[idx] == 1)
        return name[idx];
    else return name[idx]+"*";
}

void STRAND::push(string s){
    if(s[s.length()-1] == '*'){
        name.push_back(s.substr(0,s.length()-1));
        complement.push_back(1);
    }
    else {
        name.push_back(s);
        complement.push_back(0);
    }
    domainSize.push_back(TOEHOLD_DEFAULT_LENGTH);
    domainSeq.push_back("");
}

string STRAND::getNewDomain(string prefix){
    if(STRAND::nextNumber.find(prefix) == STRAND::nextNumber.end())
        STRAND::nextNumber[prefix] = 0;

    ostringstream oss1;
    oss1 << prefix << STRAND::nextNumber[prefix];
    STRAND::nextNumber[prefix]++;

    return oss1.str();
}

bool STRAND::compareDomains(STRAND to){
    vector<string> biDomain1, biDomain2;
    biDomain1 = getConcatenatedDomains(2,2);
    biDomain2 = to.getConcatenatedDomains(2,2);
    for(int i=0;i<biDomain1.size();i++)
        for(int j=0;j<biDomain2.size();j++)
            if(biDomain1[i] == biDomain2[j])
                return true;
    return false;
}

vector<STRAND> STRAND::getStrands(int bit){
    vector<STRAND> ret;
    ret.push_back(*this);
    return ret;
}

string STRAND::getID(void){
    return id[0];
}

void STRAND::createFromVector(vector<string> s){
    name.resize(s.size());
    complement.resize(s.size());
    domainSize.resize(s.size());
    domainSeq.resize(s.size());

    for(int i=0;i<s.size();i++){
        complement[i] = 0;
        name[i] = s[i];
        domainSize[i] = TOEHOLD_DEFAULT_LENGTH;
        domainSeq[i] = "";

        int last = s[i].length() - 1;
        if(s[i][last] == '*'){
            complement[i] = 1;
            name[i].erase(name[i].end()-1);
        }
    }
}

vector<string> STRAND::getConcatenatedDomains(int len, int which){
    vector<string> ret;
    int num = name.size();
    for(int i=0,j=num-1;i<num-1;i++,j--){
        string twoDomains;

        if(which == 0 || which == 2){
            twoDomains.clear();
            //i,i+1
            twoDomains += name[i];
            if(complement[i]) twoDomains += "*";
            twoDomains += name[i+1];
            if(complement[i+1]) twoDomains += "*";
            ret.push_back(twoDomains);
        }

        if(which == 1 || which == 2){
            twoDomains.clear();
            //revcomp(j,j-1)
            twoDomains += name[j];
            if(!complement[j]) twoDomains += "*";
            twoDomains += name[j-1];
            if(!complement[j-1]) twoDomains += "*";
            ret.push_back(twoDomains);
        }
    }
    return ret;
}

int STRAND::getNumberOfDomains(void){
    return name.size();
}

string reverseComplement(string seq){
    cerr << "in reverse complement for string:"
        << seq
        << endl;
    transform(seq.begin(), seq.end(),
        seq.begin(), ::toupper);
    replace(seq.begin(), seq.end(), 'A','X');
    replace(seq.begin(), seq.end(), 'T','A');
    replace(seq.begin(), seq.end(), 'X','T');
    replace(seq.begin(), seq.end(), 'C','X');
    replace(seq.begin(), seq.end(), 'G','C');
    replace(seq.begin(), seq.end(), 'X','G');
    reverse(seq.begin(), seq.end());
    cerr << "obtained:"
        << seq
        << endl;
    return seq;
}

set<pair<string, string> > STRAND::getUniqueDomains(void){
    cerr << "in the STRAND getUniqueDomains" << endl;
    cerr << "sizes: "
        << name.size() << ","
        << complement.size() << ","
        << domainSize.size() << ","
        << domainSeq.size() << ","
        << endl;
    int len = name.size();
    set<pair<string, string> > ret;
    pair<string, string> p;
    for (int i=0;i<len;i++){
        p.first = name[i];
        if(complement[i])
            p.second = reverseComplement(domainSeq[i]);
        else
            p.second = domainSeq[i];
        ret.insert(p);
    }
    cerr << "exiting the STRAND getUniqueDomains" << endl;
    return ret;
}

int STRAND::getDomainLength(int pos){
    return domainSize[pos];
}

string STRAND::getName(void){
}

int STRAND::parseDomains(const string &seq){
    int len = domainSize.size();
    int pos = 0;
    for(int i=0;i<len;i++){
        domainSeq[i] = seq.substr(pos,domainSize[i]);
        pos += domainSize[i];
    }
    if(pos != seq.size())
        throw "Error in parseDomains, sequence and domain lengths do not match";
    return len;
}

//END of STRAND member functions.
//START OF BIT CLASS

typedef class BIT : public DNAMotif{
    public:
        BIT(){}

        BIT(const string &s0, const string &s1, int num){
            DNAMotif::setType(BIT_MOTIF);
            id[0] = "A";
            id[0][0] += num;
            idNum = num;
            andStrand[0] = new STRAND(s0,num);
            andStrand[1] = new STRAND(s1,num);
        }

        void print(void){
            for(int i=0;i<2;i++)
                andStrand[i]->print();
        }

        void printForDSD(string motif = "BIT"){
            for(int i=0;i<2;i++){
                motif = getDSDName(i);
                andStrand[i]->printForDSD(motif);
            }
        }

        void printForDSDWithConc(void){
            string motif;
            for(int i=0;i<2;i++){
                motif = getDSDName(i);
                cout << "| " 
                    << DNA::conc * DNAMotif::getConcMultiplier()
                    << " * "
                    << motif
                    << endl;
            }
        }

        string getDomain(int bit, int idx, TOEHOLD_TYPE type = NORMAL){
            if(bit != 0 && bit != 1)
                throw "REQUEST FOR INPUT INVALID ";

            cerr << "in bit getdomain" << endl;
            return andStrand[bit]->getDomain(bit, idx, type);
        }

        string getComplementDomain (int bit, int idx){
            if(bit != 0 && bit != 1)
                throw "REQUEST FOR INPUT INVALID ";

            return andStrand[bit]->getComplementDomain(bit, idx);
        }

        vector<STRAND> getStrands(int bit){
            vector<STRAND> ret;
            if(bit == -1)
                for(int i=0;i<2;i++)
                    ret.push_back(*andStrand[i]);
            else
                ret.push_back(*andStrand[bit]);
            return ret;
        }

        set<pair<string, string> > getUniqueDomains(void){
            cerr << "in the BIT getUniqueDomains" << endl;
            set<pair<string, string> > ret, current;
            set<pair<string, string> >::iterator it;
            for(int i=0; i<2; i++){
                current = andStrand[i]->getUniqueDomains();
                for(it = current.begin(); it != current.end(); it++)
                    ret.insert(*it);
            }
            return ret;
        }

        vector<string> printNupackStructureAndSequence(
                map<int, DNAMotif*> &m, 
                map<int, int> &names,
                vector<vector<int> >&g){
            vector<string> ret;
            string seq;
            int numDomains;
            int len;
            for(int bit=0;bit<2;bit++){
                ostringstream oss1;
                numDomains = andStrand[bit]->getNumberOfDomains();
                len = 0;
                for(int i=0;i<numDomains;i++){
                    len += andStrand[bit]->getDomainLength(i);
                    seq += getDomain(bit,i) + " ";
                }

                oss1 << "structure "
                    << getDSDName(bit)
                    << " = U"
                    << len
                    << endl;

                oss1 << getDSDName(bit)
                    << ".seq = "
                    << seq
                    << endl;
                ret.push_back(oss1.str());
                seq.clear();
            }
            //cout << oss1.str() << endl;
            return ret;
        }

        vector<STRAND> getOutputStrandList(int bit){
            vector<STRAND> ret;
            cerr << "having trouble here it seems" << endl;
            ret.push_back(*andStrand[bit]);
            return ret;
        }

        int parseStrand(const string &name, const string &seq){
            //Assuming some standard format.
            //char 0 -> name of the strand.
            //char 1 gives out which bit 
            //this corresponds to.
            //char 2 is always _.
            //char 3 is always Z 
            int bit = name[1]-'0';
            switch(name[3]){
                case 'Z':
                    return andStrand[bit]->parseDomains(seq);
                default:
                    string err = name + "invalid strand name in BIT";
                    throw err.c_str();
                    return -1;
            }
        }
}BIT;

map<string, int> STRAND::nextNumber;

typedef class BI_INPUT : public DNAMotif {
    public:
        BI_INPUT(){}
        BI_INPUT(int num, MOTIF_TYPE type){
            DNAMotif::setType(type);
            id[0] = "A";
            id[0][0] += num;
            idNum = num;
            for(int bit2=0;bit2<4;bit2++){
                andStrand[bit2] = new STRAND(num);
                gateStrand[bit2] = new STRAND(num);
            }
        }


        vector<string> lInputGateNupackComplex(int bit2, vector<STRAND> &output, STRAND &gate){
            int count = output.size();
            string seq,subname;
            int len1, len2, numDomains;
            vector<string> ret;
            cerr << "in lInputGate function, count=" << count << endl;
            for(int i=0;i<count;i++){
                ostringstream oss1;
                //bind output[i] and the gate strand together.
                numDomains = output[i].getNumberOfDomains();
                if(numDomains != 5)
                    throw "Incorrect number of domains in left output strand";

                len1 = 0;
                for(int j=1;j<4;j++){
                    len1 += output[i].getDomainLength(j);
                }

                numDomains = gate.getNumberOfDomains();
                if(numDomains != 7)
                    throw "Incorrect number of domains in gate";

                len2 = 0;
                for(int j=0;j<4;j++){
                    len2 += gate.getDomainLength(j);
                }

                subname = "_XG";
                if(count > 1){
                    subname += "_a";
                    subname[subname.length()-1]+=i;
                }

                oss1 << "structure "
                    << getDSDName(bit2)
                    << subname << " = U"
                    << output[i].getDomainLength(0)
                    << " D"
                    << len1
                    << " (U"
                    << output[i].getDomainLength(4)
                    << " + U"
                    << len2
                    << ")"
                    << endl;

                seq.clear();
                for(int j=0;j<5;j++)
                    seq += output[i].getDomain(-1,j) + " ";
                for(int j=0;j<7;j++)
                    seq += gate.getDomain(-1,j) + " ";

                oss1 << getDSDName(bit2)
                    << subname << ".seq = "
                    << seq
                    << endl;
                ret.push_back(oss1.str());
            }
            //cout << oss1.str() << endl;
            cerr << "exiting lInputGate function" << endl;
            return ret;
        }

        vector<string> rInputGateNupackComplex(int bit2, vector<STRAND> &output, STRAND &gate){
            int count = output.size();
            string seq,subname;
            int len1, len2, numDomains;
            vector<string> ret;
            cerr << "in rInputGate function, count=" << count << endl;
            for(int i=0;i<count;i++){
                ostringstream oss1;
                //bind output[i] and the gate strand together.
                numDomains = output[i].getNumberOfDomains();
                if(numDomains != 5)
                    throw "Incorrect number of domains in right output strand";

                len1 = 0;
                for(int j=1;j<4;j++){
                    len1 += output[i].getDomainLength(j);
                }

                numDomains = gate.getNumberOfDomains();
                if(numDomains != 7)
                    throw "Incorrect number of domains in gate";

                len2 = 0;
                for(int j=3;j<7;j++){
                    len2 += gate.getDomainLength(j);
                }

                subname = "_YG";
                if(count > 1){
                    subname += "_a";
                    subname[subname.length()-1]+=i;
                }

                oss1 << "structure "
                    << getDSDName(bit2)
                    << subname << " = U"
                    << output[i].getDomainLength(0)
                    << " D"
                    << len1
                    << " (U"
                    << output[i].getDomainLength(4)
                    << " +) U"
                    << len2
                    << endl;

                seq.clear();
                for(int j=0;j<5;j++)
                    seq += output[i].getDomain(-1,j) + " ";
                for(int j=0;j<7;j++)
                    seq += gate.getDomain(-1,j) + " ";

                oss1 << getDSDName(bit2)
                    << subname << ".seq = "
                    << seq
                    << endl;
                ret.push_back(oss1.str());
            }
            //cout << oss1.str() << endl;
            cerr << "exiting rInputGate function" << endl;
            return ret;
        }

        vector<string> bothInputGateNupackComplex(int bit2, vector<STRAND> &loutput, vector<STRAND> routput, STRAND &gate){
            int lcount = loutput.size();
            int rcount = routput.size();
            string seq,subname;
            int len1, len2, lnumDomains, rnumDomains;
            vector<string> ret;
            cerr << "in bothInputGate function, lcount=" 
                << lcount 
                << ", rcount = " << rcount
                << endl;
            for(int i=0;i<lcount;i++){
            for(int k=0;k<rcount;k++){
                ostringstream oss1;
                //bind output[i] and the gate strand together.
                lnumDomains = loutput[i].getNumberOfDomains();
                if(lnumDomains != 5)
                    throw "Incorrect number of domains in both output strand (left)";

                len1 = 0;
                for(int j=1;j<4;j++){
                    len1 += loutput[i].getDomainLength(j);
                }

                rnumDomains = routput[k].getNumberOfDomains();
                if(rnumDomains != 5)
                    throw "Incorrect number of domains in both output strand (right)";

                if(gate.getNumberOfDomains() != 7)
                    throw "Incorrect number of domains in both output strand (gate)";
                len2 = 0;
                for(int j=1;j<4;j++){
                    len2 += routput[k].getDomainLength(j);
                }

                subname = "_XYG";
                if(lcount > 1 || rcount > 1){
                    subname += "_aa";
                    subname[subname.length()-2]+=i;
                    subname[subname.length()-1]+=k;
                }

                oss1 << "structure "
                    << getDSDName(bit2)
                    << subname << " = U"
                    << loutput[i].getDomainLength(0)
                    << " D"
                    << len1
                    << " (U"
                    << loutput[i].getDomainLength(4)
                    << " + U"
                    << routput[k].getDomainLength(0)
                    << " D"
                    << len2
                    << " (U"
                    << routput[k].getDomainLength(4)
                    << "+) U"
                    << gate.getDomainLength(3)
                    << ")"
                    << endl;

                seq.clear();
                for(int j=0;j<5;j++)
                    seq += loutput[i].getDomain(-1,j) + " ";
                for(int j=0;j<5;j++)
                    seq += routput[k].getDomain(-1,j) + " ";
                for(int j=0;j<7;j++)
                    seq += gate.getDomain(-1,j) + " ";

                oss1 << getDSDName(bit2)
                    << subname << ".seq = "
                    << seq
                    << endl;
                ret.push_back(oss1.str());
            }
            }
            //cout << oss1.str() << endl;
            cerr << "exiting bothInputGate function" << endl;
            return ret;
        }

        vector<string> outputGateNupackComplex(int bit2, vector<STRAND> &output, STRAND &gate){
            int count = output.size();
            string seq,subname;
            int len1, numDomains;
            vector<string> ret;
            cerr << "in outputGate function, count=" << count << endl;
            for(int i=0;i<count;i++){
                ostringstream oss1;
                //bind output[i] and the gate strand together.
                numDomains = output[i].getNumberOfDomains();
                if(numDomains != 5)
                    throw "Incorrect number of domains in left output strand";

                len1 = 0;
                for(int j=0;j<5;j++){
                    len1 += output[i].getDomainLength(j);
                }

                numDomains = gate.getNumberOfDomains();
                if(numDomains != 7)
                    throw "Incorrect number of domains in gate";

                subname = "_ZG";
                if(count > 1){
                    subname += "_a";
                    subname[subname.length()-1]+=i;
                }

                oss1 << "structure "
                    << getDSDName(bit2)
                    << subname << " = D"
                    << len1
                    << " (+U"
                    << gate.getDomainLength(0)
                    << ") U"
                    << gate.getDomainLength(6)
                    << endl;

                seq.clear();
                for(int j=0;j<5;j++)
                    seq += output[i].getDomain(-1,j) + " ";
                for(int j=0;j<7;j++)
                    seq += gate.getDomain(-1,j) + " ";

                oss1 << getDSDName(bit2)
                    << subname << ".seq = "
                    << seq
                    << endl;
                ret.push_back(oss1.str());
            }
            //cout << oss1.str() << endl;
            cerr << "exiting outputGate function" << endl;
            return ret;
        }

        vector<string> printNupackStructureAndSequence(
                map<int, DNAMotif*> &m, 
                map<int, int> &names,
                vector<vector<int> >&g){
            vector<string> ret;
            //ostringstream oss1;
            //string seq;
            int numDomains;
            int len;

            //try and get the two other left and right input strands that feed 
            //into this bi_input gate.
            int nodes = g.size();
            int lInput=-1, rInput=-1;
            for(int j=0;j<nodes;j++){
                if(g[j][names[idNum]] == 1)
                    lInput = j;
                if(g[j][names[idNum]] == 2)
                    rInput = j;
            }
            cerr << "lInput=" << lInput << ", rInput=" << rInput << endl;
            if(lInput==-1 || rInput==-1){
                for(int j=0;j<nodes;j++){
                    for(int k=0;k<nodes;k++){
                        cerr << g[j][k] << " ";
                    }
                    cerr << endl;
                }
                throw "Incorrect input given, the gate should have two inputs";
            }

            vector<STRAND> lOutputStrands, rOutputStrands;
            vector<string> collect;
            ret.clear();

            for(int bit2=0;bit2<4;bit2++){

                //get the output strands from the left input.
                int lbit = bit2/2;
                int rbit = bit2%2;
                cerr << "functioning right here." << endl;
                lOutputStrands = m[lInput]->getOutputStrandList(lbit);
                rOutputStrands = m[rInput]->getOutputStrandList(rbit);

                collect = lInputGateNupackComplex(bit2, lOutputStrands, *gateStrand[bit2]);
                for(int i=0;i<collect.size();i++) {
                    ret.push_back(collect[i]); 
                    cout << "IN:" << collect[i] << endl;
                }
                collect.clear();
                collect = rInputGateNupackComplex(bit2, rOutputStrands, *gateStrand[bit2]);
                for(int i=0;i<collect.size();i++) ret.push_back(collect[i]); collect.clear();
                collect = bothInputGateNupackComplex(bit2, lOutputStrands, rOutputStrands, *gateStrand[bit2]);
                for(int i=0;i<collect.size();i++) ret.push_back(collect[i]); collect.clear();
                vector<STRAND> outputStrs;
                outputStrs.push_back(*andStrand[bit2]);
                collect = outputGateNupackComplex(bit2, outputStrs, *gateStrand[bit2]);
                for(int i=0;i<collect.size();i++) ret.push_back(collect[i]); collect.clear();
                cerr << "passed the output strand " << bit2 << endl;

            }
            //cout << oss1.str() << endl;
            return ret;
        }

        int parseStrand(const string &name, const string &seq){
            //Assuming some standard format.
            //char 0 -> name of the strand.
            //char 1,2 give out which bit 
            //this corresponds to.
            //char 3 is always _.
            //char 4 is either Z or G based on 
            //and or gate strand.
            int bit2 = (name[1]-'0')*2+(name[2]-'0');
            switch(name[4]){
                case 'Z':
                    return andStrand[bit2]->parseDomains(seq);
                case 'G':
                    return gateStrand[bit2]->parseDomains(seq);
                default:
                    string err;
                    err = name + "invalid strand name in BI_INPUT";
                    throw err.c_str();
                    return -1;
            }
        }

}BI_INPUT;

typedef class NOT : public DNAMotif {
    public:
        NOT(){}
        NOT(int num){
            DNAMotif::setType(NOT_MOTIF);
            id[0] = "A";
            id[0][0] += num;
            idNum = num;
            for(int i=0;i<2;i++){
                andStrand[i] = new STRAND(num);
                gateStrand[i] = new STRAND(num);
            }
        }

        void constructFromInput(DNAMotif *I1, DNAMotif *I2, int f1=0, int f2=0){
            cout << "for gates " << I1->getID() <<
                ", the nums are : " << f1 << endl;
            string ids[2];
            ids[0] = id[0] + "0";
            ids[1] = id[0] + "1";
            //andStrand
            //note that when getDomain or getComplementDomain is called on I1,
            //then the "bit" should be passed to it. Whenever these two methods
            //are called instead on a strand directly, the bit argument does not
            //matter, and to denote this, a -1 is being passed.
            //
            //Also note, the 0 input creates the not motif 1, and vice versa.
            for(int bit=0;bit<2;bit++){
                andStrand[!bit]->push(I1->getDomain(bit,2));
                andStrand[!bit]->push(I1->getDomain(bit,3));
                andStrand[!bit]->push(ids[!bit]);

                if(bit == 0){
                    andStrand[!bit]->push(STRAND::getNewDomain(DOMAINPREFIX));
                    andStrand[!bit]->push(STRAND::getNewDomain(DOMAINPREFIX));
                }
                else{
                    andStrand[!bit]->push(andStrand[1]->getDomain(-1,3));
                    andStrand[!bit]->push(andStrand[1]->getDomain(-1,4));
                }

                gateStrand[!bit]->push(andStrand[!bit]->getComplementDomain(-1,2));
                gateStrand[!bit]->push(andStrand[!bit]->getComplementDomain(-1,1));
                gateStrand[!bit]->push(andStrand[!bit]->getComplementDomain(-1,0));
                gateStrand[!bit]->push(I1->getComplementDomain(bit,1));
            }

            constructCRN(I1->getID(), "", f1, f2);
        }

        string getDomain(int bit, int idx, TOEHOLD_TYPE type = NORMAL){
            return andStrand[bit]->getDomain(-1, idx, type);
        }

        string getComplementDomain (int bit, int idx){
            return andStrand[bit]->getComplementDomain(-1, idx);
        }

        void print(void){
            cout << "andStrand:: " << endl;
            for(int i=0;i<2;i++){
                andStrand[i]->print();
                gateStrand[i]->print();
            }
        }
        void printForDSD(string motif){

            for(int i=0;i<2;i++){
                motif = DNAMotif::getDSDName(i);
                cout << "def " 
                    << motif
                    << "() = {"
                    << gateStrand[i]->getDomain(-1,3, DSD) 
                    << "}["
                    << andStrand[i]->getDomain(-1,0,DSD) << " "
                    << andStrand[i]->getDomain(-1,1,DSD) << " "
                    << andStrand[i]->getDomain(-1,2,DSD)
                    << "]{"
                    << andStrand[i]->getDomain(-1,3, DSD) << " "
                    << andStrand[i]->getDomain(-1,4, DSD) 
                    << "}"
                    << endl;
            }
        }

        void printForDSDWithConc(void){
            string motif;
            for(int i=0;i<2;i++){
                motif = DNAMotif::getDSDName(i);
                cout << "| " 
                    << DNA::conc * DNAMotif::getConcMultiplier()
                    << " * "
                    << motif
                    << endl;
            }
        }

        vector<STRAND> getStrands(int bit){
            vector<STRAND> ret;
            if(bit == -1)
                for(int i=0;i<2;i++){
                    ret.push_back(*andStrand[i]);
                    ret.push_back(*gateStrand[i]);
                }
            else{
                ret.push_back(*andStrand[bit]);
                ret.push_back(*gateStrand[bit]);
            }
            return ret;
        }

        set<pair<string, string> > getUniqueDomains(void){
            cerr << "in the NOT getUniqueDomains" << endl;
            set<pair<string, string> > ret, current;
            set<pair<string, string> >::iterator it;
            for(int i=0; i<2; i++){
                current = andStrand[i]->getUniqueDomains();
                for(it = current.begin(); it != current.end(); it++)
                    ret.insert(*it);
                current = gateStrand[i]->getUniqueDomains();
                for(it = current.begin(); it != current.end(); it++)
                    ret.insert(*it);
            }
            return ret;
        }

        vector<string> printNupackStructureAndSequence(
                map<int, DNAMotif*> &m, 
                map<int, int> &names,
                vector<vector<int> >&g){
            vector<string> ret;
            string seq;
            int numDomains;
            int len;
            for(int bit=0;bit<2;bit++){
                ostringstream oss1;
                numDomains = andStrand[bit]->getNumberOfDomains();
                //Assuming numDomains is always 5.
                if(numDomains != 5)
                    throw "number of domains in NOT gate output strand is not equal to 5";
                len = 0;
                for(int i=0;i<3;i++){
                    len += andStrand[bit]->getDomainLength(i);
                    seq += andStrand[bit]->getDomain(bit,i) + " ";
                }
                oss1 << "structure "
                    << getDSDName(bit)
                    << " = D"
                    << len;

                len = 0;
                for (int i=3; i<5; i++){
                    len += andStrand[bit]->getDomainLength(i);
                    seq += andStrand[bit]->getDomain(bit,i) + " ";
                }

                oss1 << " (U"
                    << len
                    << " +) U"
                    << gateStrand[bit]->getDomainLength(3)
                    << endl;

                for(int i=0;i<4;i++){
                    seq += gateStrand[bit]->getDomain(bit,i) + " ";
                }

                oss1 << getDSDName(bit)
                    << ".seq = "
                    << seq
                    << endl;
                ret.push_back(oss1.str());
                seq.clear();
            }
            //cout << oss1.str() << endl;
            return ret;
        }

        vector<STRAND> getOutputStrandList(int bit){
            vector<STRAND> ret;
            ret.push_back(*andStrand[bit]);
            return ret;
        }

        int parseStrand(const string &name, const string &seq){
            //Assuming some standard format.
            //char 0 -> name of the strand.
            //char 1 gives out which bit 
            //this corresponds to.
            //char 2 is always _.
            //char 3,4 together are either Z or NG based on 
            //and or gate strand.
            int bit = name[1]-'0';
            switch(name[3]){
                case 'Z':
                    return andStrand[bit]->parseDomains(seq);
                case 'N':
                    return gateStrand[bit]->parseDomains(seq);
                default:
                    throw "invalid strand name in NOT";
                    return -1;
            }
        }

}NOT;

void deleteInput(int &n, map<int, DNAMotif *> &m){
    for(int i=0;i<n;i++)
        delete m[i];
}

//Returns 0 if an input starts with a # symbol.
int getUncommentedInput(ifstream &file, string &line){
    if(std::getline(file, line)){
        if(line[0] == '#')
            return 0;
    }
    else return -1;
    return 1;
}

int createInputFile(int &nodes, map<int, int> &names, vector<vector <int> > &g, map<int, DNAMotif*> &m, char* filename){
    ifstream inputfile;
    inputfile.open(filename);
    if(!inputfile) {
        cerr << "Error in file:" << filename << endl;
        return 0;
    }
    string nodesStr;
    while(!getUncommentedInput(inputfile, nodesStr));
    istringstream issN(nodesStr);
    issN >> nodes;
    string gate, line;
    int num, e1, e2, e3;

    //nand gates, and gates, input strands, output strands
    g.resize(nodes);
    for(int i=0;i<nodes;i++){
        g[i].resize(nodes);
        for(int j=0;j<nodes;j++){
            g[i][j] = 0;
        }
    }
    cout << "nodes = " << nodes << endl;

    int ct = 0;
    while(getUncommentedInput(inputfile, gate) != -1){
        cout << "gate = " 
            << gate
            << endl;

        if(gate == "nand"){
            while(!getUncommentedInput(inputfile, line));
            istringstream iss(line);
            while(iss >> num){
                names[num] = ct;
                m[ct] = new BI_INPUT(num, NAND_MOTIF);
                ct++;
            }
        }
        else if(gate == "and"){
            while(!getUncommentedInput(inputfile, line));
            istringstream iss(line);
            while(iss >> num){
                names[num] = ct;
                m[ct] = new BI_INPUT(num, AND_MOTIF);
                ct++;
            }
        }
        else if(gate == "or"){
            while(!getUncommentedInput(inputfile, line));
            istringstream iss(line);
            while(iss >> num){
                names[num] = ct;
                m[ct] = new BI_INPUT(num, OR_MOTIF);
                ct++;
            }
        }
        else if(gate == "not"){
            while(!getUncommentedInput(inputfile, line));
            istringstream iss(line);
            while(iss >> num){
                names[num] = ct;
                m[ct] = new NOT(num);
                ct++;
            }
        }
        else if(gate == "nor"){
            while(!getUncommentedInput(inputfile, line));
            istringstream iss(line);
            while(iss >> num){
                names[num] = ct;
                m[ct] = new BI_INPUT(num, NOR_MOTIF);
                ct++;
            }
        }
        else if(gate == "input"){
            while(!getUncommentedInput(inputfile, line));
            istringstream iss(line);
            string bit[2];
            while(iss >> e1){
                for(int i=0;i<2;i++){
                    while(!getUncommentedInput(inputfile, line));
                    bit[i] = line;
                }

                names[e1] = ct;
                m[ct] = new BIT(bit[0],bit[1],e1);
                ct++;
            }
        }
        else if(gate == "output"){
            while(!getUncommentedInput(inputfile, line));
            istringstream iss(line);
            while(iss >> num){
                names[num] = ct;
                m[ct] = new STRAND(num);
                ct++;
            }
        }
        else if(gate == "edges"){
            while(getUncommentedInput(inputfile, line) == 1){
                istringstream edgeIss(line);
                edgeIss >> e1 >> e2 >> e3;
                //e3 = 1 refers to left input, 2 refers to right input.
                //by default, it is left input, just for convention.
                //0 refers to no edge
                g[names[e1]][names[e2]] = e3;
            }
        }
    }
    
    //setting fanout multiplier for each motif.
    //if a motif does not have an out edge, its 
    //fan out by default is 1, and is left unchanged.
    int count = 0;
    for(int i=0;i<nodes;i++){
        count = 0;
        for(int j=0;j<nodes;j++)
            if(g[i][j]) count++;
        if(count)
            m[i]->setFanOutMultiplier(count);
    }
    inputfile.close();
    return 1;
}

//g is the edges matrix. By default, it contains only 0s, 1s and 2s.
//for the purposes of the sort, we increment each non-zero value by 2,
//and at the end, we decrement it back.
void topologicalSort (vector<vector <int > > &g, vector<int> &ret){

    int nodes = g.size();
    bool visited[nodes];
    int i,j;
    for(i=0;i<nodes;i++)
        visited[i] = false;

    int done = nodes;
    while(done>0){
        //find node with no incoming edges, add it to queue.
        for(i=0;i<nodes;i++){
            if(visited[i]) continue;
            
            for(j=0;j<nodes;j++){
                if(g[j][i] == 1 || g[j][i] == 2) break;
            }
            if(j == nodes) break;
        }
        ret.push_back(i);
        visited[i] = true;
        done--;

        for(j=0;j<nodes;j++)
            if(g[i][j] != 0)
                g[i][j] += 2;

    }

    for(i=0;i<nodes;i++){
        for(j=0;j<nodes;j++)
            if(g[i][j] != 0)
                g[i][j] -= 2;
    }
    return;
}

void readNupackMapFile(string &npMap){
    ifstream file;
    file.open(npMap.c_str());
    string line;
    string s1,s2;
    while(std::getline(file, line)){
        istringstream iss(line);
        iss >> s1;
        iss >> s2;
        if(s1.size() == 0 || s2.size() == 0){
            cerr << "nupack map encountered wrong input line" << endl;
            continue;
        }
        nupackMap[s1] = s2;
        line.clear();
        s1.clear(); s2.clear();
    }
    file.close();
}

int readNupackOutputFile(string &npOfile){
    ifstream file;
    file.open(npOfile.c_str());
    string line;
    string s1,s2;
    int input=-1;
    while(std::getline(file, line)){
        istringstream iss(line);
        iss >> s1;
        iss >> s2;
        if(s1 == "Normalized"){
            cerr << "encountered next input" << endl;
            input++;
            continue;
        }
        else if (s1.size() == 0 or s2.size() == 0){
            cerr << "nupack output encountered wrong input line" << endl;
            continue;
        }
        inputStrands[input][s1] = s2;
        line.clear();
        s1.clear(); s2.clear();
    }
    file.close();
    map<int, map<string, string> >::iterator it;
    map<string, string>::iterator nt;
    for(it = inputStrands.begin(); it!= inputStrands.end();
    it++){
        cout << "Input: " << it->first << endl;
        for(nt = (it->second).begin(); nt!=(it->second).end();
        nt++)
            cout << nupackMap[nt->first] << "," << nt->second << endl;
    }
    return input+1;
}


void printForNupackDesign(int &n, 
        map<int, DNAMotif*> &m, 
        map<int, int> &names,
        vector< vector< int> > &g,
        string npDfile
        ){

    ofstream designFile;
    designFile.open(npDfile.c_str());
    if(!designFile)
        throw "Could not open design file";

    designFile
        << "material = dna" << endl    
        << "temperature = 25" << endl
        << "trials = 10" << endl
        << "magnesium = 0.0125" << endl    
        << "sodium = 0.5" << endl   
        << "dangles = all" << endl   
        << endl
        ;

    set<pair<string, string> > allDomains, currDomains;
    set<pair<string, string> >::iterator it;

    //Get Domains
    for(int i=0;i<n;i++){
        currDomains = m[i]->getUniqueDomains();
        for(it = currDomains.begin(); it!=currDomains.end(); it++)
            allDomains.insert(*it);
    }
    for(it = allDomains.begin(); it!=allDomains.end(); it++){
        designFile
            << "domain " << it->first
            << " = ";
        //check energy thresholds condition only for non-excluded domains.
        if((domainEnergy[it->second] < LOW_DELG_THRESHOLD or
           domainEnergy[it->second] > HIGH_DELG_THRESHOLD) and
           (excludedDomain.find(it->first) == excludedDomain.end()))
            designFile << "N" << (it->second).size()
                << endl;
        else
            designFile << it->second
                << endl;
    }
    designFile << endl;

    //Get structure and seq for each dna motif
    vector<string> collect;
    for(int i=0;i<n;i++){
        collect = m[i]->printNupackStructureAndSequence(m,names,g);

        for(int j=0;j<collect.size();j++)
            designFile << collect[j];
        //designFile.close(); return;
    }

    designFile
    << endl
    << "tube AG = BIT_A_0 BIT_A_1 BIT_G_0 BIT_G_1" << endl
    << "AG.maxsize = 3" << endl
    << "tube DG = OR_D_0_ZG OR_D_1_ZG OR_D_2_ZG OR_D_3_ZG BIT_G_0 BIT_G_1" << endl
    << "DG.maxsize = 3" << endl
    << "tube AD = BIT_A_0 BIT_A_1 OR_D_0_ZG OR_D_1_ZG OR_D_2_ZG OR_D_3_ZG" << endl
    << "AD.maxsize = 3" << endl
    << "tube BC = AND_B_0_ZG AND_B_1_ZG AND_B_2_ZG AND_B_3_ZG AND_C_0_ZG AND_C_1_ZG AND_C_2_ZG AND_C_3_ZG" << endl
    << "BC.maxsize = 3" << endl
    << endl
    << "domain stem = CGGTG" << endl
    << "domain loop = GAGAAGAGAAAG" << endl
    << "structure MB = D5 U12" << endl
    << "MB.seq = stem loop stem*" << endl
    << "prevent = GGGGG" << endl
    << endl;

    designFile.close();
    cerr << "and that is a wrap for printForNupackDesign" << endl;
}

bool isBiMotif(MOTIF_TYPE type){
    switch(type){
        case NAND_MOTIF:
        case AND_MOTIF:
        case OR_MOTIF:
        case NOR_MOTIF:
            return true;
        default:
            return false;
    }
}

string int2string(int n){
   //ostringstream oss;
   //oss << n;
   //return oss.str();
   //or 
   return static_cast<ostringstream*>(&(ostringstream() << n))->str();
}

string int2binaryString(int n){
   //Assumption, only does numbers from 0-3.
   ostringstream oss;
   oss << n/2;
   oss << n%2;
   return oss.str();
}

void createCheckMapForPythonScript(
        int &n, 
        map<int, DNAMotif*> &m, 
        map<int, int> &names,
        vector< vector< int> > &g
        ){

        vector<STRAND> v1, v2;
        string suffix[2];
        suffix[0] = "Z";
        suffix[1] = "G";
        string str_lbits, str_rbits;
        vector<int>lbits;
        vector<int>rbits;
        
        for(int i=0;i<n;i++)
        for(int j=0;j<n;j++){
            //if direct edge in or out of gate
            if(g[i][j] == 1){
                
                //cout << "i=" << m[i]->getID() << ",j=" << m[j]->getID() << endl;
                //left input
                //need the following functions
                //getStrandsThatAcceptBit(bit)
                //getStrandsThatDontAcceptBit(bit)
                //getStrands()
                if(m[i]->getType() == BIT_MOTIF and 
                   m[j]->getType() == NOT_MOTIF){
                   lbits.clear(); rbits.clear();
                   lbits.push_back(0);
                   rbits.push_back(0);
                   lbits.push_back(1);
                   rbits.push_back(1);

                   for(int k=0;k<lbits.size();k++){
                      //reset the string bits.
                      str_lbits = int2string(lbits[k]);
                      str_rbits = int2string(rbits[k]);

                      v1 = m[i]->getStrands(lbits[k]);
                      v2 = m[j]->getStrands(rbits[k]);
                      for(int a=0;a<2;a++)
                         cout
                            << m[i]->getID() << str_lbits + "_" + suffix[0]+ "," 
                            << m[j]->getID() << str_rbits + "_" + suffix[a] 
                            << endl;
                   }
                }

                if(m[i]->getType() == BIT_MOTIF and 
                   isBiMotif(m[j]->getType())){
                   lbits.clear(); rbits.clear();
                   lbits.push_back(0);
                   rbits.push_back(2);
                   lbits.push_back(0);
                   rbits.push_back(3);
                   lbits.push_back(1);
                   rbits.push_back(0);
                   lbits.push_back(1);
                   rbits.push_back(1);

                   for(int k=0;k<lbits.size();k++){
                      //reset the string bits.
                      str_lbits = int2string(lbits[k]);
                      str_rbits = int2binaryString(rbits[k]);

                      v1 = m[i]->getStrands(lbits[k]);
                      v2 = m[j]->getStrands(rbits[k]);
                      for(int a=0;a<2;a++)
                         cout
                            << m[i]->getID() << str_lbits + "_" + suffix[0]+ "," 
                            << m[j]->getID() << str_rbits + "_" + suffix[a] 
                            << endl;
                   }
                }

                if(m[i]->getType() == NOT_MOTIF and m[j]->getType() == NOT_MOTIF){
                   lbits.clear(); rbits.clear();
                   lbits.push_back(0);
                   rbits.push_back(0);
                   lbits.push_back(1);
                   rbits.push_back(1);

                   for(int k=0;k<lbits.size();k++){
                      //reset the string bits.
                      str_lbits = int2string(lbits[k]);
                      str_rbits = int2string(rbits[k]);

                      v1 = m[i]->getStrands(lbits[k]);
                      v2 = m[j]->getStrands(rbits[k]);
                      for(int a=0;a<2;a++)
                        for(int b=0;b<2;b++)
                         cout
                            << m[i]->getID() << str_lbits + "_" + suffix[a]+ "," 
                            << m[j]->getID() << str_rbits + "_" + suffix[b] 
                            << endl;
                   }
                }

                if(m[i]->getType() == NOT_MOTIF and 
                   isBiMotif(m[j]->getType()) ){
                   lbits.clear(); rbits.clear();
                   lbits.push_back(0);
                   rbits.push_back(2);
                   lbits.push_back(0);
                   rbits.push_back(3);
                   lbits.push_back(1);
                   rbits.push_back(0);
                   lbits.push_back(1);
                   rbits.push_back(1);

                   for(int k=0;k<lbits.size();k++){
                      //reset the string bits.
                      str_lbits = int2string(lbits[k]);
                      str_rbits = int2binaryString(rbits[k]);

                      v1 = m[i]->getStrands(lbits[k]);
                      v2 = m[j]->getStrands(rbits[k]);
                      for(int a=0;a<2;a++)
                      for(int b=0;b<2;b++)
                         cout
                            << m[i]->getID() << str_lbits + "_" + suffix[a]+ "," 
                            << m[j]->getID() << str_rbits + "_" + suffix[b] 
                            << endl;
                   }
                }


                if(isBiMotif(m[i]->getType()) and 
                   m[j]->getType() == NOT_MOTIF){

                   const int *TT_gate;
                   TT_gate = getTruthTable[m[i]->getType()];
                   
                   lbits.clear(); rbits.clear();
                   
                   for(int bit2=0;bit2<4;bit2++){
                      lbits.push_back(bit2);
                      rbits.push_back(TT_gate[3*bit2+2]);
                   }

                   for(int k=0;k<lbits.size();k++){
                      //reset the string bits.
                      str_lbits = int2binaryString(lbits[k]);
                      str_rbits = int2string(rbits[k]);

                      v1 = m[i]->getStrands(lbits[k]);
                      v2 = m[j]->getStrands(rbits[k]);
                      for(int a=0;a<2;a++)
                      for(int b=0;b<2;b++)
                         cout
                            << m[i]->getID() << str_lbits + "_" + suffix[a]+ "," 
                            << m[j]->getID() << str_rbits + "_" + suffix[b] 
                            << endl;
                   }
                }

                if(isBiMotif(m[i]->getType()) and 
                   isBiMotif(m[j]->getType()) ){
                   const int *TT_gate;
                   TT_gate = getTruthTable[m[i]->getType()];
                   
                   lbits.clear(); rbits.clear();
                   
                   for(int lbit2=0;lbit2<4;lbit2++)
                   for(int rbit2=0;rbit2<4;rbit2++)
                      if(TT_gate[3*lbit2+2] !=
                         rbit2/2) {
                         lbits.push_back(lbit2);
                         rbits.push_back(rbit2);
                      }

                   for(int k=0;k<lbits.size();k++){
                      //reset the string bits.
                      str_lbits = int2binaryString(lbits[k]);
                      str_rbits = int2binaryString(rbits[k]);

                      v1 = m[i]->getStrands(lbits[k]);
                      v2 = m[j]->getStrands(rbits[k]);
                      for(int a=0;a<2;a++)
                      for(int b=0;b<2;b++)
                         cout
                            << m[i]->getID() << str_lbits + "_" + suffix[a]+ "," 
                            << m[j]->getID() << str_rbits + "_" + suffix[b] 
                            << endl;
                   }
                }
            }//end if g[i][j] is 1.
            else if (g[i][j] == 2){
               //cout << "right, i=" << m[i]->getID() << ",j=" << m[j]->getID() << endl;

                //right input
                //need the following functions
                //getStrandsThatAcceptBit(bit)
                //getStrandsThatDontAcceptBit(bit)
                //getStrands()

                if(m[i]->getType() == BIT_MOTIF and 
                   isBiMotif(m[j]->getType()) ){
                   lbits.clear(); rbits.clear();
                   lbits.push_back(0);
                   rbits.push_back(1);
                   lbits.push_back(0);
                   rbits.push_back(3);
                   lbits.push_back(1);
                   rbits.push_back(0);
                   lbits.push_back(1);
                   rbits.push_back(2);

                   for(int k=0;k<lbits.size();k++){
                      //reset the string bits.
                      str_lbits = int2string(lbits[k]);
                      str_rbits = int2binaryString(rbits[k]);

                      v1 = m[i]->getStrands(lbits[k]);
                      v2 = m[j]->getStrands(rbits[k]);
                      for(int a=0;a<2;a++)
                         cout
                            << m[i]->getID() << str_lbits + "_" + suffix[0]+ "," 
                            << m[j]->getID() << str_rbits + "_" + suffix[a] 
                            << endl;
                   }
                }

                if(m[i]->getType() == NOT_MOTIF and 
                   isBiMotif(m[j]->getType()) ){
                   lbits.clear(); rbits.clear();
                   lbits.push_back(0);
                   rbits.push_back(1);
                   lbits.push_back(0);
                   rbits.push_back(3);
                   lbits.push_back(1);
                   rbits.push_back(0);
                   lbits.push_back(1);
                   rbits.push_back(2);

                   for(int k=0;k<lbits.size();k++){
                      //reset the string bits.
                      str_lbits = int2string(lbits[k]);
                      str_rbits = int2binaryString(rbits[k]);

                      v1 = m[i]->getStrands(lbits[k]);
                      v2 = m[j]->getStrands(rbits[k]);
                      for(int a=0;a<2;a++)
                      for(int b=0;b<2;b++)
                         cout
                            << m[i]->getID() << str_lbits + "_" + suffix[a]+ "," 
                            << m[j]->getID() << str_rbits + "_" + suffix[b] 
                            << endl;
                   }
                }

                if(isBiMotif(m[i]->getType()) and 
                   isBiMotif(m[j]->getType()) ){
                   const int *TT_gate;
                   TT_gate = getTruthTable[m[i]->getType()];
                   
                   lbits.clear(); rbits.clear();
                   
                   for(int lbit2=0;lbit2<4;lbit2++)
                   for(int rbit2=0;rbit2<4;rbit2++)
                      if(TT_gate[3*lbit2+2] !=
                         rbit2%2) {
                         lbits.push_back(lbit2);
                         rbits.push_back(rbit2);
                      }

                   for(int k=0;k<lbits.size();k++){
                      //reset the string bits.
                      str_lbits = int2binaryString(lbits[k]);
                      str_rbits = int2binaryString(rbits[k]);

                      v1 = m[i]->getStrands(lbits[k]);
                      v2 = m[j]->getStrands(rbits[k]);
                      for(int a=0;a<2;a++)
                      for(int b=0;b<2;b++)
                         cout
                            << m[i]->getID() << str_lbits + "_" + suffix[a]+ "," 
                            << m[j]->getID() << str_rbits + "_" + suffix[b] 
                            << endl;
                   }
                }
            }//end if g[i][j] is 2.
        }//end for loop
        return;
}

void readDomainEnergy(string &deFile){
    ifstream file;
    file.open(deFile.c_str());
    string line;
    string s1;
    float energy;
    while(std::getline(file, line)){
        istringstream iss(line);
        iss >> s1;
        iss >> energy;
        if(s1.size() == 0){
            cerr << "energy encountered wrong input line" << endl;
        }
        domainEnergy[s1] = energy;
    }
    return;
}

void associateStrandWithObject(int &input,
        map<int, DNAMotif*> &m, 
        map<int, int> &names
    ){
    map<string, string>::iterator it;
    string s;
    for(it = inputStrands[input].begin();
        it!= inputStrands[input].end();
        it++){
        s = nupackMap[it->first];
        int num = s[0]-'A';
        //it->first is the name of the strand.
        //it->second is the sequence of the strand.
        if(s != "MB")
            m[names[num]]->parseStrand(s, it->second);
    }
    return;
}


void print_usage(void){
    cout << "Usage: ./a.out -f <enable fanout> -r <enable reusable> [file] [mapfile] [readFromNupackfile] [writeDesign.np] [domainenergy.txt]"
        << endl;
}

int main(int argc, char *argv[])
{
    try{

        int n;
        vector <vector<int> >g;
        map<int, DNAMotif* > m;
        map<int, int> names;


        if(argc < 6 || argc > 9){
            print_usage();
            exit(-1);
        }
    

        cout << argv[1] << endl;
        initialize();

        char c;
        while ((c = getopt (argc, argv, "fr")) != -1)
            switch (c)
            {
                case 'f':
                    DNA::fanOutEnabled = true;
                    break;
                case 'r':
                    DNA::reusable = true;
                    break;
                default:
                    print_usage();
                    return -1;
            }
        
        if(!createInputFile(n, names, g, m, argv[optind])){
            return -1;
        }
        cerr << "Input file read" << endl;

        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++)
                cout << g[i][j] << " ";
            cout << endl;
        }

        vector<int> sorted;
        sorted.clear();
        topologicalSort(g, sorted);

        map<int,int> fanOut;
        for(int i=0;i<n;i++)
            fanOut[i] = 0;

        //Construction of DNA Domains
        //input[0] refers to left input, and 
        //input[1] refers to right input.
        DNAMotif *input[2];
        for(int i=0;i<n;i++){
            //cout << "starting with vertex: " << sorted[i] << endl;
            int ct = 0;
            int index[2];
            input[0] = input[1] = NULL;
            index[0] = index[1] = -1;
            for(int j=0;j<n;j++){
                //Getting the two gates that input into gate j.
                if(g[j][sorted[i]] == 1){
                    input[0] = m[j];
                    index[0] = j;
                    ct++;
                    fanOut[j]++;
                }
                else if(g[j][sorted[i]] == 2){
                    input[1] = m[j];
                    index[1] = j;
                    ct++;
                    fanOut[j]++;
                }
            }
            if(ct == 2){
                if(!DNA::fanOutEnabled)
                    m[sorted[i]]->constructFromInput(input[0], input[1]);
                else {
                    m[sorted[i]]->constructFromInput(input[0], input[1], 
                        m[index[0]]->getFanOutMultiplier()>1?fanOut[index[0]]:0, 
                        m[index[1]]->getFanOutMultiplier()>1?fanOut[index[1]]:0);
                }
            }
            else if (ct == 1){
                //copy over values into index[0]
                if(index[0] == -1){
                    index[0] = index[1];
                    input[0] = input[1];
                }

                if(!DNA::fanOutEnabled)
                    m[sorted[i]]->constructFromInput(input[0], input[0]);
                else{
                    m[sorted[i]]->constructFromInput(input[0], input[0], 
                        m[index[0]]->getFanOutMultiplier()>1?fanOut[index[0]]:0);
                }
            }
        }
        cerr << "Domains constructed" << endl;

        //Figuring out the concentration multiplier per motif.
        //algorithm: set the concentration for each edge g[i][j].
        vector<vector<int> > edgeWeights;
        edgeWeights.resize(n);
        for(int a=0;a<n;a++){
            edgeWeights[a].resize(n);
            for(int b=0;b<n;b++)
                edgeWeights[a][b] = 0;
        }
        cout << "all edge weights have been assigned" << endl;
        int concMult = 0;
        for(int i=n-1;i>=0;i--){
            //sorted[i] is 1x concentration.
            concMult = 0;
            for(int j=0;j<n;j++)
                if(g[sorted[i]][j])
                    concMult += m[j]->getConcMultiplier();
            if(concMult)
                m[sorted[i]]->setConcMultiplier(concMult);
        }

        for(int i=0;i<n;i++){
            cout << "Motif no.: "<< m[sorted[i]]->getID() << " #";
            m[sorted[i]]->printConcentration();
            cout << endl;
        }

        cout << "directive sample 500.0 100\n"
            << "directive scale 100.0\n\n"
            << "def TMP() = <tmp>\n";


        for(int i=0;i<n;i++){
            m[sorted[i]]->printForDSD();
        }

        cout << "( 1*TMP()" << endl;;
        for(int i=0;i<n;i++){
            m[sorted[i]]->printForDSDWithConc();
        }
        cout << ")" << endl;

        cout << "Printing in CRN - for checking with LBS" << endl;
        cout << "directive sample 500.0 100\n\n"
            << "rate uni_forward = 5e-2;\n"
            << "rate bi_forward = 1e-3;\n"
            << "rate tri_forward = 3e-5;\n";


        //based on arguments.
        CIRCUIT_TYPE circ = ONE_TIME;
        if(DNA::reusable){
            circ = REUSABLE;
        }

        for(int i=0;i<n;i++){
            m[sorted[i]]->printForCRN(circ);
        }

        //Cross Validation to ensure no two gates will
        //interact with one another.
        vector<STRAND> v1, v2;
        int it[4];
        for(int i=0;i<n;i++)
        for(int j=i+1;j<n;j++){
            //if direct edge in or out of gate
            if(g[i][j] || g[j][i]) continue;

            //if both gates receive an input from
            //the same gate
            bool flag = true;
            for(int k=0;k<n;k++){
                if(g[k][i] && g[k][j]) {
                    flag = false;
                    break;
                }
            }
            if(!flag) continue;

            v1 = m[i]->getStrands();
            v2 = m[j]->getStrands();
            it[0] = v1.size();
            it[1] = v2.size();
            for(it[2]=0;it[2]<it[0];it[2]++)
            for(it[3]=0;it[3]<it[1];it[3]++){
                if(v1[it[2]].compareDomains(v2[it[3]])){
                    cout << "INTERACTIONS FOUND between: " 
                        << i << ","
                        << j
                        << endl;
                    switch(m[i]->getType()){
                        case AND_MOTIF: cout << "AND Gate: "; break;
                        case NAND_MOTIF: cout << "NAND Gate: "; break;
                        case STRAND_MOTIF: cout << "Strand: "; break;
                        default: cout << "DUNNO Gate: "; break;
                    }
                    cout << i << ":";
                    v1[it[2]].print();
                    switch(m[j]->getType()){
                        case AND_MOTIF: cout << "AND Gate: "; break;
                        case NAND_MOTIF: cout << "NAND Gate: "; break;
                        case STRAND_MOTIF: cout << "Strand: "; break;
                        default: cout << "DUNNO Gate: "; break;
                    }
                    cout << j << ":";
                    v2[it[3]].print();
                }   
            }
        }


        //print for NUPACK
        //Domains that should not be 
        excludedDomain.insert("x3");
        excludedDomain.insert("y2");
        excludedDomain.insert("D1");
        string npMap = string(argv[2]);
        readNupackMapFile(npMap);
        string npOfile = string(argv[3]);
        int numDesigns = readNupackOutputFile(npOfile);
        string nupackDesignFile = string(argv[4]);
        cout << numDesigns << endl;
        string deFile = string(argv[5]);
        readDomainEnergy(deFile);

        cerr << numDesigns << ":no. of designs" << endl;
        //associateStrandWithObject
        for(int i=0;i<numDesigns;i++){
            associateStrandWithObject(i,m,names);
            ostringstream oss;
            oss << nupackDesignFile;
            oss << "_" << i << ".np";
            printForNupackDesign(n,m,names,g, 
            oss.str());
            cerr << "Next Design" << endl;
        }

        createCheckMapForPythonScript(n,m,names,g);
        deleteInput(n, m);
    }

    catch (const char  *err){
        cout << "Caught ERROR: " << err << endl;
        return -1;
    }
}
