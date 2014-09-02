#include<iostream>
#include<string>
#include<map>
#include<vector>
#include<string.h>
#include<sstream>
#include<fstream>
using namespace std;

enum MOTIF_TYPE { NAND_MOTIF, STRAND_MOTIF, DNA_MOTIF, 
    AND_MOTIF, OR_MOTIF, NOT_MOTIF, NOR_MOTIF};
enum TOEHOLD_TYPE {NORMAL, DSD};
const string DOMAINPREFIX = "c";
const string PREFIX2 = "k";

class STRAND;

typedef class DNAMotif{
    public:
        DNAMotif(){
            istype = DNA_MOTIF;
            //concentration multiplier - default 1
            concX = 1;
            return;
        }
        MOTIF_TYPE getType(void){
            return istype;
        }
        void setType(MOTIF_TYPE type){
            istype = type;
        }

        void setConcentrationMultiplier(int m){
            concX = m;
        }

        int getConcentrationMultiplier(void){
            return concX;
        }

        virtual string getID(void){
           cout << "in the dnamotif getID" << endl;
        }

        virtual void print(void){}
        virtual void printConcentration(void){}
        virtual void printForDSD(string motif = "GATE"){
            cout << "in the dnamotif printfordsd" << endl;
        }
        virtual void printForCRN(void){
            cerr << "in the dnamotif printforcrn" << endl;
        }
        virtual void constructFromInput(DNAMotif *I1, DNAMotif *I2, int version=1){
            cout << "in the dnamotif construct from input " << endl;
        }
        virtual string getDomain(int idx, TOEHOLD_TYPE type = NORMAL){
            cout << "in the dnamotif getDomain " << endl;
        }
        virtual string getComplementDomain(int idx){
            cout << "in the dnamotif getComplementDomain " << endl;
        }
        virtual vector<STRAND> getStrands(void){
            cout << "in the dnamotif getStrands " << endl;
        }

    private:
        MOTIF_TYPE istype;
        int concX;
}DNAMotif;

typedef class STRAND : public DNAMotif{
    public:
        STRAND(){}
        STRAND(int num){
            DNAMotif::setType(STRAND_MOTIF);
            id = "A";
            id[0] += num;
        }
        
        STRAND(const string &s, int num){
            init(s);
            id = "A";
            id[0] += num;
        }

        void init(const string &s){
            DNAMotif::setType(STRAND_MOTIF);
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

        void print(void){
            for(int i=0;i<name.size();i++){
                cout << name[i];
                if(complement[i]) cout << "'";
                cout << " ";
            }
            cout << endl;
        }
        void printConcentration(void){
            cout << "Concentration: " << DNAMotif::getConcentrationMultiplier()
                << "X" << endl;
        }
        void printForDSD(string motif = "STRAND"){
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

        void printForCRN(void){
           cout << "init " << id << "0 100 |" << endl;
           cout << id[0] << "0 + " << id[0] << "1 -> {bi_forward} " << "S_" << id[0] << " |" << endl;
           cout << "S_" << id[0] << " + " << id[0] << "0 -> {bi_forward} " << "3" << id[0] << "0 |" << endl;
           cout << "S_" << id[0] << " + " << id[0] << "1 -> {bi_forward} " << "3" << id[0] << "1 |" << endl;
        }
        
        string getDomain(int idx, TOEHOLD_TYPE type = NORMAL){
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

        string getComplementDomain (int idx){
            if(name.size() <= idx){
                throw "REQUEST FOR INDEX NOT PRESENT IN STRAND";
            }

            if(complement[idx] == 1)
                return name[idx];
            else return name[idx]+"*";
        }

        void push(string s){
            if(s[s.length()-1] == '*'){
                name.push_back(s.substr(0,s.length()-1));
                complement.push_back(1);
            }
            else {
                name.push_back(s);
                complement.push_back(0);
            }
        }

        static string getNewDomain(string prefix = DOMAINPREFIX){
            if(STRAND::nextNumber.find(prefix) == STRAND::nextNumber.end())
                STRAND::nextNumber[prefix] = 0;

            ostringstream oss1;
            oss1 << prefix << STRAND::nextNumber[prefix];
            STRAND::nextNumber[prefix]++;

            return oss1.str();
        }

        bool compareDomains(STRAND to){
            vector<string> biDomain1, biDomain2;
            biDomain1 = getConcatenatedDomains(2,2);
            biDomain2 = to.getConcatenatedDomains(2,2);
            for(int i=0;i<biDomain1.size();i++)
            for(int j=0;j<biDomain2.size();j++)
                if(biDomain1[i] == biDomain2[j])
                    return true;
            return false;
        }

        vector<STRAND> getStrands(void){
            vector<STRAND> ret;
            ret.push_back(*this);
            return ret;
        }

        string getID(void){
           return id;
        }

    private:
        vector<string> name;
        vector<bool> complement;
        static map<string, int> nextNumber;
        string id;

        void createFromVector(vector<string> s){
            name.resize(s.size());
            complement.resize(s.size());

            for(int i=0;i<s.size();i++){
                complement[i] = 0;
                name[i] = s[i];

                int last = s[i].length() - 1;
                if(s[i][last] == '*'){
                    complement[i] = 1;
                    name[i].erase(name[i].end()-1);
                }
            }
        }

        vector<string> getConcatenatedDomains(int len, int which=2){
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
}STRAND;

map<string, int> STRAND::nextNumber;

typedef class NAND : public DNAMotif {
    public:
        NAND(){}
        NAND(int num){
            DNAMotif::setType(NAND_MOTIF);
            id[0] = "A";
            id[0][0] += num;
        }
        STRAND getAnd(void){
            return andStrand;
        }
        STRAND getGate(void){
            return gateStrand;
        }
        STRAND getNot(void){
            return notStrand;
        }
        void constructFromInput(DNAMotif *I1, DNAMotif *I2, int version=1){
            //andStrand
            andStrand.push(I1->getDomain(2));
            andStrand.push(I1->getDomain(3));
            andStrand.push(STRAND::getNewDomain());
            andStrand.push(I2->getDomain(1));
            andStrand.push(I2->getDomain(2));

            //gateStrand
            gateStrand.push(I2->getComplementDomain(3));
            gateStrand.push(I2->getComplementDomain(2));
            gateStrand.push(I2->getComplementDomain(1));
            gateStrand.push(andStrand.getComplementDomain(2));
            gateStrand.push(I1->getComplementDomain(3));
            gateStrand.push(I1->getComplementDomain(2));
            gateStrand.push(I1->getComplementDomain(1));

            //notStrand
            notStrand.push(andStrand.getComplementDomain(3));
            notStrand.push(andStrand.getComplementDomain(2));
            notStrand.push(andStrand.getComplementDomain(1));
            notStrand.push(STRAND::getNewDomain(PREFIX2));

            constructCRN(I1->getID(), I2->getID());
        }

        string getDomain(int idx, TOEHOLD_TYPE type = NORMAL){
            return notStrand.getDomain(idx, type);
        }

        string getComplementDomain (int idx){
            return notStrand.getComplementDomain(idx);
        }

        void print(void){
            cout << "andStrand:: ";
            andStrand.print();
            cout << "gateStrand:: ";
            gateStrand.print();
            cout << "notStrand:: ";
            notStrand.print();
        }
        void printConcentration(void){
            cout << "Conc: " << DNAMotif::getConcentrationMultiplier()
                << "X" << endl;
        }
        void printForDSD(string motif="AND"){
            cout << "def " 
                << motif
                << "() = {"
                << gateStrand.getDomain(6, DSD) 
                << "}["
                << andStrand.getDomain(0,DSD)
                << andStrand.getDomain(1,DSD)
                << andStrand.getDomain(2,DSD)
                << andStrand.getDomain(3,DSD)
                << andStrand.getDomain(4,DSD)
                << "]{"
                << gateStrand.getDomain(0, DSD) 
                << "}"
                << endl;

            motif = "NOT_" + motif;
            notStrand.printForDSD(motif);
        }

        void printForCRN(void){
           cout << "init " << id[0] << "0 100 |" << endl;
           cout << id[0] << "0 + " << id[0] << "1 -> {bi_forward} " << "S_" << id[0] << " |" << endl;
           cout << "S_" << id[0] << " + " << id[0] << "0 -> {bi_forward} " << "3" << id[0] << "0 |" << endl;
           cout << "S_" << id[0] << " + " << id[0] << "1 -> {bi_forward} " << "3" << id[0] << "1 |" << endl;
           cout << id[1] << "0 + " << id[2] << "0 + " << id[0] << "0 -> {tri_forward} " << id[0] << "1 |" << endl;
           cout << id[1] << "0 + " << id[2] << "1 + " << id[0] << "0 -> {tri_forward} " << id[0] << "1 |" << endl;
           cout << id[1] << "1 + " << id[2] << "0 + " << id[0] << "0 -> {tri_forward} " << id[0] << "1 |" << endl;
           cout << id[1] << "1 + " << id[2] << "1 + " << id[0] << "1 -> {tri_forward} " << id[0] << "0 |" << endl;
        }

        vector<STRAND> getStrands(void){
            vector<STRAND> ret;
            ret.push_back(andStrand);
            ret.push_back(gateStrand);
            ret.push_back(notStrand);
            return ret;
        }

        string getID(void){
           return id[0];
        }

    private:
        STRAND andStrand;
        STRAND gateStrand;
        STRAND notStrand;
        string id[3];
        void constructCRN(const string &id1, const string &id2){
           id[1] = id1;
           id[2] = id2;
        }
}NAND;

typedef class AND : public DNAMotif {
    public:
        AND(){}
        AND(int num){
            DNAMotif::setType(AND_MOTIF);
            id[0] = "A";
            id[0][0] += num;
        }
        STRAND getAnd(void){
            return andStrand;
        }
        STRAND getGate(void){
            return gateStrand;
        }
        void constructFromInput(DNAMotif *I1, DNAMotif *I2, int version=1){
            //andStrand
            andStrand.push(I1->getDomain(2));
            andStrand.push(I1->getDomain(3));
            andStrand.push(STRAND::getNewDomain());
            andStrand.push(I2->getDomain(1));
            andStrand.push(I2->getDomain(2));

            //gateStrand
            gateStrand.push(I2->getComplementDomain(3));
            gateStrand.push(I2->getComplementDomain(2));
            gateStrand.push(I2->getComplementDomain(1));
            gateStrand.push(andStrand.getComplementDomain(2));
            gateStrand.push(I1->getComplementDomain(3));
            gateStrand.push(I1->getComplementDomain(2));
            gateStrand.push(I1->getComplementDomain(1));

            constructCRN(I1->getID(), I2->getID());
        }

        string getDomain(int idx, TOEHOLD_TYPE type = NORMAL){
            return andStrand.getDomain(idx, type);
        }

        string getComplementDomain (int idx){
            return andStrand.getComplementDomain(idx);
        }

        void print(void){
            cout << "andStrand:: ";
            andStrand.print();
            cout << "gateStrand:: ";
            gateStrand.print();
        }
        void printConcentration(void){
            cout << "Conc: " << DNAMotif::getConcentrationMultiplier()
                << "X" << endl;
        }
        void printForDSD(string motif="AND"){
            cout << "def " 
                << motif
                << "() = {"
                << gateStrand.getDomain(6, DSD) 
                << "}["
                << andStrand.getDomain(0,DSD)
                << andStrand.getDomain(1,DSD)
                << andStrand.getDomain(2,DSD)
                << andStrand.getDomain(3,DSD)
                << andStrand.getDomain(4,DSD)
                << "]{"
                << gateStrand.getDomain(0, DSD) 
                << "}"
                << endl;

        }

        void printForCRN(void){
           cout << "init " << id[0] << "0 100 |" << endl;
           cout << id[0] << "0 + " << id[0] << "1 -> {bi_forward} " << "S_" << id[0] << " |" << endl;
           cout << "S_" << id[0] << " + " << id[0] << "0 -> {bi_forward} " << "3" << id[0] << "0 |" << endl;
           cout << "S_" << id[0] << " + " << id[0] << "1 -> {bi_forward} " << "3" << id[0] << "1 |" << endl;
           cout << id[1] << "0 + " << id[2] << "0 + " << id[0] << "1 -> {tri_forward} " << id[0] << "0 |" << endl;
           cout << id[1] << "0 + " << id[2] << "1 + " << id[0] << "1 -> {tri_forward} " << id[0] << "0 |" << endl;
           cout << id[1] << "1 + " << id[2] << "0 + " << id[0] << "1 -> {tri_forward} " << id[0] << "0 |" << endl;
           cout << id[1] << "1 + " << id[2] << "1 + " << id[0] << "0 -> {tri_forward} " << id[0] << "1 |" << endl;
        }

        vector<STRAND> getStrands(void){
            vector<STRAND> ret;
            ret.push_back(andStrand);
            ret.push_back(gateStrand);
            return ret;
        }

        string getID(void){
           return id[0];
        }

    private:
        STRAND andStrand;
        STRAND gateStrand;
        string id[3];
        void constructCRN(const string &id1, const string &id2){
           id[1] = id1;
           id[2] = id2;
        }
}AND;

typedef class NOT : public DNAMotif {
    public:
        NOT(){}
        NOT(int num){
            DNAMotif::setType(NOT_MOTIF);
            id[0] = "A";
            id[0][0] += num;
        }

        STRAND getNot(void){
            return notStrand;
        }

        void constructFromInput(DNAMotif *I1, DNAMotif *I2, int version=1){
            //notStrand
            notStrand.push(I1->getComplementDomain(3));
            notStrand.push(I1->getComplementDomain(2));
            notStrand.push(I1->getComplementDomain(1));
            notStrand.push(STRAND::getNewDomain(PREFIX2));

            constructCRN(I1->getID(), "");
        }

        string getDomain(int idx, TOEHOLD_TYPE type = NORMAL){
            return notStrand.getDomain(idx, type);
        }

        string getComplementDomain (int idx){
            return notStrand.getComplementDomain(idx);
        }

        void print(void){
            cout << "notStrand:: ";
            notStrand.print();
        }
        void printConcentration(void){
            cout << "Conc: " << DNAMotif::getConcentrationMultiplier()
                << "X" << endl;
        }
        void printForDSD(string motif="NOT"){
            motif = "NOT_" + motif;
            notStrand.printForDSD(motif);
        }

        void printForCRN(void){
           cout << "init " << id[0] << "0 100 |" << endl;
           cout << id[0] << "0 + " << id[0] << "1 -> {bi_forward} " << "S_" << id[0] << " |" << endl;
           cout << "S_" << id[0] << " + " << id[0] << "0 -> {bi_forward} " << "3" << id[0] << "0 |" << endl;
           cout << "S_" << id[0] << " + " << id[0] << "1 -> {bi_forward} " << "3" << id[0] << "1 |" << endl;
           cout << id[1] << "0 + " << id[0] << "0 -> {bi_forward} " << id[0] << "1 |" << endl;
           cout << id[1] << "1 + " << id[0] << "1 -> {bi_forward} " << id[0] << "0 |" << endl;
        }

        vector<STRAND> getStrands(void){
            vector<STRAND> ret;
            ret.push_back(notStrand);
            return ret;
        }

        string getID(void){
           return id[0];
        }

    private:
        STRAND notStrand;
        string id[3];
        void constructCRN(const string &id1, const string &id2){
           id[1] = id1;
           id[2] = id2;
        }
}NOT;

typedef class OR : public DNAMotif {
    public:
        OR(){}
        OR(int num){
            DNAMotif::setType(OR_MOTIF);
            id[0] = "A";
            id[0][0] += num;
        }
        STRAND getAnd(void){
            return tmpStrand;
        }
        STRAND getGate(void){
            return tmpStrand;
        }
        STRAND getNot(void){
            return tmpStrand;
        }
        void constructFromInput(DNAMotif *I1, DNAMotif *I2, int version=1){
            constructCRN(I1->getID(), I2->getID());
        }

        string getDomain(int idx, TOEHOLD_TYPE type = NORMAL){
            return "";
        }

        string getComplementDomain (int idx){
            return "";
        }

        void print(void){
        }
        void printConcentration(void){
            cout << "Conc: " << DNAMotif::getConcentrationMultiplier()
                << "X" << endl;
        }
        void printForDSD(string motif="OR"){
        }

        void printForCRN(void){
           cout << "init " << id[0] << "0 100 |" << endl;
           cout << id[0] << "0 + " << id[0] << "1 -> {bi_forward} " << "S_" << id[0] << " |" << endl;
           cout << "S_" << id[0] << " + " << id[0] << "0 -> {bi_forward} " << "3" << id[0] << "0 |" << endl;
           cout << "S_" << id[0] << " + " << id[0] << "1 -> {bi_forward} " << "3" << id[0] << "1 |" << endl;
           cout << id[1] << "0 + " << id[2] << "0 + " << id[0] << "1 -> {tri_forward} " << id[0] << "0 |" << endl;
           cout << id[1] << "0 + " << id[2] << "1 + " << id[0] << "0 -> {tri_forward} " << id[0] << "1 |" << endl;
           cout << id[1] << "1 + " << id[2] << "0 + " << id[0] << "0 -> {tri_forward} " << id[0] << "1 |" << endl;
           cout << id[1] << "1 + " << id[2] << "1 + " << id[0] << "0 -> {tri_forward} " << id[0] << "1 |" << endl;
        }

        vector<STRAND> getStrands(void){
            vector<STRAND> ret;
            return ret;
        }

        string getID(void){
           return id[0];
        }

    private:
        STRAND tmpStrand;
        string id[3];
        void constructCRN(const string &id1, const string &id2){
           id[1] = id1;
           id[2] = id2;
        }
}OR;

typedef class NOR : public DNAMotif {
    public:
        NOR(){}
        NOR(int num){
            DNAMotif::setType(NOR_MOTIF);
            id[0] = "A";
            id[0][0] += num;
        }
        STRAND getAnd(void){
            return tmpStrand;
        }
        STRAND getGate(void){
            return tmpStrand;
        }
        STRAND getNot(void){
            return tmpStrand;
        }
        void constructFromInput(DNAMotif *I1, DNAMotif *I2, int version=1){
            constructCRN(I1->getID(), I2->getID());
        }

        string getDomain(int idx, TOEHOLD_TYPE type = NORMAL){
            return "";
        }

        string getComplementDomain (int idx){
            return "";
        }

        void print(void){
        }
        void printConcentration(void){
            cout << "Conc: " << DNAMotif::getConcentrationMultiplier()
                << "X" << endl;
        }
        void printForDSD(string motif="NOR"){
        }

        void printForCRN(void){
           cout << "init " << id[0] << "0 100 |" << endl;
           cout << id[0] << "0 + " << id[0] << "1 -> {bi_forward} " << "S_" << id[0] << " |" << endl;
           cout << "S_" << id[0] << " + " << id[0] << "0 -> {bi_forward} " << "3" << id[0] << "0 |" << endl;
           cout << "S_" << id[0] << " + " << id[0] << "1 -> {bi_forward} " << "3" << id[0] << "1 |" << endl;
           cout << id[1] << "0 + " << id[2] << "0 + " << id[0] << "0 -> {tri_forward} " << id[0] << "1 |" << endl;
           cout << id[1] << "0 + " << id[2] << "1 + " << id[0] << "1 -> {tri_forward} " << id[0] << "0 |" << endl;
           cout << id[1] << "1 + " << id[2] << "0 + " << id[0] << "1 -> {tri_forward} " << id[0] << "0 |" << endl;
           cout << id[1] << "1 + " << id[2] << "1 + " << id[0] << "1 -> {tri_forward} " << id[0] << "0 |" << endl;
        }

        vector<STRAND> getStrands(void){
            vector<STRAND> ret;
            return ret;
        }

        string getID(void){
           return id[0];
        }

    private:
        STRAND tmpStrand;
        string id[3];
        void constructCRN(const string &id1, const string &id2){
           id[1] = id1;
           id[2] = id2;
        }
}NOR;

/*
void createInputTurberfield(int &nodes, vector<vector <int> > &g, map<int, DNAMotif*> &m){
    //F = (X and Y) or (Y and not Z)
    nodes = 8;

    //4 nand gates, 3 inputs, 1 output
    //0,1,2 are inputs
    //3,4,5,6 are nands
    //7 is output
    g.resize(nodes);
    for(int i=0;i<nodes;i++){
        g[i].resize(nodes);
        for(int j=0;j<nodes;j++){
            g[i][j] = 0;
        }
    }

    m[0] = new STRAND ("p1 x1 x2 x3 n1");
    m[1] = new STRAND ("p2 y1 y2 y3 n2");
    m[2] = new STRAND ("p3 z1 z2 z3 n3");
    m[3] = new NAND;
    m[4] = new NAND;
    m[5] = new NAND;
    m[6] = new NAND;
    m[7] = new STRAND;

    //m[0]->print();

    g[0][3] = 1;
    g[1][3] = 1;
    g[1][5] = 1;
    g[2][4] = 1;
    g[3][6] = 1;
    g[4][5] = 1;
    g[5][6] = 1;
    g[6][7] = 1;
}
*/

void deleteInput(int &n, map<int, DNAMotif *> &m){
    for(int i=0;i<n;i++)
        delete m[i];
}

void createInputFile(int &nodes, map<int, int> &names, vector<vector <int> > &g, map<int, DNAMotif*> &m, char* filename){
    ifstream inputfile;
    inputfile.open(filename);
    if(!inputfile) cout << "Error in file:" << filename << endl;
    string nodesStr;
    std::getline(inputfile, nodesStr);
    istringstream issN(nodesStr);
    issN >> nodes;
    string gate, line;
    int num, e1, e2;

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
    while(std::getline(inputfile, gate)){
        cout << "gate = " 
            << gate
            << endl;

        if(gate == "nand"){
            std::getline(inputfile, line);
            istringstream iss(line);
            while(iss >> num){
                names[num] = ct;
                m[ct] = new NAND(num);
                ct++;
            }
        }
        else if(gate == "and"){
            std::getline(inputfile, line);
            istringstream iss(line);
            while(iss >> num){
                names[num] = ct;
                m[ct] = new AND(num);
                ct++;
            }
        }
        else if(gate == "or"){
            std::getline(inputfile, line);
            istringstream iss(line);
            while(iss >> num){
                names[num] = ct;
                m[ct] = new OR(num);
                ct++;
            }
        }
        else if(gate == "not"){
            std::getline(inputfile, line);
            istringstream iss(line);
            while(iss >> num){
                names[num] = ct;
                m[ct] = new NOT(num);
                ct++;
            }
        }
        else if(gate == "nor"){
            std::getline(inputfile, line);
            istringstream iss(line);
            while(iss >> num){
                names[num] = ct;
                m[num] = new NOR(num);
                ct++;
            }
        }
        else if(gate == "input"){
            std::getline(inputfile, line);
            istringstream iss(line);
            while(iss >> e1){
                std::getline(inputfile, line);
                names[e1] = ct;
                m[ct] = new STRAND(line,e1);
                ct++;
            }
        }
        else if(gate == "output"){
            std::getline(inputfile, line);
            istringstream iss(line);
            while(iss >> num){
                names[num] = ct;
                m[ct] = new STRAND(num);
                ct++;
            }
        }
        else if(gate == "edges"){
            while(std::getline(inputfile, line)){
                istringstream edgeIss(line);
                edgeIss >> e1 >> e2;
                g[names[e1]][names[e2]] = 1;
            }
        }
    }
    inputfile.close();
}

/*
void createInputEg1(int &nodes, vector<vector <int> > &g, map<int, DNAMotif*> &m){
    nodes = 11;

    //7 nand gates, 3 inputs, 1 output
    g.resize(nodes);
    for(int i=0;i<nodes;i++){
        g[i].resize(nodes);
        for(int j=0;j<nodes;j++){
            g[i][j] = 0;
        }
    }

    m[6] = new STRAND ("p1 x1 x2 x3 n1");
    m[9] = new STRAND ("p2 y1 y2 y3 n2");
    m[4] = new STRAND ("p3 z1 z2 z3 n3");
    m[0] = new NAND;
    m[8] = new NAND;
    m[2] = new NAND;
    m[10] = new NAND;
    m[3] = new NAND;
    m[1] = new NAND;
    m[5] = new NAND;
    m[7] = new STRAND;

    g[6][0] = 1;
    g[9][2] = 1;
    g[9][0] = 1;
    g[9][10] = 1;
    g[9][8] = 1;
    g[4][8] = 1;
    g[0][2] = 1;
    g[0][10] = 1;
    g[8][1] = 1;
    g[2][3] = 1;
    g[2][1] = 1;
    g[10][3] = 1;
    g[3][5] = 1;
    g[1][5] = 1;
    g[5][7] = 1;
}
*/

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
                if(g[j][i] == 1) break;
            }
            if(j == nodes) break;
        }
        ret.push_back(i);
        visited[i] = true;
        done--;

        for(j=0;j<nodes;j++)
            if(g[i][j] == 1)
                g[i][j] = 2;

    }
    return;
}


int main(int argc, char *argv[])
{
    try{

        int n;
        vector <vector<int> >g;
        map<int, DNAMotif* > m;
        map<int, int> names;
        cout << argv[1] << endl;
        createInputFile(n, names, g, m, argv[1]);

        vector<int> sorted;
        sorted.clear();
        topologicalSort(g, sorted);

        //Construction of DNA Domains
        DNAMotif *input[2];
        for(int i=0;i<n;i++){
            //cout << "starting with vertex: " << sorted[i] << endl;
            int ct = 0;
            input[0] = input[1] = NULL;
            for(int j=0;j<n;j++){
                if(g[j][sorted[i]] != 0){
                    //cout << "trying: " << j << "," << sorted[i] << ":" << g[j][sorted[i]] << endl;
                    input[ct++] = m[j];
                }
            }
            if(ct == 2)
                m[sorted[i]]->constructFromInput(input[0], input[1]);
            else if (ct == 1)
                m[sorted[i]]->constructFromInput(input[0], input[0]);
            //m[sorted[i]]->print();
            //cout << endl << endl;
        }

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

        for(int i=n-1;i>=0;i--){
            //sorted[i] is the node which is being considered. 
            //all nodes topologically after it already have their
            //concentration multipliers set.
            //


            //concentration multiplier of this node
            //is the sum of edge weights of all outgoing
            //nodes.
            int multiplier = 0;
            for(int j=0;j<n;j++)
                if(g[sorted[i]][j] != 0)
                    multiplier += edgeWeights[sorted[i]][j];
            

            //incase the node doesn't have any output edges,
            //set it to 1 by default.
            if(multiplier < 1) multiplier = 1;

            m[sorted[i]]->setConcentrationMultiplier(multiplier);

            //set the reqd. conc. of all incoming edges to sorted[i], 
            //as this multiplier
            for(int j=0;j<n;j++)
                if(g[j][sorted[i]] != 0)
                    edgeWeights[j][sorted[i]] = multiplier;
            cout << "looping around.." << endl;
        }

        for(int i=0;i<n;i++){
            cout << "Motif no.: "<< sorted[i] << " #";
            m[sorted[i]]->printConcentration();
            cout << endl;
        }

        cout << "directive sample 50000.0 100\n"
            << "directive scale 100.0\n\n"
            << "def TMP() = <tmp>\n";


        ostringstream oss_conc;
        oss_conc << "( 1*TMP()\n";
        for(int i=0;i<n;i++){
            ostringstream oss;
            if(m[sorted[i]]->getType() == STRAND_MOTIF){
                oss << "INPUT_" << sorted[i];
                oss_conc << "| "
                    << m[sorted[i]]->getConcentrationMultiplier()*100
                    << "*"
                    << "INPUT_"
                    << sorted[i]
                    << "()\n";
            }
            else if(m[sorted[i]]->getType() == NAND_MOTIF){
                oss << "NAND_" << sorted[i];
                oss_conc << "| "
                    << m[sorted[i]]->getConcentrationMultiplier()*100
                    << "*"
                    << "NAND_"
                    << sorted[i]
                    << "()\n";
                oss_conc << "| "
                    << m[sorted[i]]->getConcentrationMultiplier()*100
                    << "*"
                    << "NOT_NAND_"
                    << sorted[i]
                    << "()\n";
            }
            m[sorted[i]]->printForDSD(oss.str());
        }
        oss_conc << ")";
        cout << oss_conc.str() << endl;

        cout << "Printing in CRN - for checking with LBS" << endl;
        cout << "directive sample 100.0 100\n"
            << "directive scale 100.0\n\n"
            << "rate bi_forward = 1e-3;\n"
            << "rate tri_forward = 1e3;\n";
        for(int i=0;i<n;i++){
            m[sorted[i]]->printForCRN();
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

        deleteInput(n, m);
    }

    catch (string err){
        cout << "Caught ERROR: " << err << endl;
        return -1;
    }
}
