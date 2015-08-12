#include<iostream>
#include<string>
#include<map>
#include<vector>
#include<string.h>
#include<sstream>
#include<fstream>
#include<unistd.h>
#include<set>
using namespace std;

typedef class DNA {
    public:
        static const string sep;
        static const string sep2;
        static const int conc;
        static bool fanOutEnabled;
        static bool reusable;
}DNA;

const string DNA::sep = "-";
const string DNA::sep2 = "_";
const int DNA::conc = 100;
bool DNA::fanOutEnabled = false;
bool DNA::reusable = false;
enum MOTIF_TYPE { NAND_MOTIF, STRAND_MOTIF, DNA_MOTIF, BIT_MOTIF,
    AND_MOTIF, OR_MOTIF, NOT_MOTIF, NOR_MOTIF};
enum TOEHOLD_TYPE {NORMAL, DSD};
enum CIRCUIT_TYPE {ONE_TIME, REUSABLE, JIANG};
const string DOMAINPREFIX = "w";
const string PREFIX2 = "k";
const int TOEHOLD_DEFAULT_LENGTH = 5;
map<string, string> nupackMap;
map<int, map<string, string> > inputStrands;
map<string, float> domainEnergy;


//TT == TRUTH_TABLE
const int TT_NAND[] = 
    {0,0,1,
    0,1,1,
    1,0,1,
    1,1,0};

const int TT_AND[] = 
    {0,0,0,
    0,1,0,
    1,0,0,
    1,1,1};

const int TT_OR[] = 
    {0,0,0,
    0,1,1,
    1,0,1,
    1,1,1};

const int TT_NOR[] = 
    {0,0,1,
    0,1,0,
    1,0,0,
    1,1,0};

const int TT_NOT[] = 
    {0,1,
    1,0};

map<MOTIF_TYPE, string> getGatePrefix;
map<MOTIF_TYPE, const int*> getTruthTable;

class STRAND;

typedef class DNAMotif{
    public:
        DNAMotif();
        MOTIF_TYPE getType(void);
        void setType(MOTIF_TYPE type);
        void setFanOutMultiplier(int m);
        int getFanOutMultiplier(void);
        void incrementMultiplier(void);
        void setConcMultiplier(int m);
        int getConcMultiplier(void);
        virtual string getID(void);
        virtual void print(void);
        virtual void printConcentration(void);
        virtual string getDSDName(int num = 0);
        virtual void printForDSD(string motif = "GATE");
        virtual void printForDSDWithConc(void);
        virtual void printForCRN(CIRCUIT_TYPE version=ONE_TIME);
        virtual void constructFromInput(DNAMotif *I1, DNAMotif *I2, int f1=0, int f2=0);
        virtual string getDomain(int bit, int idx, TOEHOLD_TYPE type = NORMAL);
        virtual string getComplementDomain(int bit, int idx);
        virtual vector<STRAND> getStrands(int bit=-1);
        virtual set<pair<string, string> > getUniqueDomains(void);
        virtual vector<string> printNupackStructureAndSequence(
                map<int, DNAMotif*> &m, 
                map<int, int> &names, vector<vector<int> >&g);
        virtual vector<STRAND> getOutputStrandList(int bit=0);
        virtual int parseStrand(const string &name, const string &seq);

    private:
        MOTIF_TYPE istype;
        int fanCount;
        int concX;

    protected:
        int idNum;
        string id[3];
        STRAND *andStrand[4];
        STRAND *gateStrand[4];
        void constructCRN(const string &id1, const string &id2, const int &f1, const int &f2);

}DNAMotif;

typedef class STRAND : public DNAMotif{
    public:
        STRAND();
        STRAND(int num);
        STRAND(const string &s, int num);
        void init(const string &s);
        void print(void);
        void printConcentration(void);
        void printForDSD(string motif = "STRAND");
        string getDomain(int bit, int idx, TOEHOLD_TYPE type = NORMAL);
        string getComplementDomain (int bit, int idx);
        void push(string s);
        static string getNewDomain(string prefix = DOMAINPREFIX);
        bool compareDomains(STRAND to);
        vector<STRAND> getStrands(int bit=-1);
        string getID(void);
        int getNumberOfDomains(void);
        set<pair<string, string> > getUniqueDomains(void);
        int getDomainLength(int pos=0);
        string getName(void);
        int parseDomains(const string &seq);

    private:
        vector<string> name;
        vector<bool> complement;
        vector<int> domainSize;
        vector<string>domainSeq;
        static map<string, int> nextNumber;

        void createFromVector(vector<string> s);
        vector<string> getConcatenatedDomains(int len, int which=2);
}STRAND;
