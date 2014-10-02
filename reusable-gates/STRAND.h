#include<iostream>
#include<string>
#include<map>
#include<vector>
#include<string.h>
#include<sstream>
#include<fstream>
#include<unistd.h>
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
const string DOMAINPREFIX = "c";
const string PREFIX2 = "k";
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
        void getGate(const int *t);
        virtual string getID(void);
        virtual void print(void);
        virtual void printConcentration(void);
        virtual void printForDSD(string motif = "GATE");
        virtual void printForCRN(CIRCUIT_TYPE version=ONE_TIME);
        virtual void constructFromInput(DNAMotif *I1, DNAMotif *I2, int f1=0, int f2=0);
        virtual string getDomain(int bit, int idx, TOEHOLD_TYPE type = NORMAL);
        virtual string getComplementDomain(int bit, int idx);
        virtual vector<STRAND> getStrands(void);

    private:
        MOTIF_TYPE istype;
        int fanCount;
        int concX;

    protected:
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
        void printForCRN_disabled(CIRCUIT_TYPE version = ONE_TIME);
        string getDomain(int bit, int idx, TOEHOLD_TYPE type = NORMAL);
        string getComplementDomain (int bit, int idx);
        void push(string s);
        static string getNewDomain(string prefix = DOMAINPREFIX);
        bool compareDomains(STRAND to);
        vector<STRAND> getStrands(void);
        string getID(void);

    private:
        vector<string> name;
        vector<bool> complement;
        static map<string, int> nextNumber;

        void createFromVector(vector<string> s);
        vector<string> getConcatenatedDomains(int len, int which=2);
}STRAND;


/*
typedef class DNAMotif{
    public:
        DNAMotif(){
            istype = DNA_MOTIF;
            // fan out multiplier - default 0
            fanCount = 0;
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

        void setFanOutMultiplier(int m){
            fanCount = m;
        }

        int getFanOutMultiplier(void){
            return fanCount;
        }

        void incrementMultiplier(void){
            fanCount++;
        }

        void setConcMultiplier(int m){
            concX = m;
        }

        int getConcMultiplier(void){
            return concX;
        }

        virtual string getID(void){
           cout << "in the dnamotif getID" << endl;
        }

        virtual void print(void){}
        virtual void printConcentration(void){}
        virtual void printForDSD(string motif = "GATE"){

            motif = motif + "_" + id[0] + "_0";
            for(int i=0;i<4;i++){
                motif[motif.length()-1] += i;
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
                motif[motif.length()-1] -= i;
            }

            cerr << "in the dnamotif printfordsd" << endl;
        }
        virtual void printForCRN(CIRCUIT_TYPE version=ONE_TIME){
           int conc;
           const int *TT_GATE;
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

               if(istype == STRAND_MOTIF){
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
               if(istype == NOT_MOTIF){
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

        void getGate(const int *t){
            switch(t){
                case AND_MOTIF: t = TT_AND;
                                break;
                case NAND_MOTIF: t = TT_NAND;
                                break;
                case OR_MOTIF: t = TT_OR;
                                break;
                case NOR_MOTIF: t = TT_NOR;
                                break;
                case NOT_MOTIF: t = TT_NOT;
                                break;
                default: t = TT_AND;
            }
        }

        virtual void constructFromInput(DNAMotif *I1, DNAMotif *I2, int f1=0, int f2=0){
            cout << "for gates " << I1->getID() << " and " << I2->getID() 
                << ", the nums are : " << f1 << "," << f2 << endl;
            //andStrand
            
            string ids[2];
            ids[0] = id[0] + "0";
            ids[1] = id[0] + "1";
            const int *TT_gate;
            getGate(TT_gate);
            I1->print();
            I2->print();

            for(int i=0;i<4;i++){
                //andStrand
                andStrand[i]->push(I1->getDomain(TT_gate[3*i],2));
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
        virtual string getDomain(int bit, int idx, TOEHOLD_TYPE type = NORMAL){
            if(bit == 1)
                return andStrand[3]->getDomain(bit, idx, type);

            return andStrand[3]->getDomain(bit, idx, type);
            cerr << "in the dnamotif getDomain " << endl;
        }
        virtual string getComplementDomain(int bit, int idx){
            if(bit == 1)
                return andStrand[3]->getComplementDomain(bit, idx);

            return andStrand[0]->getComplementDomain(bit, idx);
            cerr << "in the dnamotif getComplementDomain " << endl;
        }
        virtual vector<STRAND> getStrands(void){
            cerr << "in the dnamotif getStrands " << endl;
        }

    private:
        MOTIF_TYPE istype;
        int fanCount;
        int concX;

    protected:
        string id[3];
        STRAND *andStrand[4];
        STRAND *gateStrand[4];

        void constructCRN(const string &id1, const string &id2, const int &f1, const int &f2){
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

}DNAMotif;
*/

