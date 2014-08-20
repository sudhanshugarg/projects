#include<iostream>
#include<string>
#include<map>
#include<vector>
#include<string.h>
#include<sstream>
using namespace std;

enum MOTIF_TYPE { NAND_MOTIF, STRAND_MOTIF, DNA_MOTIF };
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

        virtual void print(void){}
        virtual void printConcentration(void){}
        virtual void constructFromInput(DNAMotif *I1, DNAMotif *I2, int version=1){
            cout << "in the dnamotif construct from input " << endl;
        }
        virtual string getDomain(int idx){
            cout << "in the dnamotif getDomain " << endl;
        }
        virtual string getComplementDomain(int idx){
            cout << "in the dnamotif getComplementDomain " << endl;
        }

    private:
        MOTIF_TYPE istype;
        int concX;
}DNAMotif;

typedef class STRAND : public DNAMotif{
    public:
        STRAND(){
            DNAMotif::setType(STRAND_MOTIF);
        }
        
        STRAND(vector<string> s){
            createFromVector(s);
        }

        STRAND(string s){
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
        
        string getDomain(int idx){
            if(name.size() <= idx){
                throw "REQUEST FOR INDEX NOT PRESENT IN STRAND";
            }

            if(complement[idx] == 0)
                return name[idx];
            else return name[idx]+"*";
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

    private:
        vector<string> name;
        vector<bool> complement;
        static map<string, int> nextNumber;

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
}STRAND;

map<string, int> STRAND::nextNumber;

typedef class NAND : public DNAMotif {
    public:
        NAND(){
            DNAMotif::setType(NAND_MOTIF);
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
        }

        string getDomain(int idx){
            return notStrand.getDomain(idx);
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
    private:
        STRAND andStrand;
        STRAND gateStrand;
        STRAND notStrand;
}NAND;

/*
 * typedef class AND : public DNAMotif{} AND;
 * typedef class OR : public DNAMotif{} OR;
 * typedef class NOT : public DNAMotif{} NOT;
 */



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

void deleteInputTurberfield(int &n, map<int, DNAMotif *> &m){
    for(int i=0;i<n;i++)
        delete m[i];
}

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

void deleteInputEg1(int &n, map<int, DNAMotif *> &m){
    for(int i=0;i<n;i++)
        delete m[i];
}
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

int main()
{
    try{

        int n;
        vector <vector<int> >g;
        map<int, DNAMotif* > m;
        createInputEg1(n, g, m);
        /*
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                cout << g[i][j]
                    << ",";
            }
            cout << endl;
        }*/

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
        }

        for(int i=0;i<n;i++){
            cout << "Motif no.: "<< sorted[i] << " #";
            m[sorted[i]]->printConcentration();

            m[sorted[i]]->print();
            cout << endl;
        }

        deleteInputEg1(n, m);
    }

    catch (string err){
        cout << "Caught ERROR: " << err << endl;
        return -1;
    }
}
