#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sys/time.h>
#include <string.h>
using namespace std;

double similarity_score(char a,char b);
double find_array_max(double array[],int length);
void insert_at(char arr[], int n, int idx, char val);
void checkfile(int open, char filename[]);
string read_sequence(ifstream& f);


double mu;
double delta;

typedef struct {
int score;
int max_i;
int max_j;
} Coordinate;

extern "C"

int SW(char *seq_a_, char *seq_b_, Coordinate *result);

int main(int argc,char *argv[]){
	char *seq_a = (char *) malloc( 1500 );
	char *seq_b = (char *) malloc( 1500 );

	memset( seq_a, 'A', 1500 );
	memset( seq_b, 'G', 1500 );
	seq_a[ 1499 ] = '\0';
	seq_b[ 1499 ] = '\0';

	Coordinate res;
	SW( seq_a, seq_b, &res );
}


int SW(char *seq_a_, char *seq_b_, Coordinate *result){

	string seq_a(seq_a_);
	string seq_b(seq_b_);
    
/*    // read info from arguments
    if(argc!=6){
        cout<<"Give me the proper number of input arguments:"<<endl<<"1 : mu"<<endl;
        cout<<"2 : delta"<<endl<<"3 : filename sequence A"<<endl<<"4 : filename sequence B"<<endl;
        cout<<"5 : maximal length N of sequences"<<endl;exit(1);
    }*/
    double mu,delta;
    mu = 0;
    delta = 0;
/*    // atof(argv[2]);
    /////////////////////////////////
    // give it the filenames
    char *nameof_seq_a = argv[3];
    char *nameof_seq_b = argv[4];
    int N_max = 1000; //atoi(argv[5]);
    string seq_a,seq_b;
    */
/*    // read the sequences into two vectors:
    ifstream stream_seq_b;                      // first define the input-streams for seq_a and seq_b
    stream_seq_b.open(nameof_seq_b);            // the same for seq_b
    checkfile(! stream_seq_b,nameof_seq_b);
    seq_b = read_sequence(stream_seq_b);
    ifstream stream_seq_a;
    stream_seq_a.open(nameof_seq_a);            // open the file for input
    checkfile(! stream_seq_a,nameof_seq_a);     // check, whether the file was opened successfully
    seq_a = read_sequence(stream_seq_a);
    
    */
    // string s_a=seq_a,s_b=seq_b;
    int N_a = seq_a.length();                     // get the actual lengths of the sequences
    int N_b = seq_b.length();
    //cout << N_b << " "<<N_a <<endl;

    ////////////////////////////////////////////////
    
    // initialize H
    double *H_cur_row = (double *) malloc( sizeof( double ) * (N_b + 1) );
    double *H_prev_row = (double *) malloc( sizeof( double ) * (N_b + 1) );

    // initialize final storage to get length of alignment
    double *Last_row = (double *) malloc( sizeof( double ) * (N_b + 1) );
    double *Last_col = (double *) malloc( sizeof( double ) * (N_a + 1) );

    for(int i=0;i<=N_b;i++){
    		H_cur_row[i]=0.;
    		H_prev_row[i]=0.;
    		Last_row[i]=0.;
    	}

    for(int i=0;i<=N_a;i++){
    		Last_col[i]=0.;
    	}

/*    double **H = (double **) malloc( sizeof( double * ) * (N_a + 1) );
    for(int i=0; i<=N_a; i++)
    {
    	H[ i ] = (double *) malloc( sizeof( double ) * (N_b + 1) );
    }

    for(int i=0;i<=N_a;i++){
        for(int j=0;j<=N_b;j++){
            H[i][j]=0.;
        }
    }*/
    
    double temp[3];
//    int I_i[N_a+1][N_b+1],I_j[N_a+1][N_b+1];     // Index matrices to remember the 'path' for backtracking

    char chr1;
    char chr2;
    for(int i=1;i<=N_a;i++){
        for(int j=1;j<=N_b;j++){
        	chr1 = seq_a.at(i-1);
        	chr2 = seq_b.at(j-1);
            temp[0] = H_prev_row[j-1]+similarity_score(chr1,chr2); //H[i-1][j-1]+similarity_score(chr1,chr2);
            temp[1] = H_prev_row[j]-delta; //H[i-1][j]-delta;
            temp[2] = H_cur_row[j-1]-delta; // H[i][j-1]-delta;
            H_cur_row[j] = find_array_max(temp,3);//H[i][j] = find_array_max(temp,3);
            if(j ==N_b){
            	Last_col[i] = H_cur_row[j];
            }

        }

        //We only keep track of two rows in the scoring matrix. When
        // we advance one row we update the current row to be the
        // new previous row. The old previous row is written over by zeros
/*        for(int k=0;k<=N_b;k++){
            cout << "previous: "<< H_prev_row[k] <<endl;
            cout << "Current: "<< H_cur_row[k] <<endl;
        	}*/

        swap( H_cur_row, H_prev_row );

        for(int j=0;j<=N_b;j++){
        	H_cur_row[j]=0.;
        	}

    }

    for(int j=0;j<=N_b;j++){
    	Last_row[j]=H_prev_row[j];
    	}




/*
    // here comes the actual algorithm
    char chr1;
    char chr2;
    for(int i=1;i<=N_a;i++){
        for(int j=1;j<=N_b;j++){
        	chr1 = seq_a.at(i-1);
        	chr2 = seq_b.at(j-1);
            temp[0] = H[i-1][j-1]+similarity_score(chr1,chr2);
            temp[1] = H[i-1][j]-delta;
            temp[2] = H[i][j-1]-delta;
            //temp[3] = 0.;
            H[i][j] = find_array_max(temp,3);
            switch(ind){
                case 0:                                  // score in (i,j) stems from a match/mismatch
                    I_i[i][j] = i-1;
                    I_j[i][j] = j-1;
                    break;
                case 1:                                  // score in (i,j) stems from a deletion in sequence A
                    I_i[i][j] = i-1;
                    I_j[i][j] = j;
                    break;
                case 2:                                  // score in (i,j) stems from a deletion in sequence B
                    I_i[i][j] = i;
                    I_j[i][j] = j-1;
                    break;
                case 3:                                  // (i,j) is the beginning of a subsequence
                    I_i[i][j] = i;
                    I_j[i][j] = j;
                    break;
            }
        }
    }

    // Print the matrix H to the console
    cout<<"**********************************************"<<endl;
    cout<<"The scoring matrix is given by  "<<endl<<endl;
    for(int i=1;i<=N_a;i++){
        for(int j=1;j<=N_b;j++){
            cout<<H[i][j]<<" ";
        }
        cout<<endl;
    }*/
    
    // search H for the maximal score TODO: Only on the last row and column (local alignment)
    double H_max = 0.;
    int i_max=0,j_max=0;
    for(int i=1;i<=N_a;i++){
    	if(Last_col[i]>H_max){
    	                H_max = Last_col[i];
    	                i_max = i;
    	                j_max = N_b;
    	            }
    }
    for(int j=1;j<=N_b;j++){
    	if(Last_row[j]>H_max){
    	                H_max = Last_row[j];
    	                i_max = N_a;
    	                j_max = j;
    	            }
    }


/*    for(int i=1;i<=N_a;i++){
        for(int j=1;j<=N_b;j++){
            if(H[i][j]>H_max){
                H_max = H[i][j];
                i_max = i;
                j_max = j;
            }
        }
    }*/
    result->score = H_max;
    result->max_i = i_max;
    result->max_j = j_max;


    //cout<<H_max<<endl;
/*
    // Backtracking from H_max
    int current_i=i_max,current_j=j_max;
    int next_i=I_i[current_i][current_j];
    int next_j=I_j[current_i][current_j];
    int tick=0;
    char consensus_a[N_a+N_b+2],consensus_b[N_a+N_b+2];

    while(((current_i!=next_i) || (current_j!=next_j)) && (next_j!=0) && (next_i!=0)){
        if(next_i==current_i)  consensus_a[tick] = '-';                  // deletion in A
        else                   consensus_a[tick] = seq_a.at(current_i-1);   // match/mismatch in A
        
        if(next_j==current_j)  consensus_b[tick] = '-';                  // deletion in B
        else                   consensus_b[tick] = seq_b.at(current_j-1);   // match/mismatch in B
        
        current_i = next_i;
        current_j = next_j;
        next_i = I_i[current_i][current_j];
        next_j = I_j[current_i][current_j];
        tick++;
    }

    
    // Output of the consensus motif to the console
    cout<<endl<<"***********************************************"<<endl;
    cout<<"The alignment of the sequences"<<endl<<endl;
    for(int i=0;i<N_a;i++){cout<<seq_a[i];}; cout<<"  and"<<endl;
    for(int i=0;i<N_b;i++){cout<<seq_b[i];}; cout<<endl<<endl;
    cout<<"is for the parameters  mu = "<<mu<<" and delta = "<<delta<<" given by"<<endl<<endl;
    for(int i=tick-1;i>=0;i--) cout<<consensus_a[i];
    cout<<endl;
    for(int j=tick-1;j>=0;j--) cout<<consensus_b[j];
    cout<<endl;*/

    free(H_cur_row);
    free(H_prev_row);
    free(Last_row);
    free(Last_col);
    return 1;
} // END of main




/////////////////////////////////////////////////////////////////////////////
// auxiliary functions used by main:
/////////////////////////////////////////////////////////////////////////////


/*void checkfile(int open, char filename[]){
    
    if (open){cout << "Error: Can't open the file "<<filename<<endl;exit(1);}
    else cout<<"Opened file "<<filename<<endl;
}*/

/////////////////////////////////////////////////////////////////////////////

double similarity_score(char a,char b){
    
    double result;
    if(a==b){
        result=1.;
    }
    else{
        result=-mu;
    }
    return result;
}

/////////////////////////////////////////////////////////////////////////////

double find_array_max(double array[],int length){
    
    double max = array[0];            // start with max = first element
    
    for(int i = 1; i<length; i++){
        if(array[i] > max){
            max = array[i];
        }
    }
    return max;                    // return highest value in array
}

/////////////////////////////////////////////////////////////////////////////
/*
string read_sequence(ifstream& f)
{
    // overflows.
    string seq;
    char line[5000];
    while( f.good() )
    {
        f.getline(line,5000);
        // 		cout << "Line:" << line << endl;
        if( line[0] == 0 || line[0]=='#' )
            continue;
        for(int i = 0; line[i] != 0; ++i)
        {
            int c = toupper(line[i]);
            if( c != 'A' && c != 'G' && c != 'C' && c != 'T' )
                continue;
            //cout << char(c);
            seq.push_back(char(c));
        }
    }
    return seq;
}*/
