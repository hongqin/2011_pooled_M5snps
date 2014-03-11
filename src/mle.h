/* Main file for prediction
 * includes functions for reading in training parameters (build_training_data)
 * and reading in duplicate probe information (build_duplicate_probes)
 */



#ifndef _MLE_H_
#define _MLE_H_


#include <map>
#include "Normalize.h"
#include "defines.h"

// training data
string train = TRAINING_DIR + "coef.txt";
string repeats = TRAINING_DIR + "dup_probes";

float probe_drops[25];
map< string, float > triplets;
float pm[25]; 
float mm[25]; 
float parameters[72][64];
//float dup_probes[3000][3000];
float **dup_probes;


//prototypes
void build_training_data( string file  );
void build_duplicate_probes( string file );
int get_col( string triplet );
int get_row( int pos );

typedef pair< string, float > entry;

void mle( vector< probe > ref_probes, vector< probe > probes, bool header, string strand_option ){
	
	// initialize variables
	int probe_num = 0; int probe_pos = 0; int forward_probe_pos = 0;
	float drop = 0; float total_prob = 0; int num_probes;
	double prob = 0; float ref_variance = 0;  float ref_mismatch;
	float intensity = 0; float ref_intensity = 0; float expected_intensity = 0;
	string trip; int poly_cnt; int nonpoly_cnt; string vote; int gc_cnt;
	int x; int y; int dup = 0;

	//initialize output file
	string out = "results_" + (probes.at(0)).chrom;
	ofstream outfile( out.c_str() );
	
	//print header for first chromosome
	if( header )
		outfile << "pos score #probes #poly_probes #non_poly_probes #duplicated_probes" << endl;
		
	//print chromosome
	outfile << (probes.at(0)).chrom << endl;
	

	//initialize full output file for additional information/debugging
	/* string full = "full_results_" + (probes.at(0)).chrom;
	*  ofstream fulloutfile( full.c_str() );
	*  fulloutfile << (probes.at(0)).chrom << endl;
	*/
	
	build_training_data( train );
	build_duplicate_probes( repeats );
	
	vector< float > variances;
	vector< float > probs;
	
	// traverse through one base pair at a time, predict liklihood of SNP at each base
	for( int i = 0; i < probes.at(probes.size()-1).position+25; i++ ){
		
		// get all probes in which basepair is included and remember the position in that probe
		probe_num -= 6;
		total_prob = 0;
		num_probes = 0;
		poly_cnt = 0;
		nonpoly_cnt = 0;
		gc_cnt = 0;
		
		// correct for first 5 probes
		if( probe_num < 6 )
			probe_num = 0;
		
		//for all the possible probes that could overlap that base pair
		while( probes.at(probe_num).position <= i && probe_num < probes.size()-1 ){
			
			x = probes.at(probe_num).xcoor;
		  	y = probes.at(probe_num).ycoor;
		  	
			probe_pos = i - probes.at(probe_num).position;
			forward_probe_pos = 24 - probe_pos;
			
			if( strand_option == "forward" )
				probe_pos = forward_probe_pos;
			
			intensity = probes.at(probe_num).mean;
			
			ref_intensity = ref_probes.at(probe_num).mean;
			ref_mismatch = ref_probes.at(probe_num).mm_mean;
			ref_variance = ref_probes.at(probe_num).variance;
			
			prob = 0;
			
			//if a probe includes the base pair
			if( i >= probes.at(probe_num).position && i < probes.at(probe_num).position + 25){ 
			
				assert( probes.at(probe_num).xcoor == ref_probes.at(probe_num).xcoor);
				assert( probes.at(probe_num).ycoor == ref_probes.at(probe_num).ycoor);
				num_probes++;
				
				//get flanking triplet
				if( probe_pos > 0 && probe_pos < 24 )
					trip = (ref_probes.at(probe_num).seq).substr( probe_pos - 1, 3 );
				else
					trip = "AAA";
				
				//get GC count
				gc_cnt = 0;
				for( int j = 0; j <= 24; j++ ){
					string nuc = (ref_probes.at(probe_num).seq).substr( j, 1 );
					if(	nuc == "G" || nuc == "C" )
						gc_cnt++;
				}
			
				int col = get_col( trip );
				int row = get_row( probe_pos );
				
				//estimate intensity drop from linear model coefficients
				drop = parameters[row][col]
				 	+(parameters[row+9][col]*gc_cnt)
					+(parameters[row+18][col]*ref_mismatch)
					+(parameters[row+27][col]*ref_intensity)
					+(parameters[row+36][col]*gc_cnt*ref_mismatch)
					+(parameters[row+45][col]*gc_cnt*ref_intensity)
					+(parameters[row+54][col]*ref_intensity*ref_mismatch)
					+(parameters[row+63][col]*ref_intensity*ref_mismatch*gc_cnt);
				

				expected_intensity = ref_intensity + drop;
				float poly = (intensity - expected_intensity) * (intensity - expected_intensity);
				float nonpoly = (intensity - ref_intensity) * (intensity - ref_intensity);
				float prob = nonpoly - poly;

				
				// compute probes standard deviation from reference mean
				float test = ( ( ref_intensity - intensity ) / sqrt( ref_variance ) );
				
				// don't allow probes far enough away from reference mean to give detrimental results					
				if(  test > 2 && prob < 0 ){
					prob = 0;
				}
				
				//if a bad spot on hybe zero out within coordinates
				//if( x > X_MIN && x < X_MAX && y > Y_MIN && y < Y_MAX )
					//prob = 0;
					
				if( dup_probes[x][y] == 1 )
					dup++;	
				
				probs.push_back( prob );
				variances.push_back( ref_variance );
				
				//print to full output file DEBUG		
				//fulloutfile << "PRE: " << drop << " " << intensity << " " << ref_intensity << " " << expected_intensity << " " << ref_variance << " " << test << endl;
			}
			
			probe_num++;
		
		}
		
		// for all sets of probes compute log ratio ( poly - nonpoly ) / 2 * ref_variance
		if ( num_probes > 0 ){
			float cutoff = var_cutoff( variances );
	        
			for( int j = 0; j < num_probes; j++ ){
				
				prob = 0;
				
				float var = variances.at(j);
				
				if( var < cutoff )
					var = cutoff;

				
				prob = probs.at(j) / ( 2 * var );
				
				//arbitrary voting scheme
				if( prob > 2 )
					poly_cnt++;
					
				if( prob < -2 )
					nonpoly_cnt++;

				total_prob += prob;
				
				// DEBUG
				//fulloutfile << "POST: " << var << " " << prob << " " << total_prob << endl;
			}
		}
			       		
		outfile << i << " " << (1/log(10.0)) * total_prob << " " << num_probes << " " << poly_cnt << " " << nonpoly_cnt << " " << dup << endl;
		dup = 0;
		
		// DEBUG
		//fulloutfile << i << " " << (1/log(10.0)) * total_prob << " " << num_probes << " " << poly_cnt << " " << nonpoly_cnt << endl;
		
		probs.clear();
		variances.clear();
		
	}

	outfile.close();
	
	// DEBUG
	//fulloutfile.close();
}

//build training data structures from linear model coefficients
void build_training_data( string file ){
	
	ifstream data( file.c_str() );
	int line_num = 0;
	int next; int start; int i;
	string first; string second; string line; 
	
	if( !data.is_open() ){
			cout << "ERROR OPENING TRAINING DATA FILE!!! " << file << endl;
			exit(1);	
	}
	
	while( !data.eof() ){
			
			//get each line of data file
			getline( data, line );
			
			//initialize variables
			next = 0; start = 0; i = 0;			
			
			//tokenize by tab
			while( line.find( "\t", start ) != string::npos ){
				
				next = line.find( "\t", start + 1 );		
				string element = line.substr( start, next - ( start ) );

				parameters[line_num][i] = atof( element.c_str() );				
				start = next + 1;							
				i++;			
				
			}
			
			//add last element
			next = line.find( "\t", start + 1 );		
			string element = line.substr( start, next - ( start ) );
			parameters[line_num][i] = atof( element.c_str() );	
			
			line_num++;
	}
}

//read in number of times each probe is duplicated on the array
void build_duplicate_probes( string file ){
	
	ifstream data( file.c_str() );
	int line_num = 0;
	int next; int start; int i; int x; int y;
	string first; string second; string line; 
	
	dup_probes = (float **) calloc(3000, sizeof(float *));
	if (dup_probes == NULL) {
	  cout << "OUT OF MEMORY!!! " << file << endl;
	  exit(1);
	}
	for( int j = 0; j < 3000; j++ ){
	        dup_probes[j] = (float *) calloc(3000, sizeof(float));
	        if (dup_probes[j] == NULL) {
		  cout << "OUT OF MEMORY!!! " << file << endl;
		  exit(1);
		}
		for( int k = 0; k < 3000; k++ ){
			dup_probes[j][k] = 0;		
		}
	}
	
	if( !data.is_open() ){
			cout << "ERROR OPENING DUPLICATE PROBE FILE!!! " << file << endl;
			exit(1);	
	}
	
	while( !data.eof() ){
			
			//get each line of data file
			getline( data, line );
			
			//initialize variables
			next = 0; start = 0; i = 0;			
			x = 0; y = 0;
			
			//tokenize by tab
			while( line.find( "\t", start ) != string::npos ){
				
				next = line.find( "\t", start + 1 );		
				string element = line.substr( start, next - ( start ) );
				if( i == 0 )
					x = atoi( element.c_str() );	
				if( i == 1 )
					y = atoi( element.c_str() );
							
				start = next + 1;							
				i++;			
				
			}
			
			//add last element
			dup_probes[x][y] = 1;
	}
}

// get corresponding training row from nucleotide position
int get_row( int pos ){
	
	int row;
	
	if( pos == 0 || pos == 24 ){ row = 0; }
	if( pos == 1 || pos == 23 ){ row = 1; }
	if( pos == 2 || pos == 22 ){ row = 2; }
	if( pos == 3 || pos == 21 ){ row = 3; }
	if( pos == 4 || pos == 20 ){ row = 4; }
	if( pos == 5 || pos == 19 ){ row = 5; }
	if( pos == 6 || pos == 18 ){ row = 6; }
	if( pos == 7 || pos == 8 || pos == 16 || pos == 17 ){ row = 7; }
	if( pos >= 9 && pos <= 15 ){ row = 8; }
	
	return row;	
	
}

// get corresponding training column from triplet sequence
int get_col( string triplet ){
	
	int col;
	if( triplet == "AAA" ){ col = 0; }
	if( triplet == "AAC" ){ col = 1; }
	if( triplet == "AAG" ){ col = 2; }
	if( triplet == "AAT" ){ col = 3; }
	if( triplet == "ACA" ){ col = 4; }
	if( triplet == "ACC" ){ col = 5; }
	if( triplet == "ACG" ){ col = 6; }
	if( triplet == "ACT" ){ col = 7; }
	if( triplet == "AGA" ){ col = 8; }
	if( triplet == "AGC" ){ col = 9; }
	if( triplet == "AGG" ){ col = 10; }
	if( triplet == "AGT" ){ col = 11; }
	if( triplet == "ATA" ){ col = 12; }
	if( triplet == "ATC" ){ col = 13; }
	if( triplet == "ATG" ){ col = 14; }
	if( triplet == "ATT" ){ col = 15; }
	if( triplet == "CAA" ){ col = 16; }
	if( triplet == "CAC" ){ col = 17; }
	if( triplet == "CAG" ){ col = 18; }
	if( triplet == "CAT" ){ col = 19; }
	if( triplet == "CCA" ){ col = 20; }
	if( triplet == "CCC" ){ col = 21; }
	if( triplet == "CCG" ){ col = 22; }
	if( triplet == "CCT" ){ col = 23; }
	if( triplet == "CGA" ){ col = 24; }
	if( triplet == "CGC" ){ col = 25; }
	if( triplet == "CGG" ){ col = 26; }
	if( triplet == "CGT" ){ col = 27; }
	if( triplet == "CTA" ){ col = 28; }
	if( triplet == "CTC" ){ col = 29; }
	if( triplet == "CTG" ){ col = 30; }
	if( triplet == "CTT" ){ col = 31; }
	if( triplet == "GAA" ){ col = 32; }
	if( triplet == "GAC" ){ col = 33; }
	if( triplet == "GAG" ){ col = 34; }
	if( triplet == "GAT" ){ col = 35; }
	if( triplet == "GCA" ){ col = 36; }
	if( triplet == "GCC" ){ col = 37; }
	if( triplet == "GCG" ){ col = 38; }
	if( triplet == "GCT" ){ col = 39; }
	if( triplet == "GGA" ){ col = 40; }
	if( triplet == "GGC" ){ col = 41; }
	if( triplet == "GGG" ){ col = 42; }
	if( triplet == "GGT" ){ col = 43; }
	if( triplet == "GTA" ){ col = 44; }
	if( triplet == "GTC" ){ col = 45; }
	if( triplet == "GTG" ){ col = 46; }
	if( triplet == "GTT" ){ col = 47; }
	if( triplet == "TAA" ){ col = 48; }
	if( triplet == "TAC" ){ col = 49; }
	if( triplet == "TAG" ){ col = 50; }
	if( triplet == "TAT" ){ col = 51; }
	if( triplet == "TCA" ){ col = 52; }
	if( triplet == "TCC" ){ col = 53; }
	if( triplet == "TCG" ){ col = 54; }
	if( triplet == "TCT" ){ col = 55; }
	if( triplet == "TGA" ){ col = 56; }
	if( triplet == "TGC" ){ col = 57; }
	if( triplet == "TGG" ){ col = 58; }
	if( triplet == "TGT" ){ col = 59; }
	if( triplet == "TTA" ){ col = 60; }
	if( triplet == "TTC" ){ col = 61; }
	if( triplet == "TTG" ){ col = 62; }
	if( triplet == "TTT" ){ col = 63; }
	
	return col;		
}

#endif //_MLE_H_
