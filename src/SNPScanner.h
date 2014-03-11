#ifndef _SNPScanner_H_
#define _SNPScanner_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include <deque>
#include <vector>


using namespace std;

struct probe{
	
	int xcoor;
	int ycoor;
	string chrom;
	int position;
	double mean;
	double mm_mean;
    string seq;
	double variance;
	double mm_variance;
	
};

struct diff_probe{
	
	double strain_mean;
	double strain_mm_mean;
	double mean_diff;
	double mm_mean_diff;
	double stddev_diff;
	
};

// read in output
vector< probe > read_probes( const char * file, bool norm ){
	
		ifstream data( file );
		string line;
		int next; int start; int i;			
		int x; int y; int mm_x; int mm_y; int pos;
		string chr; string seq; 
		double mean; double variance; double mm_mean; double mm_variance;
		
		//store affy information in vector by probe (perfect match and mismatch)
		vector< probe > all_probes;
		if( !data.is_open() ){
			cout << "ERROR OPENING DATA FILE!!! " << file << endl;
			exit(1);	
		}	
		
		while( !data.eof() ){
			
			//get each line of data file
			getline( data, line );

			//initialize variables
			next = 0; start = 0; i = 0;			
		
			//process only Perfect Match results for now
			if( line.substr( 0, 2 ) == "PM" ){
				
				x = 0; y = 0; pos = 0;
				chr = ""; seq = "";
				mean = 0; variance = 0;
				
				//tokenize by tab
				while( line.find( "\t", start ) != string::npos ){
					
					next = line.find( "\t", start + 1 );		
					string element = line.substr( start, next - ( start ) );

					if( i == 1 )
						x = atoi( element.c_str() );	
					if( i == 2 )
						y = atoi( element.c_str() );
					if( i == 3 )
						chr = element;
					if( i == 4 )
						pos = atoi( element.c_str() );
					if( i == 5 )
						seq = element;	
					if( i == 6 )
						mean = atof( element.c_str() );	
					if( i == 7 )
						variance = atof( element.c_str() );
							
					start = next + 1;							
					i++;			
				}
		}		

		// get mismatch intensity 
		else{
			
			mm_x = 0; mm_y = 0;
			mm_mean = 0; 
			
			//tokenize by tab
			while( line.find( "\t", start ) != string::npos ){
				
				next = line.find( "\t", start + 1 );		
				string element = line.substr( start, next - ( start ) );
				if( i == 1 )
					mm_x = atoi( element.c_str() );
				if( i == 2 )
					mm_y = atoi( element.c_str() );
				if( i == 6 )
					mm_mean = atof( element.c_str() );	
				if( i == 7 )
					mm_variance = atof( element.c_str() );
	
				start = next + 1;		
				
				i++;
			
			}	
			//account for possible blank lines, usually at the end
			if( line.length() > 0 ){
				
				//filter out controls
				if( pos < 1700000 && chr.substr(0,3) == "chr" ){ //|| chr == "mitochondrion" ){
					
					//ensure correct match and mismatch pairing
					assert( x == mm_x - 1 );
					assert( y == mm_y );
					
					//create new probe, store information and add to all_probes array
					probe p;
					p.xcoor = x;
					p.ycoor = y;
					p.chrom = chr;
					p.position = pos;
					
					//if probes are not already normalized
					if( !norm ){
					  	mean = log2(mean);
					 	mm_mean = log2(mm_mean);
					}
					 
				    p.mean = mean;
					p.mm_mean = mm_mean;
					p.seq = seq;
					p.variance = variance; 
					p.mm_variance = 0;		
					all_probes.push_back(p);
				}
			}
			//cout << mean << " " << mm_mean << " " << mean/mm_mean << endl;
		}
	}
	
	data.close();
	return all_probes;
}

vector< vector< probe > > separate_chrom( vector< probe > probes ){	
	
	vector< probe > chrom_probes;
	vector< vector< probe > > all_probes;
	string prev_chr;
	for( unsigned int i = 0; i < probes.size(); i++ ){
	
		probe p = probes.at(i);
		
		//if probe is from a new chromosome, push last chromosome and start fresh
		if( p.chrom != prev_chr && chrom_probes.size() > 0){
				all_probes.push_back(chrom_probes);
				chrom_probes.clear();	
		}
		
		prev_chr = p.chrom;
		chrom_probes.push_back(p);
	}
	
	//add last chromosome
	all_probes.push_back(chrom_probes);
	
	return all_probes;
}

bool probe_cmp( diff_probe a, diff_probe b ){
	double temp_a = a.stddev_diff; 
	double temp_b = b.stddev_diff;
	
	if( temp_a < 0 )
		temp_a = a.stddev_diff * -1;
	if( temp_b < 0 )
		temp_b = b.stddev_diff * -1;	
	
	return temp_a < temp_b;	
}

//compute the differences between strain probes and reference for lowess normalization
void lowess_diff( vector< probe > ref, vector< probe > strain, string norm, ofstream& debug ){
	
	float num_stddev = 1.5;
	assert( ref.size() == strain.size() );
	ofstream out( "norm_data" );
	ofstream invar_set( norm.c_str() );
	string test_string = norm + "_stddev";
	ofstream test( test_string.c_str() );
	invar_set << "index strain ref xcoor ycoor" << endl;

	
	/*float ref_dist[200];
	float strain_dist[200];
	for( unsigned int i = 0; i < 200; i++ ){
		ref_dist[i] = 0;
		strain_dist[i] = 0;	
	}
	
	for( unsigned int i = 0; i < ref.size(); i++ ){
		probe r = ref.at(i);
		probe s = strain.at(i);
		ref_dist[int(r.mean*10)]++;
		strain_dist[int(s.mean*10)]++;
	}
	
	float ref_mode = 0;
	float ref_mode_max = 0;
	float strain_mode = 0;
	float strain_mode_max = 0;
	
	for( unsigned int i = 0; i < 200; i++ ){
		//cout << float(i)/10 << " " << ref_dist[i] << " " << strain_dist[i] << endl;
		if( ref_dist[i] > ref_mode_max ){
			ref_mode = float(i)/10;
			ref_mode_max = ref_dist[i];	
		}
		if( strain_dist[i] > strain_mode_max ){
			strain_mode = float(i)/10;
			strain_mode_max = strain_dist[i];	
		}
	}
	
	float mode_adjust = strain_mode - ref_mode;
	cout << "MODE DIFFERENCE: " << ref_mode << " " << strain_mode << " " << mode_adjust << endl;
	cout << "MODE ADJUST OFF!!!!!" << endl;
	
	vector< diff_probe > sorted;
	*/
	int x = 0;
	for( unsigned int i = 0; i < ref.size(); i++ ){
		
		probe r = ref.at(i);
		probe s = strain.at(i);
		//diff_probe d;
		
		assert( r.xcoor == s.xcoor );
		assert( r.ycoor == s.ycoor );
		assert( r.chrom == s.chrom );
		assert( r.position == s.position );
		
		//(strain.at(i)).mean = (strain.at(i)).mean - mode_adjust;
		//(strain.at(i)).mm_mean = (strain.at(i)).mm_mean - mode_adjust;
		//s.mean -= mode_adjust;
		//s.mm_mean -= mode_adjust;
		
		
		float stddev = sqrt( r.variance );
		float mm_stddev = sqrt( r.mm_variance );
		
		float max = r.mean + ( stddev * num_stddev );
		float min = r.mean - ( stddev * num_stddev );
		float mm_max = r.mm_mean + ( mm_stddev * num_stddev );
		float mm_min = r.mm_mean - ( mm_stddev * num_stddev );
		
		//cout << r.mean << " " << s.mean << " " << min << " " << max << " " << stddev << endl;
		//cout << r.mm_mean << " " << s.mm_mean << " " << mm_min << " " << mm_max << " " << stddev << endl;
		
		double test_stddev = ( s.mean - r.mean ) / stddev;
		double test_stddev_mm = ( s.mm_mean - r.mm_mean ) / stddev;
		test << test_stddev << endl;
		test << test_stddev_mm << endl;
		
		/*d.strain_mean = s.mean;
		d.strain_mm_mean = s.mm_mean;
		d.mean_diff = r.mean - s.mean;
		d.mm_mean_diff = r.mm_mean - s.mm_mean;
		d.stddev_diff = (r.mean-s.mean)/stddev;
		sorted.push_back(d);
		*/
		
		//include both perfect match and mismatch probes
		if( s.mean < max && s.mean >  min ){
			out << s.mean << " " << r.mean - s.mean << endl;
			out << s.mm_mean << " " << r.mm_mean - s.mm_mean << endl;
			invar_set << x << " " << s.mean << " " << r.mean << " " << s.mean - r.mean << " " << r.xcoor << " " << r.ycoor << endl;
			x++;	
		}
		
		//if( s.mm_mean < max && s.mm_mean >  min )
			//out << s.mm_mean << " " << r.mm_mean - s.mm_mean << endl;	
	}
	
	/*sort( sorted.begin(), sorted.end(), probe_cmp );	
	
	for( unsigned int i = 0; i < 1500000; i++ ){
		out << (sorted.at(i)).strain_mean << " " << (sorted.at(i)).mean_diff << endl; 
		out << (sorted.at(i)).strain_mm_mean << " " << (sorted.at(i)).mm_mean_diff << endl; 
	
	}
	*/
	
	debug << "Number of probes used in normalization: " << x << endl;

	out.close();
	invar_set.close();
	test.close();
}



#endif //_SNPScanner_H_
