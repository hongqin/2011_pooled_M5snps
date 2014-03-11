/* File for performing normalization on data
 * generalizable to use several methods
 */

#ifndef _NORMALIZE_H_
#define _NORMALIZE_H_

#include "defines.h"

//prototypes
float get_norm( vector< probe > probes, string type, string match );
float var_cutoff( vector< float > variances );
vector< probe > variance_correction( vector< probe > probes, vector< float > var, char * type );
vector<  probe > average( vector< vector< probe > > ref_avg );
void print_probes( vector< probe > probes, string file );
map< int, double > read_norm_adj();
vector< vector< probe > > RMAnormalize( vector< vector< probe > > ref_probes  );

//perform normalization procedures
vector< probe > normalize( vector< probe > probes, ofstream& debug ){
	
	//naive normalization procedure, divide by median or 90%
	float median = get_norm( probes, "med", "PM" );
	
	map< int, double > adj = read_norm_adj();
	
	int min_mean = (adj.begin())->first;
	float min = (adj.begin())->second;
	int max_mean = (++adj.end())->first;
	float max = (++adj.end())->second;
	int norm_probes_cnt = 0;
		
	for( unsigned int i = 0; i < probes.size(); i++ ){
		
		double mean = ( probes.at(i) ).mean;
		double mm_mean = ( probes.at(i) ).mm_mean;
		int adj_mean = int( mean * 1000 );
					  
		if( adj_mean < min_mean )
			adj_mean = int( min_mean );	
		else
			if( adj_mean > max_mean )
		  		adj_mean = int( max_mean );
			else{
				map< int, double >::const_iterator it = adj.lower_bound(adj_mean);
				adj_mean = it->first;
				//cout << i << " " << mean << " " << adj_mean << " " << it->first << " " << adj[it->first] << endl;
			}
			
		//cout << "OVER/UNDER: " << mean << " " << adj_mean << endl;
		//cout << "END " << mean << " " << adj_mean << " " << adj[adj_mean] << endl; 
		( probes.at(i) ).mean = mean + adj[adj_mean];
		( probes.at(i) ).mm_mean = mm_mean + adj[adj_mean];		
		
		// naive normalization
		//( probes.at(i) ).mean = ( probes.at(i) ).mean / norm;
		//( probes.at(i) ).mm_mean = ( probes.at(i) ).mm_mean  / norm;	
		
	}
	debug << "MIN: " << min_mean/1000.0 << " Normalized Difference: " << min << endl;
	debug << "MAX: " << max_mean/1000.0 << " Normalized Difference: " << max << endl;
	debug << "Number of probes total: " << probes.size() << endl;
	debug << "Median perfect match intensity: " << median << endl;
	
	
	//output normalized probes
	//print_probes( probes, "test_out" );
	return probes;
}


// Use Robust multi-array normalization on all chips
vector< vector< probe > > RMAnormalize( vector< vector< probe > > ref_probes  ){
		
	typedef pair< float, int > entry;
	typedef pair< int, float > norm_entry;
	int num_probes = (ref_probes.at(0)).size();
	int num_arrays = ref_probes.size();
	
	vector< multimap< float, int > > probe_map;
	
	for( int i = 0; i < num_arrays; i++ ){
		multimap< float, int > mm;
		probe_map.push_back( mm );
	}
	
	cout << "ADDING MEANS..." << endl;
	for( int i = 0; i < num_probes; i++ )
		for( int j = 0; j < num_arrays; j++ )
			(probe_map.at(j)).insert( entry( ( (ref_probes.at(j)).at(i)).mean, i )  );	
	
	
	cout << "CREATING ITERATORS..." << endl;
	vector< map< float, int >::const_iterator > it;
	
	for( int i = 0; i < num_arrays; i++ ){
		map< float, int >::const_iterator iter;
		it.push_back( iter );
	}
	
	for( int i = 0; i < num_arrays; i++ )
		it.at(i) = (probe_map.at(i)).begin();
	
	float sum = 0;
	float avg = 0;
	
	vector< map< int, float > > norm_map;
	
	for( int i = 0; i < num_arrays; i++ ){
		map< int, float > map;
		norm_map.push_back( map );
	}
	

	cout << "QUANTILE NORMALIZING" << endl;
	for( int i = 0; i < num_probes; i++ ){		
		sum = 0;	
		for( int j = 0; j < num_arrays; j++ ){
			sum += it.at(j)->first;
		}
		
		avg = sum/probe_map.size();
		
		for( int j = 0; j < num_arrays; j++ ){
			(norm_map.at(j)).insert( norm_entry( it.at(j)->second, avg ) );
		//	cout << j << " " << it.at(j)->second << " " << avg << endl;
			it.at(j)++;
		}

	}
	
	
	cout << "CREATING ITERATORS..." << endl;
	vector< map< int, float >::const_iterator > itr;
	
	for( int i = 0; i < num_arrays; i++ ){
		map< int, float >::const_iterator iter;
		itr.push_back( iter );
	}	
	
	for( int i = 0; i < num_arrays; i++ )
		itr.at(i) = (norm_map.at(i)).begin();
		
	
	cout << "FINISHING... " << num_probes << endl;
	for( int i = 0; i < num_probes; i++ ){	
		for( int j = 0; j < num_arrays; j++ ){
			((ref_probes.at(j)).at(i)).mean = itr.at(j)->second;
			itr.at(j)++;
		}	
	}
	
	return ref_probes;
	
}


//read in previously calculated normalizing constants
map< int, double > read_norm_adj( ){
	
		typedef pair< int, double > entry;
		map< int, double > adj;
		string file = HOME_DIR + "norm_out";
	
		ifstream data( file.c_str() );
		string line; string key; string element;
		int next; int start; int i;			
		int line_num = 0; 
		int first; double second; double prev;
		
		if( !data.is_open() ){
			cout << "ERROR OPENING NORMALIZED DATA FILE!!! " << file << endl;
			exit(1);	
		}	
		
		while( !data.eof() ){
			
			//get each line of data file
			getline( data, line );
			
			//initialize variables
			next = 0; start = 0; i = 0;			
			
			//tokenize by tab
			while( line.find( " ", start ) != string::npos ){
				
				next = line.find( " ", start + 1 );		
				element = line.substr( start, next - ( start ) );

				if( i == 1 )
					key = element;
							
				start = next + 1;							
				i++;			
				
			}
			
			//add last element
			next = line.find( " ", start + 1 );		
			element = line.substr( start, next - ( start ) );
			first = int ( atof( key.c_str() ) * 1000 );
			second = atof( element.c_str() );
			
			if( prev != first && line_num > 0 ){
				//cout << first << " " << second << " " << adj.count( first ) <<  endl;
				adj.insert( entry( first, second ) );
				//cout << first << " " << adj[first] << endl;
			}
			prev = first;
			line_num++;
		}
	
	return adj;
}

//print out probes for debugging or further examination
void print_probes( vector< probe > probes, string file ){

	string output = HOME_DIR + file;
	
	ofstream outfile( output.c_str() );
	cout << "WRITING TO.... " << output << endl;
	for( unsigned int i = 0; i < probes.size(); i++ ){
		probe p = probes.at(i);
		outfile << "PM:\t" << p.xcoor << "\t" << p.ycoor << "\t" << p.chrom << "\t" << p.position << "\t" << p.seq << "\t" << p.mean << "\t" << p.variance << "\t9" << endl;
		outfile << "MM:\t" << p.xcoor + 1 << "\t" << p.ycoor << "\t" << p.chrom << "\t" << p.position << "\t" << p.seq << "\t" << p.mm_mean << "\t" << p.mm_variance << "\t9" << endl;		
	}

	outfile.close();
}

//get the variable to normalize the arrays by either median or 90%
float get_norm( vector< probe > probes, string type, string match ){
	
	vector< float > means;
	for( unsigned int i = 0; i < probes.size(); i++ ){
		if( match == "PM" )
			means.push_back( (probes.at(i)).mean );
		if( match == "MM" )
			means.push_back( (probes.at(i)).mm_mean );
	}
	
	sort( means.begin(), means.end() );
	
	float norm;
	
	if( type == "med" ){
		
		int odd = (int) ( probes.size() / 2 );
		int even = ( (int) ( probes.size() / 2 ) ) - 1;
		
		if( means.size() % 2 == 0 )
			norm = ( means.at( even ) + means.at( odd ) ) / 2;
		else
			norm = means.at( odd );
	}
	else{
		norm = means.at( (int) ( probes.size() * .9 ) );
	}
	
	return norm;
}


// variance correction, take the smallest variance and replace it with the second smallest to eliminate
// a single dominant probe
float var_cutoff( vector< float > variances ){
	
	sort( variances.begin(), variances.end() );
	
	float var;
	
	int size = variances.size();
	
	if( size > 1 )
		var = variances.at(1);
	else
		var = variances.at(0);
		
	return var;
}

//take all reference arrays, average them and compute the variance between the same probes
vector<  probe > average( vector< vector< probe > > ref_avg, char * type ){
	
	float sum; float mm_sum;
	float var; float mm_var;
	vector< probe > ref;
	probe p;
	int size = ( ref_avg.at( 0 ) ).size();
	int num_ref = ref_avg.size();
	float mean; float mm_mean;
	float median; float mm_median;
	float intensity; float mm_intensity;
	float min; float mm_min;
	vector< float > variances;
	vector< probe > probes;
	
	typedef pair< float, float > entry;
	map< float, float > min_intensity;
	
	cout << "SIZE: " << size << endl;
	
	//build pooled variances map
	for( int j = 0; j < size; j++ ){
		sum = mm_sum = 0;
		mean = mm_mean = 0;
		var = mm_var = 0;
		min = mm_min = 100000;
		
		//compute mean for each probe
		for( int i = 0; i < num_ref; i++ ){
			probe x = (ref_avg.at(i)).at(j);
			probes.push_back( x );
			
			sum += x.mean;
			mm_sum += x.mm_mean;
			
			if( x.mean < min ) 
				min = x.mean;
			if( x.mm_mean < mm_min ) 
				mm_min = x.mm_mean;		
		}
		
		mean = sum / num_ref;		
		mm_mean = mm_sum / num_ref;

		//compute variance for each probe
		for( int i = 0; i < num_ref; i++ ){
			intensity = ( probes.at(i) ).mean;
			mm_intensity = ( probes.at(i) ).mm_mean;
			
			var += ( intensity - mean ) * ( intensity - mean );
			mm_var += ( mm_intensity - mm_mean ) * ( mm_intensity - mm_mean );
		}
		
		probes.clear();
		
		//build sorted map of minimum intensities and variances;
		min_intensity.insert( entry( min, ( var / ( num_ref - 1 ) ) ) );
		
	}
	
	cout << "SIZE: " << size << endl;
	for( int j = 0; j < size; j++ ){
		sum = mm_sum = 0;
		mean = mm_mean = 0;
		median = mm_median = 0;
		var = mm_var = 0;
		min = mm_min = 100000;
	
		//compute mean for each probe
		for( int i = 0; i < num_ref; i++ ){
			probe x = ((ref_avg.at(i)).at(j));
			probes.push_back( x );
			
			intensity = x.mean;
			mm_intensity = x.mm_mean;
			
			if( x.mean < min ) 
				min = x.mean;
			if( x.mm_mean < mm_min ) 
				mm_min = x.mm_mean;
		}	
		
		//get median values to use for replicates
		median = get_norm( probes, "med", "PM" );
		mm_median = get_norm( probes, "med", "MM" );	
		probes.clear();
		
		//cout << min << " " << min_intensity[min] << " " << min_intensity.count(min) << endl;
		map< float, float >::const_iterator it = min_intensity.find(min);

		//perform population bases variance estimate
		int cnt = 0;
		while( it != min_intensity.begin() && cnt < 250 ){
			it--;
			sum += it->second;
			cnt ++;
		}
		cnt = 0;
		it = min_intensity.find(min);
		while( it != min_intensity.end() && cnt < 250 ){
			sum += it->second;
			it++;
			cnt++;
		}

		var = sum / 500;
		mm_var = var;
		variances.push_back( var );
		
		
		//create new probe for the reference probes
		p = (ref_avg.at(0)).at(j);
		p.mean = median;
		p.mm_mean = mm_median; 
		p.variance = var;
		p.mm_variance = mm_var;
		ref.push_back(p);
	
	}	
	
	//if( type != "none" )
		//ref = variance_correction( ref, variances, type );
	
	return ref;	
} 


//correct variance to remove really small or really large values, NOT USED IN CURRENT VERSION
vector< probe > variance_correction( vector< probe > probes, vector< float > var, char * type ){
	
	sort( var.begin(), var.end() );
	
	cout << "VARIANCE: " << endl;
	
	if( type == "both" || type == "low" )
		cout << "LOW CUTOFF (5%): " << var.size()*.05 << " " << var.at( int(var.size() * .05) ) << endl;
		
	if( type == "both" || type == "high" )
		cout << "HIGH CUTOFF (10%): " << var.size() * .95 << " " << var.at(int(var.size() * .95)) << endl;
	
	
	float low_cutoff = var.at( int( var.size() * .05 ) );
	float high_cutoff = var.at( int( var.size() * .95 ) );
	
	for( unsigned int i = 0; i < probes.size(); i++ ){
		
		if( type == "both" || type == "low" ){
			//perform small variance correction
			if((probes.at(i)).variance < low_cutoff )
				(probes.at(i)).variance = low_cutoff;
		}
		if( type == "both" || type == "high" ){
			//perform large variance correction
			if((probes.at(i)).variance > high_cutoff )
				(probes.at(i)).variance = high_cutoff;
		}
	}
	
	return probes;
	
}

#endif //_NORMALIZE_H_
