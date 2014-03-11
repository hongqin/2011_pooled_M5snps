

#ifndef _RMA_NORM_H_
#define _RMA_NORM_H_

#include <map>


//prototype
vector<  probe > average( vector< vector< probe > > ref_avg );


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


vector<  probe > average( vector< vector< probe > > ref_avg ){
	
	float sum; float mm_sum;
	float avg;
	vector< probe > ref;
	probe p;
	int size = (ref_avg.at(0)).size();
	//cout << "SIZE: " << size << endl;
	
	for( int j = 0; j < size; j++ ){
		sum = 0;
		mm_sum = 0;
		avg = 0;
		//cout << "size: " << ref_avg.size() << endl;
		for( int i = 0; i < ref_avg.size(); i++ ){
		  //cout << ((ref_avg.at(i)).at(j)).mean << endl;
		  sum += ((ref_avg.at(i)).at(j)).mean;
		  mm_sum += ((ref_avg.at(i)).at(j)).mm_mean;
		}				
		
		p = (ref_avg.at(0)).at(j);
		p.mean = sum/ref_avg.size();
		p.mm_mean = mm_sum/ref_avg.size();
		ref.push_back(p);
		//cout << p.xcoor << " " << p.ycoor << " " << p.chrom << " " << p.position << " " << p.mean << endl;
		//cout << endl << endl;
	}	
	
	return ref;	
} 

#endif //_RMA_NORM_H_
