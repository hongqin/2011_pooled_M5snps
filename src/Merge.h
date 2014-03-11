/* Combines CEL file and bpmap file based on array grid coordinates
 */

#ifndef _MERGE_H_
#define _MERGE_H_


#include <iostream>
#include <fstream>

using namespace std;

//merge CEL file and bpmap file
vector< probe > merge( string bpmap_file, string cel_file ){

	cout << "MERGE FILES: " << bpmap_file << " " << cel_file << endl;

	ifstream bpmap( bpmap_file.c_str() );
	ifstream cel( cel_file.c_str() );
	string line;
	int x; int y; int pos; int line_num = 0; int data_line = 0;
	string chr; string seq; 
	int next; int start; int i;	float mean;
	probe ** probes = new probe*[2600];
	vector< probe > all_probes;
	string version = "Version=3";
	
	for( int i = 0; i < 2600; i++ ){
		probes[i] = new probe[2600];	
	}
	if( !cel.is_open() ){
		cout << "ERROR OPENING DATA FILE FOR MERGE!!! " << cel_file << endl;
		exit(1);	
	}
	
	while( !cel.eof() ){
			
			//get each line of data file
			getline( cel, line );
			

			//initialize variables
			next = 0; start = 0; i = 0;			
			x = 0; y = 0; pos = 0; mean = 0;
			chr = ""; seq = "";
			
			if( line_num == 1 ){
				string test = line.substr( 0, 9 );
				if( test != version ){
					cout << "ERROR IN VERSION OF CEL FILE, SNPSCANNER REQUIRES VERSION 3 " << test << endl;
					exit(1);
				}	
				
			}
			
			//tokenize by tab
			while( line.find( "\t", start ) != string::npos ){
					

				next = line.find( "\t", start + 1 );		
				string element = line.substr( start, next - ( start ) );
				if( i == 0 )
					x = atoi( element.c_str() );	
				if( i == 1 )
					y = atoi( element.c_str() );
				if( i == 2 )
					mean = atof( element.c_str() );
									
				start = next + 1;							
				i++;			

			}
			next = line.find( "\t", start + 1 );		
			string element = line.substr( start, next - ( start ) );
	
			if( i == 4 ){
				if( data_line > 0 ){
					//cout << x << " " << y << " " << mean << endl;
					probe p;
					p.xcoor = x;
					p.ycoor = y;
					p.mean = log2(mean);
					probes[x][y] = p;
				}
				data_line++;
			}
			line_num++;
	}
	
	if( !bpmap.is_open() ){
		cout << "ERROR OPENING DATA FILE FOR MERGE!!! " << bpmap << endl;
		exit(1);	
	}
	
	line_num = 0;
	while( !bpmap.eof() ){
			
			//get each line of data file
			getline( bpmap, line );
			
			//initialize variables
			next = 0; start = 0; i = 0;			
			x = 0; y = 0; pos = 0;
			chr = ""; seq = "";
			
			//tokenize by tab
			while( line.find( "\t", start ) != string::npos ){
					
				next = line.find( "\t", start + 1 );		
				string element = line.substr( start, next - ( start ) );
				if( i == 0 )
					x = atoi( element.c_str() );	
				if( i == 1 )
					y = atoi( element.c_str() );
				if( i == 2 )
					chr = element;
				if( i == 3 )
					pos = atoi( element.c_str() );
				
				start = next + 1;							
				i++;			

			}
			next = line.find( "\t", start + 1 );		
			string element = line.substr( start, next - ( start ) );
			if( i == 4 )
				seq = element;	
			
			if( i == 4 ){
					if( line_num > 0 ){
						probe p = probes[x][y];
						assert( x == p.xcoor );
						assert( y == p.ycoor );
						p.chrom = chr;
						p.position = pos;
						p.seq = seq;
						//cout << "PM: " << p.xcoor << "\t" << p.ycoor << "\t" << chr << "\t" << pos << "\t" << seq << "\t" << p.mean << "\t0\t0" << endl; 
						//cout << "MM: " << p.xcoor << "\t" << p.ycoor << "\t" << chr << "\t" << pos << "\t" << seq << "\t" << p.mean << "\t0\t0" << endl; 				
						p.mm_mean = probes[x+1][y].mean;
						p.variance = 0;
						p.mm_variance = 0;
							//filter out controls
						if( pos < 1700000 && chr.substr(0,3) == "chr" ) //|| chr == "mitochondrion" )
							all_probes.push_back(p);
					}
				line_num++;
			}
	}

	bpmap.close();
	cel.close();
	for( int i = 0; i < 2600; i++ ){
		delete [] probes[i];	
	}
	delete [] probes;
	
	return all_probes;
}



#endif //_MERGE_H_
