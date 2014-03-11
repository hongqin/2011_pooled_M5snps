/* File that runs SNPScanner
 * Reads in all CEL files for both training and prediction
 * All training files are currently hard coded
 */


#include <iostream>			
#include <assert.h>
#include "SNPScanner.h"
#include "mle.h"
#include "Normalize.h"
#include "defines.h"
#include "Merge.h"


using namespace std;

int main( int argc, char *argv[] ){
   
	//inputs
	if( argc < 2 ){
		cout << "INCORRECT NUMBER OF ARGUMENTS!!!" << endl;
		cout << "USAGE: ./tiling [cel] [yes/no] [SI/RMA (optional)]" << endl;	
		exit(1);
	}
	
	
	bool ref = false;
	string arg = argv[1];
	string cel = CEL_DIR + arg;
	string norm = NORM_DIR + "norm_" + argv[1];
	string bpmap_file = TRAINING_DIR + "bpmap.txt";
	string debug_file = arg + "_debug";
	ofstream debug( debug_file.c_str() );
	
	string norm_opt = "SI";
	string strand_opt = "reverse";
	
	if( argc > 0 ){
		string ref_opt = argv[2];
		norm_opt = argv[3];
		strand_opt = argv[4];
		cout << "Normalization option: " << norm_opt << endl;
		cout << "Strand option: " << strand_opt << endl;
		if( ref_opt == "yes" )
			ref = true;
	}
	
	if( strand_opt == "forward" )
		bpmap_file = TRAINING_DIR + "bpmap_forward.txt";
	
	//if computing reference and training parameters
	if( ref ){
		
		//reference vectors
		vector< probe > ref_norm;
	
		//training vectors
		vector< probe > train_norm;    	
    	
		{//memory management
			vector< string > ref_files;
			
			//store reference files HARD-CODED!! change these back to retrain with DNAse prepared hybes
			ref_files.push_back( REF_DIR + "Klenow_FY3_2ug.CEL");
			ref_files.push_back( REF_DIR + "Klenow_FY3_2ug_2.CEL");
			ref_files.push_back( REF_DIR + "Klenow_FY3_2ug_3.CEL");
			ref_files.push_back( REF_DIR + "Klenow_FY3can1-2.CEL");
			ref_files.push_back( REF_DIR + "Klenow_FY3can1-4.CEL");
			ref_files.push_back( REF_DIR + "Klenow_FY3fcy1-4.CEL");
			ref_files.push_back( REF_DIR + "Klenow_FY3gap1-2.CEL");
			//ref_files.push_back( REF_DIR + "FY3.CEL" );
			//ref_files.push_back( REF_DIR + "FY3-2.CEL" );
			//ref_files.push_back( REF_DIR + "FY3-3.CEL" );
			//ref_files.push_back( REF_DIR + "FY3-4.CEL" );
			//ref_files.push_back( REF_DIR + "FY3-5.CEL" );
			//ref_files.push_back( REF_DIR + "FY3-6.CEL" );
			//ref_files.push_back( REF_DIR + "FY3Can1-1.CEL" );
			//ref_files.push_back( REF_DIR + "FY3Can1-2.CEL" );
			//ref_files.push_back( REF_DIR + "FY3Can1-3.CEL" );
			//ref_files.push_back( REF_DIR + "FY3Can1-4.CEL" );
			//ref_files.push_back( REF_DIR + "FY3Fcy1-7.CEL" );
			//ref_files.push_back( REF_DIR + "FY3Gap1-2.CEL" );
			//ref_files.push_back( REF_DIR + "FY3_10ug.CEL" );
			
			int num_ref = ref_files.size();
		
			vector< vector< probe > > ref_probes;
			//read in probes from files
			for( int i = 0; i < num_ref; i++ ){
				cout << "READING REFERENCE PROBES.... " << i << " " << ref_files.at(i) << endl;
				ref_probes.push_back( merge( bpmap_file, (ref_files.at(i)).c_str()) );
			}
			
			//normalize and average reference files (option to re-normalize references)
			cout << "NORMALIZING REFERENCE ARRAYS..." << endl; 
       		cout << "NUM REF: " << num_ref << endl;    	
			cout << "AVERAGING REFERENCE ARRAYS AND COMPUTING VARIANCES..." << endl;
			ref_norm = average( ref_probes, "none" );
		}
		
		//output normalized/averaged reference probes
		print_probes( ref_norm, "REF_NORM" );	
		
		{//memory management
			
			vector< string > train_files;
			vector< vector< probe > > train_probes;
			vector< probe > train;
					
			//store training files HARD-CODED!! change these to retrain using DNase prepared hybes!!
			train_files.push_back( CEL_DIR + "Klenow_RM11_2ug.CEL");
			train_files.push_back( CEL_DIR + "Klenow_RM11_hybe1.CEL");
			train_files.push_back( CEL_DIR + "Klenow_RM11_hybe2.CEL");
			//train_files.push_back( CEL_DIR + "RM11.CEL" );
			//train_files.push_back( CEL_DIR + "RM11-2.CEL" );
			//train_files.push_back( CEL_DIR + "RM11-3.CEL" );
			//train_files.push_back( CEL_DIR + "RM11-4.CEL" );
			//train_files.push_back( CEL_DIR + "RM11-5.CEL" );
			//train_files.push_back( CEL_DIR + "RM11-6.CEL" );
			//train_files.push_back( CEL_DIR + "RM11-7.CEL" );
		
	   		int num_train = train_files.size();
		
			//read in probes from files
     		for( int i = 0; i < num_train; i++ ){
				cout << "READING TRAIN PROBES.... " << i << " " << train_files.at(i) << endl;
				train_probes.push_back( merge( bpmap_file, (train_files.at(i)).c_str() ) );
       		}
       		
       		
			//normalize and average training files  
       		cout << "NORMALIZING TRAINING ARRAYS..." << endl;        
       		cout << "NUM TRAIN: " << num_train << endl;      	
       		cout << "AVERAGING TRAINING ARRAYS AND COMPUTING VARIANCES..." << endl;
       		train = average( train_probes, "none" );
		
			lowess_diff( ref_norm, train, norm, debug );
		 
			system( (" R --no-save --no-restore-data < " + TRAINING_DIR + "lowess_norm.R ").c_str() );
		
			train_norm = normalize( train, debug );
		
		}

		//output normalized/averaged strain probes
		print_probes( train_norm, "RM_NORM" );	
		
		//create RM probe file for R
		system( ("perl " + BIN_DIR + "RM_probe_data.pl rm > " + TRAINING_DIR + "rm_snps").c_str() );

		cout << "COMPUTING TRAINING PARAMETERS..." << endl;
		system( ("R --no-save --no-restore-data < " + TRAINING_DIR + "Training_log2.R").c_str() );
		
		
	}

	vector< vector< probe > > ref_chrom_probes;
	vector< vector< probe > > strain_probes;
	
	{//memory management
		
		vector< probe > ref_normed;
		vector< probe > probes;
		vector< probe > probes_norm;
		
		
		/*******************************************/

		string ref_file = HOME_DIR + "REF_NORM";
		cout << "READING NORMALIZED REFERENCE PROBES.... " << ref_file << endl;
		ref_normed = read_probes( ref_file.c_str(), true );
		cout << "REF NORM SIZE: " << ref_normed.size() << endl;
		
		cout << "READING STRAIN PROBES.... " << cel << endl;
		probes = merge( bpmap_file, cel.c_str() );
		cout << "STRAIN SIZE: " << probes.size() << endl;
		
		//normalize strain file
		cout << "NORMALIZING STRAIN..." << endl;
		
		/*******RMA1***********************/
		/* Uses the median of all the reference arrays and performs
		 * RMA over just the median and the strain to be analyzed
		 */
		if( norm_opt == "RMA1" ){

			vector< vector< probe > > rma_probes;
			rma_probes.push_back( probes );
			rma_probes.push_back( ref_normed );
			vector< vector< probe > > rma_norm = RMAnormalize( rma_probes );
			probes_norm = average( rma_norm , "none" );
		}
		/*********************************/
		
		
		/*******RMA2***********************/
		/* Uses the 5 best FY arrays and performs RMA using those
		 * 5 and the strain to analyzed...REQUIRES ALL 5FY HYBRIDIZATION FILES!!!!!
		 */
		if( norm_opt == "RMA2" ){
			vector< string > rma_files;
			rma_files.push_back( REF_DIR + "FY3-2.CEL" );
			rma_files.push_back( REF_DIR + "FY3-3.CEL" );
			rma_files.push_back( REF_DIR + "FY3-4.CEL" );
			rma_files.push_back( REF_DIR + "FY3-5.CEL" );
			rma_files.push_back( REF_DIR + "FY3-6.CEL" );
			
			vector< vector< probe > > rma_probes;
			
			for( int i = 0; i < rma_files.size(); i++ ){
				rma_probes.push_back( merge( bpmap_file, (rma_files.at(i)).c_str() ) );
			}
			rma_probes.push_back( probes );
			vector< vector< probe > > rma_norm = RMAnormalize( rma_probes );
			probes_norm = average( rma_norm , "none" );
		}
		/*********************************/
		
		
		/*******Set Invariant**************/
		/* default normalization as describe in Gresham et al. 2006 */
		
		if( norm_opt == "SI" ){
		
			lowess_diff( ref_normed, probes, norm, debug );
			system( ("R --no-save --no-restore-data < " + TRAINING_DIR + "lowess_norm.R").c_str() );
			probes_norm = normalize( probes, debug );
		}
		/*********************************/
		
		//output normalized/averaged strain probes
		print_probes( probes_norm, "STRAIN_NORM" );

		//separate probes by chromosome
		cout << "SEPARATING REFERENCE PROBES...." << endl;
		ref_chrom_probes = separate_chrom( ref_normed );
		
		cout << "SEPARATING STRAIN PROBES...." << endl;
		strain_probes = separate_chrom( probes_norm );
		
		/*******************************************/
		
		
		/*
		//if not computing reference and training parameters
		string ref_file = HOME_DIR + "REF_NORM";
		cout << "READING NORMALIZED REFERENCE PROBES.... " << ref_file << endl;
		ref_normed = read_probes( ref_file.c_str(), true );
		
		cout << "REF NORM SIZE: " << ref_normed.size() << endl;
	
		cout << "READING STRAIN PROBES.... " << cel << endl;
		probes = merge( bpmap_file, cel.c_str() );
		cout << "STRAIN SIZE: " << probes.size() << endl;
	
		//normalize strain file
		cout << "NORMALIZING STRAIN..." << endl;
		lowess_diff( ref_normed, probes, norm );
	
		system( ("R --no-save --no-restore-data < " + TRAINING_DIR + "lowess_norm.R").c_str() );
		probes_norm = normalize( probes );
	
		//output normalized/averaged strain probes
		print_probes( probes_norm, "STRAIN_NORM" );

		//separate probes by chromosome
		cout << "SEPARATING REFERENCE PROBES...." << endl;
		ref_chrom_probes = separate_chrom( ref_normed );
		
		cout << "SEPARATING STRAIN PROBES...." << endl;
		strain_probes = separate_chrom( probes_norm );
		*/
	}
	
	
	cout << "RUNNING PREDICTION ALGORITHM..." << endl;
	for( int i = 0; i < 16; i++ ){	
		vector< probe > a;
		vector< probe > b;
		
		a = ref_chrom_probes.at(i);
		b = strain_probes.at(i);

		cout << "PREDICTING... " << (a.at(0)).chrom << endl;
		assert( (a.at(0)).chrom == (b.at(0)).chrom );
		
		//print header for first chromosome
		if( i == 0 )
			mle( a, b, true, strand_opt );
		else
			mle( a, b, false, strand_opt );

		(ref_chrom_probes.at(i)).clear();
		(strain_probes.at(i)).clear();
	}
	
	debug.close();
}

